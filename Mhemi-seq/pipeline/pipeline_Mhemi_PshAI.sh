#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/15
##-------------------------

echo -e '\nAttention: This script is used only for MspJI cut and 1bp UMI + 3bp PshAI motif (also works for 4bp UMI) adjacent to adapter!\n'

show_usage="Usage: bash Mhemi_MspJI.sh -g reference_genome\n\n
			-g|--reference-genome \t reference genome\n
			-p|--thread \t\t number of thread\n
			-h|--help \t\t get help document\n"

GETOPT_ARGS=`getopt -o g:p:h -al reference-genome:,thread:,help -- "$@"`
eval set -- "$GETOPT_ARGS"

while [ -n "$1" ]
do
        case "$1" in
                -g|--reference-genome) reference_genome=$2; shift 2;;
                -p|--thread) nthread=$2; shift 2;;
                -h|--help) echo -e $show_usage; exit 0;;
				--) break;;
                *) echo -e $1,$2,$show_usage; break;;
        esac
done

if [[ -z ${reference_genome} ]]; then
        echo -e $show_usage
        exit 0
else
	if [[ -z ${nthread} ]]; then
		nthread=8
	fi

	for prefix in `ls ./fastq/*_R1.fastq.gz | sed "s/_R1.fastq.gz//g; s/\.\/fastq\///g" | sort --parallel ${nthread} -u `; do
		if [ -d ${prefix} ]
		then continue
		fi
		mkdir ./${prefix}
		cd ./${prefix}

		echo -e "============ ${prefix} ============"  | tee ${prefix}_Report.txt
		raw=`zcat -c ../fastq/${prefix}_R1.fastq.gz | wc -l`
		nRaw=`echo "$raw/4" | bc`
		echo "Raw: $nRaw (100%)" | tee -a ${prefix}_Report.txt

		# step1. Quality control: remove low quality reads and adapter sequence
		fastp --umi --umi_loc per_read --umi_len 4 --json tmp.${prefix}.json --html tmp.${prefix}.QC.html -q 20 -e 20 -l 30 -w ${nthread} --in1 ../fastq/${prefix}_R1.fastq.gz --out1 tmp.${prefix}.QC.1.fq.gz --in2 ../fastq/${prefix}_R2.fastq.gz --out2 tmp.${prefix}.QC.2.fq.gz 2> /dev/null
		trim_galore --paired --no_report_file --cores 6 --gzip tmp.${prefix}.QC.1.fq.gz tmp.${prefix}.QC.2.fq.gz 2> /dev/null
		mv tmp.${prefix}.QC.1_val_1.fq.gz ${prefix}_trimmed_R1.fastq.gz
		mv tmp.${prefix}.QC.2_val_2.fq.gz ${prefix}_trimmed_R2.fastq.gz

		QC=`zcat -c ${prefix}_trimmed_R1.fastq.gz | wc -l`
		nQC=`echo "$QC/4" | bc`
		pQC=`echo "scale=1; 100*$nQC/$nRaw" | bc`
		echo "Quality control: $nQC ($pQC%)" | tee -a ${prefix}_Report.txt

		# step2. Align clean reads to reference genome

		bowtie2 -p ${nthread} -X 1000 --no-mixed --no-discordant -x /mnt/disk1/5/share/index/bowtie2/${reference_genome} -1 ${prefix}_trimmed_R1.fastq.gz -2 ${prefix}_trimmed_R2.fastq.gz | samtools view -bh -@ ${nthread} > tmp.bam
		samtools view -@ ${nthread} -bh -F 260 tmp.bam > ${prefix}_mapped.bam
		samtools view -@ ${nthread} -bh -f 4 tmp.bam > ${prefix}_unmapped.bam

		mapping=`samtools view -@ ${nthread} ${prefix}_mapped.bam | wc -l`
		nMapping=`echo "$mapping/2" | bc`
		pMap=`echo "scale=1; 100*$nMapping/$nQC" | bc`
		echo "Mapping: $nMapping ($pMap%)" | tee -a ${prefix}_Report.txt

		# step3. filter reads by MAPQ value according to fragment length
		samtools view -h -@ ${nthread} ${prefix}_mapped.bam | awk '/^@/ || sqrt($9*$9)>=100' OFS='\t' | samtools view -hbq 10 -@ ${nthread} | samtools sort -@ ${nthread} -m 2G > tmp.${prefix}.long.bam
		samtools index tmp.${prefix}.long.bam
		samtools view -h -@ ${nthread} ${prefix}_mapped.bam | awk '/^@/ || sqrt($9*$9)<100' OFS='\t' | samtools view -hb -@ ${nthread} | samtools sort -@ ${nthread} -m 2G > tmp.${prefix}.short.bam
		samtools index tmp.${prefix}.short.bam
		samtools merge -f -@ ${nthread} tmp.${prefix}.sorted.bam tmp.${prefix}.long.bam tmp.${prefix}.short.bam
		samtools sort -@ ${nthread} -m 2G tmp.${prefix}.sorted.bam > ${prefix}_q10.bam

		q10=`samtools view -@ ${nthread} ${prefix}_q10.bam | wc -l`
		nQ10=`echo "$q10/2" | bc`
		pQ10=`echo "scale=1; 100*$nQ10/$nMapping" | bc`
		echo "MAPQ10: $nQ10 ($pQ10%)" | tee -a ${prefix}_Report.txt

		samtools view -h -@ ${nthread} ${prefix}_q10.bam | awk '$1~/^@/ || $3!~/[_,L,M]/' | samtools view -hb -@ ${nthread} > tmp.dup.bam
		Contig=`samtools view -@ ${nthread} ${prefix}_q10.bam | awk '$3~/[_,L,M]/' | wc -l `
		nContig=`echo "$Contig/2" | bc`
		pContig=`echo "scale=1; 100*$nContig/$nQ10" | bc`
		echo "Contigs | chrM: $nContig ($pContig%)" | tee -a ${prefix}_Report.txt

		# step4. deduplication by considering fragment position and UMI
		samtools view -H tmp.dup.bam > tmp.header
		samtools view -@ ${nthread} tmp.dup.bam | awk '{print substr($1,1,length($1)-10), $0}' OFS='\t' | cut -f 1,3- | cat tmp.header - | samtools view -hb -@ ${nthread} | samtools sort -m 2G -@ ${nthread} > tmp.dup2.bam
		allDup=`samtools view -@ ${nthread} tmp.dup2.bam | wc -l `
		nDup=`echo "${allDup}/2" | bc`

		picard MarkDuplicates REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT I=tmp.dup2.bam O=tmp.mkdup.bam M=tmp.output.metrics OPTICAL_DUPLICATE_PIXEL_DISTANCE=10000 MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 TAGGING_POLICY=All

		dupPCR=`samtools view -@ ${nthread} tmp.mkdup.bam | grep DT:Z:LB | wc -l `
		nDupPCR=`echo "${dupPCR}/2" | bc`
		pDupPCR=`echo "scale=1; 100*${nDupPCR}/${nDup}" | bc`
		echo "PCR duplicate (picard): ${nDupPCR} (${pDupPCR}%)" | tee -a ${prefix}_Report.txt

		samtools sort -n -m 2G -@ ${nthread} tmp.dup.bam | bedtools bamtobed -bedpe -i - | awk '{umi=substr($7,length($7)-8,9)} {if(($9=="+") && ($10=="-")) print $1, $2, $6, umi, $7; else if(($9=="-") && ($10=="+")) print $4, $5, $3, umi, $7}' OFS='\t' > tmp.unsort.bed
		sort --parallel ${nthread} -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 tmp.unsort.bed > tmp.sorted.bed
		bedtools groupby -i tmp.sorted.bed -g 1,2,3,4 -c 5 -o collapse | sed 's/,.*//g' | cut -f 5 > tmp.dedup.readnames

		grep -w -F -f tmp.dedup.readnames tmp.sorted.bed > tmp.${prefix}_fragment.bed
		samtools view -@ ${nthread} tmp.dup.bam | grep -w -F -f tmp.dedup.readnames - | cat tmp.header - | samtools view -@ ${nthread} -bh | samtools sort -m 2G -@ ${nthread} > ${prefix}_sorted_dedup.bam

		nDedup=`cat tmp.dedup.readnames | wc -l `
		nDup2=`echo "scale=1; ${nDup}-${nDedup}" | bc`
		pDup=`echo "scale=1; 100*$nDup2/$nDup" | bc`
		echo "Duplication (UMI): $nDup2 ($pDup%)" | tee -a ${prefix}_Report.txt


		pDedup=`echo "scale=1; 100*$nDedup/$nRaw" | bc`
		echo "Deduplicated: $nDedup ($pDedup%)" | tee -a ${prefix}_Report.txt

		# step5. Extract locations of reads
		forward1="[A,T,C,G]CG[A,T,C,G][A,G]"
		reverse1="[C,T][A,T,C,G]CG[A,T,C,G]"
		awk '{print $1, $2-15, $3+15}' OFS='\t' tmp.${prefix}_fragment.bed | bedtools getfasta -fi /mnt/disk1/5/share/Reference/fasta/${reference_genome}.fa -bed - | awk '{if(!/^>/) print toupper($0); else print $0}' | sed -n '{N;s/\n/\t/p}' | sed 's/:/\t/g; s/(+)//g; s/>//g; s/-/\t/g;' | awk '{headSeq1=substr($4, 1, 7); headSeq2=substr($4, 28, 7); tailSeq1=substr($4, length($4)-34, 7); tailSeq2=substr($4, length($4)-6, 7)} {if(headSeq1~/'$forward1'/ && headSeq2!~/'$reverse1'/ && tailSeq1~/'$forward1'/ && tailSeq2!~/'$reverse1'/) print "H1T1"; else if(headSeq1!~/'$forward1'/ && headSeq2~/'$reverse1'/ && tailSeq1~/'$forward1'/ && tailSeq2!~/'$reverse1'/) print "H2T1"; else if(headSeq1~/'$forward1'/ && headSeq2!~/'$reverse1'/ && tailSeq1!~/'$forward1'/ && tailSeq2~/'$reverse1'/) print "H1T2"; else if(headSeq1!~/'$forward1'/ && headSeq2~/'$reverse1'/ && tailSeq1!~/'$forward1'/ && tailSeq2~/'$reverse1'/) print "H2T2"; else print "Others"}' > tmp.${prefix}.motif_frag.txt
		awk 'BEGIN{nH1T1=0; nH1T2=0; nH2T1=0; nH2T2=0; nOthers=0} {if($1=="H1T1") nH1T1+=1; else if($1=="H1T2") nH1T2+=1; else if($1=="H2T1") nH2T1+=1; else if($1=="H2T2") nH2T2+=1; else if($1=="Others") nOthers+=1} END{print "End\tCount\nH1T1\t"nH1T1"\nH1T2\t"nH1T2"\nH2T1\t"nH2T1"\nH2T2\t"nH2T2"\nOthers\t"nOthers}' tmp.${prefix}.motif_frag.txt > ${prefix}.EndCount.txt
		paste tmp.${prefix}_fragment.bed tmp.${prefix}.motif_frag.txt | awk '$6!="Others"' > ${prefix}_cut_fragment.bed

		nMotif=`cat ${prefix}_cut_fragment.bed | wc -l`
		pMotif=`echo "scale=1; 100*$nMotif/$nDedup" | bc`
		echo "Motif(MspJI, 1bp wobble): $nMotif ($pMotif%)" | tee -a ${prefix}_Report.txt

		n32bp=`awk '($3-$2)==32' ${prefix}_cut_fragment.bed | wc -l`
		p32bp=`echo "scale=1; 100*$n32bp/$nMotif" | bc`
		echo "32bp: $n32bp ($p32bp%)" | tee -a ${prefix}_Report.txt

		# step6. Calculate methylation frequency according to cytosine location of fragment
		## methylation count of Watson
		awk '{if($6~/H1/) print $1, $2-14, $2-11}' OFS='\t' ${prefix}_cut_fragment.bed > tmp.Wat.bed
		awk '{if($6~/T1/) print $1, $3-18, $3-15}' OFS='\t' ${prefix}_cut_fragment.bed >> tmp.Wat.bed
		sort --parallel ${nthread} -k1,1 -k2,2n tmp.Wat.bed | intersectBed -c -sorted -f 1 -a /mnt/disk1/5/share/Reference/bed/${reference_genome}_CGNR_Wat.bed -b - > tmp.${prefix}.Wat.Meth-Count.bed

		## unmethylation count of Watson
		awk '{print $1, $2, $3-18}' OFS='\t' ${prefix}_cut_fragment.bed > tmp.Wat.bed
		sort --parallel ${nthread} -k1,1 -k2,2n tmp.Wat.bed | intersectBed -c -sorted -f 1 -a /mnt/disk1/5/share/Reference/bed/${reference_genome}_CGNR_Wat.bed -b - > tmp.${prefix}.Wat.unmeth-Count.bed

		paste tmp.${prefix}.Wat.Meth-Count.bed tmp.${prefix}.Wat.unmeth-Count.bed | cut -f 1-4,8 | awk '{if($4+$5>0) print $1, $2, $2+1, $4, $4+$5, "+"}' OFS='\t' > tmp.${prefix}.Wat.Mhemi.bed

		## methylation count of Crick
		awk '{if($6~/H2/) print $1, $2+15, $2+18}' OFS='\t' ${prefix}_cut_fragment.bed > tmp.Cri.bed
		awk '{if($6~/T2/) print $1, $3+11, $3+14}' OFS='\t' ${prefix}_cut_fragment.bed >> tmp.Cri.bed
		sort --parallel ${nthread} -k1,1 -k2,2n tmp.Cri.bed | intersectBed -c -sorted -f 1 -a /mnt/disk1/5/share/Reference/bed/${reference_genome}_CGNR_Cri.bed -b - > tmp.${prefix}.Cri.Meth-Count.bed

		## unmethylation count of Crick
		awk '{print $1, $2+18, $3}' OFS='\t' ${prefix}_cut_fragment.bed > tmp.Cri.bed
		sort --parallel ${nthread} -k1,1 -k2,2n tmp.Cri.bed | intersectBed -c -sorted -f 1 -a /mnt/disk1/5/share/Reference/bed/${reference_genome}_CGNR_Cri.bed -b - > tmp.${prefix}.Cri.unmeth-Count.bed

		paste tmp.${prefix}.Cri.Meth-Count.bed tmp.${prefix}.Cri.unmeth-Count.bed | cut -f 1-4,8 | awk '{if($4+$5>0) print $1, $3-1, $3, $4, $4+$5, "-"}' OFS='\t' > tmp.${prefix}.Cri.Mhemi.bed

		cat tmp.${prefix}.Wat.Mhemi.bed tmp.${prefix}.Cri.Mhemi.bed > tmp.unsort.bed
		sort --parallel ${nthread} -k1,1 -k2,2n tmp.unsort.bed > ${prefix}_CpG.bed

		nCG=`awk '{nSum+=$5} END{print nSum}' ${prefix}_CpG.bed `
		echo "CpG: $nCG" | tee -a ${prefix}_Report.txt

		# step7. Calculate dyad methylation frequency
		awk '{fragLen=$3-$2} {if(fragLen>=40) print $1, $2+18, $3-18}' OFS='\t' ${prefix}_cut_fragment.bed | sort --parallel ${nthread} -k1,1 -k2,2n | intersectBed -f 1 -sorted -c -a /mnt/disk1/5/share/Reference/bed/${reference_genome}_YNCGNR.bed -b - > tmp.unme.bed
		awk '{fragLen=$3-$2} {if(fragLen>=40 && $6=="H1T1") print $1, $3-18, $3-14}' OFS='\t' ${prefix}_cut_fragment.bed | sort --parallel ${nthread} -k1,1 -k2,2n | intersectBed -f 1 -sorted -c -a /mnt/disk1/5/share/Reference/bed/${reference_genome}_YNCGNR.bed -b - > tmp.hmei_W.bed
		awk '{fragLen=$3-$2} {if(fragLen>=40 && $6=="H2T2") print $1, $2+14, $2+18}' OFS='\t' ${prefix}_cut_fragment.bed | sort --parallel ${nthread} -k1,1 -k2,2n | intersectBed -f 1 -sorted -c -a /mnt/disk1/5/share/Reference/bed/${reference_genome}_YNCGNR.bed -b - > tmp.hmei_C.bed
		awk '{fragLen=$3-$2} {if(fragLen<=33 && fragLen>=31 && $6=="H2T1") print $1, $2+14, $2+18}' OFS='\t' ${prefix}_cut_fragment.bed | sort --parallel ${nthread} -k1,1 -k2,2n | intersectBed -f 1 -sorted -c -a /mnt/disk1/5/share/Reference/bed/${reference_genome}_YNCGNR.bed -b - > tmp.me.bed
		paste tmp.unme.bed tmp.hmei_W.bed tmp.hmei_C.bed tmp.me.bed | awk '{if($4+$8+$12+$16>0) print $1, $2, $3, $4, $8, $12, $16}' OFS='\t' > ${prefix}_CpG_dyad.bed

		nDyad=`awk '{nSum=nSum+($4+$5+$6+$7)} END{print nSum}' ${prefix}_CpG_dyad.bed `
		echo "CpG dyad: $nDyad " | tee -a ${prefix}_Report.txt
		awk '{To+=$4+$5+$6+$7; Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "\nThere are "To" intraCpGs, of which:\n\t"Un" are unmethylated,\n\t"HW" are hemi-Watson,\n\t"HC" are hemi-Crick,\n\t"Me" are methylated."}' ${prefix}_CpG_dyad.bed | tee -a ${prefix}_Report.txt
		echo -e 'Dyad\tCount' > ${prefix}_CpG_dyad_counts.txt
		awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "Unme\t"Un"\nHemiW\t"HW"\nHemiC\t"HC"\nMe\t"Me""}' ${prefix}_CpG_dyad.bed >> ${prefix}_CpG_dyad_counts.txt

		pUseful=`echo "scale=1; 100*$nMotif/$nRaw" | bc`
		echo -e "\nUseful Alignment: $nMotif ($pUseful%)" | tee -a ${prefix}_Report.txt

		echo -e "\n========== ${prefix} done. =========="  | tee -a ${prefix}_Report.txt
		mkdir -p ./process ./bam ./bed
		mv ${prefix}_cut_fragment.bed ${prefix}_CpG.bed ${prefix}_CpG_dyad.bed ./bed
		mv ${prefix}_sorted_dedup.bam ./bam
		mv ${prefix}_mapped.bam ${prefix}_q10.bam ${prefix}_unmapped.bam ${prefix}_trimmed_R1.fastq.gz ${prefix}_trimmed_R2.fastq.gz ./process
		rm tmp*
		cd ..
	done
fi
