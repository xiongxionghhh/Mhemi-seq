#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2022/12/6
##-------------------------

show_usage="Usage: bash pipeline_hpbs.sh -g reference_genome\n\n
			-g|--reference_genome \t reference genome, should be one of hg38 and mm10\n
			-h|--help \t\t get help document\n"

GETOPT_ARGS=`getopt -o g:h -al reference_genome:,help -- "$@"`
eval set -- "$GETOPT_ARGS"

while [ -n "$1" ]
do
        case "$1" in
                -g|--reference_genome) reference_genome=$2; shift 2;;
                -h|--help) echo -e $show_usage; exit 0;;
		--) break;;
                *) echo -e $1,$2,$show_usage; break;;
        esac
done

if [[ -z $reference_genome ]]; then
        echo -e $show_usage
        exit 0
else

	for prefix in `ls tmp.*.R1.fq.gz | sed 's/tmp.//g; s/.R1.fq.gz//g' | sort -u `; do
		zcat tmp.${prefix}.R1.fq.gz | awk '{if($1~/^@/) print $1; else print $0}' > tmp.${prefix}.R1.fq
		zcat tmp.${prefix}.R2.fq.gz | awk '{if($1~/^@/) print $1; else print $0}' > tmp.${prefix}.R2.fq

		grep -A 1 '@' tmp.${prefix}.R1.fq | grep -v '@' | grep -v '-' > tmp.R1

		echo -e "================= Download completed, analysis begin ============" | tee ${prefix}_Report.txt
		raw=`cat tmp.${prefix}.R1.fq | wc -l `
		nRaw=`echo "${raw}/4" | bc`
		pRaw=`echo "scale=2; 100*${nRaw}/${nRaw}" | bc`
		echo "Raw reads: $nRaw (${pRaw}%)" | tee -a ${prefix}_Report.txt
		echo -e "=================================================================\n" | tee -a ${prefix}_Report.txt

		echo -e "================= Calculate conversion rate by hairpin adapter ============" | tee -a ${prefix}_Report.txt
		egrep -o "[CT]GA[CT]GT[CT]GT[CT]GAAGTGAACGA[CT]GA[CT]GT[CT]G" tmp.R1 | awk '{mC+=gsub(/C/, "&")-1} END{print "Conversion Rate:"100-100*mC/(NR*7)"%"}' | tee -a ${prefix}_Report.txt
		echo -e "=================================================================\n" | tee -a ${prefix}_Report.txt

		echo ">hairpin_C2T" > tmp.hairpin_ada.fa
		awk 'NR==2' /mnt/disk1/5/xx/private/Mhmei/scripts/hpbs-seq/hairpin_adapter.fa | sed 'y/C/T/' | sed 's/mT/C/g' >> tmp.hairpin_ada.fa
		echo ">hairpin_G2A" >> tmp.hairpin_ada.fa
		awk 'NR==2' /mnt/disk1/5/xx/private/Mhmei/scripts/hpbs-seq/hairpin_adapter.fa | sed 'y/ATCG/TAGC/' | rev | sed 'y/G/A/' | sed 's/Am/G/' >> tmp.hairpin_ada.fa
		fastp --adapter_fasta tmp.hairpin_ada.fa -q 20 -l 20 --thread 16 --json tmp.hpBS.json --html tmp.hpBS.html --in1 tmp.${prefix}.R1.fq --out1 tmp.${prefix}.R1.QC.fq.gz --in2 tmp.${prefix}.R2.fq --out2 tmp.${prefix}.R2.QC.fq.gz

		echo -e "================= Quality control (R1) ============" | tee -a ${prefix}_Report.txt
		hpAdaC2T=`awk 'NR==2' /mnt/disk1/5/xx/private/Mhmei/scripts/hpbs-seq/hairpin_adapter.fa | sed 'y/C/T/' | sed 's/mT/C/g' `
		cutadapt -j 8 --overlap 1 -g ${hpAdaC2T} -o tmp.${prefix}.Cutadapt.R1.fq.gz tmp.${prefix}.R1.QC.fq.gz
		trim_galore --length 20 --clip_R1 5 --stringency 3 --cores 8 --gzip tmp.${prefix}.Cutadapt.R1.fq.gz
		QC=`gzip -cd tmp.${prefix}.Cutadapt.R1_trimmed.fq.gz | wc -l `
		nQC=`echo "${QC}/4" | bc`
		pQC=`echo "scale=2; 100*${nQC}/${nRaw}" | bc`
		echo "Clean reads(R1): $nQC (${pQC}%)" | tee -a ${prefix}_Report.txt
		echo -e "=================================================================\n" | tee -a ${prefix}_Report.txt

		echo -e "================= Quality control (R2) ============" | tee -a ${prefix}_Report.txt
		hpAdaG2A=`awk 'NR==2' /mnt/disk1/5/xx/private/Mhmei/scripts/hpbs-seq/hairpin_adapter.fa | sed 'y/ATCG/TAGC/' | rev | sed 'y/G/A/' | sed 's/Am/G/' `
		cutadapt -j 8 --overlap 1 -g ${hpAdaG2A} -o tmp.${prefix}.Cutadapt.R2.fq.gz tmp.${prefix}.R2.QC.fq.gz
		trim_galore --length 20 --clip_R1 5 --stringency 3 --cores 8 --gzip tmp.${prefix}.Cutadapt.R2.fq.gz
		QC=`gzip -cd tmp.${prefix}.Cutadapt.R2_trimmed.fq.gz | wc -l `
		nQC=`echo "${QC}/4" | bc`
		pQC=`echo "scale=2; 100*${nQC}/${nRaw}" | bc`
		echo "Clean reads(R2): $nQC (${pQC}%)" | tee -a ${prefix}_Report.txt
		echo -e "=================================================================\n" | tee -a ${prefix}_Report.txt

		echo -e "================= Align to reference genome (R1) ============" | tee -a ${prefix}_Report.txt
		bismark --multicore 10 --gzip --un --ambiguous -q -L 20 -D 20 -R 3 --score_min L,0,-0.4 --genome /mnt/disk1/5/share/index/bismark/${reference_genome} --se tmp.${prefix}.Cutadapt.R1_trimmed.fq.gz
		nMapping1=`samtools view -@ 20 tmp.${prefix}.Cutadapt.R1_trimmed_bismark_bt2.bam | wc -l `
		pMapping=`echo "scale=2; 100*${nMapping1}/${nQC}" | bc`
		echo "Mapped reads(R1): ${nMapping1} (${pMapping}%)" | tee -a ${prefix}_Report.txt
		echo -e "=================================================================\n" | tee -a ${prefix}_Report.txt
		samtools view -H tmp.${prefix}.Cutadapt.R1_trimmed_bismark_bt2.bam > tmp.Header
		samtools view -@ 20 tmp.${prefix}.Cutadapt.R1_trimmed_bismark_bt2.bam | awk '$3!~/[_,L,M,K,Y]/' | cat tmp.Header - | samtools view -@ 20 -bhS > tmp.R1.bam

		echo -e "================= Align to reference genome (R2) ============" | tee -a ${prefix}_Report.txt
		bismark --multicore 10 --pbat --un --ambiguous -q -L 20 -D 20 -R 3 --score_min L,0,-0.4 --genome /mnt/disk1/5/share/index/bismark/${reference_genome} --se tmp.${prefix}.Cutadapt.R2_trimmed.fq.gz
		nMapping2=`samtools view -@ 20 tmp.${prefix}.Cutadapt.R2_trimmed_bismark_bt2.bam | wc -l `
		pMapping=`echo "scale=2; 100*${nMapping2}/${nQC}" | bc`
		echo "Mapped reads(R2): ${nMapping2} (${pMapping}%)" | tee -a ${prefix}_Report.txt
		echo -e "=================================================================\n" | tee -a ${prefix}_Report.txt
		samtools view -H tmp.${prefix}.Cutadapt.R2_trimmed_bismark_bt2.bam > tmp.Header
		samtools view -@ 20 tmp.${prefix}.Cutadapt.R2_trimmed_bismark_bt2.bam | awk '$3!~/[_,L,M,K,Y]/' | cat tmp.Header - | samtools view -@ 20 -bhS > tmp.R2.bam

		echo -e "================= Extract paired reads ============" | tee -a ${prefix}_Report.txt
		samtools view -@ 20 tmp.R1.bam | cut -f 16 > tmp.strand.txt
		bamToBed -i tmp.R1.bam | awk '{print $4, $1, $2, $3, $5, $6}' OFS='\t' | paste - tmp.strand.txt | sort -k1,1 > tmp.R1.bed
		samtools view -@ 20 tmp.R2.bam | cut -f 16 > tmp.strand.txt
		bamToBed -i tmp.R2.bam | awk '{print $4, $1, $2, $3, $5, $6}' OFS='\t' | paste - tmp.strand.txt | sort -k1,1 > tmp.R2.bed

		join -j 1 tmp.R1.bed tmp.R2.bed | sed 's/ /\t/g' > tmp.all.pairs
		nAllPaired=`cat tmp.all.pairs | wc -l `
		nSingleton1=`echo "scale=2; ${nMapping1}-${nAllPaired}" | bc`
		pSingleton1=`echo "scale=2; 100*${nSingleton1}/${nMapping1}" | bc`
		echo "Singleton (R1): ${nSingleton1} (${pSingleton1}%)" | tee -a ${prefix}_Report.txt

		nSingleton2=`echo "scale=2; ${nMapping2}-${nAllPaired}" | bc`
		pSingleton2=`echo "scale=2; 100*${nSingleton2}/${nMapping2}" | bc`
		echo "Singleton (R2): ${nSingleton2} (${pSingleton2}%)" | tee -a ${prefix}_Report.txt

		awk 'BEGIN{print "Dyad\tCount"} {if($2==$8 && $7==$13) n00+=1; else if($2==$8 && $7!=$13) n01+=1; else if($2!=$8 && $7==$13) n10+=1; else if($2!=$8 && $7!=$13) n11+=1} END{print "C0S0\t"n00"\nC0S1\t"n01"\nC1S0\t"n10"\nC1S1\t"n11""}' tmp.all.pairs > tmp.Count.txt
		Rscript /mnt/disk1/5/xx/private/Mhmei/scripts/hpbs-seq/Donut_ggplot2.R
		mv CELL_METHOD.pdf ${prefix}_Chr+Strand.pdf

		awk '$2==$8 && $7!=$13' tmp.all.pairs | awk '{dis1=sqrt(($3-$9)*($3-$9)); dis2=sqrt(($4-$10)*($4-$10))} {if(dis1<=3 && dis2<=3) print $2, $3, $4, $5, $8, $9, $10, $11, $1}' OFS='\t' | sort -k1,1 -k2,2n -k3,3n > tmp.pairs
		nPaired=`cat tmp.pairs | wc -l `
		pPaired=`echo "scale=2; 100*${nPaired}/${nAllPaired}" | bc`
		echo -e "\nValid Pairs: ${nPaired} (${pPaired}%)" | tee -a ${prefix}_Report.txt
		echo -e "=================================================================\n" | tee -a ${prefix}_Report.txt

		groupBy -i tmp.pairs -g 1,2,3,4,5,6,7,8 -c 9 -o distinct | sed 's/,.*//g' > tmp.dedup.pairs
		nDedup=`cat tmp.dedup.pairs | wc -l `
		nDup=`echo "scale=2; ${nPaired}-${nDedup}" | bc`
		pDup=`echo "scale=2; 100*${nDup}/${nPaired}" | bc`
		echo "Duplication: ${nDup} (${pDup}%)" | tee -a ${prefix}_Report.txt
		echo -e "=================================================================\n" | tee -a ${prefix}_Report.txt

		awk '$4>9 || $8>9' tmp.dedup.pairs | cut -f 9 > tmp.ReadNames

		samtools view -H tmp.R1.bam > tmp.Header
		samtools view -@ 20 tmp.R1.bam | grep -w -F -f tmp.ReadNames - | cat tmp.Header - | samtools view -@ 20 -hbS > ${prefix}.R1.Paired.bam

		samtools view -H tmp.R2.bam > tmp.Header
		samtools view -@ 20 tmp.R2.bam | grep -w -F -f tmp.ReadNames - | cat tmp.Header - | samtools view -@ 20 -hbS > ${prefix}.R2.Paired.bam

		bedtools bamtobed -i ${prefix}.R1.Paired.bam | awk '{if($6=="+"){print $1, $2-46, $2+154} else{print $1,$3-154, $3+46}}' OFS='\t' | awk '$1!~/chr[C,M,T,L]/ && $2>-1' > tmp.Reads.bed
		nLine=`cat tmp.Reads.bed | wc -l`
		samtools view -H ${prefix}.R1.Paired.bam | awk '$1=="@SQ" {OFS="\t"; print $2,$3}' - | sed 's/.N://g' > tmp.size
		sort -k1,1 -k2,2n -S 9G tmp.Reads.bed | genomeCoverageBed -bg -i - -g tmp.size | awk '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4/"'$nLine'"*1000000}' | awk '{$4/=1;print}' OFS='\t' > tmp.Reads.bdg
		bedGraphToBigWig tmp.Reads.bdg tmp.size ${prefix}.bw

		nMotifReads=`bedtools intersect -a tmp.Reads.bed -b /mnt/disk4/public/RefBed/CTCF/Hs_CTCF.motif | wc -l `
		pMotifReads=`echo "scale=2; 100*${nMotifReads}/${nLine}" | bc`
		echo "Reads cover motif: ${nMotifReads} (${pMotifReads}%)" | tee -a ${prefix}_Report.txt

		bismark_methylation_extractor -s --multicore 20 --no_header --gzip ${prefix}.R1.Paired.bam
		bismark_methylation_extractor -s --multicore 20 --no_header --gzip ${prefix}.R2.Paired.bam

		gzip -cd CpG_OT_${prefix}.R1.Paired.txt.gz | awk '{print $1, $2, $3, $4}' OFS='\t' | sort -k1,1 -k2,2n > tmp.R1.txt
		gzip -cd CpG_CTOB_${prefix}.R2.Paired.txt.gz | awk '{print $1, $2, $3, $4}' OFS='\t' | sort -k1,1 -k2,2n > tmp.R2.txt
		join -j 1 tmp.R1.txt tmp.R2.txt | awk '{if($4-$7==-1) print $3, $4-1, $7, $2, $5; else if($4-$7==1) print $3, $7-1, $4, $5, $2}' OFS='\t' > tmp.pairs.bed

		gzip -cd CpG_OB_${prefix}.R1.Paired.txt.gz | awk '{print $1, $2, $3, $4}' OFS='\t' | sort -k1,1 -k2,2n > tmp.R1.txt
		gzip -cd CpG_CTOT_${prefix}.R2.Paired.txt.gz | awk '{print $1, $2, $3, $4}' OFS='\t' | sort -k1,1 -k2,2n > tmp.R2.txt
		join -j 1 tmp.R1.txt tmp.R2.txt | awk '{if($4-$7==-1) print $3, $4-1, $7, $2, $5; else if($4-$7==1) print $3, $7-1, $4, $5, $2}' OFS='\t' >> tmp.pairs.bed
		sort -k1,1 -k2,2n tmp.pairs.bed > tmp.pairs.sorted.bed
		cut -f 1-3 tmp.pairs.sorted.bed | sort -k1,1 -k2,2n -u > tmp.dyadAll.bed

		awk '$4=="-" && $5=="-"' tmp.pairs.sorted.bed | intersectBed -c -sorted -a tmp.dyadAll.bed -b - > tmp.unme.bed
		awk '$4=="+" && $5=="-"' tmp.pairs.sorted.bed | intersectBed -c -sorted -a tmp.dyadAll.bed -b - > tmp.hemiW.bed
		awk '$4=="-" && $5=="+"' tmp.pairs.sorted.bed | intersectBed -c -sorted -a tmp.dyadAll.bed -b - > tmp.hemiC.bed
		awk '$4=="+" && $5=="+"' tmp.pairs.sorted.bed | intersectBed -c -sorted -a tmp.dyadAll.bed -b - > tmp.me.bed
		paste tmp.unme.bed tmp.hemiW.bed tmp.hemiC.bed tmp.me.bed | cut -f 1-4,8,12,16 | intersectBed -u -sorted -F 1 -f 1 -a - -b /mnt/disk1/5/share/Reference/${reference_genome}_CpG.bed > ${prefix}.intraCpG.bdg

		awk '{mSum+=$5+$6+2*$7; nSum+=2*($4+$5+$6+$7)} END{print "CpG methylation(%):"100*mSum/nSum"\n"}' ${prefix}.intraCpG.bdg | tee -a ${prefix}_Report.txt
		awk '{To+=$4+$5+$6+$7; Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "There are "To" intraCpGs, of which:\n\t"Un" are unmethylated,\n\t"HW" are hemi-Watson,\n\t"HC" are hemi-Crick,\n\t"Me" are methylated."}' ${prefix}.intraCpG.bdg | tee -a ${prefix}_Report.txt

		echo -e "=================================================================" | tee -a ${prefix}_Report.txt
		echo -e 'Dyad\tCount' > tmp.Count.txt
		awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "-     -\t"Un"\n+     -\t"HW"\n-     +\t"HC"\n+     +\t"Me""}' ${prefix}.intraCpG.bdg >> tmp.Count.txt
		Rscript /mnt/disk1/5/xx/private/Mhmei/scripts/hpbs-seq/Donut_ggplot2.R
		mv CELL_METHOD.pdf ${prefix}_intraCpG.pdf

		echo -e "=================================================================" | tee -a ${prefix}_Report.txt
		gzip -cd CHG_OT_${prefix}.R1.Paired.txt.gz | awk '{print $1, $2, $3, $4}' OFS='\t' | sort -k1,1 -k2,2n > tmp.R1.txt
		gzip -cd CHG_CTOB_${prefix}.R2.Paired.txt.gz | awk '{print $1, $2, $3, $4}' OFS='\t' | sort -k1,1 -k2,2n > tmp.R2.txt
		join -j 1 tmp.R1.txt tmp.R2.txt | awk '{if($4-$7==-2) print $3, $4-1, $7, $2, $5; else if($4-$7==2) print $3, $7-1, $4, $5, $2}' OFS='\t' > tmp.pairs.bed

		gzip -cd CHG_OB_${prefix}.R1.Paired.txt.gz | awk '{print $1, $2, $3, $4}' OFS='\t' | sort -k1,1 -k2,2n > tmp.R1.txt
		gzip -cd CHG_CTOT_${prefix}.R2.Paired.txt.gz | awk '{print $1, $2, $3, $4}' OFS='\t' | sort -k1,1 -k2,2n > tmp.R2.txt
		join -j 1 tmp.R1.txt tmp.R2.txt | awk '{if($4-$7==-2) print $3, $4-1, $7, $2, $5; else if($4-$7==2) print $3, $7-1, $4, $5, $2}' OFS='\t' >> tmp.pairs.bed
		sort -k1,1 -k2,2n tmp.pairs.bed > tmp.pairs.sorted.bed
		cut -f 1-3 tmp.pairs.sorted.bed | sort -k1,1 -k2,2n -u > tmp.dyadAll.bed

		awk '$4=="-" && $5=="-"' tmp.pairs.sorted.bed | intersectBed -c -sorted -a tmp.dyadAll.bed -b - > tmp.unme.bed
		awk '$4=="+" && $5=="-"' tmp.pairs.sorted.bed | intersectBed -c -sorted -a tmp.dyadAll.bed -b - > tmp.hemiW.bed
		awk '$4=="-" && $5=="+"' tmp.pairs.sorted.bed | intersectBed -c -sorted -a tmp.dyadAll.bed -b - > tmp.hemiC.bed
		awk '$4=="+" && $5=="+"' tmp.pairs.sorted.bed | intersectBed -c -sorted -a tmp.dyadAll.bed -b - > tmp.me.bed
		paste tmp.unme.bed tmp.hemiW.bed tmp.hemiC.bed tmp.me.bed | cut -f 1-4,8,12,16 | intersectBed -u -sorted -F 1 -f 1 -a - -b /mnt/disk1/5/share/Reference/${reference_genome}_CWG.bed > ${prefix}.intraCWG.bdg

		awk '{mSum+=$5+$6+2*$7; nSum+=2*($4+$5+$6+$7)} END{print "CWG methylation(%):"100*mSum/nSum"\n"}' ${prefix}.intraCWG.bdg | tee -a ${prefix}_Report.txt
		awk '{To+=$4+$5+$6+$7; Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "There are "To" intraCWGs, of which:\n\t"Un" are unmethylated,\n\t"HW" are hemi-Watson,\n\t"HC" are hemi-Crick,\n\t"Me" are methylated."}' ${prefix}.intraCWG.bdg | tee -a ${prefix}_Report.txt

		echo -e 'Dyad\tCount' > tmp.Count.txt
		awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "-     -\t"Un"\n+     -\t"HW"\n-     +\t"HC"\n+     +\t"Me""}' ${prefix}.intraCWG.bdg >> tmp.Count.txt
		Rscript /mnt/disk1/5/xx/private/Mhmei/scripts/hpbs-seq/Donut_ggplot2.R
		mv CELL_METHOD.pdf ${prefix}_intraCWG.pdf

		echo -e "=================================================================" | tee -a ${prefix}_Report.txt
		echo -e "Job Completed." | tee -a ${prefix}_Report.txt
		echo -e "=================================================================" | tee -a ${prefix}_Report.txt

	done

fi

