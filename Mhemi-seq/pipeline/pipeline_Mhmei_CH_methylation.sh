#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/10/3
##-------------------------

for i in GM12878_Mhemi_rep123 H1_Mhemi_rep123 H9_Mhemi_rep123; do

	bamToBed -i ${i}.bam | awk '$2>20 && $1!~/[_,M,L,T,Y]/' > tmp.bamReads.bed

	awk '{if($6=="+") print $1, $2-14, $2+18, "A0", 1, $6; else print $1, $3-18, $3+14, "A0", 1, $6}' OFS='\t' tmp.bamReads.bed > tmp.unsort.bed
	sort -S 20% -k1,1 -k2,2n tmp.unsort.bed > tmp.End.bed
	awk '{if($3-$2>40) print $1, $2+18, $3-18}' OFS='\t' tmp.bamReads.bed > tmp.unsort.bed
	sort -S 20% -k1,1 -k2,2n tmp.unsort.bed > tmp.inter.bed

	fastaFromBed -s -tab -fi /mnt/disk1/5/share/Reference/fasta/hg38.fa -bed tmp.End.bed | awk '{seqs=toupper($2)} {seq1=substr(seqs, 1, 5); seq2=substr(seqs, 28, 5)} {if(seq1~/[A,T,C,G]C[A,T,C][A,T,C,G][A,G]/ && seq2!~/[C,T][A,T,C,G][A,T,C,G]G[A,T,C,G]/) print "+"; else if(seq1!~/[A,T,C,G]C[A,T,C,G][A,T,C,G][A,G]/ && seq2~/[C,T][A,T,C,G][A,T,G]G[A,T,C,G]/) print "-"; else print "Others"}' > tmp.motif.txt
	paste tmp.End.bed tmp.motif.txt | grep -v 'Others' > tmp.EndMotif.bed

	awk '{if($6==$7) print $1, $2, $2+3}' OFS='\t' tmp.EndMotif.bed | intersectBed -c -sorted -a /mnt/disk1/5/share/Reference/bed/hg38_CHGR_Wat.bed -b - > tmp.mWat.bed
	intersectBed -c -sorted -a /mnt/disk1/5/share/Reference/bed/hg38_CHGR_Wat.bed -b tmp.inter.bed > tmp.uWat.bed
	paste tmp.mWat.bed tmp.uWat.bed | awk '{if($5+$10>0) print $1, $2, $3, $5, $5+$10, "+", $4}' OFS='\t' > tmp.${i}_merged.CHG.bdg

	awk '{if($6!=$7) print $1, $3-3, $3}' OFS='\t' tmp.EndMotif.bed | intersectBed -c -sorted -a /mnt/disk1/5/share/Reference/bed/hg38_CHGR_Cri.bed -b - > tmp.mCri.bed
	intersectBed -c -sorted -a /mnt/disk1/5/share/Reference/bed/hg38_CHGR_Cri.bed -b tmp.inter.bed > tmp.uCri.bed
	paste tmp.mCri.bed tmp.uCri.bed | awk '{if($5+$10>0) print $1, $2, $3, $5, $5+$10, "-", $4}' OFS='\t' | cat - tmp.${i}_merged.CHG.bdg > tmp.unsort.bdg
	sort -S 20% -k1,1 -k2,2n tmp.unsort.bdg > ${i}_merged_CHG.bdg


	awk '{if($6==$7) print $1, $2, $2+3}' OFS='\t' tmp.EndMotif.bed | intersectBed -c -sorted -a /mnt/disk1/5/share/Reference/bed/hg38_CHHR_Wat.bed -b - > tmp.mWat.bed
	intersectBed -c -sorted -a /mnt/disk1/5/share/Reference/bed/hg38_CHHR_Wat.bed -b tmp.inter.bed > tmp.uWat.bed
	paste tmp.mWat.bed tmp.uWat.bed | awk '{if($5+$10>0) print $1, $2, $3, $5, $5+$10, "+", $4}' OFS='\t' > tmp.${i}_merged.CHH.bdg

	awk '{if($6!=$7) print $1, $3-3, $3}' OFS='\t' tmp.EndMotif.bed | intersectBed -c -sorted -a /mnt/disk1/5/share/Reference/bed/hg38_CHHR_Cri.bed -b - > tmp.mCri.bed
	intersectBed -c -sorted -a /mnt/disk1/5/share/Reference/bed/hg38_CHHR_Cri.bed -b tmp.inter.bed > tmp.uCri.bed
	paste tmp.mCri.bed tmp.uCri.bed | awk '{if($5+$10>0) print $1, $2, $3, $5, $5+$10, "-", $4}' OFS='\t' | cat - tmp.${i}_merged.CHH.bdg > tmp.unsort.bdg
	sort -S 20% -k1,1 -k2,2n tmp.unsort.bdg > ${i}_merged_CHH.bdg

	awk '{print $1, $2, $3, 100*$4/$5, $5, $6, $7}' OFS='\t' ${i}_merged_CHG.bdg > tmp.unsort.bdg
	sort -S 20% -k4,4nr -k5,5nr tmp.unsort.bdg > ${i}_merged_CHG_mSorted.bdg
	awk '{print $1, $2, $3, 100*$4/$5, $5, $6, $7}' OFS='\t' ${i}_merged_CHH.bdg > tmp.unsort.bdg
	sort -S 20% -k4,4nr -k5,5nr tmp.unsort.bdg > ${i}_merged_CHH_mSorted.bdg

	# dyad methylation
	samtools sort -n -m 2G -@ 10 ${i}.bam | bedtools bamtobed -bedpe -i - | awk '{umi=substr($7,length($7)-4,5)} {if(($9=="+") && ($10=="-")) print $1, $2, $6, umi, $7; else if(($9=="-") && ($10=="+")) print $4, $5, $3, umi, $7}' OFS='\t' > tmp.unsort.bed
	sort -S 20% -k1,1 -k2,2n tmp.unsort.bed > tmp.sorted.bed

	forward1="C[A,T,C,G][A,T,C,G][A,G]"
	reverse1="[C,T][A,T,C,G][A,T,C,G]G"

	awk '{print $1, $2-13, $3+13}' OFS='\t' tmp.sorted.bed | bedtools getfasta -tab -fi /mnt/disk1/5/share/Reference/fasta/hg38.fa -bed - | awk '{seqs=toupper($2)} {headSeq1=substr(seqs, 1, 4); headSeq2=substr(seqs, 27, 4); tailSeq1=substr(seqs, length(seqs)-29, 4); tailSeq2=substr(seqs, length(seqs)-3, 4)} {if(headSeq1~/'$forward1'/ && headSeq2!~/'$reverse1'/ && tailSeq1~/'$forward1'/ && tailSeq2!~/'$reverse1'/) print "H1T1"; else if(headSeq1!~/'$forward1'/ && headSeq2~/'$reverse1'/ && tailSeq1~/'$forward1'/ && tailSeq2!~/'$reverse1'/) print "H2T1"; else if(headSeq1~/'$forward1'/ && headSeq2!~/'$reverse1'/ && tailSeq1!~/'$forward1'/ && tailSeq2~/'$reverse1'/) print "H1T2"; else if(headSeq1!~/'$forward1'/ && headSeq2~/'$reverse1'/ && tailSeq1!~/'$forward1'/ && tailSeq2~/'$reverse1'/) print "H2T2"; else print "Others"}' > tmp.${i}.motif_frag.txt
	awk 'BEGIN{nH1T1=0; nH1T2=0; nH2T1=0; nH2T2=0; nOthers=0} {if($1=="H1T1") nH1T1+=1; else if($1=="H1T2") nH1T2+=1; else if($1=="H2T1") nH2T1+=1; else if($1=="H2T2") nH2T2+=1; else if($1=="Others") nOthers+=1} END{print "End\tCount\nH1T1\t"nH1T1"\nH1T2\t"nH1T2"\nH2T1\t"nH2T1"\nH2T2\t"nH2T2"\nOthers\t"nOthers}' tmp.${i}.motif_frag.txt > ${i}_EndCount.txt
	paste tmp.sorted.bed tmp.${i}.motif_frag.txt | awk '$6!="Others"' > ${i}_CNNR_cut_fragment.bed

	awk '{if(($3-$2)>40) print $1, $2+18, $3-18}' OFS='\t' ${i}_CNNR_cut_fragment.bed | bedtools getfasta -tab -fi /mnt/disk1/5/share/Reference/fasta/hg38.fa -bed - | awk '{if(toupper($2)~/[C,T]C[A,T]G[A,G]/) print $1}'| sed 's/:/\t/g; s/-/\t/g' | sort -S 20% -k1,1 -k2,2n | bedtools intersect -f 1 -sorted -c -a /mnt/disk1/5/share/Reference/bed/hg38_YCWGR.bed -b - > tmp.unme.bed
	awk '{if(($3-$2)>40 && $6~/T1/) print $1, $3-18, $3-13}' OFS='\t' ${i}_CNNR_cut_fragment.bed | bedtools getfasta -tab -fi /mnt/disk1/5/share/Reference/fasta/hg38.fa -bed - | awk '{if(toupper($2)~/[C,T]C[A,T]G[A,G]/) print $1}'| sed 's/:/\t/g; s/-/\t/g' | sort -S 20% -k1,1 -k2,2n | bedtools intersect -f 1 -sorted -c -a /mnt/disk1/5/share/Reference/bed/hg38_YCWGR.bed -b - > tmp.hmei_W.bed
	awk '{if(($3-$2)>40 && $6~/H2/) print $1, $2+13, $2+18}' OFS='\t' ${i}_CNNR_cut_fragment.bed | bedtools getfasta -tab -fi /mnt/disk1/5/share/Reference/fasta/hg38.fa -bed - | awk '{if(toupper($2)~/[C,T]C[A,T]G[A,G]/) print $1}'| sed 's/:/\t/g; s/-/\t/g' | sort -S 20% -k1,1 -k2,2n | bedtools intersect -f 1 -sorted -c -a /mnt/disk1/5/share/Reference/bed/hg38_YCWGR.bed -b - > tmp.hmei_C.bed
	awk '{if(($3-$2)<=40 && $6=="H2T1") print $1, $2+13, $2+18}' OFS='\t' ${i}_CNNR_cut_fragment.bed | bedtools getfasta -tab -fi /mnt/disk1/5/share/Reference/fasta/hg38.fa -bed - | awk '{if(toupper($2)~/[C,T]C[A,T]G[A,G]/) print $1}'| sed 's/:/\t/g; s/-/\t/g' | sort -S 20% -k1,1 -k2,2n | bedtools intersect -f 1 -sorted -c -a /mnt/disk1/5/share/Reference/bed/hg38_YCWGR.bed -b - > tmp.me.bed
	paste tmp.unme.bed tmp.hmei_W.bed tmp.hmei_C.bed tmp.me.bed | awk '{if($4+$8+$12+$16>0) print $1, $2, $3, $4, $8, $12, $16}' OFS='\t' > ${i}_CWG_dyad.bed

	nDyad=`awk '{nSum=nSum+($4+$5+$6+$7)} END{print nSum}' ${i}_CWG_dyad.bed `
	echo "CWG dyad: $nDyad " | tee ${i}_CHG_Report.txt
	awk '{To+=$4+$5+$6+$7; Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "\nThere are "To" intraCWGs, of which:\n\t"Un" are unmethylated,\n\t"HW" are hemi-Watson,\n\t"HC" are hemi-Crick,\n\t"Me" are methylated."}' ${i}_CWG_dyad.bed | tee -a ${i}_CHG_Report.txt
	echo -e 'Dyad\tCount' > ${i}_CWG_dyad_counts.txt
	awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "Unme\t"Un"\nHemiW\t"HW"\nHemiC\t"HC"\nMe\t"Me""}' ${i}_CWG_dyad.bed >> ${i}_CWG_dyad_counts.txt
done

