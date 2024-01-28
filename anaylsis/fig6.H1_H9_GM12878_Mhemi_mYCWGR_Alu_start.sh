#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/23
##-------------------------

awk '$5>=0' /mnt/disk1/5/xx/private/Mhmei/H1/RNA-seq/TPM/H1_SINE-AluOnly.bdg | sort -k1,1 -k2,2n -S 5G -u  > tmp.Alu.motif

awk '{print $1, $2, $2+2, $5, $4+$5+$6+$7, "+"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/GM12878_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_W.bed
awk '{print $1, $3-2, $3, $6, $4+$5+$6+$7, "-"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/GM12878_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_C.bed
cat tmp.Dyad_W.bed tmp.Dyad_C.bed | sort -k1,1 -k2,2n -S 5G -u  > tmp.Dyad.bed

awk '{if($2>180) print $1, $2-180, $3+180}' OFS='\t' tmp.Alu.motif | bedtools intersect -u -sorted -a tmp.Dyad.bed -b - > tmp.mC.sorted.bed

echo -e "Distance\tMethylation\tStrand" > tmp.mC.txt
for((pos=-100; pos<=300; pos+=1)); do
	awk '$2>310' tmp.Alu.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount, "same"}' OFS='\t' tmp.bed >> tmp.mC.txt
	else
		echo -e "$pos\tNA\tsame" >> tmp.mC.txt
	fi

	awk '$2>310' tmp.Alu.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount, "oppo"}' OFS='\t' tmp.bed >> tmp.mC.txt
	else
		echo -e "$pos\tNA\toppo" >> tmp.mC.txt
	fi

done
mv tmp.mC.txt GM12878_Mhemi_mYCWGR_Alu_start.txt
echo "done."



awk '{print $1, $2, $2+2, $5, $4+$5+$6+$7, "+"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/H1_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_W.bed
awk '{print $1, $3-2, $3, $6, $4+$5+$6+$7, "-"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/H1_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_C.bed
cat tmp.Dyad_W.bed tmp.Dyad_C.bed | sort -k1,1 -k2,2n -S 5G -u  > tmp.Dyad.bed

awk '{if($2>180) print $1, $2-180, $3+180}' OFS='\t' tmp.Alu.motif | bedtools intersect -u -sorted -a tmp.Dyad.bed -b - > tmp.mC.sorted.bed

echo -e "Distance\tMethylation\tStrand" > tmp.mC.txt
for((pos=-100; pos<=300; pos+=1)); do
	awk '$2>310' tmp.Alu.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount, "same"}' OFS='\t' tmp.bed >> tmp.mC.txt
	else
		echo -e "$pos\tNA\tsame" >> tmp.mC.txt
	fi

	awk '$2>310' tmp.Alu.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount, "oppo"}' OFS='\t' tmp.bed >> tmp.mC.txt
	else
		echo -e "$pos\tNA\toppo" >> tmp.mC.txt
	fi

done
mv tmp.mC.txt H1_Mhemi_mYCWGR_Alu_start.txt
echo "done."






awk '{print $1, $2, $2+2, $5, $4+$5+$6+$7, "+"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/H9_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_W.bed
awk '{print $1, $3-2, $3, $6, $4+$5+$6+$7, "-"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/H9_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_C.bed
cat tmp.Dyad_W.bed tmp.Dyad_C.bed | sort -k1,1 -k2,2n -S 5G -u  > tmp.Dyad.bed

awk '{if($2>180) print $1, $2-180, $3+180}' OFS='\t' tmp.Alu.motif | bedtools intersect -u -sorted -a tmp.Dyad.bed -b - > tmp.mC.sorted.bed

echo -e "Distance\tMethylation\tStrand" > tmp.mC.txt
for((pos=-100; pos<=300; pos+=1)); do
	awk '$2>310' tmp.Alu.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount, "same"}' OFS='\t' tmp.bed >> tmp.mC.txt
	else
		echo -e "$pos\tNA\tsame" >> tmp.mC.txt
	fi

	awk '$2>310' tmp.Alu.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount, "oppo"}' OFS='\t' tmp.bed >> tmp.mC.txt
	else
		echo -e "$pos\tNA\toppo" >> tmp.mC.txt
	fi

done
mv tmp.mC.txt H9_Mhemi_mYCWGR_Alu_start.txt
echo "done."
