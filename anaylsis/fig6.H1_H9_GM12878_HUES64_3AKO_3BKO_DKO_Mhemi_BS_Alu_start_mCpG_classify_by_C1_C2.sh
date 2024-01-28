#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/10
##-------------------------

awk '$8~/Group1/' H1_Alu_group_classify_by_C1_C2_delta_mCWG_oppo-same.bed | head -5000 | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif
awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bedGraph/H1_BS_CpG.bdg -b - > tmp.mC.sorted.bed
echo -e "Distance\tMethylation\tStrand" > tmp.mC.txt
echo -e "Distance\tMethylation" > tmp.DeltaC.txt
for((pos=-100; pos<=300; pos+=1)); do
	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		oppoM=`awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print mSum/nCount}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${oppoM}\toppo" >> tmp.mC.txt
	fi

	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		sameM=`awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print mSum/nCount}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${sameM}\tsame" >> tmp.mC.txt
	fi
	deltaM=`echo "${oppoM}-${sameM}" | bc `
	echo -e "${pos}\t${deltaM}" >> tmp.DeltaC.txt
done
mv tmp.DeltaC.txt H1_BS_DeltaCpG_Group1_Alu_Start.txt
mv tmp.mC.txt H1_BS_mCpG_Group1_Alu_Start.txt



awk '$8~/Group2/' H1_Alu_group_classify_by_C1_C2_delta_mCWG_oppo-same.bed | head -5000 | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif
awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bedGraph/H1_BS_CpG.bdg -b - > tmp.mC.sorted.bed
echo -e "Distance\tMethylation\tStrand" > tmp.mC.txt
echo -e "Distance\tMethylation" > tmp.DeltaC.txt
for((pos=-100; pos<=300; pos+=1)); do
	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		oppoM=`awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print mSum/nCount}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${oppoM}\toppo" >> tmp.mC.txt
	fi

	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		sameM=`awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print mSum/nCount}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${sameM}\tsame" >> tmp.mC.txt
	fi
	deltaM=`echo "${oppoM}-${sameM}" | bc `
	echo -e "${pos}\t${deltaM}" >> tmp.DeltaC.txt
done
mv tmp.DeltaC.txt H1_BS_DeltaCpG_Group2_Alu_Start.txt
mv tmp.mC.txt H1_BS_mCpG_Group2_Alu_Start.txt



awk '$8~/Group3/' H1_Alu_group_classify_by_C1_C2_delta_mCWG_oppo-same.bed | head -5000 | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif
awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bedGraph/H1_BS_CpG.bdg -b - > tmp.mC.sorted.bed
echo -e "Distance\tMethylation\tStrand" > tmp.mC.txt
echo -e "Distance\tMethylation" > tmp.DeltaC.txt
for((pos=-100; pos<=300; pos+=1)); do
	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		oppoM=`awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print mSum/nCount}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${oppoM}\toppo" >> tmp.mC.txt
	fi

	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		sameM=`awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print mSum/nCount}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${sameM}\tsame" >> tmp.mC.txt
	fi
	deltaM=`echo "${oppoM}-${sameM}" | bc `
	echo -e "${pos}\t${deltaM}" >> tmp.DeltaC.txt
done
mv tmp.DeltaC.txt H1_BS_DeltaCpG_Group3_Alu_Start.txt
mv tmp.mC.txt H1_BS_mCpG_Group3_Alu_Start.txt


awk '$8~/Group4/' H1_Alu_group_classify_by_C1_C2_delta_mCWG_oppo-same.bed | head -5000 | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif
awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bedGraph/H1_BS_CpG.bdg -b - > tmp.mC.sorted.bed
echo -e "Distance\tMethylation\tStrand" > tmp.mC.txt
echo -e "Distance\tMethylation" > tmp.DeltaC.txt
for((pos=-100; pos<=300; pos+=1)); do
	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		oppoM=`awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print mSum/nCount}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${oppoM}\toppo" >> tmp.mC.txt
	fi

	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		sameM=`awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print mSum/nCount}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${sameM}\tsame" >> tmp.mC.txt
	fi
	deltaM=`echo "${oppoM}-${sameM}" | bc `
	echo -e "${pos}\t${deltaM}" >> tmp.DeltaC.txt
done
mv tmp.DeltaC.txt H1_BS_DeltaCpG_Group4_Alu_Start.txt
mv tmp.mC.txt H1_BS_mCpG_Group4_Alu_Start.txt
