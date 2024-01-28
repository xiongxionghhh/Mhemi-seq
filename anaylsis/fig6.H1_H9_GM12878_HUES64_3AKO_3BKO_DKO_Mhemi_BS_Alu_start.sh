#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/10
##-------------------------

awk '$5>=0' /mnt/disk1/5/xx/private/Mhmei/H1/RNA-seq/TPM/H1_SINE-AluOnly.bdg | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif
awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/share/Mhemi-seq/CH/H1_Mhemi_rep123_CHG.bdg -b - | awk '{print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' > tmp.mC.sorted.bed
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
mv tmp.DeltaC.txt H1_Mhemi_DeltaCHG_Alu_Start.txt
mv tmp.mC.txt H1_Mhemi_mCHG_Alu_Start.txt

awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bedGraph/H1_BS_CHG.bdg -b - > tmp.mC.sorted.bed
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
mv tmp.DeltaC.txt H1_BS_DeltaCHG_Alu_Start.txt
mv tmp.mC.txt H1_BS_mCHG_Alu_Start.txt



awk '$5>=0' /mnt/disk1/5/xx/private/Mhmei/H9/RNA-seq/TPM/SINE-Alu.bdg | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif
awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/H9/Mhemi-seq/merged/CH/H9_merged.CHG.bdg -b - | awk '{print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' > tmp.mC.sorted.bed
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
mv tmp.DeltaC.txt H9_Mhemi_DeltaCHG_Alu_Start.txt
mv tmp.mC.txt H9_Mhemi_mCHG_Alu_Start.txt

awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/H9/BS-seq/bedGraph/H9_BS_CHG.bdg -b - > tmp.mC.sorted.bed
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
mv tmp.DeltaC.txt H9_BS_DeltaCHG_Alu_Start.txt
mv tmp.mC.txt H9_BS_mCHG_Alu_Start.txt



awk '$5>=0' /mnt/disk1/5/xx/private/Mhmei/H9/RNA-seq/TPM/SINE-Alu.bdg | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif
awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/Mhemi-seq/merge/CH/GM12878_merged.CHG.bdg -b - | awk '{print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' > tmp.mC.sorted.bed
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
mv tmp.DeltaC.txt GM12878_Mhemi_DeltaCHG_Alu_Start.txt
mv tmp.mC.txt GM12878_Mhemi_mCHG_Alu_Start.txt

awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_BS_CHG.bdg -b - > tmp.mC.sorted.bed
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
mv tmp.DeltaC.txt GM12878_BS_DeltaCHG_Alu_Start.txt
mv tmp.mC.txt GM12878_BS_mCHG_Alu_Start.txt






awk '$5>=0' /mnt/disk1/5/xx/private/Mhmei/H1/RNA-seq/TPM/H1_SINE-AluOnly.bdg | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif

awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/HUES64/BS-seq/HUES64_BS_CHG.bdg -b - > tmp.mC.sorted.bed
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
mv tmp.DeltaC.txt HUES64_BS_DeltaCHG_Alu_Start.txt
mv tmp.mC.txt HUES64_BS_mCHG_Alu_Start.txt

awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/HUES64/BS-seq/HUES64-3AKOL_BS_CHG.bdg -b - > tmp.mC.sorted.bed
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
mv tmp.DeltaC.txt HUES64-3AKOL_BS_DeltaCHG_Alu_Start.txt
mv tmp.mC.txt HUES64-3AKOL_BS_mCHG_Alu_Start.txt

awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/HUES64/BS-seq/HUES64-3BKOL_BS_CHG.bdg -b - > tmp.mC.sorted.bed
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
mv tmp.DeltaC.txt HUES64-3BKOL_BS_DeltaCHG_Alu_Start.txt
mv tmp.mC.txt HUES64-3BKOL_BS_mCHG_Alu_Start.txt


awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/HUES64/BS-seq/HUES64-DKOL_BS_CHG.bdg -b - > tmp.mC.sorted.bed
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
mv tmp.DeltaC.txt HUES64-DKOL_BS_DeltaCHG_Alu_Start.txt
mv tmp.mC.txt HUES64-DKOL_BS_mCHG_Alu_Start.txt
