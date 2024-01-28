#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/10
##-------------------------

awk '$5>=0' /mnt/disk1/5/xx/private/Mhmei/H1/RNA-seq/TPM/H1_SINE-AluOnly.bdg | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif
awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/H1/MNase-seq/bedGraph/H1_MNase-seq_134-160.bdg -b - > tmp.mC.sorted.bed
echo -e "Distance\tSignal" > tmp.mC.txt
for((pos=-100; pos<=300; pos+=1)); do
	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' > tmp.txt
	sort -k1,1 -k2,2n -S 5G tmp.txt | bedtools intersect -sorted -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		oppoM=`awk -v pos=$pos '{mSum+=$4} END{print mSum/NR}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${oppoM}" >> tmp.mC.txt
	fi
done
mv tmp.mC.txt H1_MNase_Alu_Start.txt

awk '$5>=0' /mnt/disk1/5/xx/private/Mhmei/H9/RNA-seq/TPM/SINE-Alu.bdg | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif
awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/H9/MNase-seq/H9_MNase.bdg -b - > tmp.mC.sorted.bed
echo -e "Distance\tSignal" > tmp.mC.txt
for((pos=-100; pos<=300; pos+=1)); do
	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' > tmp.txt
	sort -k1,1 -k2,2n -S 5G tmp.txt | bedtools intersect -sorted -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		oppoM=`awk -v pos=$pos '{mSum+=$4} END{print mSum/NR}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${oppoM}" >> tmp.mC.txt
	fi
done
mv tmp.mC.txt H9_MNase_Alu_Start.txt


awk '$5>=0' /mnt/disk1/5/xx/private/Mhmei/H9/RNA-seq/TPM/SINE-Alu.bdg | awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6; else print $1, $3-1, $3, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n > tmp.motif
awk '{if($2>310) print $1, $2-310, $3+310}' OFS='\t' tmp.motif | sort -k1,1 -k2,2n | bedtools intersect -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/MNase-seq/bw/GM12878_MNase.bdg -b - > tmp.mC.sorted.bed
echo -e "Distance\tSignal" > tmp.mC.txt
for((pos=-100; pos<=300; pos+=1)); do
	awk '$2>310' tmp.motif | awk -v pos=$pos '{if($6=="+") print $1, $2+pos, $3+pos, $4, $5, $6; else print $1, $2-pos, $3-pos, $4, $5, $6}' OFS='\t' > tmp.txt
	sort -k1,1 -k2,2n -S 5G tmp.txt | bedtools intersect -sorted -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		oppoM=`awk -v pos=$pos '{mSum+=$4} END{print mSum/NR}' OFS='\t' tmp.bed `
		echo -e "${pos}\t${oppoM}" >> tmp.mC.txt
	fi
done
mv tmp.mC.txt GM12878_MNase_Alu_Start.txt
