#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/21
##-------------------------


awk '$5>1' /mnt/disk4/public/RefBed/CTCF/GM12878_CTCF_RAD21_DNase.motif | sort -k1,1 -k2,2n -S 5G -u > tmp.CTCF.motif
awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | bedtools intersect -sorted -a /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG.bed -b - | awk '{print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' > tmp.unsort.bed
sort -k1,1 -k2,2n -S 5G -u tmp.unsort.bed > tmp.mC.sorted.bed
awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | bedtools intersect -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/MNase-seq/bw/GM12878_MNase.bdg -b - > tmp.unsort.bed
sort -k1,1 -k2,2n -S 5G -u tmp.unsort.bed > tmp.MNase.bed

echo -e "Distance\tMethylation" > GM12878_Mhemi_CpG_CTCF_phasing_SAME.txt
echo -e "Distance\tMethylation" > GM12878_Mhemi_CpG_CTCF_phasing_OPPO.txt
echo -e "Distance\tSignal" > GM12878_MNase_CTCF_phasing.txt

for((pos=-1000; pos<=1000; pos+=5)); do
	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - | sort -k1,1 -k2,2n -S 5G -u > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> GM12878_Mhemi_CpG_CTCF_phasing_SAME.txt
	else
		echo -e "$pos\tNA" >> GM12878_Mhemi_CpG_CTCF_phasing_SAME.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u  | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - | sort -k1,1 -k2,2n -S 5G -u > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> GM12878_Mhemi_CpG_CTCF_phasing_OPPO.txt
	else
		echo -e "$pos\tNA" >> GM12878_Mhemi_CpG_CTCF_phasing_OPPO.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos; else print $1, $2-2-pos, $3+2-pos}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u | bedtools intersect -sorted -a tmp.MNase.bed -b - | sort -k1,1 -k2,2n -S 5G -u > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0} {mSum+=$4} END{print pos, mSum/NR}' OFS='\t' tmp.bed >> GM12878_MNase_CTCF_phasing.txt
	else
		echo -e "$pos\tNA" >> GM12878_MNase_CTCF_phasing.txt
	fi
done


awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | bedtools intersect -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_BS_CGNR.bdg -b - > tmp.unsort.bed
sort -k1,1 -k2,2n -S 5G -u tmp.unsort.bed > tmp.mC.sorted.bed

echo -e "Distance\tMethylation" > GM12878_BS_CGNR_CTCF_phasing_SAME.txt
echo -e "Distance\tMethylation" > GM12878_BS_CGNR_CTCF_phasing_OPPO.txt

for((pos=-1000; pos<=1000; pos+=5)); do
	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> GM12878_BS_CGNR_CTCF_phasing_SAME.txt
	else
		echo -e "$pos\tNA" >> GM12878_BS_CGNR_CTCF_phasing_SAME.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u  | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> GM12878_BS_CGNR_CTCF_phasing_OPPO.txt
	else
		echo -e "$pos\tNA" >> GM12878_BS_CGNR_CTCF_phasing_OPPO.txt
	fi

done


awk '$5>1' /mnt/disk4/public/RefBed/CTCF/H1_CTCF_RAD21_DNase.motif | sort -k1,1 -k2,2n -S 5G -u > tmp.CTCF.motif
awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | bedtools intersect -sorted -a /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bed -b - | awk '{print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' > tmp.unsort.bed
sort -k1,1 -k2,2n -S 5G -u tmp.unsort.bed > tmp.mC.sorted.bed
awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | bedtools intersect -sorted -a /mnt/disk1/5/xx/private/Mhmei/H1/MNase-seq/bedGraph/H1_MNase-seq_134-160.bdg -b - > tmp.unsort.bed
sort -k1,1 -k2,2n -S 5G -u tmp.unsort.bed > tmp.MNase.bed

echo -e "Distance\tMethylation" > H1_Mhemi_CpG_CTCF_phasing_SAME.txt
echo -e "Distance\tMethylation" > H1_Mhemi_CpG_CTCF_phasing_OPPO.txt
echo -e "Distance\tSignal" > H1_MNase_CTCF_phasing.txt

for((pos=-1000; pos<=1000; pos+=5)); do
	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - | sort -k1,1 -k2,2n -S 5G -u > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> H1_Mhemi_CpG_CTCF_phasing_SAME.txt
	else
		echo -e "$pos\tNA" >> H1_Mhemi_CpG_CTCF_phasing_SAME.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u  | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - | sort -k1,1 -k2,2n -S 5G -u > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> H1_Mhemi_CpG_CTCF_phasing_OPPO.txt
	else
		echo -e "$pos\tNA" >> H1_Mhemi_CpG_CTCF_phasing_OPPO.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos; else print $1, $2-2-pos, $3+2-pos}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u | bedtools intersect -sorted -a tmp.MNase.bed -b - | sort -k1,1 -k2,2n -S 5G -u > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0} {mSum+=$4} END{print pos, mSum/NR}' OFS='\t' tmp.bed >> H1_MNase_CTCF_phasing.txt
	else
		echo -e "$pos\tNA" >> H1_MNase_CTCF_phasing.txt
	fi
done


awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | bedtools intersect -sorted -a /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bedGraph/H1_BS_CGNR.bdg -b - > tmp.unsort.bed
sort -k1,1 -k2,2n -S 5G -u tmp.unsort.bed > tmp.mC.sorted.bed

echo -e "Distance\tMethylation" > H1_BS_CGNR_CTCF_phasing_SAME.txt
echo -e "Distance\tMethylation" > H1_BS_CGNR_CTCF_phasing_OPPO.txt

for((pos=-1000; pos<=1000; pos+=5)); do
	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> H1_BS_CGNR_CTCF_phasing_SAME.txt
	else
		echo -e "$pos\tNA" >> H1_BS_CGNR_CTCF_phasing_SAME.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u  | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> H1_BS_CGNR_CTCF_phasing_OPPO.txt
	else
		echo -e "$pos\tNA" >> H1_BS_CGNR_CTCF_phasing_OPPO.txt
	fi

done



awk '$5>1' /mnt/disk4/public/RefBed/CTCF/H9_CTCF_DNase.motif | sort -k1,1 -k2,2n -S 5G -u > tmp.CTCF.motif
awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | bedtools intersect -sorted -a /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG.bed -b - | awk '{print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' > tmp.unsort.bed
sort -k1,1 -k2,2n -S 5G -u tmp.unsort.bed > tmp.mC.sorted.bed
awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | bedtools intersect -sorted -a /mnt/disk1/5/xx/private/Mhmei/H9/MNase-seq/H9_MNase.bdg -b - > tmp.unsort.bed
sort -k1,1 -k2,2n -S 5G -u tmp.unsort.bed > tmp.MNase.bed

echo -e "Distance\tMethylation" > H9_Mhemi_CpG_CTCF_phasing_SAME.txt
echo -e "Distance\tMethylation" > H9_Mhemi_CpG_CTCF_phasing_OPPO.txt
echo -e "Distance\tSignal" > H9_MNase_CTCF_phasing.txt

for((pos=-1000; pos<=1000; pos+=5)); do
	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - | sort -k1,1 -k2,2n -S 5G -u > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> H9_Mhemi_CpG_CTCF_phasing_SAME.txt
	else
		echo -e "$pos\tNA" >> H9_Mhemi_CpG_CTCF_phasing_SAME.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u  | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - | sort -k1,1 -k2,2n -S 5G -u > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> H9_Mhemi_CpG_CTCF_phasing_OPPO.txt
	else
		echo -e "$pos\tNA" >> H9_Mhemi_CpG_CTCF_phasing_OPPO.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos; else print $1, $2-2-pos, $3+2-pos}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u | bedtools intersect -sorted -a tmp.MNase.bed -b - | sort -k1,1 -k2,2n -S 5G -u > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0} {mSum+=$4} END{print pos, mSum/NR}' OFS='\t' tmp.bed >> H9_MNase_CTCF_phasing.txt
	else
		echo -e "$pos\tNA" >> H9_MNase_CTCF_phasing.txt
	fi
done


awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | bedtools intersect -sorted -a /mnt/disk1/5/xx/private/Mhmei/H9/BS-seq/bedGraph/H9_BS_CGNR.bdg -b - > tmp.unsort.bed
sort -k1,1 -k2,2n -S 5G -u tmp.unsort.bed > tmp.mC.sorted.bed

echo -e "Distance\tMethylation" > H9_BS_CGNR_CTCF_phasing_SAME.txt
echo -e "Distance\tMethylation" > H9_BS_CGNR_CTCF_phasing_OPPO.txt

for((pos=-1000; pos<=1000; pos+=5)); do
	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u | bedtools intersect -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> H9_BS_CGNR_CTCF_phasing_SAME.txt
	else
		echo -e "$pos\tNA" >> H9_BS_CGNR_CTCF_phasing_SAME.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G -u  | bedtools intersect -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4*$5; nCount+=$5} END{print pos, mSum/nCount}' OFS='\t' tmp.bed >> H9_BS_CGNR_CTCF_phasing_OPPO.txt
	else
		echo -e "$pos\tNA" >> H9_BS_CGNR_CTCF_phasing_OPPO.txt
	fi

done
