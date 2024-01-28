#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/23
##-------------------------

awk '$5>1' /mnt/disk4/public/RefBed/CTCF/GM12878_CTCF_RAD21_DNase.motif | sort -k1,1 -k2,2n -S 5G -u  > tmp.CTCF.motif

awk '{print $1, $2, $2+1, $5, $4+$5+$6+$7, "+"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/GM12878_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_W.bed
awk '{print $1, $3-1, $3, $6, $4+$5+$6+$7, "-"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/GM12878_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_C.bed
cat tmp.Dyad_W.bed tmp.Dyad_C.bed | sort -k1,1 -k2,2n -S 5G -u  > tmp.Dyad.bed

awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | intersectBed -u -sorted -a tmp.Dyad.bed -b - > tmp.mC.sorted.bed

echo -e "Distance\tMethylation" > GM12878_Mhemi_mYCWGR_CTCF_phasing_SAME.txt
echo -e "Distance\tMethylation" > GM12878_Mhemi_mYCWGR_CTCF_phasing_OPPO.txt
for((pos=-1000; pos<=1000; pos+=5)); do
	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | intersectBed -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount}' OFS='\t' tmp.bed >> GM12878_Mhemi_mYCWGR_CTCF_phasing_SAME.txt
	else
		echo -e "$pos\tNA" >> GM12878_Mhemi_mYCWGR_CTCF_phasing_SAME.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | intersectBed -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount}' OFS='\t' tmp.bed >> GM12878_Mhemi_mYCWGR_CTCF_phasing_OPPO.txt
	else
		echo -e "$pos\tNA" >> GM12878_Mhemi_mYCWGR_CTCF_phasing_OPPO.txt
	fi

done
echo "done."





awk '$5>1' /mnt/disk4/public/RefBed/CTCF/H1_CTCF_RAD21_DNase.motif | sort -k1,1 -k2,2n -S 5G -u > tmp.CTCF.motif
awk '{print $1, $2, $2+1, $5, $4+$5+$6+$7, "+"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/H1_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_W.bed
awk '{print $1, $3-1, $3, $6, $4+$5+$6+$7, "-"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/H1_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_C.bed
cat tmp.Dyad_W.bed tmp.Dyad_C.bed | sort -k1,1 -k2,2n -S 5G -u  > tmp.Dyad.bed

awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | intersectBed -u -sorted -a tmp.Dyad.bed -b - > tmp.mC.sorted.bed

echo -e "Distance\tMethylation" > H1_Mhemi_mYCWGR_CTCF_phasing_SAME.txt
echo -e "Distance\tMethylation" > H1_Mhemi_mYCWGR_CTCF_phasing_OPPO.txt
for((pos=-1000; pos<=1000; pos+=5)); do
	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | intersectBed -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount}' OFS='\t' tmp.bed >> H1_Mhemi_mYCWGR_CTCF_phasing_SAME.txt
	else
		echo -e "$pos\tNA" >> H1_Mhemi_mYCWGR_CTCF_phasing_SAME.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | intersectBed -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount}' OFS='\t' tmp.bed >> H1_Mhemi_mYCWGR_CTCF_phasing_OPPO.txt
	else
		echo -e "$pos\tNA" >> H1_Mhemi_mYCWGR_CTCF_phasing_OPPO.txt
	fi

done
echo "done."




awk '$5>1' /mnt/disk4/public/RefBed/CTCF/H9_CTCF_DNase.motif | sort -k1,1 -k2,2n -S 5G -u > tmp.CTCF.motif
awk '{print $1, $2, $2+1, $5, $4+$5+$6+$7, "+"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/H9_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_W.bed
awk '{print $1, $3-1, $3, $6, $4+$5+$6+$7, "-"}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/CH/H9_Mhemi_rep123_CWG_dyad.bed > tmp.Dyad_C.bed
cat tmp.Dyad_W.bed tmp.Dyad_C.bed | sort -k1,1 -k2,2n -S 5G -u  > tmp.Dyad.bed

awk '{if($2>1010) print $1, $2-1010, $3+1010}' OFS='\t' tmp.CTCF.motif | intersectBed -u -sorted -a tmp.Dyad.bed -b - > tmp.mC.sorted.bed

echo -e "Distance\tMethylation" > H9_Mhemi_mYCWGR_CTCF_phasing_SAME.txt
echo -e "Distance\tMethylation" > H9_Mhemi_mYCWGR_CTCF_phasing_OPPO.txt
for((pos=-1000; pos<=1000; pos+=5)); do
	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | intersectBed -sorted -s -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount}' OFS='\t' tmp.bed >> H9_Mhemi_mYCWGR_CTCF_phasing_SAME.txt
	else
		echo -e "$pos\tNA" >> H9_Mhemi_mYCWGR_CTCF_phasing_SAME.txt
	fi

	awk '$2>1010' tmp.CTCF.motif | awk -v pos=$pos '{if($6=="+") print $1, $2-2+pos, $3+2+pos, $4, $5, $6; else print $1, $2-2-pos, $3+2-pos, $4, $5, $6}' OFS='\t' | sort -k1,1 -k2,2n -S 5G | intersectBed -sorted -S -a tmp.mC.sorted.bed -b - > tmp.bed
	if [ -s tmp.bed ]; then
		awk -v pos=$pos 'BEGIN{mSum=0; nCount=0} {mSum+=$4; nCount+=$5} END{print pos, 100*mSum/nCount}' OFS='\t' tmp.bed >> H9_Mhemi_mYCWGR_CTCF_phasing_OPPO.txt
	else
		echo -e "$pos\tNA" >> H9_Mhemi_mYCWGR_CTCF_phasing_OPPO.txt
	fi

done
echo "done."
