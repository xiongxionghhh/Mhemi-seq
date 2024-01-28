#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/16
##-------------------------

ln -f -s /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep1/bed/GM12878_Mhemi_rep1_CpG.bed tmp.rep1.bed
ln -f -s /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep23/bed/GM12878_Mhemi_rep23_CpG.bed tmp.rep23.bed

awk '{if($5>4) print $1, $2, $3, 100*$4/$5, 0, $6}' OFS='\t' tmp.rep1.bed > tmp.C5.rep1.bed
awk '{if($5>4) print $1, $2, $3, 100*$4/$5, 0, $6}' OFS='\t' tmp.rep23.bed > tmp.C5.rep23.bed

bedtools intersect -s -wb -sorted -a tmp.C5.rep1.bed -b tmp.C5.rep23.bed | cut -f 4,10 > GM12878_Mhemi_reproduce_CpG_rep1_VS_rep23.txt

ln -f -s /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep1/bed/GM12878_Mhemi_rep1_CpG_dyad.bed tmp.C4.rep1.dyad.bed
ln -f -s /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep23/bed/GM12878_Mhemi_rep23_CpG_dyad.bed tmp.C4.rep23.dyad.bed

awk '{if(($4+$5+$6+$7)>3) print $1, $2, $3, 100*$4/($4+$5+$6+$7), 100*$5/($4+$5+$6+$7), 100*$6/($4+$5+$6+$7), 100*$7/($4+$5+$6+$7)}' OFS='\t' tmp.C4.rep1.dyad.bed > tmp.rep1.bed
awk '{if(($4+$5+$6+$7)>3) print $1, $2, $3, 100*$4/($4+$5+$6+$7), 100*$5/($4+$5+$6+$7), 100*$6/($4+$5+$6+$7), 100*$7/($4+$5+$6+$7)}' OFS='\t' tmp.C4.rep23.dyad.bed > tmp.rep23.bed

bedtools intersect -wb -sorted -a tmp.rep1.bed -b tmp.rep23.bed > tmp.intersect.bed
cut -f 4,11 tmp.intersect.bed > GM12878_Mhemi_reproduce_CpG_dyad_unme_rep1_VS_rep23.txt
cut -f 5,12 tmp.intersect.bed > GM12878_Mhemi_reproduce_CpG_dyad_hemiW_rep1_VS_rep23.txt
cut -f 6,13 tmp.intersect.bed > GM12878_Mhemi_reproduce_CpG_dyad_hemiC_rep1_VS_rep23.txt
cut -f 7,14 tmp.intersect.bed > GM12878_Mhemi_reproduce_CpG_dyad_me_rep1_VS_rep23.txt

rm tmp*
