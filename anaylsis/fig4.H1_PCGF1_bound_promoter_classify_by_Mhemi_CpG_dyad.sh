#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/30
##-------------------------

awk '{if($5>2 && $6=="+") print $1, $2-1000, $3+500, "m"NR, $5, $6, $2, $3; else if($5>2 && $6=="-") print $1, $2-500, $3+1000, "m"NR, $5, $6, $2, $3}' OFS='\t' /mnt/disk1/5/share/plot/data/H1_PCGF1_TSS_enrichment.bed > tmp.pcgf1.tss.bed
bedtools intersect -wb -wa -a tmp.pcgf1.tss.bed -b /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG_dyad.bed > tmp.dyad.bed

sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,7,8,4,5,6 -c 12,13,14,15 -o sum > H1_PCGF1_bound_promoter_Mhemi_CpG_dyad.bed

awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' H1_PCGF1_bound_promoter_Mhemi_CpG_dyad.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed

cut -f 4,7-11 H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > H1_PCGF1_bound_promoter_dyad_freq_classify_by_Mhemi_CpG_dyad_1.2_fold.txt
