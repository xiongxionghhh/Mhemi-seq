#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/10
##-------------------------

cut -f 1-6 /mnt/disk4/public/RefBed/CTCF/GM12878_CTCF_RAD21_DNase.motif > tmp.GM12878_CTCF.motif
awk '{print $0, "m"NR}' OFS='\t' tmp.GM12878_CTCF.motif > tmp.CTCF.motif

awk '{if($6=="+") print $1, $2-5, $2-3, $4, $5, $6, $7, $2, $3; else print $1, $3+3, $3+5, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk1/5/share/hpBS-seq/GM12878/bedGraph/GM12878_hpBS_rep123_CpG_dyad.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > GM12878_CTCF_motif_hpBS_5th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' GM12878_CTCF_motif_hpBS_5th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > GM12878_CTCF_motif_classify_by_hpBS_5th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 GM12878_CTCF_motif_classify_by_hpBS_5th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > GM12878_CTCF_motif_classify_by_hpBS_5th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 GM12878_CTCF_motif_classify_by_hpBS_5th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > GM12878_CTCF_motif_signal_classify_by_hpBS_5th_CpG_dyad_count_1.2_fold.txt


awk '{if($6=="+") print $1, $2-5, $2-3, $4, $5, $6, $7, $2, $3; else print $1, $3+3, $3+5, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk4/public/GM12878/BS-seq/iSA/GM12878_intraCpG.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > GM12878_CTCF_motif_iSA_5th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' GM12878_CTCF_motif_iSA_5th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > GM12878_CTCF_motif_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 GM12878_CTCF_motif_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > GM12878_CTCF_motif_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 GM12878_CTCF_motif_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > GM12878_CTCF_motif_signal_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.txt


awk '{if($6=="+") print $1, $2-5, $2-3, $4, $5, $6, $7, $2, $3; else print $1, $3+3, $3+5, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG_dyad.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > GM12878_CTCF_motif_Mhemi_5th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' GM12878_CTCF_motif_Mhemi_5th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > GM12878_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 GM12878_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > GM12878_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 GM12878_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > GM12878_CTCF_motif_signal_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt


awk '{if($6=="+") print $1, $2-3, $2-1, $4, $5, $6, $7, $2, $3; else print $1, $3+1, $3+3, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk1/5/share/hpBS-seq/GM12878/bedGraph/GM12878_hpBS_rep123_CpG_dyad.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > GM12878_CTCF_motif_hpBS_7th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' GM12878_CTCF_motif_hpBS_7th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > GM12878_CTCF_motif_classify_by_hpBS_7th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 GM12878_CTCF_motif_classify_by_hpBS_7th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > GM12878_CTCF_motif_classify_by_hpBS_7th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 GM12878_CTCF_motif_classify_by_hpBS_7th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > GM12878_CTCF_motif_signal_classify_by_hpBS_7th_CpG_dyad_count_1.2_fold.txt


awk '{if($6=="+") print $1, $2-3, $2-1, $4, $5, $6, $7, $2, $3; else print $1, $3+1, $3+3, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk4/public/GM12878/BS-seq/iSA/GM12878_intraCpG.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > GM12878_CTCF_motif_iSA_7th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' GM12878_CTCF_motif_iSA_7th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > GM12878_CTCF_motif_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 GM12878_CTCF_motif_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > GM12878_CTCF_motif_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 GM12878_CTCF_motif_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > GM12878_CTCF_motif_signal_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.txt


awk '{if($6=="+") print $1, $2-3, $2-1, $4, $5, $6, $7, $2, $3; else print $1, $3+1, $3+3, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG_dyad.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > GM12878_CTCF_motif_Mhemi_7th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' GM12878_CTCF_motif_Mhemi_7th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > GM12878_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 GM12878_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > GM12878_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 GM12878_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > GM12878_CTCF_motif_signal_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt



cut -f 1-6 /mnt/disk4/public/RefBed/CTCF/H1_CTCF_RAD21_DNase.motif > tmp.H1_CTCF.motif
awk '{print $0, "m"NR}' OFS='\t' tmp.H1_CTCF.motif > tmp.CTCF.motif

awk '{if($6=="+") print $1, $2-5, $2-3, $4, $5, $6, $7, $2, $3; else print $1, $3+3, $3+5, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk4/public/H1/BS-seq/iSA/H1_intraCpG.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > H1_CTCF_motif_iSA_5th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' H1_CTCF_motif_iSA_5th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > H1_CTCF_motif_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 H1_CTCF_motif_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > H1_CTCF_motif_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 H1_CTCF_motif_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > H1_CTCF_motif_signal_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.txt


awk '{if($6=="+") print $1, $2-5, $2-3, $4, $5, $6, $7, $2, $3; else print $1, $3+3, $3+5, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG_dyad.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > H1_CTCF_motif_Mhemi_5th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' H1_CTCF_motif_Mhemi_5th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > H1_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 H1_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > H1_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 H1_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > H1_CTCF_motif_signal_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt


awk '{if($6=="+") print $1, $2-3, $2-1, $4, $5, $6, $7, $2, $3; else print $1, $3+1, $3+3, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk4/public/H1/BS-seq/iSA/H1_intraCpG.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > H1_CTCF_motif_iSA_7th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' H1_CTCF_motif_iSA_7th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > H1_CTCF_motif_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 H1_CTCF_motif_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > H1_CTCF_motif_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 H1_CTCF_motif_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > H1_CTCF_motif_signal_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.txt


awk '{if($6=="+") print $1, $2-3, $2-1, $4, $5, $6, $7, $2, $3; else print $1, $3+1, $3+3, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG_dyad.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > H1_CTCF_motif_Mhemi_7th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' H1_CTCF_motif_Mhemi_7th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > H1_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 H1_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > H1_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 H1_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > H1_CTCF_motif_signal_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt



cut -f 1-6 /mnt/disk4/public/RefBed/CTCF/H9_CTCF_DNase.motif > tmp.H9_CTCF.motif
awk '{print $0, "m"NR}' OFS='\t' tmp.H9_CTCF.motif > tmp.CTCF.motif

awk '{if($6=="+") print $1, $2-5, $2-3, $4, $5, $6, $7, $2, $3; else print $1, $3+3, $3+5, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk1/5/share/ChIP-hpBS-seq/H9/bed/H9_CTCF_hpBS_CpG_dyad.bdg > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > H9_CTCF_motif_ChIP-hpBS_5th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' H9_CTCF_motif_ChIP-hpBS_5th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > H9_CTCF_motif_classify_by_ChIP-hpBS_5th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 H9_CTCF_motif_classify_by_ChIP-hpBS_5th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > H9_CTCF_motif_classify_by_ChIP-hpBS_5th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 H9_CTCF_motif_classify_by_ChIP-hpBS_5th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > H9_CTCF_motif_signal_classify_by_ChIP-hpBS_5th_CpG_dyad_count_1.2_fold.txt

awk '{if($6=="+") print $1, $2-5, $2-3, $4, $5, $6, $7, $2, $3; else print $1, $3+3, $3+5, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG_dyad.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > H9_CTCF_motif_Mhemi_5th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' H9_CTCF_motif_Mhemi_5th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > H9_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 H9_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > H9_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 H9_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > H9_CTCF_motif_signal_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt

awk '{if($6=="+") print $1, $2-3, $2-1, $4, $5, $6, $7, $2, $3; else print $1, $3+1, $3+3, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk1/5/share/ChIP-hpBS-seq/H9/bed/H9_CTCF_hpBS_CpG_dyad.bdg > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > H9_CTCF_motif_ChIP-hpBS_7th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' H9_CTCF_motif_ChIP-hpBS_7th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > H9_CTCF_motif_classify_by_ChIP-hpBS_7th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 H9_CTCF_motif_classify_by_ChIP-hpBS_7th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > H9_CTCF_motif_classify_by_ChIP-hpBS_7th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 H9_CTCF_motif_classify_by_ChIP-hpBS_7th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > H9_CTCF_motif_signal_classify_by_ChIP-hpBS_7th_CpG_dyad_count_1.2_fold.txt

awk '{if($6=="+") print $1, $2-3, $2-1, $4, $5, $6, $7, $2, $3; else print $1, $3+1, $3+3, $4, $5, $6, $7, $2, $3}' OFS='\t' tmp.CTCF.motif | bedtools intersect -wb -f 1 -F 1 -a - -b /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG_dyad.bed > tmp.dyad.bed
sort -k1,1 -k2,2n tmp.dyad.bed | groupBy -i - -g 1,8,9,7,5,6 -c 13,14,15,16 -o sum > H9_CTCF_motif_Mhemi_7th_CpG_dyad_count.bed
awk '{if($6=="+") print $0; else print $1, $2, $3, $4, $5, $6, $7, $9, $8, $10}' OFS='\t' H9_CTCF_motif_Mhemi_7th_CpG_dyad_count.bed > tmp.bed
awk '{if($7>$8*1.2 && $7>$9*1.2 && $7>$10*1.2) print $0, "cluster1"; if($8>$7*1.2 && $8>$9*1.2 && $8>$10*1.2) print $0, "cluster2"; else if($9>$7*1.2 && $9>$8*1.2 && $9>$10*1.2) print $0, "cluster3"; else if($10>$7*1.2 && $10>$8*1.2 && $10>$9*1.2) print $0, "cluster4"}' OFS='\t' tmp.bed > H9_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.bed
cut -f 4,7-11 H9_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.bed | awk '{mSum=$2+$3+$4+$5} {print $1, 100*$2/mSum, 100*$3/mSum, 100*$4/mSum, 100*$5/mSum, $6}' OFS='\t' > H9_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt
awk '{print $7, $0}' OFS='\t' tmp.CTCF.motif | sort -k1,1 > tmp.CTCF.sorted.motif
sort -k1,1 H9_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt | join -j 1 tmp.CTCF.sorted.motif - | awk '{print $6, $13}' OFS='\t' > H9_CTCF_motif_signal_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt



for i in CTCF RAD21_MS TRIM22 TRIM28 ZNF143; do
	awk '$11=="cluster1"' OFS='\t' GM12878_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk4/public/GM12878/ChIP-seq/GM12878_${i}.bw -R tmp.bed -a 2000 -b 2000 -bs 20 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster1"}' OFS='\t' tmp.Coverage.matrix > GM12878_CTCF_motif_${i}_enrichment_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt
	awk '$11=="cluster2"' OFS='\t' GM12878_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk4/public/GM12878/ChIP-seq/GM12878_${i}.bw -R tmp.bed -a 2000 -b 2000 -bs 20 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster2"}' OFS='\t' tmp.Coverage.matrix >> GM12878_CTCF_motif_${i}_enrichment_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt
	awk '$11=="cluster3"' OFS='\t' GM12878_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk4/public/GM12878/ChIP-seq/GM12878_${i}.bw -R tmp.bed -a 2000 -b 2000 -bs 20 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster3"}' OFS='\t' tmp.Coverage.matrix >> GM12878_CTCF_motif_${i}_enrichment_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt
	awk '$11=="cluster4"' OFS='\t' GM12878_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk4/public/GM12878/ChIP-seq/GM12878_${i}.bw -R tmp.bed -a 2000 -b 2000 -bs 20 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster4"}' OFS='\t' tmp.Coverage.matrix >> GM12878_CTCF_motif_${i}_enrichment_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt

	awk '$11=="cluster1"' OFS='\t' GM12878_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk4/public/GM12878/ChIP-seq/GM12878_${i}.bw -R tmp.bed -a 2000 -b 2000 -bs 20 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster1"}' OFS='\t' tmp.Coverage.matrix > GM12878_CTCF_motif_${i}_enrichment_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt
	awk '$11=="cluster2"' OFS='\t' GM12878_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk4/public/GM12878/ChIP-seq/GM12878_${i}.bw -R tmp.bed -a 2000 -b 2000 -bs 20 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster2"}' OFS='\t' tmp.Coverage.matrix >> GM12878_CTCF_motif_${i}_enrichment_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt
	awk '$11=="cluster3"' OFS='\t' GM12878_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk4/public/GM12878/ChIP-seq/GM12878_${i}.bw -R tmp.bed -a 2000 -b 2000 -bs 20 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster3"}' OFS='\t' tmp.Coverage.matrix >> GM12878_CTCF_motif_${i}_enrichment_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt
	awk '$11=="cluster4"' OFS='\t' GM12878_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk4/public/GM12878/ChIP-seq/GM12878_${i}.bw -R tmp.bed -a 2000 -b 2000 -bs 20 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster4"}' OFS='\t' tmp.Coverage.matrix >> GM12878_CTCF_motif_${i}_enrichment_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt
done
