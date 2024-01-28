#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/13
##-------------------------


awk 'BEGIN{print "Cell\tDyad\tCount"} {Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "GM12878\tunme\t"Un"\nGM12878\themi-Watson\t"HW"\nGM12878\themi-Crick\t"HC"\nGM12878\tme\t"Me""}' /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG_dyad.bed > H1_H9_GM12878_Mhemi_dyad_count.txt
awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "H1\tunme\t"Un"\nH1\themi-Watson\t"HW"\nH1\themi-Crick\t"HC"\nH1\tme\t"Me""}' /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG_dyad.bed >> H1_H9_GM12878_Mhemi_dyad_count.txt
awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "H9\tunme\t"Un"\nH9\themi-Watson\t"HW"\nH9\themi-Crick\t"HC"\nH9\tme\t"Me""}' /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG_dyad.bed >> H1_H9_GM12878_Mhemi_dyad_count.txt

awk 'BEGIN{print "Cell\tDyad\tCount"} {Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "GM12878\tunme\t"Un"\nGM12878\themi-Watson\t"HW"\nGM12878\themi-Crick\t"HC"\nGM12878\tme\t"Me""}' /mnt/disk1/5/share/hpBS-seq/GM12878/bedGraph/GM12878_hpBS_rep123_YNCGNR_dyad.bed > H1_H9_GM12878_hpBS_dyad_count.txt
awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "H1\tunme\t"Un"\nH1\themi-Watson\t"HW"\nH1\themi-Crick\t"HC"\nH1\tme\t"Me""}' /mnt/disk1/5/share/hpBS-seq/H1/bedGraph/H1_hpBS_YNCGNR_dyad.bed >> H1_H9_GM12878_hpBS_dyad_count.txt
awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "H9\tunme\t"Un"\nH9\themi-Watson\t"HW"\nH9\themi-Crick\t"HC"\nH9\tme\t"Me""}' /mnt/disk1/5/share/hpBS-seq/H9/bedGraph/20220704_H9_hpBS_R1.intraCpG.bdg >> H1_H9_GM12878_hpBS_dyad_count.txt











awk '{nSum+=$4+$5+$5+$7} {print $1, $2, $3, 100*($5+$6)/nSum}' OFS='\t' /mnt/disk1/5/share/hpBS-seq/GM12878/bedGraph/GM12878_hpBS_rep123_YNCGNR_dyad.bed > /mnt/disk1/5/share/upload_data/bw/GM12878_hpBS_rep123_hemi_YNCGNR.bed
awk '{nSum+=$4+$5+$5+$7} {print $1, $2, $3, 100*($5+$6)/nSum}' OFS='\t' /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG_dyad.bed > /mnt/disk1/5/share/upload_data/bw/GM12878_Mhemi_rep123_hemi_YNCGNR.bed


bedGraphToBigWig /mnt/disk1/5/share/upload_data/bw/GM12878_hpBS_rep123_hemi_YNCGNR.bed /mnt/disk1/5/share/Reference/chromosizes/hg38_chromosome_sizes.txt /mnt/disk1/5/share/upload_data/bw/GM12878_hpBS_rep123_hemi_YNCGNR.bw
bedGraphToBigWig /mnt/disk1/5/share/upload_data/bw/GM12878_Mhemi_rep123_hemi_YNCGNR.bed /mnt/disk1/5/share/Reference/chromosizes/hg38_chromosome_sizes.txt /mnt/disk1/5/share/upload_data/bw/GM12878_Mhemi_rep123_hemi_YNCGNR.bw















