#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/10/25
##-------------------------


bedtools intersect -sorted -u -a /mnt/disk1/5/share/hpBS-seq/GM12878/bedGraph/GM12878_hpBS_rep123_CpG_dyad.bed -b /mnt/disk1/5/share/Reference/bed/hg38_CGNR_Cri.bed /mnt/disk1/5/share/Reference/bed/hg38_CGNR_Wat.bed > GM12878_hpBS_CGNR_dyad.bdg
bedtools intersect -sorted -u -a /mnt/disk1/5/share/hpBS-seq/H1/bedGraph/20220704_H1_hpBS_R1.intraCpG.bdg -b /mnt/disk1/5/share/Reference/bed/hg38_CGNR_Cri.bed /mnt/disk1/5/share/Reference/bed/hg38_CGNR_Wat.bed > H1_hpBS_CGNR_dyad.bdg
bedtools intersect -sorted -u -a /mnt/disk1/5/share/hpBS-seq/H9/bedGraph/20220704_H9_hpBS_R1.intraCpG.bdg -b /mnt/disk1/5/share/Reference/bed/hg38_CGNR_Cri.bed /mnt/disk1/5/share/Reference/bed/hg38_CGNR_Wat.bed > H9_hpBS_CGNR_dyad.bdg



bedtools intersect -sorted -u -a /mnt/disk1/5/share/hpBS-seq/GM12878/bedGraph/20220402_GM12878_hpBS_R1.intraCWG.bdg -b /mnt/disk1/5/share/Reference/bed/hg38_CHGR_Cri.bed /mnt/disk1/5/share/Reference/bed/hg38_CHGR_Wat.bed > GM12878_hpBS_CWGR_dyad.bdg
bedtools intersect -sorted -u -a /mnt/disk1/5/share/hpBS-seq/H1/bedGraph/20220704_H1_hpBS_R1.intraCWG.bdg -b /mnt/disk1/5/share/Reference/bed/hg38_CHGR_Cri.bed /mnt/disk1/5/share/Reference/bed/hg38_CHGR_Wat.bed > H1_hpBS_CWGR_dyad.bdg
bedtools intersect -sorted -u -a /mnt/disk1/5/share/hpBS-seq/H9/bedGraph/20220704_H9_hpBS_R1.intraCWG.bdg -b /mnt/disk1/5/share/Reference/bed/hg38_CHGR_Cri.bed /mnt/disk1/5/share/Reference/bed/hg38_CHGR_Wat.bed > H9_hpBS_CWGR_dyad.bdg



echo -e "Cell\tDyad\tCount" > H1_H9_GM12878_hpBS_CGNR_dyad_count.txt
awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "GM12878\tunme\t"Un"\nGM12878\themi-Watson\t"HW"\nGM12878\themi-Crick\t"HC"\nGM12878\tme\t"Me""}' GM12878_hpBS_CGNR_dyad.bdg >> H1_H9_GM12878_hpBS_CGNR_dyad_count.txt
awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "H1\tunme\t"Un"\nH1\themi-Watson\t"HW"\nH1\themi-Crick\t"HC"\nH1\tme\t"Me""}' H1_hpBS_CGNR_dyad.bdg >> H1_H9_GM12878_hpBS_CGNR_dyad_count.txt
awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "H9\tunme\t"Un"\nH9\themi-Watson\t"HW"\nH9\themi-Crick\t"HC"\nH9\tme\t"Me""}' H9_hpBS_CGNR_dyad.bdg >> H1_H9_GM12878_hpBS_CGNR_dyad_count.txt



echo -e "Cell\tDyad\tCount" > H1_H9_GM12878_hpBS_CWGR_dyad_count.txt
awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "GM12878\tunme\t"Un"\nGM12878\themi-Watson\t"HW"\nGM12878\themi-Crick\t"HC"\nGM12878\tme\t"Me""}' GM12878_hpBS_CWGR_dyad.bdg >> H1_H9_GM12878_hpBS_CWGR_dyad_count.txt
awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "H1\tunme\t"Un"\nH1\themi-Watson\t"HW"\nH1\themi-Crick\t"HC"\nH1\tme\t"Me""}' H1_hpBS_CWGR_dyad.bdg >> H1_H9_GM12878_hpBS_CWGR_dyad_count.txt
awk '{Un+=$4; HC+=$6; HW+=$5; Me+=$7} END{print "H9\tunme\t"Un"\nH9\themi-Watson\t"HW"\nH9\themi-Crick\t"HC"\nH9\tme\t"Me""}' H9_hpBS_CWGR_dyad.bdg >> H1_H9_GM12878_hpBS_CWGR_dyad_count.txt


