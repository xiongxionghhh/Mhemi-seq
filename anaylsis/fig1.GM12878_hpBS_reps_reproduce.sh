#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/16
##-------------------------

ln -f -s /mnt/disk1/5/share/hpBS-seq/GM12878/bedGraph/20220402_GM12878_hpBS_R1.intraCpG.bdg tmp.rep1.bed
ln -f -s /mnt/disk1/5/share/hpBS-seq/GM12878/bedGraph/20220602_GM12878_hpBS_R1.intraCpG.bdg tmp.rep23.bed

awk '{nSum=$4+$5+$6+$7} {if(nSum>5) print $1, $2, $3, 100*($5+$6)/nSum}' OFS='\t' tmp.rep1.bed > tmp.C5.rep1.bed
awk '{nSum=$4+$5+$6+$7} {if(nSum>5) print $1, $2, $3, 100*($5+$6)/nSum}' OFS='\t' tmp.rep23.bed > tmp.C5.rep23.bed

bedtools intersect -wb -sorted -a tmp.C5.rep1.bed -b tmp.C5.rep23.bed | cut -f 4,8 > GM12878_hpBS_reproduce_hemi_CpG_dyad_rep1_VS_rep2.txt


rm tmp*
