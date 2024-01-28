#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/30
##-------------------------

awk '$1!~/[_,Y,M,L]/' /mnt/disk4/public/RefBed/refGene_hg38/hg38_TSS.bed > tmp.TSS.bed
awk '{if($6=="+") print $1, $2-1000, $3+500, NR; else print $1, $2-500, $3+1000, NR}' OFS='\t' tmp.TSS.bed > tmp.bed
bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H1/ChIP-seq/bw/H1_PCGF1.bw tmp.bed tmp.tab
paste tmp.TSS.bed tmp.tab | awk '{print $1, $2, $3, $4, $11, $6}' OFS='\t' > H1_PCGF1_TSS_enrichment.bed

echo -e 'Pcgf\tDyad\tFreq' > H1_promoter_Mhemi_CpG_dyad_count_classify_by_PCGF1.txt
awk '$5>5' H1_PCGF1_TSS_enrichment.bed | awk '{if($6=="+") print $1, $2-1000, $3+500, $4, $5, $6; else print $1, $2-500, $3+1000, $4, $5, $6}' OFS='\t' | bedtools intersect -wb -a /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG_dyad.bed -b - | cut -f 1-7,13 > tmp.dyad.bed
awk '{nSum+=$4+$5+$6+$7; Un+=$4; Me+=$7} {if($8=="+") {HOPPO+=$6; HSAME+=$5} else{HOPPO+=$5; HSAME+=$6}} END{print "strong\tunme\t"100*Un/nSum"\nstrong\tsame\t"100*HSAME/nSum"\nstrong\toppo\t"100*HOPPO/nSum"\nstrong\tme\t"100*Me/nSum""}' tmp.dyad.bed >> H1_promoter_Mhemi_CpG_dyad_count_classify_by_PCGF1.txt

awk '$5>1 && $5<2' H1_PCGF1_TSS_enrichment.bed | awk '{if($6=="+") print $1, $2-1000, $3+500, $4, $5, $6; else print $1, $2-500, $3+1000, $4, $5, $6}' OFS='\t' | bedtools intersect -wb -a /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG_dyad.bed -b - | cut -f 1-7,13 > tmp.dyad.bed
awk '{nSum+=$4+$5+$6+$7; Un+=$4; Me+=$7} {if($8=="+") {HOPPO+=$6; HSAME+=$5} else{HOPPO+=$5; HSAME+=$6}} END{print "weak\tunme\t"100*Un/nSum"\nweak\tsame\t"100*HSAME/nSum"\nweak\toppo\t"100*HOPPO/nSum"\nweak\tme\t"100*Me/nSum""}' tmp.dyad.bed >> H1_promoter_Mhemi_CpG_dyad_count_classify_by_PCGF1.txt


awk '$5<0.1' H1_PCGF1_TSS_enrichment.bed | awk '{if($6=="+") print $1, $2-1000, $3+500, $4, $5, $6; else print $1, $2-500, $3+1000, $4, $5, $6}' OFS='\t' | bedtools intersect -wb -a /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG_dyad.bed -b - | cut -f 1-7,13 > tmp.dyad.bed
awk '{nSum+=$4+$5+$6+$7; Un+=$4; Me+=$7} {if($8=="+") {HOPPO+=$6; HSAME+=$5} else{HOPPO+=$5; HSAME+=$6}} END{print "unbound\tunme\t"100*Un/nSum"\nunbound\tsame\t"100*HSAME/nSum"\nunbound\toppo\t"100*HOPPO/nSum"\nunbound\tme\t"100*Me/nSum""}' tmp.dyad.bed >> H1_promoter_Mhemi_CpG_dyad_count_classify_by_PCGF1.txt


echo -e 'Pcgf\tDyad\tFreq' > H1_promoter_hpBS_CpG_dyad_count_classify_by_PCGF1.txt
awk '$5>2' H1_PCGF1_TSS_enrichment.bed | awk '{if($6=="+") print $1, $2-1000, $3+500, $4, $5, $6; else print $1, $2-500, $3+1000, $4, $5, $6}' OFS='\t' | bedtools intersect -wb -a /mnt/disk1/5/share/hpBS-seq/H1/bedGraph/20220704_H1_hpBS_R1.intraCpG.bdg -b - | cut -f 1-7,13 > tmp.dyad.bed
awk '{nSum+=$4+$5+$6+$7; Un+=$4; Me+=$7} {if($8=="+") {HOPPO+=$6; HSAME+=$5} else{HOPPO+=$5; HSAME+=$6}} END{print "bound\tunme\t"100*Un/nSum"\nbound\tsame\t"100*HSAME/nSum"\nbound\toppo\t"100*HOPPO/nSum"\nbound\tme\t"100*Me/nSum""}' tmp.dyad.bed >> H1_promoter_hpBS_CpG_dyad_count_classify_by_PCGF1.txt

awk '$5<0.1' H1_PCGF1_TSS_enrichment.bed | awk '{if($6=="+") print $1, $2-1000, $3+500, $4, $5, $6; else print $1, $2-500, $3+1000, $4, $5, $6}' OFS='\t' | bedtools intersect -wb -a /mnt/disk1/5/share/hpBS-seq/H1/bedGraph/20220704_H1_hpBS_R1.intraCpG.bdg -b - | cut -f 1-7,13 > tmp.dyad.bed
awk '{nSum+=$4+$5+$6+$7; Un+=$4; Me+=$7} {if($8=="+") {HOPPO+=$6; HSAME+=$5} else{HOPPO+=$5; HSAME+=$6}} END{print "unbound\tunme\t"100*Un/nSum"\nunbound\tsame\t"100*HSAME/nSum"\nunbound\toppo\t"100*HOPPO/nSum"\nunbound\tme\t"100*Me/nSum""}' tmp.dyad.bed >> H1_promoter_hpBS_CpG_dyad_count_classify_by_PCGF1.txt

