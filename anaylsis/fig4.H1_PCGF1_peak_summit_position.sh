#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/28
##-------------------------

bedtools intersect -sorted -f 1 -wb -a /mnt/disk1/5/xx/private/Mhmei/H1/ChIP-seq/peaks/PCGF1_summits.bed -b /mnt/disk1/5/share/Reference/bed/hg38_GenomicPositon.bed | cut -f 4,9 > tmp.Dyad.Dis.txt
awk '{if($2=="Intron") nIntron+=1; else if($2=="Promoter") nPromoter+=1; else if($2=="Exon") nExon+=1; else if($2=="Intergenic") nIntergenic+=1} END{print "intron\t"nIntron"\nexon\t"nExon"\nintergenic\t"nIntergenic"\npromoter\t"nPromoter}' tmp.Dyad.Dis.txt > H1_PCGF1_peak_summit_hg38_annotation.txt
sed -i '1i\Pos\tCount' H1_PCGF1_peak_summit_hg38_annotation.txt
