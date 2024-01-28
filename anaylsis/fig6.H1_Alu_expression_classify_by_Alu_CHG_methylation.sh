#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/8/9
##-------------------------

cut -f 4,8 /mnt/disk2-1/xx/private/Mhmei/H1/BS-seq/bedGraph/H1_SINE-Alu_Groups.bdg | sort -k1,1 > tmp.group.txt
intersectBed -s -sorted -F 1 -F 1 -wb -a /mnt/disk2-1/xx/private/Mhmei/H1/BS-seq/bedGraph/H1_SINE-Alu_mCTG.bdg -b /mnt/disk2-1/xx/private/Mhmei/H1/RNA-seq/TPM/H1_SINE-AluOnly.bdg | awk '{print $4, $1, $2, $3, $4, $23, $6}' OFS='\t' | sort -k1,1 | join -j 1 tmp.group.txt - | awk '{print $3, $4, $5, $6, $7, $8, $2}' OFS='\t' | sort -k1,1 -k2,2n > tmp.TPM.bed

awk '{if($7~/Group1/) print $5}' tmp.TPM.bed > tmp.Group1.TPM
awk '{if($7~/Group2/) print $5}' tmp.TPM.bed > tmp.Group2.TPM
awk '{if($7~/Group3/) print $5}' tmp.TPM.bed > tmp.Group3.TPM
awk '{if($7~/Group4/) print $5}' tmp.TPM.bed > tmp.Group4.TPM

Rscript Fig26-9.Boxplot_Expr.R

