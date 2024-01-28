#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/10/25
##-------------------------

samtools view -@ 30 /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bam/GM12878_Mhemi_rep123.bam | awk '$9>0 {print $9}' > GM12878_Mhemi_insert_size.txt
samtools view -@ 30 /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bam/H1_Mhemi_rep123.bam | awk '$9>0 {print $9}' > H1_Mhemi_insert_size.txt
samtools view -@ 30 /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bam/H9_Mhemi_rep123.bam | awk '$9>0 {print $9}' > H9_Mhemi_insert_size.txt
