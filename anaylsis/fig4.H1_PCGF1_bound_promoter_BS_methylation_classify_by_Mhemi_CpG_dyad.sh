#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/30
##-------------------------

awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6, $2, $3; else print $1, $3-1, $3, $4, $5, $6, $2, $3}' OFS='\t' /mnt/disk4/public/RefBed/hg38_refGene.bed | sort -k1,1 -k2,2n > tmp.Ref.bed

awk '$11=="cluster1"' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed | bedtools intersect -u -s -wa -a tmp.Ref.bed -b - | awk '{print $1, $7, $8, $4, $5, $6}' OFS='\t' > tmp.bed
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CpG.bw -R tmp.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "cluster 1", $0}' OFS='\t' tmp.Coverage.matrix > H1_PCGF1_bound_promoter_BS_methylation.txt
awk '$11=="cluster2"' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed | bedtools intersect -u -s -wa -a tmp.Ref.bed -b - | awk '{print $1, $7, $8, $4, $5, $6}' OFS='\t' > tmp.bed
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CpG.bw -R tmp.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "cluster 2", $0}' OFS='\t' tmp.Coverage.matrix >> H1_PCGF1_bound_promoter_BS_methylation.txt
awk '$11=="cluster3"' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed | bedtools intersect -u -s -wa -a tmp.Ref.bed -b - | awk '{print $1, $7, $8, $4, $5, $6}' OFS='\t' > tmp.bed
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CpG.bw -R tmp.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "cluster 3", $0}' OFS='\t' tmp.Coverage.matrix >> H1_PCGF1_bound_promoter_BS_methylation.txt
awk '$11=="cluster4"' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed | bedtools intersect -u -s -wa -a tmp.Ref.bed -b - | awk '{print $1, $7, $8, $4, $5, $6}' OFS='\t' > tmp.bed
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CpG.bw -R tmp.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "cluster 4", $0}' OFS='\t' tmp.Coverage.matrix >> H1_PCGF1_bound_promoter_BS_methylation.txt

