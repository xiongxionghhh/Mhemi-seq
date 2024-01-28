#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/30
##-------------------------

awk '{if($6=="+") print $1, $2, $2+1, $4, $5, $6, $2, $3; else print $1, $3-1, $3, $4, $5, $6, $2, $3}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/H1/RNA-seq/TPM/H1_RNA_R1_RefGene.bed > tmp.Ref.bed
awk '$11=="cluster1"' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed | bedtools intersect -u -s -wa -a tmp.Ref.bed -b - | awk '{print $5, "cluster 1"}' OFS='\t' > H1_PCGF1_bound_promoter_gene_expression_classify_by_Mhemi_CpG_dyad.txt
awk '$11=="cluster2"' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed | bedtools intersect -u -s -wa -a tmp.Ref.bed -b - | awk '{print $5, "cluster 2"}' OFS='\t' >> H1_PCGF1_bound_promoter_gene_expression_classify_by_Mhemi_CpG_dyad.txt
awk '$11=="cluster3"' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed | bedtools intersect -u -s -wa -a tmp.Ref.bed -b - | awk '{print $5, "cluster 3"}' OFS='\t' >> H1_PCGF1_bound_promoter_gene_expression_classify_by_Mhemi_CpG_dyad.txt
awk '$11=="cluster4"' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed | bedtools intersect -u -s -wa -a tmp.Ref.bed -b - | awk '{print $5, "cluster 4"}' OFS='\t' >> H1_PCGF1_bound_promoter_gene_expression_classify_by_Mhemi_CpG_dyad.txt
