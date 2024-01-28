#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/30
##-------------------------

for i in PCGF1 KDM2B H2Aub RING1B BCOR RYBP; do
	awk '$11=="cluster1"' OFS='\t' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk1/5/xx/private/Mhmei/H1/ChIP-seq/bw/H1_${i}.bw -R tmp.bed -a 5000 -b 5000 -bs 50 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster1"}' OFS='\t' tmp.Coverage.matrix > H1_PCGF1_bound_promoter_${i}_enrichment_classify_by_Mhemi_CpG_dyad.txt
	awk '$11=="cluster2"' OFS='\t' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk1/5/xx/private/Mhmei/H1/ChIP-seq/bw/H1_${i}.bw -R tmp.bed -a 5000 -b 5000 -bs 50 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster2"}' OFS='\t' tmp.Coverage.matrix >> H1_PCGF1_bound_promoter_${i}_enrichment_classify_by_Mhemi_CpG_dyad.txt
	awk '$11=="cluster3"' OFS='\t' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk1/5/xx/private/Mhmei/H1/ChIP-seq/bw/H1_${i}.bw -R tmp.bed -a 5000 -b 5000 -bs 50 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster3"}' OFS='\t' tmp.Coverage.matrix >> H1_PCGF1_bound_promoter_${i}_enrichment_classify_by_Mhemi_CpG_dyad.txt
	awk '$11=="cluster4"' OFS='\t' H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_1.2_fold.bed > tmp.bed
	computeMatrix reference-point -S /mnt/disk1/5/xx/private/Mhmei/H1/ChIP-seq/bw/H1_${i}.bw -R tmp.bed -a 5000 -b 5000 -bs 50 -p 30 -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
	awk '{if(NR>3) print $0, "cluster4"}' OFS='\t' tmp.Coverage.matrix >> H1_PCGF1_bound_promoter_${i}_enrichment_classify_by_Mhemi_CpG_dyad.txt
done

