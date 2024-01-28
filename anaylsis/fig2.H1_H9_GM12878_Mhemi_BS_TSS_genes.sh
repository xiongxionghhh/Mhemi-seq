#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/22
##-------------------------

awk '{if($3-$2>1500) print $1, $2, $3, $4, "gene"NR, $6}' OFS='\t' /mnt/disk4/public/RefBed/refGene_hg38/hg38_refGene_noOL.bed > tmp.hg38.genes.bed
awk '{if($6=="+") print $1, $2-1500, $2, $4, $5, $6; else print $1, $3, $3+1500, $4, $5, $6}' OFS='\t' tmp.hg38.genes.bed > tmp.hg38.tss.bed


bedtools makewindows -b tmp.hg38.genes.bed -n 300 > hg38_genes_300bins.bed
bedtools makewindows -b tmp.hg38.tss.bed -n 300 > hg38_TSS_300bins.bed


awk '{print $1, $2, $3, "GenesBin"NR}' OFS='\t' hg38_genes_300bins.bed > tmp.genes.300bins.bed
awk '{print $1, $2, $3, "TSSBin"NR}' OFS='\t' hg38_TSS_300bins.bed > tmp.tss.300bins.bed


ln -s -f /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG.bw tmp.GM12878.Mhemi.bw
ln -s -f /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bw/GM12878_BS_CGNR.bw tmp.GM12878.BS.bw

ln -s -f /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bw tmp.H1.Mhemi.bw
ln -s -f /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CGNR.bw tmp.H1.BS.bw

ln -s -f /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG.bw tmp.H9.Mhemi.bw
ln -s -f /mnt/disk1/5/xx/private/Mhmei/H9/BS-seq/bw/H9_BS_CGNR.bw tmp.H9.BS.bw

for Cell in GM12878 H1 H9; do
	for Method in Mhemi BS; do
		bigWigAverageOverBed tmp.${Cell}.${Method}.bw tmp.genes.300bins.bed tmp.${Cell}.${Method}.genes.tab
		bigWigAverageOverBed tmp.${Cell}.${Method}.bw tmp.tss.300bins.bed tmp.${Cell}.${Method}.tss.tab
		cut -f 6 tmp.${Cell}.${Method}.genes.tab > tmp.${Cell}.${Method}.genes.signal
		cut -f 6 tmp.${Cell}.${Method}.tss.tab > tmp.${Cell}.${Method}.tss.signal
	done
done



echo -e "TSS_BS\tTSS_Mhemi\tgene_BS\tgene_Mhemi" > GM12878_BS_Mhemi_TSS_gene_body_300bins.txt
paste tmp.GM12878.BS.tss.signal tmp.GM12878.Mhemi.tss.signal tmp.GM12878.BS.genes.signal tmp.GM12878.Mhemi.genes.signal >> GM12878_BS_Mhemi_TSS_gene_body_300bins.txt

echo -e "TSS_BS\tTSS_Mhemi\tgene_BS\tgene_Mhemi" > H1_BS_Mhemi_TSS_gene_body_300bins.txt
paste tmp.H1.BS.tss.signal tmp.H1.Mhemi.tss.signal tmp.H1.BS.genes.signal tmp.H1.Mhemi.genes.signal >> H1_BS_Mhemi_TSS_gene_body_300bins.txt

echo -e "TSS_BS\tTSS_Mhemi\tgene_BS\tgene_Mhemi" > H9_BS_Mhemi_TSS_gene_body_300bins.txt
paste tmp.H9.BS.tss.signal tmp.H9.Mhemi.tss.signal tmp.H9.BS.genes.signal tmp.H9.Mhemi.genes.signal >> H9_BS_Mhemi_TSS_gene_body_300bins.txt


