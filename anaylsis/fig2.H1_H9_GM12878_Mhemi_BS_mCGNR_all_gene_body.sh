#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/22
##-------------------------

ln -s -f /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_refGene_TPM.bed tmp.GM.tpm.bed
ln -s -f /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG.bw tmp.bw

computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.GM.tpm.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "Mhemi-seq", $0}' OFS='\t' tmp.Coverage.matrix > GM12878_Mhemi_BS_CGNR_all_gene_body_1kb_up_down.txt

ln -s -f /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bw/GM12878_BS_CGNR.bw tmp.bw
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.GM.tpm.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "BS-seq", $0}' OFS='\t' tmp.Coverage.matrix >> GM12878_Mhemi_BS_CGNR_all_gene_body_1kb_up_down.txt


ln -s -f /mnt/disk1/5/xx/private/Mhmei/H1/RNA-seq/TPM/H1_RNA_R1_RefGene.bed tmp.H1.tpm.bed
ln -s -f /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bw tmp.bw

computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.H1.tpm.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "Mhemi-seq", $0}' OFS='\t' tmp.Coverage.matrix > H1_Mhemi_BS_CGNR_all_gene_body_1kb_up_down.txt

ln -s -f /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CGNR.bw tmp.bw
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.H1.tpm.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "BS-seq", $0}' OFS='\t' tmp.Coverage.matrix >> H1_Mhemi_BS_CGNR_all_gene_body_1kb_up_down.txt


ln -s -f /mnt/disk1/5/xx/private/Mhmei/H9/RNA-seq/TPM/H9_R1r1_RefGene.bed tmp.H9.tpm.bed
ln -s -f /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG.bw tmp.bw

computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.H9.tpm.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "Mhemi-seq", $0}' OFS='\t' tmp.Coverage.matrix > H9_Mhemi_BS_CGNR_all_gene_body_1kb_up_down.txt

ln -s -f /mnt/disk1/5/xx/private/Mhmei/H9/BS-seq/bw/H9_BS_CGNR.bw tmp.bw
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.H9.tpm.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "BS-seq", $0}' OFS='\t' tmp.Coverage.matrix >> H9_Mhemi_BS_CGNR_all_gene_body_1kb_up_down.txt
