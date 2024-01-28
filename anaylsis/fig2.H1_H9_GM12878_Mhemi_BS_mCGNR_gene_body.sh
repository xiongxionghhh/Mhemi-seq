#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/13
##-------------------------

awk '{if($2>1000 && $6=="+") print $1, $2-1000, $2+500, $4, $5, $6; else if($2>1000 && $6=="-") print $1, $3-500, $2+1000, $4, $5, $6}' OFS='\t' /mnt/disk4/public/RefBed/refGene_hg38/hg38_TSS.bed > tmp.promoter.bed

ln -s -f /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG.bed tmp.GM.Mhemi.bed
ln -s -f /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_BS_CpG.bdg tmp.GM.BS.bed

awk '{if($5>4) print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' tmp.GM.Mhemi.bed | bedtools intersect -sorted -s -wb -a - -b tmp.GM.BS.bed > tmp.OL.bed
cut -f 4,10 tmp.OL.bed > GM12878_Mhemi_BS_total_CpG_methylation.txt
sort -k1,1 -k2,2n tmp.promoter.bed | bedtools intersect -sorted -wa -u -a tmp.OL.bed -b - | cut -f 4,10 > GM12878_Mhemi_BS_promoter_CpG_methylation.txt
bedtools intersect -sorted -wa -u -a tmp.OL.bed -b /mnt/disk4/public/RefBed/hg38_CGI.bed | cut -f 4,10 > GM12878_Mhemi_BS_CGIs_CpG_methylation.txt

ln -f -s /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bed tmp.GM.Mhemi.bed
ln -f -s /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bedGraph/H1_BS_CpG.bdg tmp.GM.BS.bed

awk '{if($5>4) print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' tmp.GM.Mhemi.bed | bedtools intersect -sorted -s -wb -a - -b tmp.GM.BS.bed > tmp.OL.bed
cut -f 4,10 tmp.OL.bed > H1_Mhemi_BS_total_CpG_methylation.txt
sort -k1,1 -k2,2n tmp.promoter.bed | bedtools intersect -sorted -wa -u -a tmp.OL.bed -b - | cut -f 4,10 > H1_Mhemi_BS_promoter_CpG_methylation.txt
bedtools intersect -sorted -wa -u -a tmp.OL.bed -b /mnt/disk4/public/RefBed/hg38_CGI.bed | cut -f 4,10 > H1_Mhemi_BS_CGIs_CpG_methylation.txt


ln -f -s /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG.bed tmp.GM.Mhemi.bed
ln -f -s /mnt/disk1/5/xx/private/Mhmei/H9/BS-seq/bedGraph/H9_BS_CpG.bdg tmp.GM.BS.bed

awk '{if($5>4) print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' tmp.GM.Mhemi.bed | bedtools intersect -sorted -s -wb -a - -b tmp.GM.BS.bed > tmp.OL.bed
cut -f 4,10 tmp.OL.bed > H9_Mhemi_BS_total_CpG_methylation.txt
sort -k1,1 -k2,2n tmp.promoter.bed | bedtools intersect -sorted -wa -u -a tmp.OL.bed -b - | cut -f 4,10 > H9_Mhemi_BS_promoter_CpG_methylation.txt
bedtools intersect -sorted -wa -u -a tmp.OL.bed -b /mnt/disk4/public/RefBed/hg38_CGI.bed | cut -f 4,10 > H9_Mhemi_BS_CGIs_CpG_methylation.txt


#-------------------------------------------
ln -s -f /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_refGene_TPM.bed tmp.GM.tpm.bed
awk '$5>20' tmp.GM.tpm.bed > tmp.high.bed
awk '$5==0' tmp.GM.tpm.bed > tmp.low.bed

ln -s -f /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG.bw tmp.bw

computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.high.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "High", $0}' OFS='\t' tmp.Coverage.matrix > GM12878_Mhemi_CpG_gene_body_1kb_up_down.txt
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.low.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "Low", $0}' OFS='\t' tmp.Coverage.matrix >> GM12878_Mhemi_CpG_gene_body_1kb_up_down.txt


ln -s -f /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bw/GM12878_BS_CpG.bw tmp.bw

computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.high.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "High", $0}' OFS='\t' tmp.Coverage.matrix > GM12878_BS_CpG_gene_body_1kb_up_down.txt
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.low.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "Low", $0}' OFS='\t' tmp.Coverage.matrix >> GM12878_BS_CpG_gene_body_1kb_up_down.txt




ln -s -f /mnt/disk1/5/xx/private/Mhmei/H1/RNA-seq/TPM/H1_RNA_R1_RefGene.bed tmp.H1.tpm.bed
awk '$5>20' tmp.H1.tpm.bed > tmp.high.bed
awk '$5==0' tmp.H1.tpm.bed > tmp.low.bed

ln -s -f /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bw tmp.bw

computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.high.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "High", $0}' OFS='\t' tmp.Coverage.matrix > H1_Mhemi_CpG_gene_body_1kb_up_down.txt
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.low.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "Low", $0}' OFS='\t' tmp.Coverage.matrix >> H1_Mhemi_CpG_gene_body_1kb_up_down.txt


ln -s -f /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CpG.bw tmp.bw

computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.high.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "High", $0}' OFS='\t' tmp.Coverage.matrix > H1_BS_CpG_gene_body_1kb_up_down.txt
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.low.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "Low", $0}' OFS='\t' tmp.Coverage.matrix >> H1_BS_CpG_gene_body_1kb_up_down.txt


ln -s -f /mnt/disk1/5/xx/private/Mhmei/H9/RNA-seq/TPM/H9_R1r1_RefGene.bed tmp.H9.tpm.bed
awk '$5>20' tmp.H9.tpm.bed > tmp.high.bed
awk '$5==0' tmp.H9.tpm.bed > tmp.low.bed

ln -s -f /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG.bw tmp.bw

computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.high.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "High", $0}' OFS='\t' tmp.Coverage.matrix > H9_Mhemi_CpG_gene_body_1kb_up_down.txt
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.low.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "Low", $0}' OFS='\t' tmp.Coverage.matrix >> H9_Mhemi_CpG_gene_body_1kb_up_down.txt


ln -s -f /mnt/disk1/5/xx/private/Mhmei/H9/BS-seq/bw/H9_BS_CpG.bw tmp.bw

computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.high.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "High", $0}' OFS='\t' tmp.Coverage.matrix > H9_BS_CpG_gene_body_1kb_up_down.txt
computeMatrix scale-regions -bs 25 --regionBodyLength 3000 --upstream 1000 --downstream 1000 -p 16 -S tmp.bw -R tmp.low.bed -o tmp.plotFile.gz --outFileNameMatrix tmp.Coverage.matrix
awk '{if(NR>3) print "Low", $0}' OFS='\t' tmp.Coverage.matrix >> H9_BS_CpG_gene_body_1kb_up_down.txt
