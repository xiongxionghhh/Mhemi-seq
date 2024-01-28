#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/25
##-------------------------

























echo -e "Gene\tAllele\tRep\tTPM" > GM12878_imprinting_genes_expression.txt
awk '{if($4=="CRELD2") print "CRELD2\tPaternal\tR1r1\t"$5"\nCRELD2\tMaternal\tR1r1\t"$7}' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r1.ImpGenes.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="CRELD2") print "CRELD2", "Paternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="CRELD2") print "CRELD2", "Maternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="CRELD2") print "CRELD2", "Paternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="CRELD2") print "CRELD2", "Maternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="CRELD2") print "CRELD2", "Paternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="CRELD2") print "CRELD2", "Maternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($4=="ACCS") print "ACCS\tPaternal\tR1r1\t"$5"\nACCS\tMaternal\tR1r1\t"$7}' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r1.ImpGenes.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="ACCS") print "ACCS", "Paternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="ACCS") print "ACCS", "Maternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="ACCS") print "ACCS", "Paternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="ACCS") print "ACCS", "Maternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="ACCS") print "ACCS", "Paternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="ACCS") print "ACCS", "Maternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($4=="SNRPN") print "SNRPN\tPaternal\tR1r1\t"$5"\nSNRPN\tMaternal\tR1r1\t"$7}' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r1.ImpGenes.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNRPN") print "SNRPN", "Paternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNRPN") print "SNRPN", "Maternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNRPN") print "SNRPN", "Paternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNRPN") print "SNRPN", "Maternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNRPN") print "SNRPN", "Paternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNRPN") print "SNRPN", "Maternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($4=="SNURF") print "SNURF\tPaternal\tR1r1\t"$5"\nSNURF\tMaternal\tR1r1\t"$7}' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r1.ImpGenes.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNURF") print "SNURF", "Paternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNURF") print "SNURF", "Maternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNURF") print "SNURF", "Paternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNURF") print "SNURF", "Maternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNURF") print "SNURF", "Paternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="SNURF") print "SNURF", "Maternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($4=="TMSB4X") print "TMSB4X\tPaternal\tR1r1\t"$5"\nTMSB4X\tMaternal\tR1r1\t"$7}' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r1.ImpGenes.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TMSB4X") print "TMSB4X", "Paternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TMSB4X") print "TMSB4X", "Maternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TMSB4X") print "TMSB4X", "Paternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TMSB4X") print "TMSB4X", "Maternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TMSB4X") print "TMSB4X", "Paternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TMSB4X") print "TMSB4X", "Maternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($4=="TUBB6") print "TUBB6\tPaternal\tR1r1\t"$5"\nTUBB6\tMaternal\tR1r1\t"$7}' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r1.ImpGenes.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TUBB6") print "TUBB6", "Paternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TUBB6") print "TUBB6", "Maternal", "R1r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TUBB6") print "TUBB6", "Paternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TUBB6") print "TUBB6", "Maternal", "R2r1", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r1_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TUBB6") print "TUBB6", "Paternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Paternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt
awk '{if($7=="TUBB6") print "TUBB6", "Maternal", "R2r2", $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R2r2_Maternal.TPM.bdg >> GM12878_imprinting_genes_expression.txt




Rscript Fig15.Boxplot_ARNTexpr.R


cat /mnt/disk1/5/xx/private/Mhmei/GM12878/ChIP-seq/AlleleSpecific/biasPeaks/*.bed > tmp.peaks.bed

awk '{if($6=="+") print $1, $2-1000, $3+500, $4, $5, $6, $7, $8; else print $1, $2-500, $3+1000, $4, $5, $6, $7, $8}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/RNA-seq/bedGraph/GM12878_RNA_R1r1.ImpGenes.bdg | bedtools intersect -wb -wa -a - -b tmp.peaks.bed > GM12878_imprinting_genes_overlapped_with_allele_specific_peaks.txt


cut -f 4,12 GM12878_imprinting_genes_overlapped_with_allele_specific_peaks.txt | sort -u



bash /mnt/disk1/5/share/Reference/Annotation/genome_tracks/snapshot_pygenometracks.sh hg38 GM12878_CTCF_Maternal.bw GM12878_CTCF_Paternal.bw GM12878_ARNT_Maternal.bw GM12878_ARNT_Paternal.bw GM12878_RNA_Maternal.bw GM12878_RNA_Paternal.bw GM12878_BS_Maternal.bw GM12878_BS_Paternal.bw
pyGenomeTracks --tracks tmp.plot.ini --region chr22:49,904,479-49,941,663 --outFileName CRELD2_Tracks.pdf



bash /mnt/disk1/5/share/Reference/Annotation/genome_tracks/snapshot_pygenometracks.sh hg38 GM12878_MTA2_Maternal.bw GM12878_MTA2_Paternal.bw GM12878_ARNT_Maternal.bw GM12878_ARNT_Paternal.bw GM12878_RNA_Maternal.bw GM12878_RNA_Paternal.bw GM12878_BS_Maternal.bw GM12878_BS_Paternal.bw
pyGenomeTracks --tracks tmp.plot.ini --region chr1:235,838,162-236,202,496 --outFileName NID1_Tracks.pdf


bash /mnt/disk1/5/share/Reference/Annotation/genome_tracks/snapshot_pygenometracks.sh hg38 GM12878_ELF1_Maternal.bw GM12878_ELF1_Paternal.bw GM12878_RNA_Maternal.bw GM12878_RNA_Paternal.bw GM12878_BS_Maternal.bw GM12878_BS_Paternal.bw
pyGenomeTracks --tracks tmp.plot.ini --region chr11:44,050,456-44,099,878 --outFileName GM12878_ACCS_ELF1_Tracks.pdf
mv tmp.plot.ini GM12878_ACCS_ELF1_plot.ini



bash /mnt/disk1/5/share/Reference/Annotation/genome_tracks/snapshot_pygenometracks.sh hg38 GM12878_IKZF1_Maternal.bw GM12878_IKZF1_Paternal.bw GM12878_RNA_Maternal.bw GM12878_RNA_Paternal.bw GM12878_BS_Maternal.bw GM12878_BS_Paternal.bw
pyGenomeTracks --tracks tmp.plot.ini --region chr15:24,602,430-25,199,464 --outFileName GM12878_SNRPN_IKZF1_Tracks.pdf
mv tmp.plot.ini GM12878_SNRPN_IKZF1_plot.ini




bash /mnt/disk1/5/share/Reference/Annotation/genome_tracks/snapshot_pygenometracks.sh hg38 GM12878_ZBTB40_Maternal.bw GM12878_ZBTB40_Paternal.bw GM12878_RNA_Maternal.bw GM12878_RNA_Paternal.bw GM12878_BS_Maternal.bw GM12878_BS_Paternal.bw
pyGenomeTracks --tracks tmp.plot.ini --region chr16:3,409,612-3,475,714 --outFileName GM12878_ZNF597_ZBTB40_Tracks.pdf
mv tmp.plot.ini GM12878_ZNF597_ZBTB40_plot.ini




bash /mnt/disk1/5/share/Reference/Annotation/genome_tracks/snapshot_pygenometracks.sh hg38 ./bw/GM12878_ATF7_Maternal.bw ./bw/GM12878_ATF7_Paternal.bw ./bw/GM12878_ELF1_Maternal.bw ./bw/GM12878_ELF1_Paternal.bw ./bw/GM12878_FOXK2_Maternal.bw ./bw/GM12878_FOXK2_Paternal.bw ./bw/GM12878_NKRF_Maternal.bw ./bw/GM12878_NKRF_Paternal.bw ./bw/GM12878_RELB_Maternal.bw ./bw/GM12878_RELB_Paternal.bw ./bw/GM12878_RNA_Maternal.bw ./bw/GM12878_RNA_Paternal.bw ./bw/GM12878_BS_Maternal.bw ./bw/GM12878_BS_Paternal.bw
pyGenomeTracks --tracks tmp.plot.ini --region chr2:113,162,978-113,340,196 --outFileName GM12878_PAX8-AS1_CTCF_Tracks.pdf
mv tmp.plot.ini GM12878_PAX8-AS1_CTCF_plot.ini



bash /mnt/disk1/5/share/Reference/Annotation/genome_tracks/snapshot_pygenometracks.sh hg38 ./bw/GM12878_BHLHE40_Maternal.bw ./bw/GM12878_BHLHE40_Paternal.bw ./bw/GM12878_IKZF1_Maternal.bw ./bw/GM12878_IKZF1_Paternal.bw ./bw/GM12878_IKZF2_Maternal.bw ./bw/GM12878_IKZF2_Paternal.bw ./bw/GM12878_RNA_Maternal.bw ./bw/GM12878_RNA_Paternal.bw ./bw/GM12878_BS_Maternal.bw ./bw/GM12878_BS_Paternal.bw
pyGenomeTracks --tracks tmp.plot.ini --region chrX:9,134,460-10,047,890 --outFileName GM12878_TBL1X_Tracks.pdf
mv tmp.plot.ini GM12878_ZNF597_ZBTB40_plot.ini



bash /mnt/disk1/5/share/Reference/Annotation/genome_tracks/snapshot_pygenometracks.sh hg38 ./bw/GM12878_BMI1_Maternal.bw ./bw/GM12878_BMI1_Paternal.bw ./bw/GM12878_MLLT1_Maternal.bw ./bw/GM12878_MLLT1_Paternal.bw ./bw/GM12878_SMAD1_Maternal.bw ./bw/GM12878_SMAD1_Paternal.bw ./bw/GM12878_RNA_Maternal.bw ./bw/GM12878_RNA_Paternal.bw ./bw/GM12878_BS_Maternal.bw ./bw/GM12878_BS_Paternal.bw
pyGenomeTracks --tracks tmp.plot.ini --region chrX:12,968,937-12,983,208 --outFileName GM12878_TMSB4X_Tracks.pdf
mv tmp.plot.ini GM12878_TMSB4X_plot.ini

