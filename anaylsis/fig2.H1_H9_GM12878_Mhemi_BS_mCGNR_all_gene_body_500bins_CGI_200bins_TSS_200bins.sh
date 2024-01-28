#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/22
##-------------------------

bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bw/GM12878_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_genes_1kb_500bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_genes_1kb_500bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > GM12878_BS_Mhemi_mCGNR_genes_1kb_500bins.txt


bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_genes_1kb_500bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_genes_1kb_500bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > H1_BS_Mhemi_mCGNR_genes_1kb_500bins.txt


bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H9/BS-seq/bw/H9_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_genes_1kb_500bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_genes_1kb_500bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > H9_BS_Mhemi_mCGNR_genes_1kb_500bins.txt





bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bw/GM12878_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_CGI_200bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_CGI_200bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > GM12878_BS_Mhemi_mCGNR_CGI_200bins.txt


bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_CGI_200bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_CGI_200bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > H1_BS_Mhemi_mCGNR_CGI_200bins.txt


bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H9/BS-seq/bw/H9_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_CGI_200bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_CGI_200bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > H9_BS_Mhemi_mCGNR_CGI_200bins.txt






bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bw/GM12878_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_TSS_+-2kb_200bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_TSS_+-2kb_200bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > GM12878_BS_Mhemi_mCGNR_TSS_200bins.txt


bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_TSS_+-2kb_200bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_TSS_+-2kb_200bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > H1_BS_Mhemi_mCGNR_TSS_200bins.txt


bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H9/BS-seq/bw/H9_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_TSS_+-2kb_200bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_TSS_+-2kb_200bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > H9_BS_Mhemi_mCGNR_TSS_200bins.txt










bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bw/GM12878_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_gene_body_300bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_gene_body_300bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > GM12878_BS_Mhemi_mCGNR_gene_body_300bins.txt


bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H1/BS-seq/bw/H1_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_gene_body_300bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_gene_body_300bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > H1_BS_Mhemi_mCGNR_gene_body_300bins.txt


bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H9/BS-seq/bw/H9_BS_CGNR.bw /mnt/disk1/5/share/Reference/bed/hg38_gene_body_300bins.bed tmp.BS.tab
bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_gene_body_300bins.bed tmp.Mhemi.tab

cut -f 6 tmp.BS.tab > tmp.BS.signal
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.BS.signal tmp.Mhemi.signal | grep -v 'e' | grep -v 'nan' | awk 'BEGIN{print "BS\tMhemi"} {print $0}' > H9_BS_Mhemi_mCGNR_gene_body_300bins.txt



