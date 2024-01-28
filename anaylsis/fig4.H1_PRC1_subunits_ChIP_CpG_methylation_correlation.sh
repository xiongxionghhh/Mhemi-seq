#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/30
##-------------------------

for i in BCOR H2Aub KDM2B PCGF1 RING1B RYBP; do
	bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H1/ChIP-seq/bw/H1_${i}.bw /mnt/disk1/5/share/Reference/bed/hg38_1kb_window.bed tmp.${i}.tab
done

for i in H3K4me3 H3K27me3; do
	bigWigAverageOverBed /mnt/disk4/public/H1/ChIP-seq_his/H1_${i}.bw /mnt/disk1/5/share/Reference/bed/hg38_1kb_window.bed tmp.${i}.tab
done

bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bw /mnt/disk1/5/share/Reference/bed/hg38_1kb_window.bed tmp.Mhemi.tab

for i in BCOR H2Aub KDM2B PCGF1 RING1B RYBP H3K4me3 H3K27me3; do
	cut -f 5 tmp.${i}.tab > tmp.${i}.signal
done
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

paste tmp.PCGF1.signal tmp.BCOR.signal tmp.KDM2B.signal tmp.RING1B.signal tmp.RYBP.signal tmp.H2Aub.signal tmp.H3K4me3.signal tmp.H3K27me3.signal tmp.Mhemi.signal | awk '$1+$2+$3+$4+$5+$6+$7+$8!=0' | grep -v 'e' | grep -v 'nan' > H1_PRC1_subunits_ChIP_CpG_methylation_correlation_signal.txt
sed -i '1i\PCGF1\tBCOR\tKDM2B\tRING1B\tRYBP\tH2Aub\tH3K4me3\tH3K27me3\tDNAme' H1_PRC1_subunits_ChIP_CpG_methylation_correlation_signal.txt










awk '{if($4=="Promoter") print $1, $2, $3, "window"NR}' OFS='\t' /mnt/disk1/5/share/Reference/bed/hg38_GenomicPositon.bed > tmp.bed

for i in BCOR H2Aub KDM2B PCGF1 RING1B RYBP; do
	bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/H1/ChIP-seq/bw/H1_${i}.bw tmp.bed tmp.${i}.tab
done

for i in H3K4me3 H3K27me3; do
	bigWigAverageOverBed /mnt/disk4/public/H1/ChIP-seq_his/H1_${i}.bw tmp.bed tmp.${i}.tab
done

for i in BCOR H2Aub KDM2B PCGF1 RING1B RYBP H3K4me3 H3K27me3; do
	cut -f 5 tmp.${i}.tab > tmp.${i}.signal
done

bigWigAverageOverBed /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG.bw tmp.bed tmp.Mhemi.tab
cut -f 6 tmp.Mhemi.tab > tmp.Mhemi.signal

bigWigAverageOverBed /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bw/GM12878_BS_CGNR.bw tmp.bed tmp.Mhemi.tab
cut -f 6 tmp.Mhemi.tab > tmp.BS.signal

paste tmp.PCGF1.signal tmp.BCOR.signal tmp.KDM2B.signal tmp.RING1B.signal tmp.RYBP.signal tmp.H2Aub.signal tmp.H3K4me3.signal tmp.H3K27me3.signal tmp.BS.signal tmp.Mhemi.signal | awk '$1+$2+$3+$4+$5+$6+$7+$8!=0' | grep -v 'e' | grep -v 'nan' > tmp.txt
awk 'BEGIN{print "PCGF1\tBCOR\tKDM2B\tRING1B\tRYBP\tH2Aub\tH3K4me3\tH3K27me3\tBS\tMhemi"} {print $0}' tmp.txt > H1_PRC1_subunits_ChIP_BS_Mhemi_promoter_correlation_signal.txt

paste tmp.PCGF1.signal tmp.BCOR.signal tmp.KDM2B.signal tmp.RING1B.signal tmp.RYBP.signal tmp.H2Aub.signal tmp.H3K4me3.signal tmp.H3K27me3.signal tmp.Mhemi.signal | awk '$1+$2+$3+$4+$5+$6+$7+$8!=0' | grep -v 'e' | grep -v 'nan' > tmp.txt
awk 'BEGIN{print "PCGF1\tBCOR\tKDM2B\tRING1B\tRYBP\tH2Aub\tH3K4me3\tH3K27me3\tDNAme"} {print $0}' tmp.txt > H1_PRC1_subunits_ChIP_Mhemi_promoter_correlation_signal.txt
