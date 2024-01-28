#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2021/12/30
##-------------------------

echo -e "TF\tBinds\tCount" > GM12878_TF_ChIP_allele_specific_peaks_count.txt
for i in ARID3A ARNT ATF2 ATF7 BACH1 BCLAF1 BHLHE40 BMI1 CTCF DPF2 E2F8 E4F1 ELF1 ETV6 FOXK2 GATAD2B HDGF IKZF1 IKZF2 IRF3 IRF5 KLF5 LARP7 MAZ MEF2B MLLT1 MTA2 NBN NFATC3 NFXL1 NKRF NR2C1 NR2F1 PKNOX1 RAD51 RB1 RELB SKIL SMAD1 SMARCA5 SRF STAT1 TARDBP TBX21 TCF12 TRIM22 YBX1 ZBTB33 ZBTB40 ZNF143 ZNF207 ZNF217 ZNF24 ZNF592 ZNF687 ZSCAN29; do
	cut -f 4,5 /mnt/disk1/5/xx/private/Mhmei/GM12878/ChIP-seq/AlleleSpecific/biasPeaks/${i}'_'summits.bed | sort -u | awk '{if($2=="P_Imprinted") nMa+=1; else nPa+=1} END{if(nMa>=10 && nPa>=10) print "'${i}'\tMaternal\t"nMa"\n'${i}'\tPaternal\t"nPa}' >> GM12878_TF_ChIP_allele_specific_peaks_count.txt
done
