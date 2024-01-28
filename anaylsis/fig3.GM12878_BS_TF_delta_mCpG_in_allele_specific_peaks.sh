#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/25
##-------------------------

for i in ARID3A ARNT ATF2 ATF7 BACH1 BCLAF1 BHLHE40 BMI1 CTCF DPF2 E2F8 E4F1 ELF1 ETV6 FOXK2 GATAD2B HDGF IKZF1 IKZF2 IRF3 IRF5 KLF5 LARP7 MAZ MEF2B MLLT1 MTA2 NBN NFATC3 NFXL1 NKRF NR2C1 NR2F1 PKNOX1 RAD51 RB1 RELB SKIL SMAD1 SMARCA5 SRF STAT1 TARDBP TBX21 TCF12 TRIM22 YBX1 ZBTB33 ZBTB40 ZNF143 ZNF207 ZNF217 ZNF24 ZNF592 ZNF687 ZSCAN29; do
	awk '{print $1, $2-100, $2+100, $4, $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/ChIP-seq/AlleleSpecific/biasPeaks/${i}'_'summits.bed > tmp.regions.bed
	awk '$5=="M_Imprinted"' tmp.regions.bed > tmp.bed
	intersectBed -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_R2_Maternal.bdg -b tmp.bed > tmp.unbound.bdg
	intersectBed -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_R1_Maternal.bdg -b tmp.bed >> tmp.unbound.bdg

	intersectBed -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_R2_Paternal.bdg -b tmp.bed > tmp.bound.bdg
	intersectBed -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_R1_Paternal.bdg -b tmp.bed >> tmp.bound.bdg

	awk '$5=="P_Imprinted"' tmp.regions.bed > tmp.bed
	intersectBed -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_R2_Maternal.bdg -b tmp.bed >> tmp.bound.bdg
	intersectBed -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_R1_Maternal.bdg -b tmp.bed >> tmp.bound.bdg

	intersectBed -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_R2_Paternal.bdg -b tmp.bed >> tmp.unbound.bdg
	intersectBed -u -sorted -a /mnt/disk1/5/xx/private/Mhmei/GM12878/BS-seq/bedGraph/GM12878_R1_Paternal.bdg -b tmp.bed >> tmp.unbound.bdg

	awk '{mSum+=$4*$5; nSum+=$5} END{print "'${i}'", "unbound", mSum/(nSum+0.001)}' OFS='\t' tmp.unbound.bdg >> tmp.TF.AS.txt
	awk '{mSum+=$4*$5; nSum+=$5} END{print "'${i}'", "bound", mSum/(nSum+0.001)}' OFS='\t' tmp.bound.bdg >> tmp.TF.AS.txt
done

awk '{if($2=="unbound") print $1, $3}' OFS='\t' tmp.TF.AS.txt | sort -k1,1 > tmp.unbound
awk '{if($2=="bound") print $1, $3}' OFS='\t' tmp.TF.AS.txt | sort -k1,1 | paste tmp.unbound - | awk '{if($1==$3) print $1, $2, $4}' OFS='\t' > GM12878_BS_delta_mCpG_in_allele_specific_peaks.txt
awk '{if($2=="bound") print $1, $3}' OFS='\t' tmp.TF.AS.txt | sort -k1,1 | paste tmp.unbound - | awk '{dev=$2-$4} {if($1==$3 && sqrt(dev*dev)>10) print $1, dev}' OFS='\t' | awk '{if($2>0) print $1}' > GM12878_BS_TF_allele_specific_delta_mCpG_over15_names.txt
