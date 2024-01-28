#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/25
##-------------------------

for i in `cat GM12878_BS_TF_allele_specific_delta_mCpG_over15_names.txt`; do
	awk '{print $1, $2-100, $2+100, $4, $5}' OFS='\t' /mnt/disk1/5/xx/private/Mhmei/GM12878/ChIP-seq/AlleleSpecific/biasPeaks/${i}'_'summits.bed > tmp.regions.bed
	awk '$5=="M_Imprinted"' tmp.regions.bed | sort -k1,1 -k2,2n -u > tmp.bed
	bedtools intersect -u -sorted -a /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_Maternal_CpG_dyad.bed -b tmp.bed > tmp.unbound.bdg
	bedtools intersect -u -sorted -a /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_Paternal_CpG_dyad.bed -b tmp.bed > tmp.bound.bdg

	awk '$5=="P_Imprinted"' tmp.regions.bed | sort -k1,1 -k2,2n -u > tmp.bed
	bedtools intersect -u -sorted -a /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_Maternal_CpG_dyad.bed -b tmp.bed >> tmp.bound.bdg
	bedtools intersect -u -sorted -a /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_Paternal_CpG_dyad.bed -b tmp.bed >> tmp.unbound.bdg
done

awk 'BEGIN{print "allele", "unme", "hemi", "me"} {hemi+=$5+$6; unme+=$4; me+=$7} END{print "unbound", unme, hemi, me}' OFS='\t' tmp.unbound.bdg > GM12878_Mhemi_TF_allele_specific_CpG_dyad_count.txt
awk '{hemi+=$5+$6; unme+=$4; me+=$7} END{print "bound", unme, hemi, me}' tmp.bound.bdg >> GM12878_Mhemi_TF_allele_specific_CpG_dyad_count.txt
