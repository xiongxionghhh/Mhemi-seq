#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2022/12/31
##-------------------------

[ $# != 2 ] && { echo "Usage: bash dyad_overlap_ctcf.sh <input> <output>"
		 echo "<input> file name of CpG dyad methylation in bedGraph format"
		 echo "<output> output file name"
		 exit 1
		}

awk '{if($2>10) print $1, $2-10, $3+10, $4, $5, $6}' OFS='\t' /mnt/disk4/public/RefBed/CTCF/Hs_CTCF.motif | \
	bedtools intersect -wa -a ${1} -b - > tmp.dyad.ctcf.bdg

echo -e "cyto\tun\tsame\toppo\tme" > ${2}
for((pos=-9; pos<=9; pos+=1)); do
	awk '{if($6=="+") print $1, $2+"'${pos}'", $3+1+"'${pos}'", $4, $5, $6; else print $1, $2-1-"'${pos}'", $3-"'${pos}'", $4, $5, $6}' OFS='\t' /mnt/disk4/public/RefBed/CTCF/Hs_CTCF.motif | \
		bedtools intersect -F 1 -f 1 -wa -wb -a tmp.dyad.ctcf.bdg -b - > tmp.bdg
		sort -u -k1,1 -k2,2n tmp.bdg | \
		awk '{if($13=="+") print $1, $2, $3, $4, $5, $6, $7; else print $1, $2, $3, $4, $6, $5, $7}' | \
		awk -v pos=${pos} '{unme+=$4; same+=$5; oppo+=$6; me+=$7} END{print "C"pos+10, unme, same, oppo, me}' OFS='\t' >> ${2}
done
