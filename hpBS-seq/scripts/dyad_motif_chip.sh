#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/2/1
##-------------------------

[ $# != 4 ] && { echo "Usage: bash dyad_motif_chip.sh <motif> <dyad> <position> <name>"
		 echo "<motif> motif position in bedGraph format"
		 echo "<dyad> dyad methylation in bedGraph format"
		 echo "<position> distance (bp) from motif center, separated by comma"
		 echo "<name> output file name"
		 exit 1
		}

awk '{print $0, "m"NR}' OFS='\t' ${1} > tmp.name.motif
touch tmp.${4}.txt

for pos in `echo ${3} | sed 's/,/ /g' `; do
	awk '{if($6=="+") print $1, $2+"'${pos}'", $3+1+"'${pos}'", $4, $5, $6, $(NF); else print $1, $2-1-"'${pos}'", $3-"'${pos}'", $4, $5, $6, $(NF)}' OFS='\t' tmp.name.motif | sort -k1,1 -k2,2n | \
	bedtools intersect -sorted -wb -f 1 -F 1 -a - -b ${2} | \
	awk '{nCount=$11+$12+$13+$14} {if($6=="+") print $7, 100*$11/nCount, 100*$12/nCount, 100*$13/nCount, 100*$14/nCount; else print $7, 100*$11/nCount, 100*$13/nCount, 100*$12/nCount, 100*$14/nCount}' OFS='\t' | \
	awk '{if($2>$3 && $2>$4 && $2>$5) print $1, "unme", "'${pos}'", "'${4}'"; else if($3>$2 && $3>$4 && $3>$5) print $1, "same", "'${pos}'", "'${4}'"; else if($4>$2 && $4>$3 && $4>$5) print $1, "oppo", "'${pos}'", "'${4}'"; else if($5>$2 && $5>$3 && $5>$4) print $1, "me", "'${pos}'", "'${4}'"}' OFS='\t' >> tmp.${4}.txt
done

sort -k1,1 -k2,2 -k3,3 tmp.${4}.txt > tmp.${4}.sorted.txt
awk '{print $(NF), $5}' OFS='\t' tmp.name.motif > tmp.unsorted.motif
sort -k1,1 tmp.unsorted.motif | join -j 1 - tmp.${4}.sorted.txt | sed 's/ /\t/g' | sort -u > ${4}_signal.txt
