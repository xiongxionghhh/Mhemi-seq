#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/1/31
##-------------------------

[ $# != 2 ] && { echo "Usage: bash dyad_mC_motif.sh <dyad_mC> <motif>"
		 echo "<dyad_mC> file name of CpG dyad methylation in bedGraph format"
		 echo "<motif> position of motif"
		 exit 1
		}

awk '{if($2>180) print $1, $2-180, $3+180, $4, $5, $6}' OFS='\t' ${2} | bedtools intersect -sorted -u -a ${1} -b - > tmp.bdg
echo -e "Distance\tMethylation" > SAME.txt
echo -e "Distance\tMethylation" > OPPO.txt
for((pos=-170; pos<=170; pos+=1)); do
	awk '$2>180' ${2} | awk -v pos=$pos '{if($6=="+") print $1, $2-5+pos, $3+5+1+pos, $4, $5, $6; else print $1, $2-5-1-pos, $3+5-pos, $4, $5, $6}' OFS='\t' | bedtools intersect -wb -f 1 -a tmp.bdg -b - > tmp.bed

	if [ -s tmp.bed ]; then
		awk -v pos=$pos '{nsum+=$4+$5+$6+$7} {if($13=="+") nsame+=$6; else nsame+=$5} END{print pos, 100*nsame/nsum}' OFS='\t' tmp.bed >> SAME.txt
		awk -v pos=$pos '{nsum+=$4+$5+$6+$7} {if($13=="+") noppo+=$5; else noppo+=$6} END{print pos, 100*noppo/nsum}' OFS='\t' tmp.bed >> OPPO.txt
	else
		echo -e "$pos\tNA" >> SAME.txt
		echo -e "$pos\tNA" >> OPPO.txt
	fi

	echo "position $pos has done!"
done
