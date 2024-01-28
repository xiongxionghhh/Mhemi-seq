#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/22
##-------------------------

ln -s -f /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bed/GM12878_Mhemi_rep123_CpG_dyad.bed tmp.GM.dyad.bed
awk '{if(($4+$5+$6+$7)>9) print $1, $2, $3, $4, $5+$6, $7}' OFS='\t' tmp.GM.dyad.bed | awk '{if($4>$5 && $4>$6) print $1, $2, $3, "unme"; else if($5>$4 && $5>$6) print $1, $2, $3, "hemi"; else if($6>$4 && $6>$5) print $1, $2, $3, "me"}' OFS='\t' > tmp.C10.dyad
awk '{print $0, "+"}' OFS='\t' tmp.C10.dyad > tmp.unsort.bed
awk '{print $0, "-"}' OFS='\t' tmp.C10.dyad >> tmp.unsort.bed
sort -k1,1 -k2,2n --parallel 10 tmp.unsort.bed | awk '{print $1, $2, $3, "dyad"NR, $4, $5}' OFS='\t' > tmp.GM.C10.dyad.bed

awk '$5=="unme"' tmp.GM.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "unme"}' OFS='\t' > tmp.anno.txt

awk '$5=="hemi"' tmp.GM.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "hemi"}' OFS='\t' >> tmp.anno.txt

awk '$5=="me"' tmp.GM.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "me"}' OFS='\t' >> tmp.anno.txt

samtools view -H /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep1/bam/GM12878_Mhemi_rep1_sorted_dedup.bam | awk '{if($1=="@SQ") print $2, $3}' OFS='\t' | sed 's/.N://g' > tmp.sizes
bedtools shuffle -i tmp.GM.C10.dyad.bed -g tmp.sizes > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "random"}' OFS='\t' >> tmp.anno.txt

sed 's/ \..*//g; s/-TSS//g; s/ //g' tmp.anno.txt | sed s#\'#-#g | awk '$1!="NA" && $2!=""' > GM12878_Mhemi_CpG_dyad_hg38_annotation.txt


ln -s -f /mnt/disk4/public/GM12878/BS-seq/iSA/GM12878_intraCpG.bed tmp.GM.dyad.bed
awk '{if(($4+$5+$6+$7)>9) print $1, $2, $3, $4, $5+$6, $7}' OFS='\t' tmp.GM.dyad.bed | awk '{if($4>$5 && $4>$6) print $1, $2, $3, "unme"; else if($5>$4 && $5>$6) print $1, $2, $3, "hemi"; else if($6>$4 && $6>$5) print $1, $2, $3, "me"}' OFS='\t' > tmp.C10.dyad
awk '{print $0, "+"}' OFS='\t' tmp.C10.dyad > tmp.unsort.bed
awk '{print $0, "-"}' OFS='\t' tmp.C10.dyad >> tmp.unsort.bed
sort -k1,1 -k2,2n --parallel 10 tmp.unsort.bed | awk '{print $1, $2, $3, "dyad"NR, $4, $5}' OFS='\t' > tmp.GM.C10.dyad.bed

awk '$5=="unme"' tmp.GM.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "unme"}' OFS='\t' > tmp.anno.txt

awk '$5=="hemi"' tmp.GM.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "hemi"}' OFS='\t' >> tmp.anno.txt

awk '$5=="me"' tmp.GM.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "me"}' OFS='\t' >> tmp.anno.txt

samtools view -H /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep1/bam/GM12878_Mhemi_rep1_sorted_dedup.bam | awk '{if($1=="@SQ") print $2, $3}' OFS='\t' | sed 's/.N://g' > tmp.sizes
bedtools shuffle -i tmp.GM.C10.dyad.bed -g tmp.sizes > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "random"}' OFS='\t' >> tmp.anno.txt

sed 's/ \..*//g; s/-TSS//g; s/ //g' tmp.anno.txt | sed s#\'#-#g | awk '$1!="NA" && $2!=""' > GM12878_iSA_CpG_dyad_hg38_annotation.txt


ln -s -f /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bed/H1_Mhemi_rep123_CpG_dyad.bed tmp.H1.dyad.bed
awk '{if(($4+$5+$6+$7)>9) print $1, $2, $3, $4, $5+$6, $7}' OFS='\t' tmp.H1.dyad.bed | awk '{if($4>$5 && $4>$6) print $1, $2, $3, "unme"; else if($5>$4 && $5>$6) print $1, $2, $3, "hemi"; else if($6>$4 && $6>$5) print $1, $2, $3, "me"}' OFS='\t' > tmp.C10.dyad
awk '{print $0, "+"}' OFS='\t' tmp.C10.dyad > tmp.unsort.bed
awk '{print $0, "-"}' OFS='\t' tmp.C10.dyad >> tmp.unsort.bed
sort -k1,1 -k2,2n --parallel 10 tmp.unsort.bed | awk '{print $1, $2, $3, "dyad"NR, $4, $5}' OFS='\t' > tmp.H1.C10.dyad.bed

awk '$5=="unme"' tmp.H1.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "unme"}' OFS='\t' > tmp.anno.txt

awk '$5=="hemi"' tmp.H1.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "hemi"}' OFS='\t' >> tmp.anno.txt

awk '$5=="me"' tmp.H1.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "me"}' OFS='\t' >> tmp.anno.txt

samtools view -H /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep1/bam/H1_Mhemi_rep1_sorted_dedup.bam | awk '{if($1=="@SQ") print $2, $3}' OFS='\t' | sed 's/.N://g' > tmp.sizes
bedtools shuffle -i tmp.H1.C10.dyad.bed -g tmp.sizes > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "random"}' OFS='\t' >> tmp.anno.txt

sed 's/ \..*//g; s/-TSS//g; s/ //g' tmp.anno.txt | sed s#\'#-#g | awk '$1!="NA" && $2!=""' > H1_Mhemi_CpG_dyad_hg38_annotation.txt


ln -s -f /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bed/H9_Mhemi_rep123_CpG_dyad.bed tmp.H9.dyad.bed
awk '{if(($4+$5+$6+$7)>9) print $1, $2, $3, $4, $5+$6, $7}' OFS='\t' tmp.H9.dyad.bed | awk '{if($4>$5 && $4>$6) print $1, $2, $3, "unme"; else if($5>$4 && $5>$6) print $1, $2, $3, "hemi"; else if($6>$4 && $6>$5) print $1, $2, $3, "me"}' OFS='\t' > tmp.C10.dyad
awk '{print $0, "+"}' OFS='\t' tmp.C10.dyad > tmp.unsort.bed
awk '{print $0, "-"}' OFS='\t' tmp.C10.dyad >> tmp.unsort.bed
sort -k1,1 -k2,2n --parallel 10 tmp.unsort.bed | awk '{print $1, $2, $3, "dyad"NR, $4, $5}' OFS='\t' > tmp.H9.C10.dyad.bed

awk '$5=="unme"' tmp.H9.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "unme"}' OFS='\t' > tmp.anno.txt

awk '$5=="hemi"' tmp.H9.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "hemi"}' OFS='\t' >> tmp.anno.txt

awk '$5=="me"' tmp.H9.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "me"}' OFS='\t' >> tmp.anno.txt

samtools view -H /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep1/bam/H9_Mhemi_rep1_sorted_dedup.bam | awk '{if($1=="@SQ") print $2, $3}' OFS='\t' | sed 's/.N://g' > tmp.sizes
bedtools shuffle -i tmp.H9.C10.dyad.bed -g tmp.sizes > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "random"}' OFS='\t' >> tmp.anno.txt

sed 's/ \..*//g; s/-TSS//g; s/ //g' tmp.anno.txt | sed s#\'#-#g | awk '$1!="NA" && $2!=""' > H9_Mhemi_CpG_dyad_hg38_annotation.txt







ln -s -f /mnt/disk1/5/share/hpBS-seq/GM12878/bedGraph/GM12878_hpBS_rep123_YNCGNR_dyad.bed tmp.GM.dyad.bed
awk '{if(($4+$5+$6+$7)>9) print $1, $2, $3, $4, $5+$6, $7}' OFS='\t' tmp.GM.dyad.bed | awk '{if($4>$5 && $4>$6) print $1, $2, $3, "unme"; else if($5>$4 && $5>$6) print $1, $2, $3, "hemi"; else if($6>$4 && $6>$5) print $1, $2, $3, "me"}' OFS='\t' > tmp.C10.dyad
awk '{print $0, "+"}' OFS='\t' tmp.C10.dyad > tmp.unsort.bed
awk '{print $0, "-"}' OFS='\t' tmp.C10.dyad >> tmp.unsort.bed
sort -k1,1 -k2,2n --parallel 10 tmp.unsort.bed | awk '{print $1, $2, $3, "dyad"NR, $4, $5}' OFS='\t' > tmp.GM.C10.dyad.bed

awk '$5=="unme"' tmp.GM.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "unme"}' OFS='\t' > tmp.anno.txt

awk '$5=="hemi"' tmp.GM.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "hemi"}' OFS='\t' >> tmp.anno.txt

awk '$5=="me"' tmp.GM.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "me"}' OFS='\t' >> tmp.anno.txt

samtools view -H /mnt/disk1/5/share/hpBS-seq/GM12878/bam/20220402_GM12878_hpBS_R1.R1.Paired.bam | awk '{if($1=="@SQ") print $2, $3}' OFS='\t' | sed 's/.N://g' > tmp.sizes
bedtools shuffle -i tmp.GM.C10.dyad.bed -g tmp.sizes > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "random"}' OFS='\t' >> tmp.anno.txt

sed 's/ \..*//g; s/-TSS//g; s/ //g' tmp.anno.txt | sed s#\'#-#g | awk '$1!="NA" && $2!=""' > GM12878_hpBS_CpG_dyad_hg38_annotation.txt



ln -s -f /mnt/disk1/5/share/hpBS-seq/H1/bedGraph/20220704_H1_hpBS_R1.intraCpG.bdg tmp.H1.dyad.bed
awk '{if(($4+$5+$6+$7)>9) print $1, $2, $3, $4, $5+$6, $7}' OFS='\t' tmp.H1.dyad.bed | awk '{if($4>$5 && $4>$6) print $1, $2, $3, "unme"; else if($5>$4 && $5>$6) print $1, $2, $3, "hemi"; else if($6>$4 && $6>$5) print $1, $2, $3, "me"}' OFS='\t' > tmp.C10.dyad
awk '{print $0, "+"}' OFS='\t' tmp.C10.dyad > tmp.unsort.bed
awk '{print $0, "-"}' OFS='\t' tmp.C10.dyad >> tmp.unsort.bed
sort -k1,1 -k2,2n --parallel 10 tmp.unsort.bed | awk '{print $1, $2, $3, "dyad"NR, $4, $5}' OFS='\t' > tmp.H1.C10.dyad.bed

awk '$5=="unme"' tmp.H1.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "unme"}' OFS='\t' > tmp.anno.txt

awk '$5=="hemi"' tmp.H1.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "hemi"}' OFS='\t' >> tmp.anno.txt

awk '$5=="me"' tmp.H1.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "me"}' OFS='\t' >> tmp.anno.txt

samtools view -H /mnt/disk1/5/share/hpBS-seq/H1/bam/20220704_H1_hpBS_R1.R1.Paired.bam | awk '{if($1=="@SQ") print $2, $3}' OFS='\t' | sed 's/.N://g' > tmp.sizes
bedtools shuffle -i tmp.H1.C10.dyad.bed -g tmp.sizes > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "random"}' OFS='\t' >> tmp.anno.txt

sed 's/ \..*//g; s/-TSS//g; s/ //g' tmp.anno.txt | sed s#\'#-#g | awk '$1!="NA" && $2!=""' > H1_hpBS_CpG_dyad_hg38_annotation.txt


ln -s -f /mnt/disk1/5/share/hpBS-seq/H9/bedGraph/20220704_H9_hpBS_R1.intraCpG.bdg tmp.H9.dyad.bed
awk '{if(($4+$5+$6+$7)>9) print $1, $2, $3, $4, $5+$6, $7}' OFS='\t' tmp.H9.dyad.bed | awk '{if($4>$5 && $4>$6) print $1, $2, $3, "unme"; else if($5>$4 && $5>$6) print $1, $2, $3, "hemi"; else if($6>$4 && $6>$5) print $1, $2, $3, "me"}' OFS='\t' > tmp.C10.dyad
awk '{print $0, "+"}' OFS='\t' tmp.C10.dyad > tmp.unsort.bed
awk '{print $0, "-"}' OFS='\t' tmp.C10.dyad >> tmp.unsort.bed
sort -k1,1 -k2,2n --parallel 10 tmp.unsort.bed | awk '{print $1, $2, $3, "dyad"NR, $4, $5}' OFS='\t' > tmp.H9.C10.dyad.bed

awk '$5=="unme"' tmp.H9.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "unme"}' OFS='\t' > tmp.anno.txt

awk '$5=="hemi"' tmp.H9.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "hemi"}' OFS='\t' >> tmp.anno.txt

awk '$5=="me"' tmp.H9.C10.dyad.bed > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "me"}' OFS='\t' >> tmp.anno.txt

samtools view -H /mnt/disk1/5/share/hpBS-seq/H9/bam/20220704_H9_hpBS_R1.R1.Paired.bam | awk '{if($1=="@SQ") print $2, $3}' OFS='\t' | sed 's/.N://g' > tmp.sizes
bedtools shuffle -i tmp.H9.C10.dyad.bed -g tmp.sizes > tmp.bed
perl /mnt/disk1/5/xx/Software/HOMER-v4.11/bin/annotatePeaks.pl tmp.bed hg38 > tmp.annotation.txt
cut -f 8 tmp.annotation.txt | sed 's/(.*)//g' | awk '{if(NR>1) print $0, "random"}' OFS='\t' >> tmp.anno.txt

sed 's/ \..*//g; s/-TSS//g; s/ //g' tmp.anno.txt | sed s#\'#-#g | awk '$1!="NA" && $2!=""' > H9_hpBS_CpG_dyad_hg38_annotation.txt
