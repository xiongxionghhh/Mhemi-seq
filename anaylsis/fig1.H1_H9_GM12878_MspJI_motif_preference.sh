#!/bin/bash
##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/13
##-------------------------

samtools view -@ 30 /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bam/H1_Mhemi_rep123.bam | awk '{if(length($10)==32 && $10~/[A,T,C,G]{32}/) print $10}' | head -10000 > H1_32bp_sequence.txt
samtools view -@ 30 /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bam/H9_Mhemi_rep123.bam | awk '{if(length($10)==32 && $10~/[A,T,C,G]{32}/) print $10}' | head -10000 > H9_32bp_sequence.txt
samtools view -@ 30 /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bam/GM12878_Mhemi_rep123.bam | awk '{if(length($10)==32 && $10~/[A,T,C,G]{32}/) print $10}' | head -10000 > GM12878_32bp_sequence.txt

bamToBed -i /mnt/disk1/5/share/Mhemi-seq/H1/H1_Mhemi_rep123/bam/H1_Mhemi_rep123.bam | awk '$1!~/[_,M,C,Y,L]/ && $2>13' > tmp.seqs
awk '{if($6=="+") print $1, $2-13, $2+17; else print $1, $3-17, $3+13}' OFS='\t' tmp.seqs | fastaFromBed -tab -fi /ssd/index/bismark/hg38XX/hg38XX.fa -bed - | awk '{headPos=toupper(substr($2, 1, 1)); tailPos=toupper(substr($2, 30, 1))} {if(headPos=="C" && tailPos!="G") print "Head"; else if(headPos!="C" && tailPos=="G") print "Tail"; else print "N"}' > tmp.ends
paste tmp.seqs tmp.ends | awk '$7!="N"' | awk '{if($6=="+" && $7=="Head") print $1, $2-13, $2-9, NR, 0, "+"; else if($6=="+" && $7=="Tail") print $1, $2+13, $2+17, NR, 0, "-"; else if($6=="-" && $7=="Head") print $1, $3-17, $3-13, NR, 0, "+"; else if($6=="-" && $7=="Tail") print $1, $3+9, $3+13, NR, 0, "-"}' OFS='\t' | fastaFromBed -s -tab -fi /ssd/index/bismark/hg38XX/hg38XX.fa -bed - | awk '{print toupper($2)}' > tmp.motif.seqs
awk 'BEGIN{print "motif\tCount"} {if($1~/CG[A,T,C,G]A/) nCGNA+=1; else if($1~/CG[A,T,C,G]G/) nCGNG+=1; else if($1~/C[A,T,C][A,T,C,G]A/) nCHNA+=1; else if($1~/C[A,T,C][A,T,C,G]G/) nCHNG+=1; else if($1~/C[A,T,C,G][A,T,C,G]C/) nCNNC+=1; else if($1~/C[A,T,C,G][A,T,C,G]T/) nCNNT+=1} END{print "CGNA\t"nCGNA"\nCGNG\t"nCGNG"\nCHNA\t"nCHNA"\nCHNG\t"nCHNG"\nCNNC\t"nCNNC"\nCNNT\t"nCNNT}' tmp.motif.seqs > H1_fragment_end_motif.Counts

bamToBed -i /mnt/disk1/5/share/Mhemi-seq/H9/H9_Mhemi_rep123/bam/H9_Mhemi_rep123.bam | awk '$1!~/[_,M,C,Y,L]/ && $2>13' > tmp.seqs
awk '{if($6=="+") print $1, $2-13, $2+17; else print $1, $3-17, $3+13}' OFS='\t' tmp.seqs | fastaFromBed -tab -fi /ssd/index/bismark/hg38XX/hg38XX.fa -bed - | awk '{headPos=toupper(substr($2, 1, 1)); tailPos=toupper(substr($2, 30, 1))} {if(headPos=="C" && tailPos!="G") print "Head"; else if(headPos!="C" && tailPos=="G") print "Tail"; else print "N"}' > tmp.ends
paste tmp.seqs tmp.ends | awk '$7!="N"' | awk '{if($6=="+" && $7=="Head") print $1, $2-13, $2-9, NR, 0, "+"; else if($6=="+" && $7=="Tail") print $1, $2+13, $2+17, NR, 0, "-"; else if($6=="-" && $7=="Head") print $1, $3-17, $3-13, NR, 0, "+"; else if($6=="-" && $7=="Tail") print $1, $3+9, $3+13, NR, 0, "-"}' OFS='\t' | fastaFromBed -s -tab -fi /ssd/index/bismark/hg38XX/hg38XX.fa -bed - | awk '{print toupper($2)}' > tmp.motif.seqs
awk 'BEGIN{print "motif\tCount"} {if($1~/CG[A,T,C,G]A/) nCGNA+=1; else if($1~/CG[A,T,C,G]G/) nCGNG+=1; else if($1~/C[A,T,C][A,T,C,G]A/) nCHNA+=1; else if($1~/C[A,T,C][A,T,C,G]G/) nCHNG+=1; else if($1~/C[A,T,C,G][A,T,C,G]C/) nCNNC+=1; else if($1~/C[A,T,C,G][A,T,C,G]T/) nCNNT+=1} END{print "CGNA\t"nCGNA"\nCGNG\t"nCGNG"\nCHNA\t"nCHNA"\nCHNG\t"nCHNG"\nCNNC\t"nCNNC"\nCNNT\t"nCNNT}' tmp.motif.seqs > H9_fragment_end_motif.Counts

bamToBed -i /mnt/disk1/5/share/Mhemi-seq/GM12878/GM12878_Mhemi_rep123/bam/GM12878_Mhemi_rep123.bam | awk '$1!~/[_,M,C,Y,L]/ && $2>13' > tmp.seqs
awk '{if($6=="+") print $1, $2-13, $2+17; else print $1, $3-17, $3+13}' OFS='\t' tmp.seqs | fastaFromBed -tab -fi /ssd/index/bismark/hg38XX/hg38XX.fa -bed - | awk '{headPos=toupper(substr($2, 1, 1)); tailPos=toupper(substr($2, 30, 1))} {if(headPos=="C" && tailPos!="G") print "Head"; else if(headPos!="C" && tailPos=="G") print "Tail"; else print "N"}' > tmp.ends
paste tmp.seqs tmp.ends | awk '$7!="N"' | awk '{if($6=="+" && $7=="Head") print $1, $2-13, $2-9, NR, 0, "+"; else if($6=="+" && $7=="Tail") print $1, $2+13, $2+17, NR, 0, "-"; else if($6=="-" && $7=="Head") print $1, $3-17, $3-13, NR, 0, "+"; else if($6=="-" && $7=="Tail") print $1, $3+9, $3+13, NR, 0, "-"}' OFS='\t' | fastaFromBed -s -tab -fi /ssd/index/bismark/hg38XX/hg38XX.fa -bed - | awk '{print toupper($2)}' > tmp.motif.seqs
awk 'BEGIN{print "motif\tCount"} {if($1~/CG[A,T,C,G]A/) nCGNA+=1; else if($1~/CG[A,T,C,G]G/) nCGNG+=1; else if($1~/C[A,T,C][A,T,C,G]A/) nCHNA+=1; else if($1~/C[A,T,C][A,T,C,G]G/) nCHNG+=1; else if($1~/C[A,T,C,G][A,T,C,G]C/) nCNNC+=1; else if($1~/C[A,T,C,G][A,T,C,G]T/) nCNNT+=1} END{print "CGNA\t"nCGNA"\nCGNG\t"nCGNG"\nCHNA\t"nCHNA"\nCHNG\t"nCHNG"\nCNNC\t"nCNNC"\nCNNT\t"nCNNT}' tmp.motif.seqs > GM12878_fragment_end_motif.Counts

