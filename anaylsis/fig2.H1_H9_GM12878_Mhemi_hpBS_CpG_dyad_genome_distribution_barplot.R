##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/22
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)
library(dplyr)

load('/mnt/disk1/5/share/others/colorset.RData')
dyad <- read.table('GM12878_Mhemi_CpG_dyad_hg38_annotation.txt', sep='\t')

dyadFreq <- dyad %>% group_by(V2, V1) %>% 
	summarise(n=n()) %>% 
	mutate(freq=100*n/sum(n))

repDyad <- as.data.frame(dyadFreq)
repDyad$V2 <- factor(repDyad$V2, levels=c('random', 'unme', 'hemi', 'me'))
repDyad$V1 <- factor(repDyad$V1, levels=c('5-UTR', '3-UTR', 'exon', 'intron', 'Intergenic', 'non-coding', 'promoter', 'TTS'))

pic <- ggplot(repDyad, aes(x=V2, y=freq, fill=V1))+geom_bar(stat='identity', width=.6)+
	theme_classic()+labs(title='GM12878', x='', y='fraction (%)')+
	scale_fill_manual(values=colorSet[1:8])+ #c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533')
	scale_y_continuous(limits=c(-.18, 100.18), breaks=seq(0, 100, 25), expand=c(0, 0))+
	guides(fill=guide_legend(nrow=2))+
	theme(plot.margin=unit(c(-.02, .1, .18, .1), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		legend.background=element_blank(),
		legend.title=element_blank(),
		legend.text=element_text(size=8.5),
		legend.key.size=unit(.12, 'in'),
		legend.spacing.x=unit(.01, 'in'),
		legend.spacing.y=unit(0, 'in'),
		legend.position=c(.4, -.22),
		axis.line.x=element_blank(),
		axis.line.y=element_line(colour='black', linewidth=.25),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(linewidth=.25),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10),
		axis.text.y=element_text(colour='black', size=10, angle=90, hjust=.5))
ggsave('fig2.GM12878_Mhemi_CpG_dyad_genome_distribution.pdf', pic, width=2.5, height=2.5)

dyad <- read.table('H1_Mhemi_CpG_dyad_hg38_annotation.txt', sep='\t')

dyadFreq <- dyad %>% group_by(V2, V1) %>% 
	summarise(n=n()) %>% 
	mutate(freq=100*n/sum(n))

repDyad <- as.data.frame(dyadFreq)
repDyad$V2 <- factor(repDyad$V2, levels=c('random', 'unme', 'hemi', 'me'))
repDyad$V1 <- factor(repDyad$V1, levels=c('5-UTR', '3-UTR', 'exon', 'intron', 'Intergenic', 'non-coding', 'promoter', 'TTS'))

pic <- ggplot(repDyad, aes(x=V2, y=freq, fill=V1))+geom_bar(stat='identity', width=.6)+
	theme_classic()+labs(title='H1', x='', y='fraction (%)')+
	scale_fill_manual(values=colorSet[1:8])+ #c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533')
	scale_y_continuous(limits=c(-.18, 100.18), breaks=seq(0, 100, 25), expand=c(0, 0))+
	guides(fill=guide_legend(nrow=2))+
	theme(plot.margin=unit(c(-.02, .1, .18, .1), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		legend.background=element_blank(),
		legend.title=element_blank(),
		legend.text=element_text(size=8.5),
		legend.key.size=unit(.12, 'in'),
		legend.spacing.x=unit(.01, 'in'),
		legend.spacing.y=unit(0, 'in'),
		legend.position=c(.4, -.22),
		axis.line.x=element_blank(),
		axis.line.y=element_line(colour='black', linewidth=.25),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(linewidth=.25),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10),
		axis.text.y=element_text(colour='black', size=10, angle=90, hjust=.5))
ggsave('fig2.H1_CpG_dyad_genome_distribution.pdf', pic, width=2.5, height=2.5)

dyad <- read.table('H9_Mhemi_CpG_dyad_hg38_annotation.txt', sep='\t')

dyadFreq <- dyad %>% group_by(V2, V1) %>% 
	summarise(n=n()) %>% 
	mutate(freq=100*n/sum(n))

repDyad <- as.data.frame(dyadFreq)
repDyad$V2 <- factor(repDyad$V2, levels=c('random', 'unme', 'hemi', 'me'))
repDyad$V1 <- factor(repDyad$V1, levels=c('5-UTR', '3-UTR', 'exon', 'intron', 'Intergenic', 'non-coding', 'promoter', 'TTS'))

pic <- ggplot(repDyad, aes(x=V2, y=freq, fill=V1))+geom_bar(stat='identity', width=.6)+
	theme_classic()+labs(title='H9', x='', y='fraction (%)')+
	scale_fill_manual(values=colorSet[1:8])+ #c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533')
	scale_y_continuous(limits=c(-.18, 100.18), breaks=seq(0, 100, 25), expand=c(0, 0))+
	guides(fill=guide_legend(nrow=2))+
	theme(plot.margin=unit(c(-.02, .1, .18, .1), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		legend.background=element_blank(),
		legend.title=element_blank(),
		legend.text=element_text(size=8.5),
		legend.key.size=unit(.12, 'in'),
		legend.spacing.x=unit(.01, 'in'),
		legend.spacing.y=unit(0, 'in'),
		legend.position=c(.4, -.22),
		axis.line.x=element_blank(),
		axis.line.y=element_line(colour='black', linewidth=.25),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(linewidth=.25),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10),
		axis.text.y=element_text(colour='black', size=10, angle=90, hjust=.5))
ggsave('fig2.H9_CpG_dyad_genome_distribution.pdf', pic, width=2.5, height=2.5)




dyad <- read.table('GM12878_hpBS_CpG_dyad_hg38_annotation.txt', sep='\t')

dyadFreq <- dyad %>% group_by(V2, V1) %>% 
	summarise(n=n()) %>% 
	mutate(freq=100*n/sum(n))

repDyad <- as.data.frame(dyadFreq)
repDyad$V2 <- factor(repDyad$V2, levels=c('random', 'unme', 'hemi', 'me'))
repDyad$V1 <- factor(repDyad$V1, levels=c('5-UTR', '3-UTR', 'exon', 'intron', 'Intergenic', 'non-coding', 'promoter', 'TTS'))

pic <- ggplot(repDyad, aes(x=V2, y=freq, fill=V1))+geom_bar(stat='identity', width=.6)+
	theme_classic()+labs(title='hpBS (GM12878)', x='', y='fraction (%)')+
	scale_fill_manual(values=colorSet[1:8])+ #c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533')
	scale_y_continuous(limits=c(-.18, 100.18), breaks=seq(0, 100, 25), expand=c(0, 0))+
	guides(fill=guide_legend(nrow=2))+
	theme(plot.margin=unit(c(-.02, .1, .18, .1), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		legend.background=element_blank(),
		legend.title=element_blank(),
		legend.text=element_text(size=8.5),
		legend.key.size=unit(.12, 'in'),
		legend.spacing.x=unit(.01, 'in'),
		legend.spacing.y=unit(0, 'in'),
		legend.position=c(.4, -.22),
		axis.line.x=element_blank(),
		axis.line.y=element_line(colour='black', linewidth=.25),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(linewidth=.25),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10),
		axis.text.y=element_text(colour='black', size=10, angle=90, hjust=.5))
ggsave('fig2.GM12878_hpBS_CpG_dyad_genome_distribution.pdf', pic, width=2.5, height=2.5)

dyad <- read.table('GM12878_iSA_CpG_dyad_hg38_annotation.txt', sep='\t')

dyadFreq <- dyad %>% group_by(V2, V1) %>% 
	summarise(n=n()) %>% 
	mutate(freq=100*n/sum(n))

repDyad <- as.data.frame(dyadFreq)
repDyad$V2 <- factor(repDyad$V2, levels=c('random', 'unme', 'hemi', 'me'))
repDyad$V1 <- factor(repDyad$V1, levels=c('5-UTR', '3-UTR', 'exon', 'intron', 'Intergenic', 'non-coding', 'promoter', 'TTS'))

pic <- ggplot(repDyad, aes(x=V2, y=freq, fill=V1))+geom_bar(stat='identity', width=.6)+
	theme_classic()+labs(title='iSA (GM12878)', x='', y='fraction (%)')+
	scale_fill_manual(values=colorSet[1:8])+ #c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533')
	scale_y_continuous(limits=c(-.18, 100.18), breaks=seq(0, 100, 25), expand=c(0, 0))+
	guides(fill=guide_legend(nrow=2))+
	theme(plot.margin=unit(c(-.02, .1, .18, .1), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		legend.background=element_blank(),
		legend.title=element_blank(),
		legend.text=element_text(size=8.5),
		legend.key.size=unit(.12, 'in'),
		legend.spacing.x=unit(.01, 'in'),
		legend.spacing.y=unit(0, 'in'),
		legend.position=c(.4, -.22),
		axis.line.x=element_blank(),
		axis.line.y=element_line(colour='black', linewidth=.25),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(linewidth=.25),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10),
		axis.text.y=element_text(colour='black', size=10, angle=90, hjust=.5))
ggsave('fig2.GM12878_iSA_CpG_dyad_genome_distribution.pdf', pic, width=2.5, height=2.5)




dyad <- read.table('H1_hpBS_CpG_dyad_hg38_annotation.txt', sep='\t')

dyadFreq <- dyad %>% group_by(V2, V1) %>% 
	summarise(n=n()) %>% 
	mutate(freq=100*n/sum(n))

repDyad <- as.data.frame(dyadFreq)
repDyad$V2 <- factor(repDyad$V2, levels=c('random', 'unme', 'hemi', 'me'))
repDyad$V1 <- factor(repDyad$V1, levels=c('5-UTR', '3-UTR', 'exon', 'intron', 'Intergenic', 'non-coding', 'promoter', 'TTS'))

pic <- ggplot(repDyad, aes(x=V2, y=freq, fill=V1))+geom_bar(stat='identity', width=.6)+
	theme_classic()+labs(title='hpBS (H1)', x='', y='fraction (%)')+
	scale_fill_manual(values=colorSet[1:8])+ #c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533')
	scale_y_continuous(limits=c(-.18, 100.18), breaks=seq(0, 100, 25), expand=c(0, 0))+
	guides(fill=guide_legend(nrow=2))+
	theme(plot.margin=unit(c(-.02, .1, .18, .1), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		legend.background=element_blank(),
		legend.title=element_blank(),
		legend.text=element_text(size=8.5),
		legend.key.size=unit(.12, 'in'),
		legend.spacing.x=unit(.01, 'in'),
		legend.spacing.y=unit(0, 'in'),
		legend.position=c(.4, -.22),
		axis.line.x=element_blank(),
		axis.line.y=element_line(colour='black', linewidth=.25),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(linewidth=.25),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10),
		axis.text.y=element_text(colour='black', size=10, angle=90, hjust=.5))
ggsave('fig2.H1_hpBS_CpG_dyad_genome_distribution.pdf', pic, width=2.5, height=2.5)

dyad <- read.table('H9_hpBS_CpG_dyad_hg38_annotation.txt', sep='\t')

dyadFreq <- dyad %>% group_by(V2, V1) %>% 
	summarise(n=n()) %>% 
	mutate(freq=100*n/sum(n))

repDyad <- as.data.frame(dyadFreq)
repDyad$V2 <- factor(repDyad$V2, levels=c('random', 'unme', 'hemi', 'me'))
repDyad$V1 <- factor(repDyad$V1, levels=c('5-UTR', '3-UTR', 'exon', 'intron', 'Intergenic', 'non-coding', 'promoter', 'TTS'))

pic <- ggplot(repDyad, aes(x=V2, y=freq, fill=V1))+geom_bar(stat='identity', width=.6)+
	theme_classic()+labs(title='hpBS (H9)', x='', y='fraction (%)')+
	scale_fill_manual(values=colorSet[1:8])+ #c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533')
	scale_y_continuous(limits=c(-.18, 100.18), breaks=seq(0, 100, 25), expand=c(0, 0))+
	guides(fill=guide_legend(nrow=2))+
	theme(plot.margin=unit(c(-.02, .1, .18, .1), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		legend.background=element_blank(),
		legend.title=element_blank(),
		legend.text=element_text(size=8.5),
		legend.key.size=unit(.12, 'in'),
		legend.spacing.x=unit(.01, 'in'),
		legend.spacing.y=unit(0, 'in'),
		legend.position=c(.4, -.22),
		axis.line.x=element_blank(),
		axis.line.y=element_line(colour='black', linewidth=.25),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(linewidth=.25),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10),
		axis.text.y=element_text(colour='black', size=10, angle=90, hjust=.5))
ggsave('fig2.H9_hpBS_CpG_dyad_genome_distribution.pdf', pic, width=2.5, height=2.5)

