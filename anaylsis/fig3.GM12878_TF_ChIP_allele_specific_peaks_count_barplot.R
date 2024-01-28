##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/25
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)

peakAS <- read.table('GM12878_TF_ChIP_allele_specific_peaks_count.txt', header=T, as.is=T)
peakAS$Count <- log2(peakAS$Count)
peakAS$totalCount <- rep(peakAS$Count[which(peakAS$Binds=='Maternal')] + peakAS$Count[which(peakAS$Binds=='Paternal')], times=1, each=2)
peakAS <- peakAS[order(peakAS$totalCount, decreasing=T), ]
peakAS$TF <- factor(peakAS$TF, levels=unique(peakAS$TF))

pic <- ggplot(peakAS, aes(x=TF, y=Count, fill=Binds))+geom_bar(stat='identity', width=.6)+
	theme_classic()+labs(title='Number of allele-specific peaks', x='', y=expression(log[2] (Count)))+
	scale_fill_manual(values=c('#708090', '#E69F00'))+
	scale_y_continuous(limits=c(-.05, 22.05), breaks=seq(0, 22, 11), expand=c(0, 0))+
	guides(fill=guide_legend(nrow=1))+
	theme(plot.margin=unit(c(.05, 0, -.18, .02), 'in'),
		plot.title=element_text(size=11, hjust=.5),
		legend.background=element_blank(),
		legend.title=element_blank(),
		legend.text=element_text(size=9),
		legend.key.size=unit(.15, 'in'),
		legend.spacing.x=unit(.05, 'in'),
		legend.spacing.y=unit(.0, 'in'),
		legend.position=c(.8, .9),
		axis.line.x=element_blank(),
		axis.line.y=element_line(colour='black', linewidth=.3),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(linewidth=.3),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10, angle=90, hjust=1, vjust=.5),
		axis.text.y=element_text(colour='black', size=10))
ggsave('fig3.GM12878_TF_ChIP_allele_specific_peaks_count_barplot.pdf', pic, width=5, height=2.3)
