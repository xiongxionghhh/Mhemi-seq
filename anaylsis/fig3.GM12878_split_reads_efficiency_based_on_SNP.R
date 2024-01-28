##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/25
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)

effAS <- read.table('GM12878_split_reads_efficiency_based_on_SNP.txt', header=T)
effAS$Method <- factor(effAS$Method, levels=c('Mhemi-seq', 'BS-seq', 'ATAC-seq', 'ChIP-seq'))

pic <- ggplot(effAS, aes(x=Method, y=Efficiency))+
	geom_col(position=position_dodge(width=.7), width=.7, size=0, colour=NA)+
	geom_text(aes(y=Efficiency+.5, label=Efficiency), size=3.5, position=position_dodge(width=.6))+
	theme_classic()+labs(title='', x='', y='% of all alignments')+
	scale_y_continuous(limits=c(-.03, 12.03), breaks=seq(0, 12, 4), expand=c(0, 0))+
	guides(fill=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(.12, 0, -.18, .02), 'in'),
		plot.title=element_blank(),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_line(colour='black', linewidth=.3),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(linewidth=.3),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10, angle=30, hjust=1),
		axis.text.y=element_text(colour='black', size=10))
ggsave('fig3.GM12878_split_reads_efficiency_based_on_SNP_barplot.pdf', pic, width=2, height=2.5)
