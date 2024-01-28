##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/30
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)

repDyad <- read.table('H1_promoter_Mhemi_CpG_dyad_count_classify_by_PCGF1.txt', header=T, as.is=T)
repDyad$Dyad <- factor(repDyad$Dyad, levels=rev(c('unme', 'same', 'oppo', 'me')))
repDyad$Pcgf <- factor(repDyad$Pcgf, levels=rev(c('strong', 'weak', 'unbound')))

pic <- ggplot(repDyad, aes(x=Pcgf, y=Freq, fill=Dyad))+geom_bar(stat='identity', width=.4)+
	theme_classic()+labs(title='', x='', y='frequency (%)')+
	scale_fill_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	scale_y_continuous(limits=c(-.18, 100.18), breaks=seq(0, 100, 25), expand=c(0, 0))+
	guides(fill=guide_legend(nrow=1))+
	theme(plot.margin=unit(c(.25, .05, -.15, .01), 'inches'),
		plot.title=element_blank(),
		panel.grid.major.y=element_line(color='#C0C0C0', linewidth=.5),
		legend.background=element_blank(),
		legend.title=element_blank(),
		legend.text=element_text(size=9),
		legend.key=element_blank(),
		legend.key.size=unit(.1, 'in'),
		legend.spacing.x=unit(.02, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		legend.position=c(.48, 1.05),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10),
		axis.text.y=element_text(colour='black', size=10, angle=90, hjust=.5))
ggsave('fig4.H1_promoter_Mhemi_CpG_dyad_freq_classify_by_PCGF1_barplot.pdf', pic, width=2.2, height=2.5)



repDyad <- read.table('H1_promoter_hpBS_CpG_dyad_count_classify_by_PCGF1.txt', header=T, as.is=T)
repDyad$Dyad <- factor(repDyad$Dyad, levels=rev(c('unme', 'same', 'oppo', 'me')))

pic <- ggplot(repDyad, aes(x=Pcgf, y=Freq, fill=Dyad))+geom_bar(stat='identity', width=.4)+
	theme_classic()+labs(title='', x='', y='frequency (%)')+
	scale_fill_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	scale_y_continuous(limits=c(-.18, 100.18), breaks=seq(0, 100, 25), expand=c(0, 0))+
	guides(fill=guide_legend(nrow=1))+
	theme(plot.margin=unit(c(.25, .05, -.2, .01), 'inches'),
		plot.title=element_blank(),
		panel.grid.major.y=element_line(color='#C0C0C0', linewidth=.5),
		legend.background=element_blank(),
		legend.title=element_blank(),
		legend.text=element_text(size=9),
		legend.key=element_blank(),
		legend.key.size=unit(.1, 'in'),
		legend.spacing.x=unit(.02, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		legend.position=c(.48, 1.05),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.y=element_text(colour='black', size=10),
		axis.text.x=element_text(colour='black', size=10),
		axis.text.y=element_text(colour='black', size=10, angle=90, hjust=.5))
ggsave('fig4.H1_promoter_hpBS_CpG_dyad_freq_classify_by_PCGF1_barplot.pdf', pic, width=2.2, height=2.5)
