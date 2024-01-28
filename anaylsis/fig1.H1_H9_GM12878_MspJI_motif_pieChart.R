##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/13
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
# windowsFonts(consolas=windowsFont("consolas"))
##-----------------------------
library(ggplot2)
library(ggsci)
library(ggrepel)

for(celine in c('H1', 'H9', 'GM12878')){
	motifs <- read.table(paste0(celine, '_fragment_end_motif.Counts'), header=T)
	motifs$percents <- paste0(as.character(round(100*motifs$Count/sum(motifs$Count), 2)), '%')

	pic <- ggplot(motifs, aes(x='', y=Count, fill=motif))+
		labs(title=celine)+
		geom_bar(stat='identity', position='stack', width=1)+scale_fill_jco()+
		geom_text_repel(aes(x=1, label=percents), vjust='inward', size=4, stat='identity', position=position_stack(vjust=.5))+
		coord_polar(theta='y', start=0)+
		guides(fill=guide_legend(nrow=2))+
		theme_void()+
		theme(plot.margin=unit(c(-.5, 0, -.2, 0), 'in'),
			plot.title=element_text(size=12, hjust=.5, vjust=-5),
			legend.background=element_blank(),
			legend.key=element_blank(),
			legend.key.size=unit(.15, units='in'),
			legend.text=element_text(colour='black', size=11),
			legend.title=element_blank(),
			legend.position=c(.5, .04))
	ggsave(paste0('fig1.', celine, '_MspJI_motif_PieChart.pdf'), width=2.5, height=2.7)
}
