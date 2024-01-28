##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/28
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2) 

oriCount <- read.table('H1_PCGF1_peak_summit_hg38_annotation.txt', header=T, sep='\t')
oriCount$CountLab <- gsub(' ', '', format(oriCount$Count, big.mark=','))
oriCount$Pos <- factor(oriCount$Pos, levels=c('intron', 'intergenic', 'exon', 'promoter'))

pic <- ggplot(oriCount, aes(x=.5, y=Count, fill=Pos)) +
	geom_col(width=.72)+
	geom_text(aes(label=CountLab), position=position_stack(vjust=.5), size=4)+
	coord_polar(theta='y')+
	labs(title='PCGF1 peak summit')+
	scale_x_continuous(limits=c(-1.2, 1.2))+
	scale_fill_manual(values=c('#2A557F', '#45BC9C', '#F05076', '#FFCD6E'))+
	theme_classic()+
	theme(plot.margin=unit(c(-.4, -.5, -.4, -.5), 'in'),
		plot.title=element_text(colour='black', size=11, hjust=.5, vjust=-12),
		legend.background=element_blank(),
		legend.key=element_blank(),
		legend.key.size=unit(.15, 'in'),
		legend.title=element_blank(),
		legend.text=element_text(colour='black', size=10),
		legend.position=c(.5, .5),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
ggsave('fig4.H1_PCGF1_peak_summit_position_donut.pdf', pic, width=2, height=2)
