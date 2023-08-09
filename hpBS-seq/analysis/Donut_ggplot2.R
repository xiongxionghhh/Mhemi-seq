##----------------------------
# @Author: Xiong Xiong
# @Date: 2021/1/13 10:45
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2) 

oriCount <- read.table('tmp.Count.txt', header=T, sep='\t')
oriCount$CountLab <- gsub(' ', '', format(oriCount$Count, big.mark=','))

ggplot(oriCount, aes(x=.5, y=Count, fill=Dyad)) +
	geom_col()+
	geom_text(aes(label=CountLab), position=position_stack(vjust=.5), size=4)+
	coord_polar(theta='y')+
	scale_x_continuous(limits=c(-1.25, 1.25))+
	scale_fill_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	theme_classic()+
	theme(plot.margin=unit(c(-.3, -1, -.35, -1), 'in'),
		plot.title=element_blank(),
		legend.background=element_blank(),
		legend.key=element_blank(),
		legend.key.size=unit(.15, 'in'),
		legend.title=element_blank(),
		legend.text=element_text(colour='black', size=11),
		legend.position=c(.5, .5),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
ggsave('CELL_METHOD.pdf', width=2, height=2)
