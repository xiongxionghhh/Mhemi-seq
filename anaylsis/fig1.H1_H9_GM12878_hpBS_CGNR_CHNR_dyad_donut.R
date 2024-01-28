##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/10/28
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2) 

oriCount <- read.table('H1_H9_GM12878_hpBS_CGNR_dyad_count.txt', header=T, sep='\t')
oriCount$Dyad <- factor(oriCount$Dyad, levels=c('unme', 'hemi-Watson', 'hemi-Crick', 'me'))
oriCount$CountLab <- gsub(' ', '', format(oriCount$Count, big.mark=','))

dyadCount <- subset(oriCount, Cell=='GM12878', select=-1)

pic <- ggplot(dyadCount, aes(x=.5, y=Count, fill=Dyad)) +
	geom_col(width=.72)+labs(title=('CGNR (GM12878)'))+
	geom_text(aes(label=CountLab), position=position_stack(vjust=.5), size=3)+
	coord_polar(theta='y')+
	scale_x_continuous(limits=c(-1.25, 1.25))+
	scale_fill_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	theme_classic()+
	theme(plot.margin=unit(c(-.42, -1, -.35, -1), 'in'),
		plot.title=element_text(colour='black', size=10, hjust=.5, vjust=-13.5),
		legend.background=element_blank(),
		legend.key=element_blank(),
		legend.key.size=unit(.11, 'in'),
		legend.spacing.x=unit(.02, 'in'),
		legend.title=element_blank(),
		legend.text=element_text(colour='black', size=8),
		legend.position=c(.5, .5),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
ggsave('fig1.GM12878_hpBS_CGNR_dyad_count_donut.pdf', width=2, height=2)


dyadCount <- subset(oriCount, Cell=='H1', select=-1)

pic <- ggplot(dyadCount, aes(x=.5, y=Count, fill=Dyad)) +
	geom_col(width=.72)+labs(title=('CGNR (H1)'))+
	geom_text(aes(label=CountLab), position=position_stack(vjust=.5), size=3)+
	coord_polar(theta='y')+
	scale_x_continuous(limits=c(-1.25, 1.25))+
	scale_fill_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	theme_classic()+
	theme(plot.margin=unit(c(-.42, -1, -.35, -1), 'in'),
		plot.title=element_text(colour='black', size=10, hjust=.5, vjust=-13.5),
		legend.background=element_blank(),
		legend.key=element_blank(),
		legend.key.size=unit(.11, 'in'),
		legend.spacing.x=unit(.02, 'in'),
		legend.title=element_blank(),
		legend.text=element_text(colour='black', size=8),
		legend.position=c(.5, .5),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
ggsave('fig1.H1_hpBS_CGNR_dyad_count_donut.pdf', width=2, height=2)

dyadCount <- subset(oriCount, Cell=='H9', select=-1)

pic <- ggplot(dyadCount, aes(x=.5, y=Count, fill=Dyad)) +
	geom_col(width=.72)+labs(title=('CGNR (H9)'))+
	geom_text(aes(label=CountLab), position=position_stack(vjust=.5), size=3)+
	coord_polar(theta='y')+
	scale_x_continuous(limits=c(-1.25, 1.25))+
	scale_fill_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	theme_classic()+
	theme(plot.margin=unit(c(-.42, -1, -.35, -1), 'in'),
		plot.title=element_text(colour='black', size=10, hjust=.5, vjust=-13.5),
		legend.background=element_blank(),
		legend.key=element_blank(),
		legend.key.size=unit(.11, 'in'),
		legend.spacing.x=unit(.02, 'in'),
		legend.title=element_blank(),
		legend.text=element_text(colour='black', size=8),
		legend.position=c(.5, .5),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
ggsave('fig1.H9_hpBS_CGNR_dyad_count_donut.pdf', width=2, height=2)





oriCount <- read.table('H1_H9_GM12878_hpBS_CWGR_dyad_count.txt', header=T, sep='\t')
oriCount$Dyad <- factor(oriCount$Dyad, levels=c('unme', 'hemi-Watson', 'hemi-Crick', 'me'))
oriCount$CountLab <- gsub(' ', '', format(oriCount$Count, big.mark=','))

dyadCount <- subset(oriCount, Cell=='GM12878', select=-1)

pic <- ggplot(dyadCount, aes(x=.5, y=Count, fill=Dyad)) +
	geom_col(width=.72)+labs(title=('CWGR (GM12878)'))+
	geom_text(aes(label=CountLab), position=position_stack(vjust=.5), size=3)+
	coord_polar(theta='y')+
	scale_x_continuous(limits=c(-1.25, 1.25))+
	scale_fill_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	theme_classic()+
	theme(plot.margin=unit(c(-.42, -1, -.35, -1), 'in'),
		plot.title=element_text(colour='black', size=10, hjust=.5, vjust=-13.5),
		legend.background=element_blank(),
		legend.key=element_blank(),
		legend.key.size=unit(.11, 'in'),
		legend.spacing.x=unit(.02, 'in'),
		legend.title=element_blank(),
		legend.text=element_text(colour='black', size=8),
		legend.position=c(.5, .5),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
ggsave('fig1.GM12878_hpBS_CWGR_dyad_count_donut.pdf', width=2, height=2)


dyadCount <- subset(oriCount, Cell=='H1', select=-1)

pic <- ggplot(dyadCount, aes(x=.5, y=Count, fill=Dyad)) +
	geom_col(width=.72)+labs(title=('CWGR (H1)'))+
	geom_text(aes(label=CountLab), position=position_stack(vjust=.5), size=3)+
	coord_polar(theta='y')+
	scale_x_continuous(limits=c(-1.25, 1.25))+
	scale_fill_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	theme_classic()+
	theme(plot.margin=unit(c(-.42, -1, -.35, -1), 'in'),
		plot.title=element_text(colour='black', size=10, hjust=.5, vjust=-13.5),
		legend.background=element_blank(),
		legend.key=element_blank(),
		legend.key.size=unit(.11, 'in'),
		legend.spacing.x=unit(.02, 'in'),
		legend.title=element_blank(),
		legend.text=element_text(colour='black', size=8),
		legend.position=c(.5, .5),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
ggsave('fig1.H1_hpBS_CWGR_dyad_count_donut.pdf', width=2, height=2)

dyadCount <- subset(oriCount, Cell=='H9', select=-1)

pic <- ggplot(dyadCount, aes(x=.5, y=Count, fill=Dyad)) +
	geom_col(width=.72)+labs(title=('CWGR (H9)'))+
	geom_text(aes(label=CountLab), position=position_stack(vjust=.5), size=3)+
	coord_polar(theta='y')+
	scale_x_continuous(limits=c(-1.25, 1.25))+
	scale_fill_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	theme_classic()+
	theme(plot.margin=unit(c(-.42, -1, -.35, -1), 'in'),
		plot.title=element_text(colour='black', size=10, hjust=.5, vjust=-13.5),
		legend.background=element_blank(),
		legend.key=element_blank(),
		legend.key.size=unit(.11, 'in'),
		legend.spacing.x=unit(.02, 'in'),
		legend.title=element_blank(),
		legend.text=element_text(colour='black', size=8),
		legend.position=c(.5, .5),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
ggsave('fig1.H9_hpBS_CWGR_dyad_count_donut.pdf', width=2, height=2)

















