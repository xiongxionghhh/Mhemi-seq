##----------------------------
# @Author: Xiong Xiong
# @Date: 2022/2/15
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)

mC <- read.table('HUES64_BS_mCHG_Alu_Start.txt', sep='\t', header=T)
mC$Methylation[which(mC$Strand=='oppo')] <- -1*mC$Methylation[which(mC$Strand=='oppo')]

p <- ggplot(mC, aes(x=Distance, y=Methylation, color=Strand))+
	geom_bar(fill=NA, stat='identity', width=.1)+
	geom_hline(yintercept=0, linewidth=.3)+scale_color_manual(values=c('#E69F00', '#708090'))+
	annotate('text', label='oppo', -50, -9.5, size=3.5, color='#E69F00')+
	annotate('text', label='same', -50, 9.5, size=3.5, color='#708090')+
	scale_y_continuous(limits=c(-13, 13), breaks=seq(-12, 12, 12), labels=c('12', '0', '12'))+
	labs(title='HUES64', x='Distance from Alu start sites (bp)', y='BS-seq mCHG (%)')+
	theme_classic()+
	theme(plot.margin=unit(c(0, .2, 0, .01), 'in'),
		plot.title=element_text(hjust=.5, colour='black', size=13, vjust=-.5),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_text(colour='black', size=12),
		axis.title.y=element_text(colour='black', size=12),
		axis.text.x=element_text(colour='black', size=12),
		axis.text.y=element_text(colour='black', size=12))
ggsave('fig6.HUES64_BS_Alu_start.pdf', p, height=2.5, width=2.7)







mC <- read.table('HUES64-3AKOL_BS_mCHG_Alu_Start.txt', sep='\t', header=T)
mC$Methylation[which(mC$Strand=='oppo')] <- -1*mC$Methylation[which(mC$Strand=='oppo')]

p2 <- ggplot(mC, aes(x=Distance, y=Methylation, color=Strand))+
	geom_bar(fill=NA, stat='identity', width=.1)+
	geom_hline(yintercept=0, linewidth=.3)+scale_color_manual(values=c('#E69F00', '#708090'))+
	annotate('text', label='oppo', -50, -9.5, size=3.5, color='#E69F00')+
	annotate('text', label='same', -50, 9.5, size=3.5, color='#708090')+
	scale_y_continuous(limits=c(-13, 13), breaks=seq(-12, 12, 12), labels=c('12', '0', '12'))+
	labs(title='Dnmt3a-', x='Distance from Alu start sites (bp)', y='BS-seq mCHG (%)')+
	theme_classic()+
	theme(plot.margin=unit(c(0, .2, 0, .01), 'in'),
		plot.title=element_text(hjust=.5, colour='black', face='italic', size=13, vjust=-.5),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_text(colour='black', size=12),
		axis.title.y=element_text(colour='black', size=12),
		axis.text.x=element_text(colour='black', size=12),
		axis.text.y=element_text(colour='black', size=12))
ggsave('fig6.HUES64-3AKOL_BS_Alu_start.pdf', p2, height=2.5, width=2.7)






mC <- read.table('HUES64-3BKOL_BS_mCHG_Alu_Start.txt', sep='\t', header=T)
mC$Methylation[which(mC$Strand=='oppo')] <- -1*mC$Methylation[which(mC$Strand=='oppo')]

p3 <- ggplot(mC, aes(x=Distance, y=Methylation, color=Strand))+
	geom_bar(fill=NA, stat='identity', width=.1)+
	geom_hline(yintercept=0, linewidth=.3)+scale_color_manual(values=c('#E69F00', '#708090'))+
	annotate('text', label='oppo', -50, -9.5, size=3.5, color='#E69F00')+
	annotate('text', label='same', -50, 9.5, size=3.5, color='#708090')+
	scale_y_continuous(limits=c(-13, 13), breaks=seq(-12, 12, 12), labels=c('12', '0', '12'))+
	labs(title='Dnmt3b-', x='Distance from Alu start sites (bp)', y='BS-seq mCHG (%)')+
	theme_classic()+
	theme(plot.margin=unit(c(0, .2, 0, .01), 'in'),
		plot.title=element_text(hjust=.5, colour='black', face='italic', size=13, vjust=-.5),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_text(colour='black', size=12),
		axis.title.y=element_text(colour='black', size=12),
		axis.text.x=element_text(colour='black', size=12),
		axis.text.y=element_text(colour='black', size=12))
ggsave('fig6.HUES64-3BKOL_BS_Alu_start.pdf', height=2.5, p3, width=2.7)




mC <- read.table('HUES64-DKOL_BS_mCHG_Alu_Start.txt', sep='\t', header=T)
mC$Methylation[which(mC$Strand=='oppo')] <- -1*mC$Methylation[which(mC$Strand=='oppo')]

p4 <- ggplot(mC, aes(x=Distance, y=Methylation, color=Strand))+
	geom_bar(fill=NA, stat='identity', width=.1)+
	geom_hline(yintercept=0, linewidth=.3)+scale_color_manual(values=c('#E69F00', '#708090'))+
	annotate('text', label='oppo', -50, -9.5, size=3.5, color='#E69F00')+
	annotate('text', label='same', -50, 9.5, size=3.5, color='#708090')+
	scale_y_continuous(limits=c(-13, 13), breaks=seq(-12, 12, 12), labels=c('12', '0', '12'))+
	labs(title='Dnmt3a-/3b-', x='Distance from Alu start sites (bp)', y='BS-seq mCHG (%)')+
	theme_classic()+
	theme(plot.margin=unit(c(0, .2, 0, .01), 'in'),
		plot.title=element_text(hjust=.5, colour='black', face='italic', size=13, vjust=-.5),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_text(colour='black', size=12),
		axis.title.y=element_text(colour='black', size=12),
		axis.text.x=element_text(colour='black', size=12),
		axis.text.y=element_text(colour='black', size=12))
ggsave('fig6.HUES64-DKOL_BS_Alu_start.pdf', p4, height=2.5, width=2.7)
