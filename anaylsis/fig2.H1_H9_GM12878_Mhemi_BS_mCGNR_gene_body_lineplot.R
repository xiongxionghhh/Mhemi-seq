options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)

meth <- read.table('GM12878_Mhemi_CpG_gene_body_1kb_up_down.txt')
methVal <- as.vector(colMeans(subset(meth, V1=='high', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- methVal

methVal <- as.vector(colMeans(subset(meth, V1=='low', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- c(mCset, methVal)

exprSet <- c(rep('High', 200), rep('Low', 200))
poSet <- rep(1:200, 2)

methData <- data.frame('Distance'=poSet, 'Methylation'=mCset, 'Expr'=exprSet)
methData$Expr <- factor(methData$Expr, levels=c('High', 'Low'))

pic <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
	geom_line(lwd=.5)+scale_color_manual(values=c('#FF9521', '#00498F'))+
	geom_vline(xintercept=c(40, 160), color='#708090', linetype='dashed')+
	theme_bw()+labs(title='Mhemi-seq (GM12878)', x='', y='Methylation (%)')+
	scale_x_continuous(limits=c(-15, 215), breaks=c(0, 40, 160, 200), expand=c(0, 0), label=c('-1kb', 'TSS', 'TTS', '1kb'))+
	scale_y_continuous(limits=c(-8, 108), breaks=seq(0, 100, 25), expand=c(0, 0))+
	theme(plot.margin=unit(c(-.05, .1, -.17, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.5, .15),
		legend.key=element_blank(),
		legend.key.height=unit(.2, 'in'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=10),
		legend.spacing.x=unit(.02, 'in'),
		legend.spacing.y=unit(-.05, 'in'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.title.x=element_text(size=10),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=.5))
ggsave('fig2.GM12878_Mhemi_CpG_gene_body_1kb_up_down.pdf', pic, width=2.5, height=2.3)


meth <- read.table('GM12878_BS_CpG_gene_body_1kb_up_down.txt')
methVal <- as.vector(colMeans(subset(meth, V1=='High', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- methVal

methVal <- as.vector(colMeans(subset(meth, V1=='Low', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- c(mCset, methVal)

exprSet <- c(rep('High', 200), rep('Low', 200))
poSet <- rep(1:200, 2)

methData <- data.frame('Distance'=poSet, 'Methylation'=mCset, 'Expr'=exprSet)
methData$Expr <- factor(methData$Expr, levels=c('High', 'Low'))

pic2 <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
	geom_line(lwd=.5)+scale_color_manual(values=c('#FF9521', '#00498F'))+
	geom_vline(xintercept=c(40, 160), color='#708090', linetype='dashed')+
	theme_bw()+labs(title='BS-seq (GM12878)', x='', y='Methylation (%)')+
	scale_x_continuous(limits=c(-15, 215), breaks=c(0, 40, 160, 200), expand=c(0, 0), label=c('-1kb', 'TSS', 'TTS', '1kb'))+
	scale_y_continuous(limits=c(-8, 108), breaks=seq(0, 100, 25), expand=c(0, 0))+
	theme(plot.margin=unit(c(-.05, .1, -.17, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.5, .15),
		legend.key=element_blank(),
		legend.key.height=unit(.2, 'in'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=10),
		legend.spacing.x=unit(.02, 'in'),
		legend.spacing.y=unit(-.05, 'in'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.title.x=element_text(size=10),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=.5))
ggsave('fig2.GM12878_BS_CpG_gene_body_1kb_up_down.pdf', pic2, width=2.5, height=2.3)



meth <- read.table('H1_Mhemi_CpG_gene_body_1kb_up_down.txt')
methVal <- as.vector(colMeans(subset(meth, V1=='High', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- methVal

methVal <- as.vector(colMeans(subset(meth, V1=='Low', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- c(mCset, methVal)

exprSet <- c(rep('High', 200), rep('Low', 200))
poSet <- rep(1:200, 2)

methData <- data.frame('Distance'=poSet, 'Methylation'=mCset, 'Expr'=exprSet)
methData$Expr <- factor(methData$Expr, levels=c('High', 'Low'))

pic <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
	geom_line(lwd=.5)+scale_color_manual(values=c('#FF9521', '#00498F'))+
	geom_vline(xintercept=c(40, 160), color='#708090', linetype='dashed')+
	theme_bw()+labs(title='Mhemi-seq (H1)', x='', y='Methylation (%)')+
	scale_x_continuous(limits=c(-15, 215), breaks=c(0, 40, 160, 200), expand=c(0, 0), label=c('-1kb', 'TSS', 'TTS', '1kb'))+
	scale_y_continuous(limits=c(-8, 108), breaks=seq(0, 100, 25), expand=c(0, 0))+
	theme(plot.margin=unit(c(-.05, .1, -.17, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.5, .15),
		legend.key=element_blank(),
		legend.key.height=unit(.2, 'in'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=10),
		legend.spacing.x=unit(.02, 'in'),
		legend.spacing.y=unit(-.05, 'in'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.title.x=element_text(size=10),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=.5))
ggsave('fig2.H1_Mhemi_CpG_gene_body_1kb_up_down.pdf', pic, width=2.5, height=2.3)


meth <- read.table('H1_BS_CpG_gene_body_1kb_up_down.txt')
methVal <- as.vector(colMeans(subset(meth, V1=='High', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- methVal

methVal <- as.vector(colMeans(subset(meth, V1=='Low', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- c(mCset, methVal)

exprSet <- c(rep('High', 200), rep('Low', 200))
poSet <- rep(1:200, 2)

methData <- data.frame('Distance'=poSet, 'Methylation'=mCset, 'Expr'=exprSet)
methData$Expr <- factor(methData$Expr, levels=c('High', 'Low'))

pic2 <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
	geom_line(lwd=.5)+scale_color_manual(values=c('#FF9521', '#00498F'))+
	geom_vline(xintercept=c(40, 160), color='#708090', linetype='dashed')+
	theme_bw()+labs(title='BS-seq (H1)', x='', y='Methylation (%)')+
	scale_x_continuous(limits=c(-15, 215), breaks=c(0, 40, 160, 200), expand=c(0, 0), label=c('-1kb', 'TSS', 'TTS', '1kb'))+
	scale_y_continuous(limits=c(-8, 108), breaks=seq(0, 100, 25), expand=c(0, 0))+
	theme(plot.margin=unit(c(-.05, .1, -.17, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.5, .15),
		legend.key=element_blank(),
		legend.key.height=unit(.2, 'in'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=10),
		legend.spacing.x=unit(.02, 'in'),
		legend.spacing.y=unit(-.05, 'in'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.title.x=element_text(size=10),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=.5))
ggsave('fig2.H1_BS_CpG_gene_body_1kb_up_down.pdf', pic2, width=2.5, height=2.3)



meth <- read.table('H9_Mhemi_CpG_gene_body_1kb_up_down.txt')
methVal <- as.vector(colMeans(subset(meth, V1=='High', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- methVal

methVal <- as.vector(colMeans(subset(meth, V1=='Low', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- c(mCset, methVal)

exprSet <- c(rep('High', 200), rep('Low', 200))
poSet <- rep(1:200, 2)

methData <- data.frame('Distance'=poSet, 'Methylation'=mCset, 'Expr'=exprSet)
methData$Expr <- factor(methData$Expr, levels=c('High', 'Low'))

pic <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
	geom_line(lwd=.5)+scale_color_manual(values=c('#FF9521', '#00498F'))+
	geom_vline(xintercept=c(40, 160), color='#708090', linetype='dashed')+
	theme_bw()+labs(title='Mhemi-seq (H9)', x='', y='Methylation (%)')+
	scale_x_continuous(limits=c(-15, 215), breaks=c(0, 40, 160, 200), expand=c(0, 0), label=c('-1kb', 'TSS', 'TTS', '1kb'))+
	scale_y_continuous(limits=c(-8, 108), breaks=seq(0, 100, 25), expand=c(0, 0))+
	theme(plot.margin=unit(c(-.05, .1, -.17, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.5, .15),
		legend.key=element_blank(),
		legend.key.height=unit(.2, 'in'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=10),
		legend.spacing.x=unit(.02, 'in'),
		legend.spacing.y=unit(-.05, 'in'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.title.x=element_text(size=10),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=.5))
ggsave('fig2.H9_Mhemi_CpG_gene_body_1kb_up_down.pdf', pic, width=2.5, height=2.3)


meth <- read.table('H9_BS_CpG_gene_body_1kb_up_down.txt')
methVal <- as.vector(colMeans(subset(meth, V1=='High', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- methVal

methVal <- as.vector(colMeans(subset(meth, V1=='Low', select=-1), na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
mCset <- c(mCset, methVal)

exprSet <- c(rep('High', 200), rep('Low', 200))
poSet <- rep(1:200, 2)

methData <- data.frame('Distance'=poSet, 'Methylation'=mCset, 'Expr'=exprSet)
methData$Expr <- factor(methData$Expr, levels=c('High', 'Low'))

pic2 <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
	geom_line(lwd=.5)+scale_color_manual(values=c('#FF9521', '#00498F'))+
	geom_vline(xintercept=c(40, 160), color='#708090', linetype='dashed')+
	theme_bw()+labs(title='BS-seq (H9)', x='', y='Methylation (%)')+
	scale_x_continuous(limits=c(-15, 215), breaks=c(0, 40, 160, 200), expand=c(0, 0), label=c('-1kb', 'TSS', 'TTS', '1kb'))+
	scale_y_continuous(limits=c(-8, 108), breaks=seq(0, 100, 25), expand=c(0, 0))+
	theme(plot.margin=unit(c(-.05, .1, -.17, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.5, .15),
		legend.key=element_blank(),
		legend.key.height=unit(.2, 'in'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=10),
		legend.spacing.x=unit(.02, 'in'),
		legend.spacing.y=unit(-.05, 'in'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.title.x=element_text(size=10),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=.5))
ggsave('fig2.H9_BS_CpG_gene_body_1kb_up_down.pdf', pic2, width=2.5, height=2.3)


