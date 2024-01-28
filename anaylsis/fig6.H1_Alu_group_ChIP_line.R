options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)

for(i in c('H3K4me2', 'H3K4me3', 'H3K27ac', 'H3K27me3', 'H3K36me2')){
	
	meth <- read.table(paste0('H1_Alu_group_', i, '_filter_by_TPM.txt'))
	tmpSignal <- subset(meth, V201=='Group1', select=1:200)
	methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
	methVal[is.na(methVal)] <- median(methVal, na.rm=T)
	methylation <- methVal
	distance <- rep(1:length(methVal))
	expr <- rep('Group1', length(methVal))

	tmpSignal <- subset(meth, V201=='Group2', select=1:200)
	methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
	methVal[is.na(methVal)] <- median(methVal, na.rm=T)
	methylation <- c(methylation, methVal)
	distance <- c(distance, rep(1:length(methVal)))
	expr <- c(expr, rep('Group2', length(methVal)))

	tmpSignal <- subset(meth, V201=='Group3', select=1:200)
	methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
	methVal[is.na(methVal)] <- median(methVal, na.rm=T)
	methylation <- c(methylation, methVal)
	distance <- c(distance, rep(1:length(methVal)))
	expr <- c(expr, rep('Group3', length(methVal)))

	tmpSignal <- subset(meth, V201=='Group4', select=1:200)
	methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
	methVal[is.na(methVal)] <- median(methVal, na.rm=T)
	methylation <- c(methylation, methVal)
	distance <- c(distance, rep(1:length(methVal)))
	expr <- c(expr, rep('Group4', length(methVal)))

	methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
	methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

	p <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
		geom_line(lwd=.5)+scale_color_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
		geom_vline(xintercept=c(50, 150), color='#000000', linetype='dashed', lwd=.5)+
		theme_bw()+labs(title=i, x='', y='RPKM')+
		scale_x_continuous(limits=c(-10, 210), breaks=c(50, 150), expand=c(0, 0), label=c('C1', 'C2'))+
		scale_y_continuous(limits=c(-.05, .65), breaks=c(0, .6), expand=c(0, 0))+
		theme(plot.margin=unit(c(-.05, .01, -.2, .01), 'inches'),
			plot.title=element_text(size=11, hjust=.5, vjust=-1),
			panel.grid=element_blank(),
			panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
			legend.position='none',
			axis.line.x=element_blank(),
			axis.line.y=element_blank(),
			axis.ticks.x=element_blank(),
			axis.ticks.y=element_line(colour='black', linewidth=.3),
			axis.title.x=element_text(size=10, colour='black'),
			axis.title.y=element_text(size=10, colour='black'),
			axis.text.x=element_text(size=10, colour='black'),
			axis.text.y=element_text(size=10, colour='black'))
	ggsave(paste0('fig6.H1_Alu_group_', i, '_filter_by_TPM_line.pdf'), p, width=1.8, height=1.3)

}

meth <- read.table(paste0('H1_Alu_group_', 'H3K36me3', '_filter_by_TPM.txt'))
tmpSignal <- subset(meth, V201=='Group1', select=1:200)
methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
methylation <- methVal
distance <- rep(1:length(methVal))
expr <- rep('Group1', length(methVal))

tmpSignal <- subset(meth, V201=='Group2', select=1:200)
methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
methylation <- c(methylation, methVal)
distance <- c(distance, rep(1:length(methVal)))
expr <- c(expr, rep('Group2', length(methVal)))

tmpSignal <- subset(meth, V201=='Group3', select=1:200)
methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
methylation <- c(methylation, methVal)
distance <- c(distance, rep(1:length(methVal)))
expr <- c(expr, rep('Group3', length(methVal)))

tmpSignal <- subset(meth, V201=='Group4', select=1:200)
methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
methylation <- c(methylation, methVal)
distance <- c(distance, rep(1:length(methVal)))
expr <- c(expr, rep('Group4', length(methVal)))

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

p <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
	geom_line(lwd=.5)+scale_color_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	geom_vline(xintercept=c(50, 150), color='#000000', linetype='dashed', lwd=.5)+
	annotate('text', label=c('Group 1', 'Group 2', 'Group 3', 'Group 4'), x=c(20,20,179,179), y=c(.5,.4,.5,.4), size=2.3, color=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	theme_bw()+labs(title='H3K36me3', x='', y='RPKM')+
	scale_x_continuous(limits=c(-10, 210), breaks=c(50, 150), expand=c(0, 0), label=c('C1', 'C2'))+
	scale_y_continuous(limits=c(-.05, .65), breaks=c(0, .6), expand=c(0, 0))+
	theme(plot.margin=unit(c(-.05, .01, -.1, .01), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.title.x=element_text(size=10, colour='black'),
		axis.title.y=element_text(size=10, colour='black'),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black'))
ggsave(paste0('fig6.H1_Alu_group_', 'H3K36me3', '_filter_by_TPM_line.pdf'), p, width=1.8, height=1.3)






for(i in c('H3K4me2', 'H3K4me3', 'H3K27ac', 'H3K27me3', 'H3K36me2')){
	
	meth <- read.table(paste0('H1_Alu_group_', i, '.txt'))
	tmpSignal <- subset(meth, V201=='Group1', select=1:200)
	methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
	methVal[is.na(methVal)] <- median(methVal, na.rm=T)
	methylation <- methVal
	distance <- rep(1:length(methVal))
	expr <- rep('Group1', length(methVal))

	tmpSignal <- subset(meth, V201=='Group2', select=1:200)
	methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
	methVal[is.na(methVal)] <- median(methVal, na.rm=T)
	methylation <- c(methylation, methVal)
	distance <- c(distance, rep(1:length(methVal)))
	expr <- c(expr, rep('Group2', length(methVal)))

	tmpSignal <- subset(meth, V201=='Group3', select=1:200)
	methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
	methVal[is.na(methVal)] <- median(methVal, na.rm=T)
	methylation <- c(methylation, methVal)
	distance <- c(distance, rep(1:length(methVal)))
	expr <- c(expr, rep('Group3', length(methVal)))

	tmpSignal <- subset(meth, V201=='Group4', select=1:200)
	methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
	methVal[is.na(methVal)] <- median(methVal, na.rm=T)
	methylation <- c(methylation, methVal)
	distance <- c(distance, rep(1:length(methVal)))
	expr <- c(expr, rep('Group4', length(methVal)))

	methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
	methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

	p <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
		geom_line(lwd=.5)+scale_color_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
		geom_vline(xintercept=c(50, 150), color='#000000', linetype='dashed', lwd=.5)+
		theme_bw()+labs(title=i, x='', y='RPKM')+
		scale_x_continuous(limits=c(-10, 210), breaks=c(50, 150), expand=c(0, 0), label=c('C1', 'C2'))+
		scale_y_continuous(limits=c(-.05, .65), breaks=c(0, .6), expand=c(0, 0))+
		theme(plot.margin=unit(c(-.05, .01, -.2, .01), 'inches'),
			plot.title=element_text(size=11, hjust=.5, vjust=-1),
			panel.grid=element_blank(),
			panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
			legend.position='none',
			axis.line.x=element_blank(),
			axis.line.y=element_blank(),
			axis.ticks.x=element_blank(),
			axis.ticks.y=element_line(colour='black', linewidth=.3),
			axis.title.x=element_text(size=10, colour='black'),
			axis.title.y=element_text(size=10, colour='black'),
			axis.text.x=element_text(size=10, colour='black'),
			axis.text.y=element_text(size=10, colour='black'))
	ggsave(paste0('fig6.H1_Alu_group_', i, '_line.pdf'), p, width=1.8, height=1.3)

}

meth <- read.table(paste0('H1_Alu_group_', 'H3K36me3', '.txt'))
tmpSignal <- subset(meth, V201=='Group1', select=1:200)
methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
methylation <- methVal
distance <- rep(1:length(methVal))
expr <- rep('Group1', length(methVal))

tmpSignal <- subset(meth, V201=='Group2', select=1:200)
methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
methylation <- c(methylation, methVal)
distance <- c(distance, rep(1:length(methVal)))
expr <- c(expr, rep('Group2', length(methVal)))

tmpSignal <- subset(meth, V201=='Group3', select=1:200)
methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
methylation <- c(methylation, methVal)
distance <- c(distance, rep(1:length(methVal)))
expr <- c(expr, rep('Group3', length(methVal)))

tmpSignal <- subset(meth, V201=='Group4', select=1:200)
methVal <- as.vector(colMeans(tmpSignal, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
methylation <- c(methylation, methVal)
distance <- c(distance, rep(1:length(methVal)))
expr <- c(expr, rep('Group4', length(methVal)))

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

p <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
	geom_line(lwd=.5)+scale_color_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	geom_vline(xintercept=c(50, 150), color='#000000', linetype='dashed', lwd=.5)+
	annotate('text', label=c('Group 1', 'Group 2', 'Group 3', 'Group 4'), x=c(20,20,179,179), y=c(.5,.4,.5,.4), size=2.3, color=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	theme_bw()+labs(title='H3K36me3', x='', y='RPKM')+
	scale_x_continuous(limits=c(-10, 210), breaks=c(50, 150), expand=c(0, 0), label=c('C1', 'C2'))+
	scale_y_continuous(limits=c(-.05, .65), breaks=c(0, .6), expand=c(0, 0))+
	theme(plot.margin=unit(c(-.05, .01, -.1, .01), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.title.x=element_text(size=10, colour='black'),
		axis.title.y=element_text(size=10, colour='black'),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black'))
ggsave(paste0('fig6.H1_Alu_group_', 'H3K36me3', '_line.pdf'), p, width=1.8, height=1.3)

