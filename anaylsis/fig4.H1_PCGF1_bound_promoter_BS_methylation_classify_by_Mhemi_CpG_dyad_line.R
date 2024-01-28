##-------------------------
# @Author: XiongXiong
# @Date: 2023/7/30
##-------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)
library(mFilter)

meth <- read.table('H1_PCGF1_bound_promoter_BS_methylation.txt', sep='\t')

meth2 <- subset(meth, V1=='cluster 1', select=-V1)
methVal <- as.vector(colMeans(meth2, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
tmp <- mFilter(as.ts(methVal), freq=10)
methylation <- as.vector(tmp$trend)
distance <- rep(1:length(methVal))
expr <- rep('cluster 1', length(methVal))

meth2 <- subset(meth, V1=='cluster 2', select=-V1)
methVal <- as.vector(colMeans(meth2, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
tmp <- mFilter(as.ts(methVal), freq=10)
methylation <- c(methylation, as.vector(tmp$trend))
distance <- c(distance, rep(1:length(methVal)))
expr <- c(expr, rep('cluster 2', length(methVal)))

meth2 <- subset(meth, V1=='cluster 3', select=-V1)
methVal <- as.vector(colMeans(meth2, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
tmp <- mFilter(as.ts(methVal), freq=10)
methylation <- c(methylation, as.vector(tmp$trend))
distance <- c(distance, rep(1:length(methVal)))
expr <- c(expr, rep('cluster 3', length(methVal)))

meth2 <- subset(meth, V1=='cluster 4', select=-V1)
methVal <- as.vector(colMeans(meth2, na.rm=T))
methVal[is.na(methVal)] <- median(methVal, na.rm=T)
tmp <- mFilter(as.ts(methVal), freq=10)
methylation <- c(methylation, as.vector(tmp$trend))
distance <- c(distance, rep(1:length(methVal)))
expr <- c(expr, rep('cluster 4', length(methVal)))

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

pic <- ggplot(methData, aes(x=Distance, y=Methylation, color=Expr))+
	geom_line(lwd=.45)+scale_color_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	geom_vline(xintercept=c(40, 160), color='#708090', linetype='dashed')+
	theme_bw()+labs(title='', x='', y='Methylation (%)')+
	scale_x_continuous(limits=c(-15, 215), breaks=c(0, 40, 160, 200), expand=c(0, 0), label=c('-1kb', 'TSS', 'TTS', '1kb'))+
	scale_y_continuous(limits=c(-8, 108), breaks=seq(0, 100, 25), expand=c(0, 0))+
	theme(plot.margin=unit(c(.02, .1, -.17, .03), 'inches'),
		plot.title=element_blank(),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.5, .15),
		legend.key=element_blank(),
		legend.key.height=unit(.12, 'in'),
		legend.title=element_blank(),
		legend.box="",
		legend.text=element_text(size=9),
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
ggsave('fig4.H1_PCGF1_bound_promoter_BS_methylation_classify_by_Mhemi_CpG_dyad_line.pdf', pic, width=2.5, height=2.3)
