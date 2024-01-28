##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/22
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
# windowsFonts(consolas=windowsFont("consolas"))
##-----------------------------
library(ggplot2)
library(mFilter)

load('/mnt/disk1/5/share/others/colorset.RData')
meth <- read.table('GM12878_Mhemi_CpG_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('GM12878_Mhemi_CpG_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

meth <- read.table('GM12878_MNase_CTCF_phasing.txt', header=T)
tmp <- meth$Signal

methData <- data.frame('Distance'=distance,
						'Methylation'=methylation,
						'Expr'=expr,
						'MNase'=rep(tmp*32-.2, 2))
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_line(lwd=.4, aes(y=MNase), color='#5F9EA0')+
	annotate('segment', x=250, xend=400, y=5, yend=5, linewidth=.4, color='#5F9EA0')+
	annotate('text', x=680, y=5, label='MNase', color='#000000', size=3)+
	geom_vline(xintercept=c(-160, 160), color='#708090', linetype='dashed', linewidth=.4)+
	theme_bw()+labs(title='Mhemi-seq (GM12878)', x='Distance from CTCF motifs (bp)')+
	scale_x_continuous(limits=c(-1100, 1100), breaks=seq(-1000, 1050, 500), expand=c(0, 0))+
	scale_y_continuous(limits=c(-5.5, 85.5), breaks=seq(0, 80, 20), expand=c(0, 0), 
		name='Methylation (%)',
		sec.axis=sec_axis(~./32+.2, breaks=seq(.2, 2.6, .8), name='MNase-seq signal'))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(-.05, 0, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.2, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
		legend.title=element_blank(),
		legend.box="",
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.GM12878_Mhemi_CpG_MNase_CTCF_phasing.pdf', pic, height=2.3, width=2.75)



meth <- read.table('GM12878_BS_CGNR_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('GM12878_BS_CGNR_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

meth <- read.table('GM12878_MNase_CTCF_phasing.txt', header=T)
tmp <- meth$Signal

methData <- data.frame('Distance'=distance,
						'Methylation'=methylation,
						'Expr'=expr,
						'MNase'=rep(tmp*32-.2, 2))
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_line(lwd=.4, aes(y=MNase), color='#5F9EA0')+
	annotate('segment', x=250, xend=400, y=5, yend=5, linewidth=.4, color='#5F9EA0')+
	annotate('text', x=680, y=5, label='MNase', color='#000000', size=3)+
	geom_vline(xintercept=c(-160, 160), color='#708090', linetype='dashed', linewidth=.4)+
	theme_bw()+labs(title='BS-seq (GM12878)', x='Distance from CTCF motifs (bp)')+
	scale_x_continuous(limits=c(-1100, 1100), breaks=seq(-1000, 1050, 500), expand=c(0, 0))+
	scale_y_continuous(limits=c(-5.5, 85.5), breaks=seq(0, 80, 20), expand=c(0, 0), 
		name='Methylation (%)',
		sec.axis=sec_axis(~./32+.2, breaks=seq(.2, 2.6, .8), name='MNase-seq signal'))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(-.05, 0, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.2, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
		legend.title=element_blank(),
		legend.box="",
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.GM12878_BS_CGNR_MNase_CTCF_phasing.pdf', pic, height=2.3, width=2.75)


meth <- read.table('H1_Mhemi_CpG_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H1_Mhemi_CpG_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

meth <- read.table('H1_MNase_CTCF_phasing.txt', header=T)
tmp <- meth$Signal

methData <- data.frame('Distance'=distance,
						'Methylation'=methylation,
						'Expr'=expr,
						'MNase'=rep(tmp*70, 2))
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_line(lwd=.4, aes(y=MNase), color='#5F9EA0')+
	annotate('segment', x=250, xend=400, y=5, yend=5, linewidth=.4, color='#5F9EA0')+
	annotate('text', x=680, y=5, label='MNase', color='#000000', size=3)+
	geom_vline(xintercept=c(-160, 160), color='#708090', linetype='dashed', linewidth=.4)+
	theme_bw()+labs(title='Mhemi-seq (H1)', x='Distance from CTCF motifs (bp)')+
	scale_x_continuous(limits=c(-1100, 1100), breaks=seq(-1000, 1050, 500), expand=c(0, 0))+
	scale_y_continuous(limits=c(-5.5, 85.5), breaks=seq(0, 80, 20), expand=c(0, 0), 
		name='Methylation (%)',
		sec.axis=sec_axis(~./70, breaks=seq(0, 1.2, .4), name='MNase-seq signal'))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(-.05, 0, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.2, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
		legend.title=element_blank(),
		legend.box="",
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H1_Mhemi_CpG_MNase_CTCF_phasing.pdf', pic, height=2.3, width=2.75)



meth <- read.table('H1_BS_CGNR_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H1_BS_CGNR_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

meth <- read.table('H1_MNase_CTCF_phasing.txt', header=T)
tmp <- meth$Signal

methData <- data.frame('Distance'=distance,
						'Methylation'=methylation,
						'Expr'=expr,
						'MNase'=rep(tmp*70, 2))
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_line(lwd=.4, aes(y=MNase), color='#5F9EA0')+
	annotate('segment', x=250, xend=400, y=5, yend=5, linewidth=.4, color='#5F9EA0')+
	annotate('text', x=680, y=5, label='MNase', color='#000000', size=3)+
	geom_vline(xintercept=c(-160, 160), color='#708090', linetype='dashed', linewidth=.4)+
	theme_bw()+labs(title='BS-seq (H1)', x='Distance from CTCF motifs (bp)')+
	scale_x_continuous(limits=c(-1100, 1100), breaks=seq(-1000, 1050, 500), expand=c(0, 0))+
	scale_y_continuous(limits=c(-5.5, 85.5), breaks=seq(0, 80, 20), expand=c(0, 0), 
		name='Methylation (%)',
		sec.axis=sec_axis(~./70, breaks=seq(0, 1.2, .4), name='MNase-seq signal'))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(-.05, 0, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.2, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
		legend.title=element_blank(),
		legend.box="",
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H1_BS_CGNR_MNase_CTCF_phasing.pdf', pic, height=2.3, width=2.75)


meth <- read.table('H9_Mhemi_CpG_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H9_Mhemi_CpG_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

meth <- read.table('H9_MNase_CTCF_phasing.txt', header=T)
tmp <- meth$Signal

methData <- data.frame('Distance'=distance,
						'Methylation'=methylation,
						'Expr'=expr,
						'MNase'=rep(tmp*133.33, 2))
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_line(lwd=.4, aes(y=MNase), color='#5F9EA0')+
	annotate('segment', x=250, xend=400, y=5, yend=5, linewidth=.4, color='#5F9EA0')+
	annotate('text', x=680, y=5, label='MNase', color='#000000', size=3)+
	geom_vline(xintercept=c(-160, 160), color='#708090', linetype='dashed', linewidth=.4)+
	theme_bw()+labs(title='Mhemi-seq (H9)', x='Distance from CTCF motifs (bp)')+
	scale_x_continuous(limits=c(-1100, 1100), breaks=seq(-1000, 1050, 500), expand=c(0, 0))+
	scale_y_continuous(limits=c(-5.5, 85.5), breaks=seq(0, 80, 20), expand=c(0, 0), 
		name='Methylation (%)',
		sec.axis=sec_axis(~./133.33, breaks=seq(0, .6, .2), name='MNase-seq signal'))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(-.05, 0, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.2, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
		legend.title=element_blank(),
		legend.box="",
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H9_Mhemi_CpG_MNase_CTCF_phasing.pdf', pic, height=2.3, width=2.75)



meth <- read.table('H9_BS_CGNR_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H9_BS_CGNR_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=50))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

meth <- read.table('H9_MNase_CTCF_phasing.txt', header=T)
tmp <- meth$Signal

methData <- data.frame('Distance'=distance,
						'Methylation'=methylation,
						'Expr'=expr,
						'MNase'=rep(tmp*133.33, 2))
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_line(lwd=.4, aes(y=MNase), color='#5F9EA0')+
	annotate('segment', x=250, xend=400, y=5, yend=5, linewidth=.4, color='#5F9EA0')+
	annotate('text', x=680, y=5, label='MNase', color='#000000', size=3)+
	geom_vline(xintercept=c(-160, 160), color='#708090', linetype='dashed', linewidth=.4)+
	theme_bw()+labs(title='BS-seq (H9)', x='Distance from CTCF motifs (bp)')+
	scale_x_continuous(limits=c(-1100, 1100), breaks=seq(-1000, 1050, 500), expand=c(0, 0))+
	scale_y_continuous(limits=c(-5.5, 85.5), breaks=seq(0, 80, 20), expand=c(0, 0), 
		name='Methylation (%)',
		sec.axis=sec_axis(~./133.33, breaks=seq(0, .6, .2), name='MNase-seq signal'))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(-.05, 0, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.2, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
		legend.title=element_blank(),
		legend.box="",
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H9_BS_CGNR_MNase_CTCF_phasing.pdf', pic, height=2.3, width=2.75)
