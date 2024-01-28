##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/23
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)
library(mFilter)

load('/mnt/disk1/5/share/others/colorset.RData')

meth <- read.table('GM12878_Mhemi_mYNCGNR_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('GM12878_Mhemi_mYNCGNR_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='Mhemi-seq (GM12878)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-1.5, 16.5), breaks=seq(0, 15, 5), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.GM12878_Mhemi_mYNCGNR_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)




meth <- read.table('GM12878_hpBS_mYNCGNR_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('GM12878_hpBS_mYNCGNR_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='hpBS-seq (GM12878)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-0.5, 4.5), breaks=seq(0, 4, 1), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.GM12878_hpBS_mYNCGNR_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)



meth <- read.table('H1_Mhemi_mYNCGNR_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H1_Mhemi_mYNCGNR_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='Mhemi-seq (H1)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-1.5, 16.5), breaks=seq(0, 15, 5), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H1_Mhemi_mYNCGNR_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)




meth <- read.table('H1_hpBS_mYNCGNR_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H1_hpBS_mYNCGNR_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='hpBS-seq (H1)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-1, 10), breaks=seq(0, 9, 3), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H1_hpBS_mYNCGNR_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)



meth <- read.table('H9_Mhemi_mYNCGNR_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H9_Mhemi_mYNCGNR_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='Mhemi-seq (H9)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-1.5, 16.5), breaks=seq(0, 15, 5), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H9_Mhemi_mYNCGNR_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)


meth <- read.table('H9_ChIP-hpBS_mYNCGNR_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H9_ChIP-hpBS_mYNCGNR_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=30))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='ChIP-hpBS-seq (H9)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-1.5, 15.5), breaks=seq(0, 15, 5), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H9_ChIP-hpBS_mYNCGNR_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)



meth <- read.table('GM12878_hpBS_mCpG_dyad_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=10))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('GM12878_hpBS_mCpG_dyad_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=10))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='hpBS-seq (GM12878)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-0.5, 4.5), breaks=seq(0, 4, 1), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.GM12878_hpBS_mCpG_dyad_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)




meth <- read.table('H1_hpBS_mCpG_dyad_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=10))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H1_hpBS_mCpG_dyad_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=10))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='hpBS-seq (H1)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-1, 10), breaks=seq(0, 9, 3), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H1_hpBS_mCpG_dyad_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)




meth <- read.table('H9_ChIP-hpBS_mCpG_dyad_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=10))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H9_ChIP-hpBS_mCpG_dyad_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=10))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='ChIP-hpBS-seq (H9)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-1.5, 15.5), breaks=seq(0, 15, 5), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H9_ChIP-hpBS_mCpG_dyad_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)


meth <- read.table('H1_iSA_mCpG_dyad_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=10))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H1_iSA_mCpG_dyad_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=10))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='BS-seq+iSA (H1)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-1.5, 15), breaks=seq(0, 15, 5), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H1_iSA_mCpG_dyad_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)




meth <- read.table('H9_nasBS-iSA_mCpG_dyad_CTCF_phasing_OPPO.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=10))
methylation <- as.vector(tmp$trend)
expr <- rep('oppo', nrow(meth))
distance <- meth$Distance

meth <- read.table('H9_nasBS-iSA_mCpG_dyad_CTCF_phasing_SAME.txt', header=T)
tmp <- mFilter(as.ts(meth$Methylation, freq=10))
methylation <- c(methylation, as.vector(tmp$trend))
expr <- c(expr, rep('same', nrow(meth)))
distance <- c(distance, meth$Distance)

methData <- data.frame('Distance'=distance, 'Methylation'=methylation, 'Expr'=expr)
methData$Expr <- factor(methData$Expr, levels=unique(methData$Expr))

pic <- ggplot(methData, aes(x=Distance))+
	geom_line(lwd=.4, aes(y=Methylation, color=Expr))+scale_color_manual(values=colorSet[1:2])+
	geom_vline(xintercept=c(-160, 160), lwd=.4, color='#708090', linetype='dashed')+
	theme_bw()+labs(title='nasBS-seq+iSA (H9)', x='Distance from CTCF motifs (bp)', y='Hemi-methylation (%)')+
	scale_x_continuous(limits=c(-200, 200), breaks=c(-160, 0, 160), expand=c(0, 0))+
	scale_y_continuous(limits=c(-.5, 6.5), breaks=seq(0, 6, 2), expand=c(0, 0))+
	guides(color=guide_legend(ncol=1))+
	theme(plot.margin=unit(c(0, .3, 0, .03), 'in'),
		plot.title=element_text(size=11, hjust=.5, vjust=-.75),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.background=element_blank(),
		legend.position=c(.29, .12),
		legend.key=element_blank(),
		legend.key.height=unit(.15, 'in'),
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
		axis.text.y=element_text(size=10, colour='black'))
ggsave('fig2.H9_nasBS-iSA_mCpG_dyad_CTCF_phasing_lineplot.pdf', pic, height=2.35, width=2.5)



