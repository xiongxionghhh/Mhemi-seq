##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/10
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)

motif <- read.table('GM12878_CTCF_motif_signal_classify_by_hpBS_5th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 5 (GM12878; hpBS-seq)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.GM12878_CTCF_motif_signal_classify_by_hpBS_5th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)


motif <- read.table('GM12878_CTCF_motif_signal_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 5 (GM12878; iSA)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.GM12878_CTCF_motif_signal_classify_by_iSA_5th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)


motif <- read.table('GM12878_CTCF_motif_signal_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 5 (GM12878; Mhemi-seq)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.GM12878_CTCF_motif_signal_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)



motif <- read.table('GM12878_CTCF_motif_signal_classify_by_hpBS_7th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 7 (GM12878; hpBS-seq)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.GM12878_CTCF_motif_signal_classify_by_hpBS_7th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)

motif <- read.table('GM12878_CTCF_motif_signal_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 7 (GM12878; iSA)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.GM12878_CTCF_motif_signal_classify_by_iSA_7th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)




motif <- read.table('GM12878_CTCF_motif_signal_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 7 (GM12878; Mhemi-seq)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.GM12878_CTCF_motif_signal_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)














motif <- read.table('H1_CTCF_motif_signal_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 5 (H1; iSA)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.H1_CTCF_motif_signal_classify_by_iSA_5th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)


motif <- read.table('H1_CTCF_motif_signal_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 5 (H1; Mhemi-seq)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.H1_CTCF_motif_signal_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)



motif <- read.table('H1_CTCF_motif_signal_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 7 (H1; iSA)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.H1_CTCF_motif_signal_classify_by_iSA_7th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)


motif <- read.table('H1_CTCF_motif_signal_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 7 (H1; Mhemi-seq)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.H1_CTCF_motif_signal_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)





motif <- read.table('H9_CTCF_motif_signal_classify_by_ChIP-hpBS_5th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 5 (H9; ChIP-hpBS-seq)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.H9_CTCF_motif_signal_classify_by_ChIP-hpBS_5th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)


motif <- read.table('H9_CTCF_motif_signal_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 5 (H9; Mhemi-seq)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.H9_CTCF_motif_signal_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)



motif <- read.table('H9_CTCF_motif_signal_classify_by_ChIP-hpBS_7th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 7 (H9; ChIP-hpBS-seq)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.H9_CTCF_motif_signal_classify_by_ChIP-hpBS_7th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)


motif <- read.table('H9_CTCF_motif_signal_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt')

pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster2')], motif$V1[which(motif$V2=='cluster3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster3')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(motif$V1[which(motif$V2=='cluster1')], motif$V1[which(motif$V2=='cluster4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- subset(motif, V1<=10.5)
nCounts1 <- paste0('(n=', format(sum(plotData$V2=='cluster1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V2=='cluster2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V2=='cluster3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V2=='cluster4'), big.mark=','), ')')


plotData$V2 <- gsub('cluster', 'cluster ', plotData$V2)
plotData$V2 <- factor(plotData$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y='CTCF ChIP-seq (RPM)', title='Dyad 7 (H9; Mhemi-seq)')+
#	geom_violin(aes(color=V2), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=-1, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=-1, size=2.75, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#0072B2', '#E69F00', '#5F9EA0', '#646464'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
		plot.title=element_text(colour='black', size=10.5, hjust=.5, vjust=-1),
		panel.grid=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.5, color='black'),
		legend.position='none',
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_line(colour='black', linewidth=.3),
		axis.ticks.y=element_line(colour='black', linewidth=.3),
		axis.text.x=element_text(colour='black', size=10, angle=15, hjust=1),
		axis.text.y=element_text(colour='black', size=10),
		axis.title.x=element_text(colour='black', size=10),
		axis.title.y=element_text(colour='black', size=10))
ggsave('fig5.H9_CTCF_motif_signal_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold_boxplot.pdf', p, width=2.5, height=2.5)









