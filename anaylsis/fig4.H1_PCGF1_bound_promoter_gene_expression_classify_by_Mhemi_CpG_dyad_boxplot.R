##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/30
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
# windowsFonts(consolas=windowsFont("consolas"))
##-----------------------------
library(ggplot2)

exprSet <- read.table('H1_PCGF1_bound_promoter_gene_expression_classify_by_Mhemi_CpG_dyad.txt', sep='\t')
exprSet$nCounts <- c(rep(paste0('(n=', format(sum(exprSet$V2=='cluster 1'), big.mark=','), ')'), sum(exprSet$V2=='cluster 1')),
					rep(paste0('(n=', format(sum(exprSet$V2=='cluster 2'), big.mark=','), ')'), sum(exprSet$V2=='cluster 2')),
					rep(paste0('(n=', format(sum(exprSet$V2=='cluster 3'), big.mark=','), ')'), sum(exprSet$V2=='cluster 3')),
					rep(paste0('(n=', format(sum(exprSet$V2=='cluster 4'), big.mark=','), ')'), sum(exprSet$V2=='cluster 4')))

exprSet$V2 <- factor(exprSet$V2, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

pValue <- signif(t.test(exprSet$V1[which(exprSet$V2=='cluster 1')], exprSet$V1[which(exprSet$V2=='cluster 2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(exprSet$V1[which(exprSet$V2=='cluster 2')], exprSet$V1[which(exprSet$V2=='cluster 3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(exprSet$V1[which(exprSet$V2=='cluster 3')], exprSet$V1[which(exprSet$V2=='cluster 4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(exprSet$V1[which(exprSet$V2=='cluster 1')], exprSet$V1[which(exprSet$V2=='cluster 4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData <- exprSet
plotData$V1 <- log2(plotData$V1+1)


p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y=expression(paste(log[2] (TPM))), title=)+
	geom_violin(aes(color=V2), fill=NA, width=.75, linewidth=.35, trim=TRUE, scale='area')+
#	geom_boxplot(aes(color=V2), fill=NA, width=.4, linewidth=.35, outlier.size=.5, notch=TRUE)+
	geom_text(aes(label=nCounts), y=-1, size=3, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(.02, .01, -.15, .02), 'inches'),
		plot.title=element_blank(),
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
ggsave('fig4.H1_PCGF1_bound_promoter_gene_expression_classify_by_Mhemi_CpG_dyad_violinplot.pdf', p, width=2.5, height=2.5)

p <- ggplot(plotData, aes(x=V2, y=V1))+
	theme_classic()+labs(x='', y=expression(paste(log[2] (TPM))), title=)+
#	geom_violin(aes(color=V2), fill=NA, width=.75, linewidth=.35, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V2), fill=NA, width=.4, linewidth=.35, outlier.size=.5, notch=TRUE)+
	geom_text(aes(label=nCounts), y=-1, size=3, color='#2F4F4F')+
	stat_summary(fun='mean', geom='point', shape=23, size=1, fill='#7CFC00')+
	annotate('text', label=anno1, x=1.5, y=11.2, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=10.7, yend=10.7, linewidth=.25)+
	annotate('text', label=anno2, x=2.5, y=10.0, size=3, parse=TRUE)+
	annotate('segment', x=2, xend=3, y=9.5, yend=9.5, linewidth=.25)+
	annotate('text', label=anno3, x=3.5, y=8.6, size=3, parse=TRUE)+
	annotate('segment', x=3, xend=4, y=8.1, yend=8.1, linewidth=.25)+
	annotate('text', label=anno4, x=2.5, y=12.6, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=4, y=12.1, yend=12.1, linewidth=.25)+
	scale_color_manual(values=c('#386CAF', '#7FC87F', '#9f79d3', '#fc9533'))+
	scale_y_continuous(limits=c(-1.8, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(.02, .01, -.15, .02), 'inches'),
		plot.title=element_blank(),
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
ggsave('fig4.H1_PCGF1_bound_promoter_gene_expression_classify_by_Mhemi_CpG_dyad_boxplot.pdf', p, width=2.5, height=2.5)

