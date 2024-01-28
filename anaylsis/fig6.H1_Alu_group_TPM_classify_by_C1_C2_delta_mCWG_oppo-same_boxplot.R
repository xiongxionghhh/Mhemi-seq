##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/8/10
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
# windowsFonts(consolas=windowsFont("consolas"))
##-----------------------------
library(ggplot2)

exprSet <- read.table('H1_Alu_group_TPM_classify_by_C1_C2_delta_mCWG_oppo-same_filter_by_TPM.bed')
plotData <- exprSet[, c(5,7)]

nCounts1 <- paste0('(n=', format(sum(plotData$V7=='Group1'), big.mark=','), ')')
nCounts2 <- paste0('(n=', format(sum(plotData$V7=='Group2'), big.mark=','), ')')
nCounts3 <- paste0('(n=', format(sum(plotData$V7=='Group3'), big.mark=','), ')')
nCounts4 <- paste0('(n=', format(sum(plotData$V7=='Group4'), big.mark=','), ')')


pValue <- signif(t.test(plotData$V5[which(plotData$V7=='Group1')], plotData$V5[which(plotData$V7=='Group2')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(plotData$V5[which(plotData$V7=='Group2')], plotData$V5[which(plotData$V7=='Group3')])$p.value, 2)
anno2 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(plotData$V5[which(plotData$V7=='Group3')], plotData$V5[which(plotData$V7=='Group4')])$p.value, 2)
anno3 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
pValue <- signif(t.test(plotData$V5[which(plotData$V7=='Group1')], plotData$V5[which(plotData$V7=='Group4')])$p.value, 2)
anno4 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))

plotData$V7 <- gsub('Group', 'Group ', plotData$V7)
plotData$V7 <- factor(plotData$V7, levels=c('Group 1', 'Group 2', 'Group 3', 'Group 4'))


plotData$V5 <- log2(plotData$V5+1)



p <- ggplot(plotData, aes(x=V7, y=V5))+
	theme_classic()+labs(x='', y='TPM', title='')+
#	geom_violin(aes(color=V7), fill=NA, width=.75, size=.45, trim=TRUE, scale='area')+
	geom_boxplot(aes(color=V7), fill=NA, width=.55, size=.35, outlier.size=.5, notch=TRUE)+
	annotate('text', label=nCounts1, x=1, y=.5, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=.5, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=.5, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=.5, size=2.75, color='#2F4F4F')+
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
	scale_y_continuous(limits=c(-.5, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
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
ggsave('fig6.H1_Alu_group_TPM_classify_by_C1_C2_delta_mCWG_oppo-same_filter_by_TPM_boxplot.pdf', p, width=2.5, height=2.5)




p2 <- ggplot(plotData, aes(x=V7, y=V5))+
	theme_classic()+labs(x='', y='TPM', title='')+
	geom_violin(aes(color=V7), fill=NA, width=.75, linewidth=.45, trim=TRUE, scale='area')+
	annotate('text', label=nCounts1, x=1, y=.5, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts2, x=2, y=.5, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts3, x=3, y=.5, size=2.75, color='#2F4F4F')+
	annotate('text', label=nCounts4, x=4, y=.5, size=2.75, color='#2F4F4F')+
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
	scale_y_continuous(limits=c(-.5, 13.5), breaks=seq(0, 12, 4), expand=c(0, 0))+
	theme(plot.margin=unit(c(0, .1, -.15, .02), 'inches'),
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
ggsave('fig6.H1_Alu_group_TPM_classify_by_C1_C2_delta_mCWG_oppo-same_filter_by_TPM_violinplot.pdf', p2, width=2.5, height=2.5)

