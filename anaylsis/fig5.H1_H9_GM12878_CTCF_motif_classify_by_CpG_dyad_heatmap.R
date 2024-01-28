##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/30
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ComplexHeatmap)

oriData <- read.table('GM12878_CTCF_motif_classify_by_hpBS_5th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.GM12878_CTCF_motif_classify_by_hpBS_5th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#FE817D'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 5 (GM12878; hpBS-seq)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()


oriData <- read.table('GM12878_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.GM12878_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#FE817D'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 5 (GM12878; Mhemi-seq)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()


oriData <- read.table('GM12878_CTCF_motif_classify_by_hpBS_7th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.GM12878_CTCF_motif_classify_by_hpBS_7th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#81B8DF'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 7 (GM12878; hpBS-seq)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()



oriData <- read.table('GM12878_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.GM12878_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#81B8DF'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 7 (GM12878; Mhemi-seq)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()

















oriData <- read.table('H1_CTCF_motif_classify_by_iSA_5th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.H1_CTCF_motif_classify_by_iSA_5th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#FE817D'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 5 (H1; iSA)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()


oriData <- read.table('H1_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.H1_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#FE817D'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 5 (H1; Mhemi-seq)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()


oriData <- read.table('H1_CTCF_motif_classify_by_iSA_7th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.H1_CTCF_motif_classify_by_iSA_7th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#81B8DF'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 7 (H1; iSA)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()



oriData <- read.table('H1_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.H1_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#81B8DF'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 7 (H1; Mhemi-seq)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()











oriData <- read.table('H9_CTCF_motif_classify_by_ChIP-hpBS_5th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.H9_CTCF_motif_classify_by_ChIP-hpBS_5th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#FE817D'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 5 (H9; ChIP-hpBS-seq)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()


oriData <- read.table('H9_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.H9_CTCF_motif_classify_by_Mhemi_5th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#FE817D'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 5 (H9; Mhemi-seq)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()


oriData <- read.table('H9_CTCF_motif_classify_by_ChIP-hpBS_7th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.H9_CTCF_motif_classify_by_ChIP-hpBS_7th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#81B8DF'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 7 (H9; ChIP-hpBS-seq)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()



oriData <- read.table('H9_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig5.H9_CTCF_motif_classify_by_Mhemi_7th_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#81B8DF'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='Dyad 7 (H9; Mhemi-seq)', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()
