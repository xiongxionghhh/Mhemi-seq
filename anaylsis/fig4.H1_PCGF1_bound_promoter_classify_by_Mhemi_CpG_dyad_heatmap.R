##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/30
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ComplexHeatmap)

oriData <- read.table('H1_PCGF1_bound_promoter_dyad_freq_classify_by_Mhemi_CpG_dyad_1.2_fold.txt', row.names=1)
oriData$V6 <- gsub('cluster', 'cluster ', oriData$V6)
oriData$V6 <- factor(oriData$V6, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

colnames(oriData) <-c('unme', 'same', 'oppo', 'me', 'cluster')
figRes <- as.matrix(oriData[, -ncol(oriData)])

pdf('fig4.H1_PCGF1_bound_promoter_classify_by_Mhemi_CpG_dyad_heatmap.pdf', height=2.5, width=3.2)
print(Heatmap(figRes, width=unit(1.85, 'in'), height=unit(2.0, 'in'),
	cluster_rows=F, cluster_columns=F, col=colorRampPalette(c('#F6F6F6', '#EB835E'))(100),
	row_split=factor(oriData$cluster), row_gap=unit(.01, 'in'),
	row_title_rot=0, row_title_gp=gpar(fontsize=12),
	column_title='PCGF1 bound promoters', column_title_gp=gpar(fontsize=13),
	column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
	heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
		legend_height=unit(.7, 'in'), at=c(0, 100), labels=c(0, 100))))
dev.off()
