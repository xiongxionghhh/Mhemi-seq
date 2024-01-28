##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/25
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
# options(expressions=50000)
##-----------------------------
library(ComplexHeatmap)
library(circlize)
library(myplot)

oriData <- read.table('GM12878_BS_delta_mCpG_in_allele_specific_peaks.txt', row.names=1)
colnames(oriData) <-c('unbound', 'bound')
oriData$deltaC <- oriData$unbound - oriData$bound
results <- oriData[order(oriData$deltaC, decreasing=T), ]

figRes <- as.matrix(results[, ncol(results)])
rownames(figRes) <- rownames(results)
colFun <- colorRamp2(c(-30, 0, 30), c('#0073c2', '#FFFFF0', '#efc000'))

p <- myplot({
	print(Heatmap(figRes, width=unit(.618, 'in'), height=unit(4.7, 'in'),
		cluster_rows=F, cluster_columns=F, col=colFun,
		column_title=NULL,
		row_names_centered=F, show_row_names=T, row_names_gp=gpar(fontsize=9),
		heatmap_legend_param=list(title=expression(paste('mC' ^ 'unbound', '-mC' ^ 'bound')),
			title_position='leftcenter-rot', legend_height=unit(1.15, 'in'), direction='vertical',
			at=c(-40, 0, 40), labels_gp=gpar(labels_rot=0, fontsize=11))))
		})
plotsave(p, file='fig3.GM12878_BS_TF_delta_mCpG_in_allele_specific_peaks_heatmap.pdf', height=4.71, width=2.05)
