##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/30
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ComplexHeatmap)
library(circlize)
library(myplot)

ht_opt$message=FALSE

tfSet <- c('CTCF', 'TRIM22', 'TRIM28', 'ZNF143', 'RAD21')
colorSet <- c('#008080', '#32CD32', '#FFD700', '#BB6C27', '#FFA500')

for(i in 1:5){
	oriData <- read.table(paste0('GM12878_CTCF_motif_', tfSet[i], '_enrichment_classify_by_Mhemi_5th_CpG_dyad_count_1.2_fold.txt'))

	tmpSum <- rowSums(oriData[, -ncol(oriData)], na.rm=T)
	oriData$SumVal <- tmpSum
	oriData <- oriData[which(!is.na(tmpSum)), ]

	groupData <- arrange(oriData, V201, -SumVal)
	groupData$V201 <- gsub('cluster', 'cluster ', groupData$V201)
	groupData$V201 <- factor(groupData$V201, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

	figRes <- as.matrix(subset(groupData, select=-c(V201, SumVal)))
	colnames(figRes) <- c('-2kb', rep('', 198), '2kb')
	colFun <- colorRamp2(quantile(figRes, c(0, 0.985), na.rm=T), c('#FFFFFF', colorSet[i]))

	p <- myplot({
		print(Heatmap(figRes, width=unit(1.1, 'in'), height=unit(3.0, 'in'),
			cluster_rows=F, cluster_columns=F, col=colFun,
			row_split=factor(groupData$V201), row_gap=unit(.05, 'in'), border=TRUE,
			row_title_rot=0, row_title_gp=gpar(fontsize=10),
			column_title=tfSet[i], column_title_gp=gpar(fontsize=13, fontface='bold'),
			show_column_names=T, column_names_rot=0, column_names_centered=T, column_names_gp=gpar(fontsize=11),
			show_row_names=F,
			heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
				legend_height=unit(.7, 'in'), at=c(0, ceiling(quantile(figRes, c(0, 0.985), na.rm=T)[2])))))
	})
	plotsave(p, file=paste0('fig5.GM12878_CTCF_motif_', tfSet[i], '_ChIP_signal_classify_by_Mhemi_5th_CpG_dyad_heatmap.pdf'), height=3.5, width=2.1)



	oriData <- read.table(paste0('GM12878_CTCF_motif_', tfSet[i], '_enrichment_classify_by_Mhemi_7th_CpG_dyad_count_1.2_fold.txt'))

	tmpSum <- rowSums(oriData[, -ncol(oriData)], na.rm=T)
	oriData$SumVal <- tmpSum
	oriData <- oriData[which(!is.na(tmpSum)), ]

	groupData <- arrange(oriData, V201, -SumVal)
	groupData$V201 <- gsub('cluster', 'cluster ', groupData$V201)
	groupData$V201 <- factor(groupData$V201, levels=c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'))

	figRes <- as.matrix(subset(groupData, select=-c(V201, SumVal)))
	colnames(figRes) <- c('-2kb', rep('', 198), '2kb')
	colFun <- colorRamp2(quantile(figRes, c(0, 0.985), na.rm=T), c('#FFFFFF', colorSet[i]))

	p <- myplot({
		print(Heatmap(figRes, width=unit(1.1, 'in'), height=unit(3.0, 'in'),
			cluster_rows=F, cluster_columns=F, col=colFun,
			row_split=factor(groupData$V201), row_gap=unit(.05, 'in'), border=TRUE,
			row_title_rot=0, row_title_gp=gpar(fontsize=10),
			column_title=tfSet[i], column_title_gp=gpar(fontsize=13, fontface='bold'),
			show_column_names=T, column_names_rot=0, column_names_centered=T, column_names_gp=gpar(fontsize=11),
			show_row_names=F,
			heatmap_legend_param=list(title=NULL, title_gp=gpar(fontsize=12), title_position='topcenter',
				legend_height=unit(.7, 'in'), at=c(0, ceiling(quantile(figRes, c(0, 0.985), na.rm=T)[2])))))
	})
	plotsave(p, file=paste0('fig5.GM12878_CTCF_motif_', tfSet[i], '_ChIP_signal_classify_by_Mhemi_7th_CpG_dyad_heatmap.pdf'), height=3.5, width=2.1)
}
