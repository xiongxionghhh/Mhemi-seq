##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/8/10
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(myplot)

ht_opt$message = FALSE

tmpSame <- read.table('H1_BS_Alu_C1_C2_same_strand.txt', row.names=1)
tmpOppo <- read.table('H1_BS_Alu_C1_C2_oppo_strand.txt', row.names=1)
deltaC <- tmpOppo-tmpSame
mGroup <- ifelse(deltaC$V7 > 10 & deltaC$V16 > 10, 'Group1', ifelse(deltaC$V7 > 10 & deltaC$V16 <= 2, 'Group2', ifelse(deltaC$V7 <= 2 & deltaC$V16 > 10, 'Group3', ifelse(deltaC$V7 <= 2 & deltaC$V16 <= 2, 'Group4', 'Others'))))
deltaC$mGroup <- mGroup
write.table(subset(deltaC, mGroup!='Others', select=mGroup), 'H1_Alu_group_classify_by_C1_C2_delta_mCWG_oppo-same.txt', col.names=F, sep='\t', quote=F)

groupDeltaC <- subset(deltaC, mGroup!='Others')
groupDeltaC$mGroup <- factor(groupDeltaC$mGroup, levels=c('Group1', 'Group2', 'Group3', 'Group4'))

groupData <- arrange(groupDeltaC, mGroup, -V7, -V16) 

figRes <- as.matrix(subset(groupData, select=-mGroup))
colnames(figRes) <- c(rep('', 5), 'C1', rep('', 8), 'C2', rep('', 5))
colFun <- colorRamp2(c(-50, 0, 50), c('#0073c2', '#FFFFF0', '#efc000')) # quantile(figRes, c(.001, .35, .999), na.rm=T)

p <- myplot({
	print(Heatmap(figRes, width=unit(1.9, 'in'), height=unit(3.02, 'in'), use_raster=TRUE,
		cluster_rows=F, cluster_columns=F, col=colFun, na_col='#DCDCDC',
		row_split=factor(groupData$mGroup), row_gap=unit(.05, 'in'),
		row_title_rot=0, row_title_gp=gpar(fontsize=12),
		column_title='Methylation bias @ Alu', column_title_gp=gpar(fontsize=13),
		column_names_rot=0, column_names_centered=T, show_row_names=F, column_names_gp=gpar(fontsize=12),
		heatmap_legend_param=list(title=expression(paste('mCHG' ^ 'oppo', '-mCHG' ^ 'same')),
			title_position='leftcenter-rot', legend_height=unit(1.5, 'in'), direction='vertical',
			title_gp=gpar(fontsize=11), at=c(-100, 0, 100))))
	})
plotsave(p, file='fig6.H1_Alu_C1_C2_delta_mCWG_oppo-same_heatmap.pdf', height=3.5, width=3.5)
