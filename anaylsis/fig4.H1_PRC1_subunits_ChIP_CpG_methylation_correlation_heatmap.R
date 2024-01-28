##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/30
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(GGally)

singalMatrix <- read.table('H1_PRC1_subunits_ChIP_CpG_methylation_correlation_signal.txt', header=T)

pic <- ggcorr(singalMatrix, hjust=.82, label=T, low='#0073c2', mid='#FFFFFF', high='#008080',
	label_color='#FFD700', label_size=2.5, layout.exp=2, label_round=2)+
	theme(plot.margin=unit(c(-.05, 0, -.05, -.25), 'in'),
		legend.position='none')
ggsave('fig4.H1_PRC1_subunits_ChIP_CpG_methylation_correlation_heatmap.pdf', pic, width=3, height=2.8)


singalMatrix <- read.table('H1_PRC1_subunits_ChIP_BS_Mhemi_promoter_correlation_signal.txt', header=T)

pic <- ggcorr(singalMatrix, hjust=.82, label=T, low='#09224E', mid='#EDF4F6', high='#84001D',
	label_color='#DCDCDC', label_size=2.5, layout.exp=2, label_round=2)+
	theme(plot.margin=unit(c(-.05, 0, -.05, -.25), 'in'),
		legend.position='none')
ggsave('fig4.H1_PRC1_subunits_ChIP_BS_Mhemi_promoter_heatmap.pdf', pic, width=3, height=2.8)



singalMatrix <- read.table('H1_PRC1_subunits_ChIP_Mhemi_promoter_correlation_signal.txt', header=T)

pic <- ggcorr(singalMatrix, hjust=.82, label=T, low='#09224E', mid='#EDF4F6', high='#84001D',
	label_color='#DCDCDC', label_size=2.5, layout.exp=2, label_round=2)+
	theme(plot.margin=unit(c(-.05, 0, -.05, -.25), 'in'),
		legend.position='none')
ggsave('fig4.H1_PRC1_subunits_ChIP_Mhemi_promoter_heatmap.pdf', pic, width=3, height=2.8)
