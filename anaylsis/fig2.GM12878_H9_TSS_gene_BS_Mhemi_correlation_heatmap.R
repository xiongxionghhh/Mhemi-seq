##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/8/10
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(tidyverse)
library(ggtree)
library(aplot)
library(reshape2)


oriData <- read.table('GM12878_BS_Mhemi_TSS_gene_body_300bins.txt', header=T)
corrMat <- cor(oriData, method='spearman')
corrDf <- melt(corrMat)
corrDf$Labels <- round(corrDf$value, digits=2)

p <- ggplot(corrDf, aes(x=Var2, y=Var1))+
	geom_raster(aes(fill=value))+
	geom_text(aes(label=Labels), color='#FFFFFF', size=4)+
	scale_fill_gradient2(low='#4B74B2', high='#DB3124', mid='#FFDF92', midpoint=.5, na.value='black')+
	scale_x_discrete(labels=c('BS-seq', 'Mhemi-seq', 'BS-seq', 'Mhemi-seq'))+
	scale_y_discrete(labels=c('BS-seq', 'Mhemi-seq', 'BS-seq', 'Mhemi-seq'), position='right')+
	theme_classic()+labs(title='', x='', y='', fill='SCC')+
	theme(plot.title=element_blank(),
		plot.margin=unit(c(.01, .12, .35, .01), 'in'),
		axis.line=element_blank(),
		legend.position=c(.35, -1.5),
		legend.direction='vertical',
		legend.title=element_text(size=12),
		legend.text=element_text(size=11),
		legend.key.height=unit(.3, 'in'),
		legend.key.width=unit(.3, 'in'),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_text(size=12, color='black', angle=90, hjust=1),
		axis.text.y=element_text(size=12, color='black'),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank())

featureSet <- data.frame('X'=c('TSS_BS', 'TSS_Mhemi', 'gene_BS', 'gene_Mhemi'),
	'Features'=c('TSS', 'TSS', 'gene body', 'gene body'),
	'Y'=rep('Feature', 4))

pX <- ggplot(featureSet, aes(x=X,y=Y, fill=Features))+
	geom_tile()+
	scale_fill_manual(values=c('#F79324', '#304872'))+
	theme_void()+
	labs(fill='Features')+
	theme(
		legend.position=c(.5, -1.05),
		legend.title=element_text(size=12),
		legend.text=element_text(size=11),
		legend.key.height=unit(.2, 'in'),
		legend.key.width=unit(.2, 'in'))

phY <- corrMat %>%
	dist(method='euclidean')%>%
	hclust(method='complete')%>%
	ggtree(layout='rectangular', branch.length='one')

phX <- corrMat %>%
	t()%>%
	dist(method='euclidean')%>%
	hclust(method='complete')%>%
	ggtree(layout='rectangular', branch.length='none')+
	layout_dendrogram()


p2 <- p %>%
	insert_top(pX, height=.05) %>%
	insert_left(phY,width=.2) %>%
	insert_top(phX,height=.2)

ggsave('fig2.GM12878_TSS_gene_BS_Mhemi_correlation_heatmap.pdf', p2, width=6, height=5)




oriData <- read.table('H1_BS_Mhemi_TSS_gene_body_300bins.txt', header=T)
corrMat <- cor(oriData, method='spearman')
corrDf <- melt(corrMat)
corrDf$Labels <- round(corrDf$value, digits=2)


p <- ggplot(corrDf, aes(x=Var2, y=Var1))+
	geom_raster(aes(fill=value))+
	geom_text(aes(label=Labels), color='#FFFFFF', size=4)+
	scale_fill_gradient2(low='#4B74B2', high='#DB3124', mid='#FFDF92', midpoint=.5, na.value='black')+
	scale_x_discrete(labels=c('BS-seq', 'Mhemi-seq', 'BS-seq', 'Mhemi-seq'))+
	scale_y_discrete(labels=c('BS-seq', 'Mhemi-seq', 'BS-seq', 'Mhemi-seq'), position='right')+
	theme_classic()+labs(title='', x='', y='', fill='SCC')+
	theme(plot.title=element_blank(),
		plot.margin=unit(c(.01, .12, .35, .01), 'in'),
		axis.line=element_blank(),
		legend.position=c(.35, -1.5),
		legend.direction='vertical',
		legend.title=element_text(size=12),
		legend.text=element_text(size=11),
		legend.key.height=unit(.3, 'in'),
		legend.key.width=unit(.3, 'in'),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_text(size=12, color='black', angle=90, hjust=1),
		axis.text.y=element_text(size=12, color='black'),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank())

featureSet <- data.frame('X'=c('TSS_BS', 'TSS_Mhemi', 'gene_BS', 'gene_Mhemi'),
	'Features'=c('TSS', 'TSS', 'gene body', 'gene body'),
	'Y'=rep('Feature', 4))

pX <- ggplot(featureSet, aes(x=X,y=Y, fill=Features))+
	geom_tile()+
	scale_fill_manual(values=c('#F79324', '#304872'))+
	theme_void()+
	labs(fill='Features')+
	theme(
		legend.position=c(.5, -1.05),
		legend.title=element_text(size=12),
		legend.text=element_text(size=11),
		legend.key.height=unit(.2, 'in'),
		legend.key.width=unit(.2, 'in'))

phY <- corrMat %>%
	dist(method='euclidean')%>%
	hclust(method='complete')%>%
	ggtree(layout='rectangular', branch.length='one')

phX <- corrMat %>%
	t()%>%
	dist(method='euclidean')%>%
	hclust(method='complete')%>%
	ggtree(layout='rectangular', branch.length='none')+
	layout_dendrogram()


p2 <- p %>%
	insert_top(pX, height=.05) %>%
	insert_left(phY,width=.2) %>%
	insert_top(phX,height=.2)

ggsave('fig2.H1_TSS_gene_BS_Mhemi_correlation_heatmap.pdf', p2, width=6, height=5)





oriData <- read.table('H9_BS_Mhemi_TSS_gene_body_300bins.txt', header=T)
corrMat <- cor(oriData, method='spearman')
corrDf <- melt(corrMat)
corrDf$Labels <- round(corrDf$value, digits=2)


p <- ggplot(corrDf, aes(x=Var2, y=Var1))+
	geom_raster(aes(fill=value))+
	geom_text(aes(label=Labels), color='#FFFFFF', size=4)+
	scale_fill_gradient2(low='#4B74B2', high='#DB3124', mid='#FFDF92', midpoint=.5, na.value='black')+
	scale_x_discrete(labels=c('BS-seq', 'Mhemi-seq', 'BS-seq', 'Mhemi-seq'))+
	scale_y_discrete(labels=c('BS-seq', 'Mhemi-seq', 'BS-seq', 'Mhemi-seq'), position='right')+
	theme_classic()+labs(title='', x='', y='', fill='SCC')+
	theme(plot.title=element_blank(),
		plot.margin=unit(c(.01, .12, .35, .01), 'in'),
		axis.line=element_blank(),
		legend.position=c(.35, -1.5),
		legend.direction='vertical',
		legend.title=element_text(size=12),
		legend.text=element_text(size=11),
		legend.key.height=unit(.3, 'in'),
		legend.key.width=unit(.3, 'in'),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.text.x=element_text(size=12, color='black', angle=90, hjust=1),
		axis.text.y=element_text(size=12, color='black'),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank())

featureSet <- data.frame('X'=c('TSS_BS', 'TSS_Mhemi', 'gene_BS', 'gene_Mhemi'),
	'Features'=c('TSS', 'TSS', 'gene body', 'gene body'),
	'Y'=rep('Feature', 4))

pX <- ggplot(featureSet, aes(x=X,y=Y, fill=Features))+
	geom_tile()+
	scale_fill_manual(values=c('#F79324', '#304872'))+
	theme_void()+
	labs(fill='Features')+
	theme(
		legend.position=c(.5, -1.05),
		legend.title=element_text(size=12),
		legend.text=element_text(size=11),
		legend.key.height=unit(.2, 'in'),
		legend.key.width=unit(.2, 'in'))

phY <- corrMat %>%
	dist(method='euclidean')%>%
	hclust(method='complete')%>%
	ggtree(layout='rectangular', branch.length='one')

phX <- corrMat %>%
	t()%>%
	dist(method='euclidean')%>%
	hclust(method='complete')%>%
	ggtree(layout='rectangular', branch.length='none')+
	layout_dendrogram()


p2 <- p %>%
	insert_top(pX, height=.05) %>%
	insert_left(phY,width=.2) %>%
	insert_top(phX,height=.2)

ggsave('fig2.H9_TSS_gene_BS_Mhemi_correlation_heatmap.pdf', p2, width=6, height=5)





