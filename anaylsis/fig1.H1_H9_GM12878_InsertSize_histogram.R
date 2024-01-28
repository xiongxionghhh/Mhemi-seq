##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/13
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggplot2)

for(celine in c('H1', 'H9', 'GM12878')){
	## in silico
	sizeFreq <- as.numeric(scan(paste0('/mnt/disk1/5/share/DigInSilico/', celine, '-MspJI.dist.txt'), what='character'))
	sizeFreq <- sizeFreq[sizeFreq<=1000]
	coverage <- 3088486552/sum(sizeFreq) # NOTE:119150527 for Arabidopsis; 2725779669 for mm10; 3088486552 for hg38
	tmpTable <- table(sizeFreq)
	inSilicoSizeFreq <- as.data.frame(prop.table(coverage*tmpTable), stringsAsFactors=F)
	inSilicoSizeFreq$Freq <- 100*inSilicoSizeFreq$Freq
	inSilicoSizeFreq$sizeFreq <- as.numeric(inSilicoSizeFreq$sizeFreq)


	sizeFreq <- as.numeric(scan(paste0(celine, '_Mhemi_insert_size.txt'), what='character'))
	sizeFreq <- sizeFreq[sizeFreq<=1000]
	coverage <- 3088486552/sum(sizeFreq) # NOTE:119150527 for Arabidopsis; 2725779669 for mm10; 3088486552 for hg38
	tmpTable <- table(sizeFreq)
	realSizeFreq <- as.data.frame(prop.table(coverage*tmpTable), stringsAsFactors=F)
	realSizeFreq$Freq <- 100*realSizeFreq$Freq
	realSizeFreq$sizeFreq <- as.numeric(realSizeFreq$sizeFreq)

	mergeSizeFreq <- merge(inSilicoSizeFreq, realSizeFreq, by='sizeFreq')
	colnames(mergeSizeFreq) <- c('size', 'insilico', 'real')

	sizeFreq40 <- subset(mergeSizeFreq, size<=1000 & size>40)

	pic1 <- ggplot(sizeFreq40, aes(x=size))+geom_bar(aes(y=real), stat='identity', color='#808080')+
		geom_line(aes(y=insilico), color='#FF8C00')+
		theme_classic()+labs(title=celine, x='Insert size (bp)', y='Percentage')+
		scale_x_continuous(breaks=seq(250, 1000, 250), labels=as.character(seq(250, 1000, 250)))+
		scale_y_continuous(limits=c(0, 1.1), breaks=seq(0, 1, .5), expand=c(0, 0))+
		theme(plot.margin=unit(c(.01, .08, 0, .01), 'inches'),
			plot.title=element_text(colour='black', size=12, hjust=.5),
			legend.title=element_blank(),
			legend.position='none',
			axis.line.x=element_blank(),
			axis.line.y=element_blank(),
			axis.ticks.x=element_blank(),
			axis.ticks.y=element_blank(),
			axis.title.x=element_text(colour='black', size=12),
			axis.title.y=element_text(colour='black', size=12, margin=margin(r=.1, unit='inches')),
			axis.text.x=element_text(colour='black', size=11),
			axis.text.y=element_text(colour='black', size=11))
	ggsave(paste0('fig1.', celine, '_InsertSize_Histogram.pdf'), pic1, width=2.5, height=2.5)

	##---------------------------- 32 bp region
	sizeFreq30 <- subset(mergeSizeFreq, size<=50 & size>=20)

	pic2 <- ggplot(sizeFreq30, aes(x=size))+geom_bar(aes(y=real), stat='identity')+
		geom_line(aes(y=insilico), color='#FF8C00')+
		annotate('segment', x=34, xend=38.5, y=3, yend=3, linewidth=.5, color='#FF8C00')+
		annotate('text', x=45, y=3, label="italic('in silico')", color='#FF8C00', size=4, parse=TRUE)+
		theme_classic()+labs(title='', x='', y='')+
		scale_x_continuous(limits=c(19.9, 50.1), breaks=seq(20, 50, 5), labels=seq(20, 50, 5))+
		scale_y_continuous(limits=c(-.07, 18.07), breaks=seq(0, 18, 6), expand=c(0, 0))+
		theme(plot.margin=unit(c(-.21, 0, -.2, -.18), 'inches'),
			plot.background=element_blank(),
			panel.background=element_blank(),
			legend.position='none',
			axis.line.x=element_blank(),
			axis.line.y=element_line(colour='black', linewidth=.3),
			axis.ticks.x=element_blank(),
			axis.text.x=element_text(colour='black', size=7),
			axis.text.y=element_text(colour='black', size=7))
	ggsave(paste0('fig1.', celine, '_InsertSize_32bp_Histogram.pdf'), pic2, width=1.5, height=1.5)
}

