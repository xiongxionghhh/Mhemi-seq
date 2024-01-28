##----------------------------
# @Author: Xiong Xiong
# @Date: 2021/1/30
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
# windowsFonts(consolas=windowsFont("consolas"))
##-----------------------------
library(ggplot2)
library(ggsci)

motifs <- scan('tmp.SNP.txt', what='character')
motifCount <- as.data.frame(table(motifs))
motifCount$percents <- paste0(as.character(round(100*motifCount$Freq/sum(motifCount$Freq), 2)), '%')
motifCount <- motifCount[order(motifCount$Freq), ]
motifCount$motifs <- factor(motifCount$motifs, levels=motifCount$motifs)

ggplot(motifCount, aes(x='', y=Freq, fill=motifs))+
	geom_bar(stat='identity', position='stack', width=2)+scale_fill_igv()+
	geom_text(aes(x=1.4, label=percents), vjust='outward', angle=20, size=4, stat='identity', position=position_stack(vjust=.5))+
	coord_polar(theta='y', start=0)+
	labs(title='')+
	guides(fill=guide_legend(ncol=1))+
	theme_void()+
	theme(plot.margin=unit(c(-.18, .5, -.2, -.1), 'in'),
		plot.title=element_blank(),
		legend.background=element_blank(),
		legend.key=element_blank(),
		legend.key.size=unit(.16, units='in'),
		legend.text=element_text(colour='black', size=11),
		legend.title=element_blank(),
		legend.position=c(1.05, .5))
ggsave('Fig11.SNP.tiff', width=3.2, height=2.5, compression='lzw')
ggsave('Fig11.SNP.pdf', width=3.2, height=2.5)
