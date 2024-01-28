##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/13
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggseqlogo)
library(ggplot2)

for(celine in c('H1', 'H9', 'GM12878')){
	tmpSequence <- scan(paste0(celine, '_32bp_sequence.txt'), what='character')
	pic <- ggseqlogo(tmpSequence, method='probability')+
		labs(title=celine)+
		scale_x_continuous(breaks=c(1, 14, 19,  32))+
		scale_y_continuous(limits=c(0, 1))+
		theme_classic()+
		theme(plot.margin=unit(c(.03, -.17, -.18, .02), 'in'),
			plot.title=element_text(size=12, hjust=.5),
			panel.border=element_blank(),
			axis.line.x=element_blank(),
			axis.line.y=element_blank(),
			axis.ticks.x=element_blank(),
			axis.ticks.y=element_blank(),
			axis.title.y=element_text(size=12),
			axis.text.x=element_text(size=11, colour='black'),
			axis.text.y=element_text(size=11, colour='black'))
	ggsave(paste0('fig1.', celine, '_32bp_sequence_logo_prob.pdf'), pic, height=2.2, width=5)
	
	pic2 <- ggseqlogo(tmpSequence, method='bits')+
		labs(title=celine)+
		scale_x_continuous(breaks=c(1, 14, 19,  32))+
		scale_y_continuous(limits=c(0, 1))+
		theme_classic()+
		theme(plot.margin=unit(c(.03, -.17, -.18, .02), 'in'),
			plot.title=element_text(size=12, hjust=.5),
			panel.border=element_blank(),
			axis.line.x=element_blank(),
			axis.line.y=element_blank(),
			axis.ticks.x=element_blank(),
			axis.ticks.y=element_blank(),
			axis.title.y=element_text(size=12),
			axis.text.x=element_text(size=11, colour='black'),
			axis.text.y=element_text(size=11, colour='black'))
	ggsave(paste0('fig1.', celine, '_32bp_sequence_logo.pdf'), pic2, height=2.2, width=5)
}

