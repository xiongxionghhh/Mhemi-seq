##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/8/9
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(ggseqlogo)
library(ggplot2)

tmpSequence <- scan('hg38_Alu_start_sites_5_35nt_denovo.seqs', what='character')
p <- ggseqlogo(tmpSequence[1:10000], method='bits')+
	geom_point(aes(x=12, y=1.8), shape=8)+
	scale_y_continuous(limits=c(0, 2))+
	labs(title='P1 (E-value = 9.7e-5859)')+
	theme_classic()+
	theme(plot.margin=unit(c(.01, 0, -.2, 0), 'in'),
		plot.title=element_text(size=14, hjust=.5),
		panel.border=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.y=element_text(size=12),
		axis.text.x=element_blank(),
		axis.text.y=element_text(size=11, colour='black'))
ggsave('fig6.hg38_Alu_start_site_5_35_denovo_motif_sequence_logo.pdf', p, height=2.2, width=4.8)


tmpSequence <- scan('hg38_Alu_start_sites_140_170nt_denovo.seqs', what='character')
p <- ggseqlogo(tmpSequence[1:10000], method='bits')+
	geom_point(aes(x=11, y=1.8), shape=8)+
	scale_y_continuous(limits=c(0, 2))+
	labs(title='P2 (E-value = 8.4e-3222)')+
	theme_classic()+
	theme(plot.margin=unit(c(.01, 0, -.2, 0), 'in'),
		plot.title=element_text(size=14, hjust=.5),
		panel.border=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.title.y=element_text(size=12),
		axis.text.x=element_blank(),
		axis.text.y=element_text(size=11, colour='black'))
ggsave('fig6.hg38_Alu_start_sites_140_170nt_denovo_motif_sequence_logo.pdf', p, height=2.2, width=4.8)
