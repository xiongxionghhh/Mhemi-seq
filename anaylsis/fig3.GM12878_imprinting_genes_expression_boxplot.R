##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/10
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
# windowsFonts(consolas=windowsFont("consolas"))
##-----------------------------
library(ggplot2)

load('/mnt/disk1/5/share/others/colorset.RData')

exprSet <- read.table('GM12878_imprinting_genes_expression.txt', header=T)

exprGene <- subset(exprSet, Gene=='CRELD2')
pValue <- signif(t.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Creld2')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Creld2_expression_boxplot.pdf', p, width=2, height=2.5)




exprGene <- subset(exprSet, Gene=='ACCS')
pValue <- signif(t.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Accs')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Accs_expression_boxplot.pdf', p, width=2, height=2.5)



exprGene <- subset(exprSet, Gene=='SNRPN')
pValue <- signif(t.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Snrpn')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Snrpn_expression_boxplot.pdf', p, width=2, height=2.5)



exprGene <- subset(exprSet, Gene=='SNURF')
pValue <- signif(t.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Snurf')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Snurf_expression_boxplot.pdf', p, width=2, height=2.5)



exprGene <- subset(exprSet, Gene=='TMSB4X')
pValue <- signif(t.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Tmsb4x')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Tmsb4x_expression_boxplot.pdf', p, width=2, height=2.5)



exprGene <- subset(exprSet, Gene=='TUBB6')
pValue <- signif(t.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Tubb6')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Tubb6_expression_boxplot.pdf', p, width=2, height=2.5)







exprGene <- subset(exprSet, Gene=='CRELD2')
pValue <- signif(wilcox.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Creld2')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Creld2_expression_wilcox_boxplot.pdf', p, width=2, height=2.5)




exprGene <- subset(exprSet, Gene=='ACCS')
pValue <- signif(wilcox.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Accs')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Accs_expression_wilcox_boxplot.pdf', p, width=2, height=2.5)



exprGene <- subset(exprSet, Gene=='SNRPN')
pValue <- signif(wilcox.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Snrpn')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Snrpn_expression_wilcox_boxplot.pdf', p, width=2, height=2.5)



exprGene <- subset(exprSet, Gene=='SNURF')
pValue <- signif(wilcox.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Snurf')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Snurf_expression_wilcox_boxplot.pdf', p, width=2, height=2.5)



exprGene <- subset(exprSet, Gene=='TMSB4X')
pValue <- signif(wilcox.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Tmsb4x')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Tmsb4x_expression_wilcox_boxplot.pdf', p, width=2, height=2.5)



exprGene <- subset(exprSet, Gene=='TUBB6')
pValue <- signif(wilcox.test(exprGene$TPM[which(exprGene$Allele=='Maternal')], exprGene$TPM[which(exprGene$Allele=='Paternal')])$p.value, 2)
anno1 <- ifelse(pValue<=2.2e-16, 'italic(P) < 2.2e-16', paste0('italic(P) == ', pValue))
maxVal <- max(exprGene$TPM)

p <- ggplot(exprGene, aes(x=Allele))+
	geom_boxplot(aes(y=TPM, color=Allele), fill=NA, width=.2, size=.25, outlier.size=.01, notch=FALSE)+
	annotate('text', label=anno1, x=1.5, y=1.15*maxVal, size=3, parse=TRUE)+
	annotate('segment', x=1, xend=2, y=1.1*maxVal, yend=1.1*maxVal, linewidth=.25)+
	theme_classic()+labs(title=expression(italic('Tubb6')), x='', y='TPM')+
	scale_color_manual(values=colorSet[1:length(unique(exprGene$Allele))])+
	scale_y_continuous(limits=c(0, 1.2*maxVal), expand=expansion(mult=c(.05, .05), add=0))+
	theme(plot.margin=unit(c(.03, .03, .03, .03), 'inches'),
		plot.title=element_text(size=11, hjust=.5, vjust=-1),
		panel.grid.major.x=element_blank(),
		panel.grid.major.y=element_blank(),
		panel.border=element_rect(fill=NA, linewidth=.75, color='black'),
		legend.background=element_blank(),
		legend.position='none',
		legend.key=element_blank(),
		legend.key.height=unit(.1, 'inches'),
		legend.title=element_blank(),
		legend.box='',
		legend.text=element_text(size=9),
		legend.spacing.x=unit(.01, 'inches'),
		legend.spacing.y=unit(0, 'inches'),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=10),
		axis.text.x=element_text(size=10, colour='black'),
		axis.text.y=element_text(size=10, colour='black', angle=90, hjust=0.5))
ggsave('fig3.GM12878_imprinting_gene_Tubb6_expression_wilcox_boxplot.pdf', p, width=2, height=2.5)
