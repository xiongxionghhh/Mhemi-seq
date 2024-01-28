##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/16
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(myplot)

meth <- read.table('GM12878_BS_Mhemi_mCGNR_genes_1kb_500bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='gene body (GM12878)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.GM12878_BS_Mhemi_mCGNR_genes_1kb_500bins_density.pdf', height=1.5, width=1.5)

meth <- read.table('H1_BS_Mhemi_mCGNR_genes_1kb_500bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='gene body (H1)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.H1_BS_Mhemi_mCGNR_genes_1kb_500bins_density.pdf', height=1.5, width=1.5)

meth <- read.table('H9_BS_Mhemi_mCGNR_genes_1kb_500bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='gene body (H9)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.H9_BS_Mhemi_mCGNR_genes_1kb_500bins_density.pdf', height=1.5, width=1.5)







meth <- read.table('GM12878_BS_Mhemi_mCGNR_CGI_200bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='CGI (GM12878)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.GM12878_BS_Mhemi_mCGNR_CGI_200bins_density.pdf', height=1.5, width=1.5)

meth <- read.table('H1_BS_Mhemi_mCGNR_CGI_200bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='CGI (H1)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.H1_BS_Mhemi_mCGNR_CGI_200bins_density.pdf', height=1.5, width=1.5)

meth <- read.table('H9_BS_Mhemi_mCGNR_CGI_200bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='CGI (H9)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.H9_BS_Mhemi_mCGNR_CGI_200bins_density.pdf', height=1.5, width=1.5)






meth <- read.table('GM12878_BS_Mhemi_mCGNR_TSS_200bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='TSS (GM12878)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.GM12878_BS_Mhemi_mCGNR_TSS_200bins_density.pdf', height=1.5, width=1.5)

meth <- read.table('H1_BS_Mhemi_mCGNR_TSS_200bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='TSS (H1)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.H1_BS_Mhemi_mCGNR_TSS_200bins_density.pdf', height=1.5, width=1.5)

meth <- read.table('H9_BS_Mhemi_mCGNR_TSS_200bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='TSS (H9)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.H9_BS_Mhemi_mCGNR_TSS_200bins_density.pdf', height=1.5, width=1.5)












meth <- read.table('GM12878_BS_Mhemi_mCGNR_gene_body_300bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='gene body (GM12878)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.GM12878_BS_Mhemi_mCGNR_gene_body_300bins_density.pdf', height=1.5, width=1.5)

meth <- read.table('H1_BS_Mhemi_mCGNR_gene_body_300bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='gene body (H1)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.H1_BS_Mhemi_mCGNR_gene_body_300bins_density.pdf', height=1.5, width=1.5)

meth <- read.table('H9_BS_Mhemi_mCGNR_gene_body_300bins.txt', header=T)
pcc <- paste0('Pearson\'s r=', round(cor(meth$BS, meth$Mhemi), 3))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='gene body (H9)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Methylation in BS-seq (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Methylation in Mhemi-seq (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
		})
plotsave(p, file='fig2.H9_BS_Mhemi_mCGNR_gene_body_300bins_density.pdf', height=1.5, width=1.5)
