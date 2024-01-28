##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/16
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(myplot)

meth <- read.table('GM12878_hpBS_reproduce_hemi_CpG_dyad_rep1_VS_rep2.txt')
pcc <- paste0('Pearson\'s r=', round(cor(meth$V1, meth$V2), 3))
nCounts <- paste0('n=', sep=gsub(' ', '', format(nrow(meth), big.mark=',')))
p <- myplot({
	par(mar=c(1.2, 1.3, 1.1, .3), lwd=.5)
	smoothScatter(meth, nrpoints=0, xlim=c(0, 100), ylim=c(0, 100), main='',
		xlab='', ylab='', mai=c(0, 0, 0, 0), xaxt='n', yaxt='n', bty='n')
	axis(side=1, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=1, cex.axis=.6, mgp=c(0, -.25, 0))
	axis(side=2, at=seq(0, 100, 25), tcl=-.05, lwd=0, lwd.ticks=.65, las=0, cex.axis=.6, mgp=c(0, .01, 0))
	title(main='hpBS-seq (hemi)', cex.main=.75, font.main=1, mgp=c(.25, 0, 0),
		xlab='Frequency in replicate 1 (%)',
		cex.lab=.6)
	title(mgp=c(.6, 0, 0),
		ylab='Frequency in replicate 2 (%)',
		cex.lab=.6)
	text(50, 3, pcc, cex=.5, col='#708090')
	text(0, 95, nCounts, cex=.5, col='#FF8C00', adj=0)
		})
plotsave(p, file='fig1.GM12878_hpBS_reproduce_hemi_CpG_dyad_rep1_VS_rep2.pdf', height=1.5, width=1.5)

