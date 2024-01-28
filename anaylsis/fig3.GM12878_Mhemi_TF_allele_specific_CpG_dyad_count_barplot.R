##----------------------------
# @Author: Xiong Xiong
# @Date: 2023/7/25
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(myplot)

dyadAS <- as.matrix(read.table('GM12878_Mhemi_TF_allele_specific_CpG_dyad_count.txt', header=T, row.names=1))

p <- myplot({
par(mar=c(1.2, 1.2, 1.5, 1.2))
y <- barplot(dyadAS, beside=T, border=NA, xlim=c(0, 300), width=.5, xaxp=c(0, 300, 3), col=c('#7FC87F', '#fc9533'),
	xlab='', ylab='', tcl=-.15, cex.axis=.9, cex.names=.9, horiz=T, mgp=c(1.2, .2, 0))
legend('topright', rownames(dyadAS), pch=15, pt.cex=1.2,
	bty='n', col=c('#7FC87F', '#fc9533'))
title(main='Number of allele-specific dyad', font.main=1, cex.main=1, mgp=c(1.2, .2, 0),
	xlab='', ylab='', cex.lab=.9)
text(dyadAS+1, y, labels=as.character(dyadAS), cex=.8, adj=0)
	})
plotsave(p, file='fig3.GM12878_Mhemi_TF_allele_specific_CpG_dyad_count_barplot.pdf', height=2.5, width=2.5)
