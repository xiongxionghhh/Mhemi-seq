##----------------------------
# @Author: Xiong X2023/7/9
##----------------------------
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
##-----------------------------
library(myplot)


mC <- read.table('H1_BS_mCpG_Group1_Alu_Start.txt', sep='\t', header=T)
mC$Methylation[which(mC$Strand=='oppo')] <- -1*mC$Methylation[which(mC$Strand=='oppo')]
tmpOppo <- subset(mC, Strand=='oppo', select=-3)
colnames(tmpOppo)[2] <- 'oppo'
tmpSame <- subset(mC, Strand=='same', select=-3)
colnames(tmpSame)[2] <- 'same'
mnase <- read.table('H1_MNase_Group1_Alu_Start.txt', sep='\t', header=T)

comSignal <- merge.data.frame(merge.data.frame(tmpOppo, tmpSame, by='Distance'), mnase, by='Distance')
coMC <- t(as.matrix(comSignal[, c(2, 3)]))

p <- myplot({
	par(mar=c(2.2, 2, 1, 2.2))
	barplot(coMC, beside=TRUE, col=c('#E69F00', '#708090'), border=c('#E69F00', '#708090'), ylim=c(-106, 106),
		axes=F, xpd=F)
	axis(2, at=c(-100,0,100), labels=c(100,0,100), tcl=-.1, lwd=1, lwd.ticks=1, las=0, cex.axis=.8, mgp=c(0, .05, 0))
	mtext('BS-seq mCpG (%)', side=2, cex=0.9, line=1)
	par(new=TRUE)
	plot(comSignal$Distance, comSignal$Signal, type='l', col='#000000', lwd=1, ylim=c(0, 1.2),,
		axes=F, xpd=F)
	axis(1, at=c(-100,0,100,200,300), labels=c(-100,0,100,200,300), tcl=-.1, lwd=1, lwd.ticks=1, las=1, cex.axis=.8, mgp=c(0, 0, 0))
	mtext('Distance from Alu start sites (bp)', side=1, cex=0.9, line=1)
	axis(4, at=c(0,.6,1.2), labels=c(0,.6,1.2), tcl=-.1, lwd=1, lwd.ticks=1, las=1, cex.axis=.8, mgp=c(0, .25, 0))
	mtext('MNase-seq signal', side=4, line=1)
	title(main='Group 1 (H1)', cex.main=.9, mgp=c(1.2, .2, 0))
#	legend(c(300, 1.2), legend=c('oppo', 'same'), ncol=1,
#		pch=15,	cex=.72, bty='n', col=c('#E69F00', '#708090'))
	legend('bottomright', legend='MNase', col='#000000',
		lty=1, lwd=1.2, cex=.72, bty='n')
})
plotsave(p, file='fig6.H1_BS_MNase_Group1_Alu_start_mCpG.pdf', height=2.3, width=3.2, units='in')



mC <- read.table('H1_BS_mCpG_Group2_Alu_Start.txt', sep='\t', header=T)
mC$Methylation[which(mC$Strand=='oppo')] <- -1*mC$Methylation[which(mC$Strand=='oppo')]
tmpOppo <- subset(mC, Strand=='oppo', select=-3)
colnames(tmpOppo)[2] <- 'oppo'
tmpSame <- subset(mC, Strand=='same', select=-3)
colnames(tmpSame)[2] <- 'same'
mnase <- read.table('H1_MNase_Group2_Alu_Start.txt', sep='\t', header=T)

comSignal <- merge.data.frame(merge.data.frame(tmpOppo, tmpSame, by='Distance'), mnase, by='Distance')
coMC <- t(as.matrix(comSignal[, c(2, 3)]))

p <- myplot({
	par(mar=c(2.2, 2, 1, 2.2))
	barplot(coMC, beside=TRUE, col=c('#E69F00', '#708090'), border=c('#E69F00', '#708090'), ylim=c(-106, 106),
		axes=F, xpd=F)
	axis(2, at=c(-100,0,100), labels=c(100,0,100), tcl=-.1, lwd=1, lwd.ticks=1, las=0, cex.axis=.8, mgp=c(0, .05, 0))
	mtext('BS-seq mCpG (%)', side=2, cex=0.9, line=1)
	par(new=TRUE)
	plot(comSignal$Distance, comSignal$Signal, type='l', col='#000000', lwd=1, ylim=c(0, 1.2),,
		axes=F, xpd=F)
	axis(1, at=c(-100,0,100,200,300), labels=c(-100,0,100,200,300), tcl=-.1, lwd=1, lwd.ticks=1, las=1, cex.axis=.8, mgp=c(0, 0, 0))
	mtext('Distance from Alu start sites (bp)', side=1, cex=0.9, line=1)
	axis(4, at=c(0,.6,1.2), labels=c(0,.6,1.2), tcl=-.1, lwd=1, lwd.ticks=1, las=1, cex.axis=.8, mgp=c(0, .25, 0))
	mtext('MNase-seq signal', side=4, line=1)
	title(main='Group 2 (H1)', cex.main=.9, mgp=c(1.2, .2, 0))
#	legend('bottomleft', legend=c('oppo', 'same'),
#		pch=15,	cex=.72, bty='n', col=c('#E69F00', '#708090'))
	legend('bottomright', legend='MNase', col='#000000',
		lty=1, lwd=1.2, cex=.72, bty='n')
})
plotsave(p, file='fig6.H1_BS_MNase_Group2_Alu_start_mCpG.pdf', height=2.3, width=3.2, units='in')


mC <- read.table('H1_BS_mCpG_Group3_Alu_Start.txt', sep='\t', header=T)
mC$Methylation[which(mC$Strand=='oppo')] <- -1*mC$Methylation[which(mC$Strand=='oppo')]
tmpOppo <- subset(mC, Strand=='oppo', select=-3)
colnames(tmpOppo)[2] <- 'oppo'
tmpSame <- subset(mC, Strand=='same', select=-3)
colnames(tmpSame)[2] <- 'same'
mnase <- read.table('H1_MNase_Group3_Alu_Start.txt', sep='\t', header=T)

comSignal <- merge.data.frame(merge.data.frame(tmpOppo, tmpSame, by='Distance'), mnase, by='Distance')
coMC <- t(as.matrix(comSignal[, c(2, 3)]))

p <- myplot({
	par(mar=c(2.2, 2, 1, 2.2))
	barplot(coMC, beside=TRUE, col=c('#E69F00', '#708090'), border=c('#E69F00', '#708090'), ylim=c(-106, 106),
		axes=F, xpd=F)
	axis(2, at=c(-100,0,100), labels=c(100,0,100), tcl=-.1, lwd=1, lwd.ticks=1, las=0, cex.axis=.8, mgp=c(0, .05, 0))
	mtext('BS-seq mCpG (%)', side=2, cex=0.9, line=1)
	par(new=TRUE)
	plot(comSignal$Distance, comSignal$Signal, type='l', col='#000000', lwd=1, ylim=c(0, 1.2),
		axes=F, xpd=F)
	axis(1, at=c(-100,0,100,200,300), labels=c(-100,0,100,200,300), tcl=-.1, lwd=1, lwd.ticks=1, las=1, cex.axis=.8, mgp=c(0, 0, 0))
	mtext('Distance from Alu start sites (bp)', side=1, cex=0.9, line=1)
	axis(4, at=c(0,.6,1.2), labels=c(0,.6,1.2), tcl=-.1, lwd=1, lwd.ticks=1, las=1, cex.axis=.8, mgp=c(0, .25, 0))
	mtext('MNase-seq signal', side=4, line=1)
	title(main='Group 3 (H1)', cex.main=.9, mgp=c(1.2, .2, 0))
#	legend('bottomleft', legend=c('oppo', 'same'),
#		pch=15,	cex=.72, bty='n', col=c('#E69F00', '#708090'))
	legend('bottomright', legend='MNase', col='#000000',
		lty=1, lwd=1.2, cex=.72, bty='n')
})
plotsave(p, file='fig6.H1_BS_MNase_Group3_Alu_start_mCpG.pdf', height=2.3, width=3.2, units='in')



mC <- read.table('H1_BS_mCpG_Group4_Alu_Start.txt', sep='\t', header=T)
mC$Methylation[which(mC$Strand=='oppo')] <- -1*mC$Methylation[which(mC$Strand=='oppo')]
tmpOppo <- subset(mC, Strand=='oppo', select=-3)
colnames(tmpOppo)[2] <- 'oppo'
tmpSame <- subset(mC, Strand=='same', select=-3)
colnames(tmpSame)[2] <- 'same'
mnase <- read.table('H1_MNase_Group4_Alu_Start.txt', sep='\t', header=T)

comSignal <- merge.data.frame(merge.data.frame(tmpOppo, tmpSame, by='Distance'), mnase, by='Distance')
coMC <- t(as.matrix(comSignal[, c(2, 3)]))

p <- myplot({
	par(mar=c(2.2, 2, 1, 2.2))
	barplot(coMC, beside=TRUE, col=c('#E69F00', '#708090'), border=c('#E69F00', '#708090'), ylim=c(-106, 106),
		axes=F, xpd=F)
	axis(2, at=c(-100,0,100), labels=c(100,0,100), tcl=-.1, lwd=1, lwd.ticks=1, las=0, cex.axis=.8, mgp=c(0, .05, 0))
	mtext('BS-seq mCpG (%)', side=2, cex=0.9, line=1)
	par(new=TRUE)
	plot(comSignal$Distance, comSignal$Signal, type='l', col='#000000', lwd=1, ylim=c(0, 1.2),
		axes=F, xpd=F)
	axis(1, at=c(-100,0,100,200,300), labels=c(-100,0,100,200,300), tcl=-.1, lwd=1, lwd.ticks=1, las=1, cex.axis=.8, mgp=c(0, 0, 0))
	mtext('Distance from Alu start sites (bp)', side=1, cex=0.9, line=1)
	axis(4, at=c(0,.6,1.2), labels=c(0,.6,1.2), tcl=-.1, lwd=1, lwd.ticks=1, las=1, cex.axis=.8, mgp=c(0, .25, 0))
	mtext('MNase-seq signal', side=4, line=1)
	title(main='Group 4 (H1)', cex.main=.9, mgp=c(1.2, .2, 0))
#	legend('bottomleft', legend=c('oppo', 'same'),
#		pch=15,	cex=.72, bty='n', col=c('#E69F00', '#708090'))
	legend('bottomright', legend='MNase', col='#000000',
		lty=1, lwd=1.2, cex=.72, bty='n')
})
plotsave(p, file='fig6.H1_BS_MNase_Group4_Alu_start_mCpG.pdf', height=2.3, width=3.2, units='in')
