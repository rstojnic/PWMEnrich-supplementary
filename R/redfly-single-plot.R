library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)
library(RColorBrewer)

data(PWMLogn.dm3.MotifDb.Dmel)

source("rocRanks.R")

load("../data/redfly/crms-v3_2014-10-13-intersectPWMEnrich.RData")
load("../data/res-redfly-single.RData")

calcRanks = function(res, expected.tfs, bg){
	ranks = list()
	for(i in 1:nrow(res$sequence.nobg)){
		r = sequenceReport(res, i, bg=bg)
		targets = r$target
		rr = r$rank
		
		ranks[[i]] = rr[as.vector(na.omit(match(expected.tfs[[i]], targets)))]
	}

	# remove zero-length ranks
	ranks.len = sapply(ranks, length)
	ranks[ranks.len != 0]
}

num.motifs = ncol(res$logn$sequence.nobg)

# I should first plot all of them, and then decide on a subset that perform the best

# everything background
ranks = list()
ranks$logn = calcRanks(res$logn, crms.tfbs, TRUE)
ranks$pval1e2 = calcRanks(res$pval1e2, crms.tfbs, TRUE)
ranks$pval1e3 = calcRanks(res$pval1e3, crms.tfbs, TRUE)
ranks$pval1e4 = calcRanks(res$pval1e4, crms.tfbs, TRUE)
ranks$cutoff4 = calcRanks(res$cutoff4, crms.tfbs, TRUE)
ranks$cutoff5 = calcRanks(res$cutoff5, crms.tfbs, TRUE)
ranks$gev = calcRanks(res$gev, crms.tfbs, TRUE)

# everything without
ranks$affinity = calcRanks(res$logn, crms.tfbs, FALSE)
ranks$count1e2 = calcRanks(res$pval1e2, crms.tfbs, FALSE)
ranks$count1e3 = calcRanks(res$pval1e3, crms.tfbs, FALSE)
ranks$count1e4 = calcRanks(res$pval1e4, crms.tfbs, FALSE)
ranks$count4 = calcRanks(res$cutoff4, crms.tfbs, FALSE)
ranks$count5 = calcRanks(res$cutoff5, crms.tfbs, FALSE)

# plot everything 
pdf("figs/redfly-single.pdf", width=12, height=6)
par(mfrow=c(1,2))
col = brewer.pal(9, "Set1")

lwd = 2
# without BG
plotPrcRanks(ranks$count5, num.motifs, xlim=c(0, 0.2), col=col[5], rank.pch=6, lwd=lwd, main="Without background correction")
plotPrcRanks(ranks$count4, num.motifs, col=col[4], add=TRUE, rank.pch=5, lwd=lwd)
plotPrcRanks(ranks$count1e4, num.motifs, col=col[3], add=TRUE, rank.pch=4, lwd=lwd)
plotPrcRanks(ranks$count1e3, num.motifs, col=col[2], add=TRUE, rank.pch=3, lwd=lwd)
plotPrcRanks(ranks$count1e2, num.motifs, col=col[1], add=TRUE, rank.pch=2, lwd=lwd)
plotPrcRanks(ranks$affinity, num.motifs, col="black", add=TRUE, rank.pch=1, lwd=lwd)
legend("topright", lty=rep(1, 6), col=c("black", col[1:5]), lwd=2, pch=1:6,
	legend=c("Affinity", "P-value=1e-2", "P-value=1e-3", "P-value=1e-4",
	         "Cuttof=4", "Cutoff=5"))

# with BG
plotPrcRanks(ranks$gev, num.motifs, col="grey", xlim=c(0, 0.2), rank.pch=7, lwd=lwd, main="With background correction")
plotPrcRanks(ranks$cutoff5, num.motifs, col=col[5], add=TRUE, rank.pch=6, lwd=lwd)
plotPrcRanks(ranks$cutoff4, num.motifs, col=col[4], add=TRUE, rank.pch=5, lwd=lwd)
plotPrcRanks(ranks$pval1e4, num.motifs, col=col[3], add=TRUE, rank.pch=4, lwd=lwd)
plotPrcRanks(ranks$pval1e3, num.motifs, col=col[2], add=TRUE, rank.pch=3, lwd=lwd)
plotPrcRanks(ranks$pval1e2, num.motifs, col=col[1], add=TRUE, rank.pch=2, lwd=lwd)
plotPrcRanks(ranks$logn, num.motifs, col="black", add=TRUE, rank.pch=1, lwd=lwd)
legend("topright", lty=rep(1, 7), col=c("black", col[1:5], "grey"), lwd=2, pch=1:7,
	legend=c("Lognormal", "P-value=1e-2", "P-value=1e-3", "P-value=1e-4",
	         "Cuttof=4", "Cutoff=5", "GEV"))
dev.off()
