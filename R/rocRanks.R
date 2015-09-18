# roc curve for ranked data

#' Draw a ROC curve for ranks
#'
#' @param ranks a list of expected TFs and their ranks 
#' @param num.motif the total number of motifs
rocRanks = function(ranks, num.motifs){	

	uranks = unlist(ranks)
	max.rank = max(uranks)
	
	roc = matrix(0, ncol=2, nrow=max.rank+2)
	colnames(roc) = c("fpr", "tpr")
	for(r in 1:max.rank){
		tp = sum(uranks <= r)
		fn = sum(uranks > r)
		
		fp = sum(sapply(ranks, function(x) r - sum(x <= r)))
		tn = sum(sapply(ranks, function(x) (num.motifs-r) - sum(x > r)))
		
		stopifnot(sum(tp, fn, fp, tn) == (num.motifs*length(ranks)))
		roc[r+1, 2] = tp / (tp+fn)
		roc[r+1, 1] = fp / (fp+tn)
	}
	roc[max.rank+2,] = 1
	roc
}

#' Draw a PRC curve for ranks
#'
#' @param ranks a list of expected TFs and their ranks 
#' @param num.motif the total number of motifs
prcRanks = function(ranks, num.motifs){	

	uranks = unlist(ranks)
	max.rank = num.motifs
	
	prc = matrix(0, ncol=4, nrow=max.rank)
	colnames(prc) = c("prec", "recall", "tp", "fp")
	for(r in 1:max.rank){
		tp = sum(uranks <= r)
		fn = sum(uranks > r)
		
		fp = sum(sapply(ranks, function(x) r - sum(x <= r)))
		tn = sum(sapply(ranks, function(x) (num.motifs-r) - sum(x > r)))
		
		stopifnot(sum(tp, fn, fp, tn) == (num.motifs*length(ranks)))
		prc[r, "prec"] = tp / (tp+fp)
		prc[r, "recall"] = tp / (tp+fn)
		prc[r, "tp"] = tp
		prc[r, "fp"] = fp
	}
	prc
}

#' Calculate the ROC curve from score
#'
#' @param scores a matrix of scores, where columns are motifs and rows sequences
#' @param expected.tfs a list of expected TF names for each of the sequences
rocScore = function(scores, expected.tfs){
	stopifnot(length(expected.tfs) == nrow(scores))
	
	r = range(scores)	
	p = seq(r[1], r[2], length.out=100)
	
	roc = matrix(0, ncol=2, nrow=length(p))
	colnames(roc) = c("fpr", "tpr")
	motif.names = colnames(scores)
	
	for(i in 1:length(p)){
		tp = 0
		fp = 0
		fn = 0
		tn = 0
		for(j in 1:nrow(scores)){
			pos = motif.names[scores[j,] > p[i]]
			tp = tp + sum(pos %in% expected.tfs[[j]])
			fp = fp + length(pos) - sum(pos %in% expected.tfs[[j]])
			
			neg = motif.names[scores[j,] <= p[i]]
			tn = tn + sum(!(neg %in% expected.tfs[[j]]))
			fn = fn + length(neg) - sum(!(neg %in% expected.tfs[[j]]))
		}
		
		stopifnot(sum(tp, fn, fp, tn) == nrow(scores)*ncol(scores))
		roc[i, "tpr"] = tp / (tp+fn)
		roc[i, "fpr"] = fp / (fp+tn)
	}
	
	roc[nrow(roc):1,]
}

#' Calculate AUC from ROC curve
#'
#' @param roc two column matrix with fpr and tpr rates
calcAuc = function(roc){
	auc = 0
	for(i in 2:nrow(roc)){
		h = (roc[i,2] + roc[i-1,2])/2
		w = roc[i,1] - roc[i-1,1]
		auc = auc + h*w
	}

	names(auc) = "AUC"
	auc
}

#' Plot the ROC curve for ranked data
#'
#' @param ranks a list of expected TFs and their ranks 
#' @param num.motif the total number of motifs
#' @param add if to add to existing plot
plotRocRanks = function(ranks, num.motifs, add=FALSE, xlim=c(0,1), ylim=c(0,1), ...){
	roc = rocRanks(ranks, num.motifs)
	if(add){
		lines(roc, ...)
	} else{
		plot(roc, type="l", xlim=xlim, ylim=ylim, ...)
	}
	abline(0, 1, col="black", lty=2)
	
	calcAuc(roc)
}

#' Plot the ROC curve for ranked data
#'
#' @param ranks a list of expected TFs and their ranks 
#' @param num.motif the total number of motifs
#' @param add if to add to existing plot
plotPrcRanks = function(ranks, num.motifs, add=FALSE, xlim=c(0,1), ylim=c(0,1), prc=NULL, rank.pch=19, rank.cex=0.6, ...){
	if(is.null(prc))
		prc = prcRanks(ranks, num.motifs)
	
	# perform correct interpolation by interpolating on true positive rate
	
	# generate tp vector
	# a list of intervals in increment of one
	seq.list = apply(cbind(prc[1:(nrow(prc)-1),"tp"], prc[2:nrow(prc),"tp"]), 1, function(x) seq(x[1], x[2]))
	# remove the last element unless the list is of length one
	tp = sapply(seq.list, function(x){
		if(length(x) == 1)
			x
		else
			x[-length(x)]
	})
	
	tp = c(unlist(tp), prc[nrow(prc), "tp"]) # we always skip the last one, so include
	fp = rep(0, length(tp))
	# do the interpolation calculation
	for(i in 1:length(tp)){
		# value already present in the table, just use it
		if(tp[i] %in% prc[,"tp"]){
			# there might be multiple values!
			if(length(which(tp == tp[i])) != length(which(prc[,"tp"] %in% tp[i])))
				browser()
			fp[which(tp == tp[i])] = prc[ prc[,"tp"] %in% tp[i] ,"fp"]
		} else {		
			# find interval
			start = max(which(tp[i] > prc[,"tp"])) 
			end = min(which(tp[i] < prc[,"tp"]))
			stopifnot( (end-start) == 1)
			
			skew = (prc[end, "fp"] - prc[start, "fp"]) / (prc[end, "tp"] - prc[start, "tp"])
			fp[i] = prc[start, "fp"] + skew * (tp[i] - prc[start, "tp"])
			
		}
	}
	
	prec = tp / (tp + fp)
	recall = tp / max(tp)
	
	if(add){
		lines(prec, recall, ...)
	} else{
		#colnames(prc) = c("precision", "recall")
		plot(prec, recall, type="l", xlim=xlim, ylim=ylim, xlab="precision", ylab="recall", ...)
	}
	points(prc[,1:2], cex=rank.cex, pch=rank.pch, ...)
	#abline(0, 1, col="black", lty=2)
	
	structure(calcAuc(cbind(sort(prec), recall[order(prec)])), names="AUC-PR")
}


