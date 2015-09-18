library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)
library(GenomicRanges)

data(PWMLogn.dm3.MotifDb.Dmel)

#' Convert 2L:10293651..10293655 into GRanges
coordToGRanges = function(coord){
	p = strsplit(coord, ":|\\.\\.")
	p.chr = sapply(p, function(x) x[1])
	p.start = as.numeric(sapply(p, function(x) x[2]))
	p.end = as.numeric(sapply(p, function(x) x[3]))
	
	GRanges(paste("chr", p.chr, sep=""), IRanges(p.start, p.end))	
	
}

crms = read.csv("../data/redfly/redfly_download_crms-v3_2014-10-13.csv", as.is=TRUE)
# select only sequence between 50 and 8000 bp
crms.len = nchar(crms$sequence)
crms = crms[which(crms.len >= 50 & crms.len <= 8000),]

tfbs = read.csv("../data/redfly/redfly_download_TFBS-v3_2014-10-13.csv", as.is=TRUE)

# add the column with the TF name
tfbs = data.frame(tfbs, tf=sapply(strsplit(tfbs$name, "_"), function(x) x[1]), stringsAsFactors=FALSE)

# intersect the tfbs with CRMs
crms.r = coordToGRanges(crms$coordinates)
tfbs.r = coordToGRanges(tfbs$coordinates)

o = as.data.frame(findOverlaps(crms.r, tfbs.r))

# group by CRM
crms.tfbs = tapply(o$subjectHits, o$queryHits, function(x) tfbs$tf[x])
crms.tfbs = lapply(crms.tfbs, function(x) names(sort(table(x), decreasing=TRUE)))

sel.crms = as.numeric(names(crms.tfbs))
names(crms.tfbs) = crms$name[sel.crms]
crms = crms[sel.crms,]

# and now the dna sequences
sequences = DNAStringSet(crms$sequence)

save(crms.tfbs, sequences, file="../data/redfly/crms-v3_2014-10-13.RData")

# intersect with PWMEnrich motifs
all.targets = sapply(PWMLogn.dm3.MotifDb.Dmel$pwms, function(x) x$name)
keep = sapply(crms.tfbs, function(x) any(x %in% all.targets))

# keep only sequences with 1+ known motifs
sel = which(keep)
crms.tfbs = crms.tfbs[sel]
sequences = sequences[sel]

# filter only those we know about
crms.tfbs = lapply(crms.tfbs, function(x) x[x %in% all.targets])
save(crms.tfbs, sequences, file="../data/redfly/crms-v3_2014-10-13-intersectPWMEnrich.RData")

#############################################
# now group by motif
all.motifs = sort(unique(unlist(crms.tfbs)))
crms.by.tf = lapply(all.motifs, function(x) names(which(sapply(crms.tfbs, function(y) any(x %in% y)))))
names(crms.by.tf) = all.motifs
crms.by.tf.all = lapply(crms.by.tf, function(x) unique(unlist(crms.tfbs[x])))
names(sequences) = names(crms.tfbs)
save(crms.by.tf, crms.by.tf.all, sequences, file="../data/redfly/crms-v3_2014-10-13-intersectPWMEnrich-by.tf.RData")

