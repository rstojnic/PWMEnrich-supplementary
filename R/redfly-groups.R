library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)

data(PWMLogn.dm3.MotifDb.Dmel)
data(PWMCutoff4.dm3.MotifDb.Dmel)
data(PWMCutoff5.dm3.MotifDb.Dmel)
data(PWMPvalueCutoff1e2.dm3.MotifDb.Dmel)
data(PWMPvalueCutoff1e3.dm3.MotifDb.Dmel)
data(PWMPvalueCutoff1e4.dm3.MotifDb.Dmel)

load("../data/PWMGEV.dm3.MotifDb.Dmel.RData")

load("../data/redfly/crms-v3_2014-10-13-intersectPWMEnrich-by.tf.RData")

useBigMemoryPWMEnrich(TRUE)
registerCoresPWMEnrich(8)

res = list()
res$logn = list()
res$cutoff4 = list()
res$cutoff5 = list()
res$pval1e2 = list()
res$pval1e3 = list()
res$pval1e4 = list()
res$clover = list()

for(i in 1:length(crms.by.tf)){
	cat("Doing", i, "/", length(crms.by.tf), "\n")
	s = sequences[crms.by.tf[[i]]]
	res$logn[[i]] = motifEnrichment(s, PWMLogn.dm3.MotifDb.Dmel)
	res$cutoff4[[i]] = motifEnrichment(s, PWMCutoff4.dm3.MotifDb.Dmel)
	res$cutoff5[[i]] = motifEnrichment(s, PWMCutoff5.dm3.MotifDb.Dmel)
	res$pval1e2[[i]] = motifEnrichment(s, PWMPvalueCutoff1e2.dm3.MotifDb.Dmel)
	res$pval1e3[[i]] = motifEnrichment(s, PWMPvalueCutoff1e3.dm3.MotifDb.Dmel)
	res$pval1e4[[i]] = motifEnrichment(s, PWMPvalueCutoff1e4.dm3.MotifDb.Dmel)
	res$clover[[i]] = motifEnrichment(s, PWMLogn.dm3.MotifDb.Dmel, score="clover")	
}

save(res, file="data-2014/res-redfly-groups.RData")

