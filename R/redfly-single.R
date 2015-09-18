library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)

data(PWMLogn.dm3.MotifDb.Dmel)
data(PWMCutoff4.dm3.MotifDb.Dmel)
data(PWMCutoff5.dm3.MotifDb.Dmel)
data(PWMPvalueCutoff1e2.dm3.MotifDb.Dmel)
data(PWMPvalueCutoff1e3.dm3.MotifDb.Dmel)
data(PWMPvalueCutoff1e4.dm3.MotifDb.Dmel)

load("../data/PWMGEV.dm3.MotifDb.Dmel.RData")

load("../data/redfly/crms-v3_2014-10-13-intersectPWMEnrich.RData")

useBigMemoryPWMEnrich(TRUE)
registerCoresPWMEnrich(8)

res = list()

res$logn = motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)
res$cutoff4 = motifEnrichment(sequences, PWMCutoff4.dm3.MotifDb.Dmel)
res$cutoff5 = motifEnrichment(sequences, PWMCutoff5.dm3.MotifDb.Dmel)
res$pval1e2 = motifEnrichment(sequences, PWMPvalueCutoff1e2.dm3.MotifDb.Dmel)
res$pval1e3 = motifEnrichment(sequences, PWMPvalueCutoff1e3.dm3.MotifDb.Dmel)
res$pval1e4 = motifEnrichment(sequences, PWMPvalueCutoff1e4.dm3.MotifDb.Dmel)
res$gev = motifEnrichment(sequences, PWMGEV.dm3.MotifDb.Dmel)

save(res, file="../data/res-redfly-single.RData")

