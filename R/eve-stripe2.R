library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)

# load the pre-compiled lognormal background
data(PWMLogn.dm3.MotifDb.Dmel)

# load the stripe2 sequences from a FASTA file for motif enrichment
sequence = readDNAStringSet(system.file(package="PWMEnrich", 
  dir="extdata", file="stripe2.fa"))

# perform motif enrichment!
res = motifEnrichment(sequence, PWMLogn.dm3.MotifDb.Dmel)

report = sequenceReport(res, 1)

pdf("figs/eve-stripe2.pdf", width=16, height=12)
plot(report[1:34])
dev.off()

