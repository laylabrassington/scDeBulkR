version

print("loading packages")
.libPaths("/filepath/R_singularity/4.4.2")
library(readr)
library(bseqsc)

setwd('/filepath/bseqsc_job/')

print("loading data")
e_this <- readRDS("e_this_7oct25_model.rds")
B <- readRDS("B_7oct25_model.rds")

print("running deconvolution")
bseq <- bseqsc:::bseqsc_proportions(e_this, B, verbose = TRUE, QN = FALSE)
out <- cbind(pData(e_this), t(coef(bseq)))

print("saving data") 
saveRDS(bseq, "bseq_7oct25.rds")
saveRDS(out, "out_7oct25.rds")

print("done :)")
