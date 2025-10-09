version

print("loading packages")
.libPaths("/filepath/R_singularity/4.4.2")
library(readr)
library(bseqsc)
library(parallel)

setwd('/filepath/bseqsc_job/')

print("loading data")
e_this <- readRDS("e_this_7oct25_model.rds")
B <- readRDS("B_7oct25_model.rds")

sample_chunks <- split(colnames(e_this), 1:detectCores())
#ignore warning 
str(sample_chunks)

print("running deconvolution")
res_list <- mclapply(sample_chunks, function(cols) {
  bseqsc:::bseqsc_proportions(e_this[, cols, drop = FALSE], B, QN = FALSE, verbose = TRUE)
}, mc.cores = 8)  # adjust to available cores

print("combining results")
prop_list <- lapply(res_list, coef)
prop_mat <- do.call(cbind, prop_list)

print("saving outputs")
saveRDS(prop_mat, "bseq_parallel_props_9oct25.rds")

out <- data.frame(Sample = colnames(prop_mat), t(prop_mat))
saveRDS(out, "out_parallel_9oct25.rds")


print("done :)")
