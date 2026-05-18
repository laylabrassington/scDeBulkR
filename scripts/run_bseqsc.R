
#install.packages("readr")
#install.packages("parallel")
#install.packages("https://cran.r-project.org/src/contrib/Archive/pkgmaker/pkgmaker_0.32.10.tar.gz", repos = NULL, type = "source")
#install.packages("remotes")
#install.packages("NMF", dependencies = c("Depends", "Imports"))
#remotes::install_github("shenorrLabTRDF/bseqsc") #say no updates

#this code runs on R 4.4.2 

library(readr)
library(bseqsc)
library(parallel)


e_this <- readRDS("e_this_model.rds")
B <- as.matrix(readRDS("B_model.rds"))

sample_chunks <- split(colnames(e_this), 1:detectCores())

res_list <- mclapply(seq_along(sample_chunks), function(i) {
  cols <- sample_chunks[[i]]
  cat("Processing chunk", i, "of", length(sample_chunks), "\n",
      file = "progress_log.txt", append = TRUE)
  flush.console()
  bseqsc:::bseqsc_proportions(e_this[, cols, drop = FALSE], B, QN = FALSE, verbose = TRUE)
}, mc.cores = 8)


saveRDS(res_list, "res_full.rds")

prop_list <- lapply(res_list, coef)
prop_mat <- do.call(cbind, prop_list)
saveRDS(prop_mat, "prop_mat_full.rds")

out <- data.frame(Sample=colnames(prop_mat), t(prop_mat))
saveRDS(out, "out_full.rds")
