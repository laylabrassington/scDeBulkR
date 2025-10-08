setwd('~/scRNA/Diane/')
library(Biobase)
dat = readRDS('17apr25_voomWithQualityWeights_normalized.rds')
e.keep = dat$E
keep.genes <- rownames(e.keep)
meta <- read.delim('18apr25_OA_merged_filt.txt', header = TRUE)
meta$tid <- ifelse(is.na(meta$tid), NA, paste0("X", meta$tid))

setdiff(colnames(e.keep), meta$tid)
setdiff(meta$tid, colnames(e.keep))

meta$tid <- as.character(meta$tid)
meta <- meta[order(match(meta$tid, colnames(e.keep))),]

matching_samples <- intersect(colnames(e.keep), meta$tid)
if (length(matching_samples) == 0) {
  stop("No matching samples between 'e.keep' and 'meta$tid'.")
}

expression_matrix <- e.keep[keep.genes, colnames(e.keep) %in% matching_samples]
if (ncol(expression_matrix) == 0) {
  stop("After subsetting, there are no samples in the expression matrix.")
}

meta_sub <- meta[meta$tid %in% matching_samples, ]
rownames(meta_sub) <- meta_sub$tid
identical(rownames(meta_sub), colnames(expression_matrix))
phenoData <- AnnotatedDataFrame(data = meta_sub)

e.this <- ExpressionSet(assayData = as.matrix(expression_matrix), phenoData = phenoData)

saveRDS(e.this, file = "e_this_7oct25.rds")