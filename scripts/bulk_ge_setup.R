

library(limma)
library(Biobase)

# Load your normalized voom matrix
v <- readRDS('voomWithQualityWeights_normalized.rds')$E
meta <- readRDS('meta_data.rds') 

# Ensure samples match
v <- v[, meta$sample_name]
stopifnot(identical(meta$sample_name, colnames(v)))

# Remove batch effects while preserving overall expression 
v_corrected <- removeBatchEffect(v, 
                                 batch = as.factor(meta$batch), 
                                 covariates = meta$Uniquely.mapped.reads)


# Build ExpressionSet for downstream bseqsc
keep.genes <- rownames(v_corrected)
meta$sample_name <- as.character(meta$sample_name)
meta <- meta[order(match(meta$sample_name, colnames(v_corrected))),]

matching_samples <- intersect(colnames(v_corrected), meta$sample_name)
expression_matrix <- v_corrected[keep.genes, colnames(v_corrected) %in% matching_samples]
meta_sub <- meta[meta$sample_name %in% matching_samples, ]
rownames(meta_sub) <- meta_sub$sample_name

phenoData <- AnnotatedDataFrame(data = meta_sub)
e.this <- ExpressionSet(assayData = as.matrix(expression_matrix), phenoData = phenoData)

saveRDS(e.this, file = "e_this_batch_corrected.rds")