setwd('/filepath/')

library(Seurat)
library(zellkonverter)
library(SummarizedExperiment)

adata <- readH5AD("/filepath/filtered_genes_adata_7oct25.h5ad") 
colnames(colData(adata))[colnames(colData(adata)) == "cell_type_low"] <- "cell_type"

# Convert to Seurat object
tmp=adata@assays@data@listData[["X"]]
seurat_obj <- CreateSeuratObject(counts = tmp)

# Extract cell_type column
cell_type <- colData(adata)$cell_type

# Add to Seurat object's metadata
seurat_obj$cell_type <- cell_type

ref_matrix <- AggregateExpression(seurat_obj, group.by = "cell_type", slot = "counts")
B <- ref_matrix$RNA

e_this <- as.matrix(readRDS("e_this_7oct25.rds"))
common_genes <- intersect(rownames(e_this), rownames(B)) 
e_this <- e_this[common_genes, ] 
B <- B[common_genes, ] 
e_this=as.matrix(e_this) 
B=as.matrix(B)

saveRDS(e_this, "e_this_7oct25_model.rds")
saveRDS(B, "B_7oct25_model.rds")
