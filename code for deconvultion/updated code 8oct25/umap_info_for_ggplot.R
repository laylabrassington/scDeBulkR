#R 4.5 

library(zellkonverter)
library(reticulate)
library(dplyr)
library(ggplot2)
library(umap)
library(SingleCellExperiment)
library(ggrepel)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(tidyr)
library(colorspace)

adata <- readH5AD("/filepath/adata_highconf_filtered_6oct25.h5ad")

umap_coords <- reducedDims(adata)$X_umap  
umap_df <- as.data.frame(umap_coords)
umap_df$cell_id <- colnames(adata)  

umap_df$cell_type_low <- adata$cell_type_low
umap_df$cell_type_high <- adata$cell_type_high
umap_df$leiden <- adata$leiden
umap_df$doublet_score <- adata$doublet_score
umap_df$total_counts <- adata$total_counts
umap_df$n_genes_by_counts <- adata$n_genes_by_counts
umap_df$age <- adata$age.
umap_df$sex <- adata$sex.
umap_df$sampling_location <- adata$sampling_location
umap_df$sample <- adata$sample
umap_df$h_sol <- adata$h_sol
umap_df$urb_score <- adata$urb_score
umap_df$acculturation_score <- adata$acculturation_score
umap_df$material_wealth <- adata$material_wealth
umap_df$traditional_diet_score <- adata$traditional_diet_score
umap_df$market_diet_index <- adata$market_diet_index

#saveRDS(umap_df, "/filepath/umap_plot_df.rds")
