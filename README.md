# scDeBulkR


## order of scripts for updated code 8jan26
- starts with the adata object from extra_cell_filtering.py from 8oct25
1. scRNA_cell_type_markers.Rmd




## order of scripts for updated code 8oct25
1. making_e_this.R
2. plot_old_umap.py
3. removing_low_confidence_cells.py
4. extra_cell_filtering.py
5. sc_markers.py
6. getting_gene_lists.R
7. filtering_sc_data.py
8. prep_data_for_deconvo.R
9. bseqsc_job_code.R
10. bseqsc_job_submission.sh





## step 1: make_bulk_expr_set
This code requires: 
1. bulk RNA-seq expression matrix file (bulk_data_19Jun24.txt)
2. metadata file (bulk_data_metadata.txt)
3. Biobase R package

This code generates:
1. cleaned expression_matrix containing only samples present in both the metadata and expression data
2. AnnotatedDataFrame for the phenotype data (phenoData)
3. e.this which combined expression and phenotype data
4. e_this.rds that stores the ExpressionSet object 

## step 2: making_sc_markers
This code requires:
1. StdEnv/2020 and python/3.11.2 modules loaded
2. python packages: anndata, scanpy, numpy, pandas, and celltypist
3. raw single-cell RNA-seq dataset in .h5ad format (raw_sc_data_24jul24.h5ad)

This code generates:
1. filtered AnnData object (adata_filtered) excluding cell types with only one cell
2. Differential expression results identifying marker genes per cell type using t-test 
3. A .csv file (sc_markers.csv) listing marker genes for each cell type

## step 3: matching_bulk_sc_genes
This code requires:
1. marker gene list from single-cell data (sc_markers_26mar25.csv)
2. bulk RNA-seq ExpressionSet object (e_this.rds)
3. dplyr R package
   
This code generates:
1. filtered list of genes present in both single-cell and bulk RNA-seq datasets
2. list of marker genes per cell type (filtered by adjusted p-value threshold)
3. markers.rds and markers.csv containing the cell-type-specific marker genes
4. unique_genes.csv containing all unique marker genes across all cell types

## step 4: filtering_sc_data_genes
This code requires:
1. python packages: anndata, scanpy, celltypist, numpy, and pandas
2. single-cell RNA-seq dataset (raw_sc_data_24jul24.h5ad)
3. CSV file of unique marker genes (unique_genes.csv)

This code generates:
1. filtered AnnData object (adata_filtered) containing only genes matching the marker list
2. filtered_adata.h5ad which contains the filtered gene expression data

## step 5: run_deconvolution
This code requires:
1. R packages: Seurat, zellkonverter, SummarizedExperiment, and bseqsc.
2. filtered_adata.h5ad and e_this.rds

This code generates:
1. normalized Seurat object (seurat_obj) with cell type metadata
2. reference matrix (B) of average expression by cell type
3. matched gene expression matrices from bulk and single-cell data
4. estimated cell type proportions in bulk samples (out)
