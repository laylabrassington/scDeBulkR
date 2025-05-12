# bulk_rnaseq_deconvo


## step 1: making_expr_set_19mar25
This code requires: 
1. A bulk RNA-seq expression matrix file (bulk_data_19Jun24.txt)
2. A metadata file (bulk_data_metadata.txt)
3. A Biobase R package

This code generates:
1. A cleaned expression_matrix containing only samples present in both the metadata and expression data.
2. An AnnotatedDataFrame for the phenotype data (phenoData).
3. e.this which combined expression and phenotype data.
4. e_this.rds that stores the ExpressionSet object 

## step 2: making_sc_markers
This code requires:
1. StdEnv/2020 and python/3.11.2 modules loaded
2. python packages: anndata, scanpy, numpy, pandas, and celltypist
3. A raw single-cell RNA-seq dataset in .h5ad format (raw_sc_data_24jul24.h5ad)

This code generates:
1. A filtered AnnData object (adata_filtered) excluding cell types with only one cell
2. Differential expression results identifying marker genes per cell type using t-test 
3. A .csv file (sc_markers.csv) listing marker genes for each cell type


