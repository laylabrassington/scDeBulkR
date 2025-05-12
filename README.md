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

