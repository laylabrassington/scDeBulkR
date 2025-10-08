ml StdEnv/2020 python/3.11.2
conda activate scRNA_data_processing_26mar25
python

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import celltypist
from collections import Counter

adata = ad.read_h5ad('/filepath/adata_highconf_filtered_6oct25.h5ad')

markers = pd.read_csv("/filepath/unique_genes_7oct25.csv")                       

marker_genes = markers['x'].unique()          
adata.var.index = adata.var['gene_name'] 

adata_filtered = adata[:, adata.var_names.isin(marker_genes)] 

adata_filtered.write_h5ad("/filepath/filtered_genes_adata_7oct25.h5ad")
