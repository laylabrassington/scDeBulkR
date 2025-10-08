ml StdEnv/2020 python/3.11.2
conda activate scRNA_data_processing_26mar25
python

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import celltypist
from collections import Counter

# Load the single-cell data
adata = ad.read_h5ad('/filepath/adata_highconf_filtered_6oct25.h5ad')


group_counts = Counter(adata.obs['cell_type_low'])
valid_groups = [group for group, count in group_counts.items() if count > 1]
adata_filtered = adata[adata.obs['cell_type_low'].isin(valid_groups)].copy()
#didn't remove any cells 

adata_filtered.var.index = adata_filtered.var.gene_name

# re-normalized
sc.pp.normalize_total(adata_filtered, target_sum=1e4)

# sc.pp.log1p(adata_filtered) 
# got a warning that it was already log transformed 

adata_filtered.raw = adata_filtered

# Identify markers by cell type
# Run rank_genes_groups on the filtered dataset
sc.tl.rank_genes_groups(adata_filtered, 'cell_type_low', method='t-test')

# Get all group names
groups = adata_filtered.obs['cell_type_low'].unique()

# Create a dataframe for all groups
markers_list = [sc.get.rank_genes_groups_df(adata_filtered, group=grp) for grp in groups]
markers = pd.concat(markers_list, keys=groups, names=['group'])

markers.reset_index().to_csv("/filepath/sc_markers_6oct25.csv", index=False)
