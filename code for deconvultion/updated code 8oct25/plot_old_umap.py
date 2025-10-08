ml StdEnv/2020 python/3.11.2
conda activate scRNA_data_processing_26mar25
python
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import celltypist
from collections import Counter
from scipy import sparse
from sklearn.cluster import DBSCAN
adata = ad.read_h5ad('/filepath/raw_sc_data_24jul24.h5ad')

import matplotlib.pyplot as plt
plt.switch_backend('Agg') 

sc.settings.autoshow = False
sc.settings.figdir = '/filepath/figures' 

sc.pl.umap(
    adata,
    color='cell_type_low',
    legend_loc='on data',
    size=20,
    legend_fontsize=10,
    legend_fontoutline=2,
    frameon=False,
    show=False)

plt.gcf().set_size_inches(10, 8)
plt.tight_layout()
plt.savefig('/filepath/figures/umap_cell_type_low_big.png', dpi=300)


sc.pl.umap(
    adata,
    color='cell_type_high',
    legend_loc='on data',
    size=20,
    legend_fontsize=10,
    legend_fontoutline=2,
    frameon=False,
    show=False)

plt.gcf().set_size_inches(10, 8)
plt.tight_layout()
plt.savefig('/filepath/figures/umap_cell_type_high_big.png', dpi=300)



