adata = adata_highconf.copy()  # work on a copy

print(f"Starting with {adata.n_obs} cells")

conn = adata.obsp['connectivities']

def neighbor_indices(row):
    return row.indices

from collections import Counter
cell_type_array = adata.obs['cell_type_low'].values
neighbor_labels = []

for i in range(conn.shape[0]):
    neighs = neighbor_indices(conn[i])
    if len(neighs) == 0:
        neighbor_labels.append(None)
    else:
        labels = cell_type_array[neighs]
        majority_label = Counter(labels).most_common(1)[0][0]
        neighbor_labels.append(majority_label)

adata.obs['nn_majority_celltype'] = neighbor_labels
adata.obs['nn_majority_disagree'] = (
    adata.obs['nn_majority_celltype'] != adata.obs['cell_type_low']
)

n_disagree = int(adata.obs['nn_majority_disagree'].sum())
print(f"{n_disagree} cells disagree with kNN majority vote "
      f"({100 * n_disagree / adata.n_obs:.2f}% of total)")

adata_knn_filtered = adata[~adata.obs['nn_majority_disagree']].copy()
print(f"Kept {adata_knn_filtered.n_obs} after kNN majority filter "
      f"(removed {adata.n_obs - adata_knn_filtered.n_obs})")


sc.pl.umap(
    adata_knn_filtered,
    color='cell_type_high',
    legend_loc='on data',
    size=20,
    legend_fontsize=10,
    legend_fontoutline=2,
    show=False)

plt.gcf().set_size_inches(10, 8)
plt.tight_layout()
plt.savefig('/home/brassl1/scRNA/Diane/figures/umap_cell_type_high_knn_filt.png', dpi=300)


sc.pl.umap(
    adata_knn_filtered,
    color='cell_type_low',
    legend_loc='on data',
    size=20,
    legend_fontsize=10,
    legend_fontoutline=2,
    show=False)

plt.gcf().set_size_inches(10, 8)
plt.tight_layout()
plt.savefig('/home/brassl1/scRNA/Diane/figures/umap_cell_type_low_knn_filt.png', dpi=300)


umap = adata_knn_filtered.obsm['X_umap']

# you can tune eps and min_samples
db = DBSCAN(eps=0.2, min_samples=4).fit(umap)
adata_knn_filtered.obs['dbscan'] = db.labels_.astype(str)

n_noise = np.sum(db.labels_ == -1)
print(f" DBSCAN flagged {n_noise} cells as noise "
      f"({100 * n_noise / adata_knn_filtered.n_obs:.2f}% of total)")

# filter out DBSCAN noise points
adata_filtered = adata_knn_filtered[db.labels_ != -1].copy()
print(f" Kept {adata_filtered.n_obs} cells after DBSCAN filter "
      f"(removed {adata_knn_filtered.n_obs - adata_filtered.n_obs})")


sc.pl.umap(
    adata_filtered,
    color='cell_type_high',
    legend_loc='on data',
    size=20,
    legend_fontsize=10,
    legend_fontoutline=2,
    show=False)

plt.gcf().set_size_inches(10, 8)
plt.tight_layout()
plt.savefig('/home/brassl1/scRNA/Diane/figures/umap_cell_type_high_knn&noise_filt.png', dpi=300)

sc.pl.umap(
    adata_filtered,
    color='cell_type_low',
    legend_loc='on data',
    size=20,
    legend_fontsize=10,
    legend_fontoutline=2,
    show=False)

plt.gcf().set_size_inches(15, 10)
plt.tight_layout()
plt.savefig('/home/brassl1/scRNA/Diane/figures/umap_cell_type_low_knn&noise_filt.png', dpi=300)

#SAVE 
adata_filtered.write('/home/brassl1/scRNA/Diane/adata_highconf_filtered_6oct25.h5ad')