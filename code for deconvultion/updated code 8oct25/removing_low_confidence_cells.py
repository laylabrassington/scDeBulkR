

#plot low confidence cells 
adata.obs['low_confidence'] = adata.obs['cell_type_low_prob'] < 0.7

sc.pl.umap(
    adata,
    color='low_confidence', 
    size=20)

plt.gcf().set_size_inches(10, 8)
plt.tight_layout()
plt.savefig('/home/brassl1/scRNA/Diane/figures/umap_confidence.png', dpi=300)

# --- 1. Filter out low-confidence cells ---
print(f"Starting with {adata.n_obs} total cells")

# keep only cells where low_confidence == False
adata_highconf = adata[~adata.obs['low_confidence']].copy()

print(f"Retained {adata_highconf.n_obs} high-confidence cells "
      f"({adata_highconf.n_obs / adata.n_obs:.2%} of total)")

# --- 2. Recompute neighbors, clusters, and UMAP ---
# Adjust parameters as needed
sc.pp.neighbors(adata_highconf, n_neighbors=20, n_pcs=50, use_rep='X_pca')

# Leiden clustering (tune resolution)
sc.tl.leiden(adata_highconf, resolution=1.0)

# UMAP embedding
sc.tl.umap(adata_highconf)


# --- 3. plot --- 
sc.pl.umap(
    adata_highconf,
    color='leiden', legend_loc='right margin', size=20)

plt.gcf().set_size_inches(10, 8)
plt.tight_layout()
plt.savefig('/home/brassl1/scRNA/Diane/figures/umap_new_clusters.png', dpi=300)

sc.pl.umap(
    adata_highconf,
    color='cell_type_high',
    legend_loc='on data',
    size=20,
    legend_fontsize=10,
    legend_fontoutline=2,
    show=False)

plt.gcf().set_size_inches(10, 8)
plt.tight_layout()
plt.savefig('/home/brassl1/scRNA/Diane/figures/umap_cell_type_high_new_clusters.png', dpi=300)


sc.pl.umap(
    adata_highconf,
    color='cell_type_low',
    legend_loc='on data',
    size=20,
    legend_fontsize=10,
    legend_fontoutline=2,
    show=False)

plt.gcf().set_size_inches(10, 8)
plt.tight_layout()
plt.savefig('/home/brassl1/scRNA/Diane/figures/umap_cell_type_low_new_clusters.png', dpi=300)


