
import numpy as np
import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
import harmonypy as hm
from tqdm import tqdm
import celltypist
import scrublet as scr
import warnings
import os
import gc
import json

adata1 = sc.read_10x_h5('sample_filtered_feature_bc_matrix1.h5')
adata2 = sc.read_10x_h5('sample_filtered_feature_bc_matrix2.h5')
adata3 = sc.read_10x_h5('sample_filtered_feature_bc_matrix3.h5')
adata4 = sc.read_10x_h5('sample_filtered_feature_bc_matrix4.h5')


metadata = pd.read_csv('scRNA_metadata.csv', sep=',')

adata1.obs['SampleID'] = '1'
metadata['SampleID'] = metadata['SampleID'].astype(str)
adata1.obs = adata1.obs.reset_index().merge(metadata, on='SampleID').set_index('index')

adata2.obs['SampleID'] = '2'
metadata['SampleID'] = metadata['SampleID'].astype(str)
adata2.obs = adata2.obs.reset_index().merge(metadata, on='SampleID').set_index('index')

adata3.obs['SampleID'] = '3'
metadata['SampleID'] = metadata['SampleID'].astype(str)
adata3.obs = adata3.obs.reset_index().merge(metadata, on='SampleID').set_index('index')

adata4.obs['SampleID'] = '4'
metadata['SampleID'] = metadata['SampleID'].astype(str)
adata4.obs = adata4.obs.reset_index().merge(metadata, on='SampleID').set_index('index')

for adata in [adata1, adata2, adata3, adata4]:
    # Check for duplicates in observations
    if adata.obs.index.duplicated().any():
        adata.obs.index = adata.obs.index + "_" + adata.obs.groupby(level=0).cumcount().astype(str)
    
    # Check for duplicates in variables
    if adata.var.index.duplicated().any():
        adata.var.index = adata.var.index + "_" + adata.var.groupby(level=0).cumcount().astype(str)



combined_obs = ad.concat([adata1, adata2, adata3, adata4], join='outer', axis=0)


combined_vars = ad.concat([adata1, adata2, adata3, adata4], join='outer', axis=1)


combined_adata = ad.AnnData(
    X=combined_obs.X,
    obs=combined_obs.obs,
    var=adata.var
)



# Remove "_0" from each entry in adata.var that we had to add to combine everything 
combined_adata.var.index = adata.var.index.str.replace('_0$', '', regex=True)


# Split by batch (SampleID)
combined_adata.obs['SampleID'] = combined_adata.obs['SampleID'].astype(str)  # Ensure SampleID is string
adata_batches = [combined_adata[combined_adata.obs['SampleID'] == b].copy() for b in combined_adata.obs['SampleID'].unique()]



# Create DataFrame for reporting the cells we are removing from each filtering step
lost_cells_df = pd.DataFrame({'SampleID': [batch.obs['SampleID'][0] for batch in adata_batches]})
lost_cells_df['initial_cells'] = [len(batch) for batch in adata_batches]

# Convert SampleID to categorical and reorder categories
lost_cells_df['SampleID'] = lost_cells_df['SampleID'].astype('category')
lost_cells_df['SampleID'] = lost_cells_df['SampleID'].cat.reorder_categories(
    lost_cells_df.sort_values('initial_cells')['SampleID'].values
)

# Function to plot lost cells
def plot_lost_cells(lost_cells_df):
    lost_cells_df = lost_cells_df.sort_values('initial_cells')
    fig = plt.figure(figsize=(15, 10))
    ax = plt.gca()
    sns.barplot(data=lost_cells_df, x='SampleID', y='initial_cells', ax=ax, color='#bbb')
    cmap = plt.get_cmap('Set1')
    for i, cat in enumerate(lost_cells_df.columns[2:][::-1]):
        sns.barplot(data=lost_cells_df, x='SampleID', y=cat, ax=ax, label=cat, color=cmap(i))
    plt.xticks(rotation=90)
    plt.xlabel('SampleID')
    plt.ylabel('Number of Cells')
    plt.legend()
    plt.show()




warnings.filterwarnings("ignore")
for adata_batch in tqdm(adata_batches):
    # Ensure there are genes and data
    if adata_batch.var_names.size == 0 or adata_batch.X.shape[1] == 0:
        print("Empty data encountered")
        continue
    
    # Annotate mitochondrial and Hbb genes
    adata_batch.var['mt'] = adata_batch.var_names.str.startswith('MT-')
    adata_batch.var['Hbb'] = adata_batch.var_names.str.startswith('HBB')

    # Calculate metrics for mitochondrial and Hbb genes
    sc.pp.calculate_qc_metrics(adata_batch, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.calculate_qc_metrics(adata_batch, qc_vars=['Hbb'], percent_top=None, log1p=False, inplace=True)

warnings.filterwarnings("default")



%matplotlib inline
#plot the mitochondrial reads versus the total number of reads
# without filtering, many warnings about data being converted to categorical
warnings.filterwarnings("ignore")
fig, ax = plt.subplots(1,1,figsize=(15,15))
for adata_batch in adata_batches:
    # plot total number of reads and percent of those reads from mitochondrial genes
    sc.pl.scatter(adata_batch, x='total_counts', y='pct_counts_mt', alpha=0.05, ax=ax, show=False)
warnings.filterwarnings("default")

ax.hlines(5, 0, ax.get_xlim()[-1], color='r', linestyle='dashed')
plt.tight_layout()
plt.show()


# now filter mitochondrial genes and hemoglobin counts

# filter mito counts
lost_cells_df['High MT counts'] = 0
for i in range(len(adata_batches)):
    # kep only observations with low percent of mitochondrial genes
    adata_batches[i] = adata_batches[i][(adata_batches[i].obs.pct_counts_mt < 5)]

    # keep track of total counts after filtering
    batch_name = adata_batches[i].obs.SampleID.values[0]
    lost_cells_df.loc[lost_cells_df.SampleID==batch_name, 'High MT counts'] = lost_cells_df['initial_cells'][lost_cells_df.SampleID==batch_name] - adata_batches[i].n_obs

# filter HBB counts
lost_cells_df['High HBB counts'] = 0
for i in range(len(adata_batches)):

    # kep only observations with low percent of mitochondrial genes and total counts that are not too high
    adata_batches[i] = adata_batches[i][(adata_batches[i].obs.total_counts_Hbb < 2)]

    # view data after filtering
    batch_name = adata_batches[i].obs.SampleID.values[0]
    lost_cells_df.loc[lost_cells_df.SampleID==batch_name, 'High HBB counts'] = lost_cells_df['initial_cells'][lost_cells_df.SampleID==batch_name] - adata_batches[i].n_obs


# filter out cells with a high cell cycle gene content
cell_cycle_genes = ["MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1","HELLS","RFC2","RPA2",
                    "NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1",
                    "BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8","HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
                    "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","JPT1","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1",
                    "KIF20B","HJURP","CDCA3","PIMREG","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA",
                    "PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"]

warnings.filterwarnings("ignore")
for adata_batch in tqdm(adata_batches):
    total_cc_counts = np.zeros_like(adata_batch.obs['SampleID'])
    total_cc_counts = total_cc_counts.astype(np.float64)
    counts = adata_batch.X.toarray()
    for g in cell_cycle_genes:
        try:
            # get index of gene
            index = np.where(adata_batch.var.index.str.upper()==g)[0][0]
            # add total number of reads from gene
            total_cc_counts += counts[:, index]
        except:
            print(f'{g} not found in gene ids')
    # annotate with frac of cell cycle counts
    adata_batch.obs['frac_cc_counts'] = 0
    adata_batch.obs.loc[:, 'frac_cc_counts'] = total_cc_counts/np.sum(counts, axis=1)
warnings.filterwarnings("default")

fig, ax = plt.subplots(1,1,figsize=(15,10))
for adata_batch in adata_batches:
    # plot fraction of genes belonging to cell cycle
    sns.histplot(adata_batch.obs, x='frac_cc_counts', ax=ax, kde=True)
ax.vlines(0.01, 0, ax.get_ylim()[1], color='r', linestyle='dashed')
plt.tight_layout()
plt.show()




lost_cells_df['High CC counts'] = 0

# keep only cells with <1% of counts from CC genes

for i in range(len(adata_batches)):

    # kep only observations with low percent cell cycle genes
    adata_batches[i] = adata_batches[i][(adata_batches[i].obs.frac_cc_counts < 0.01), ]

    # view data after filtering
    batch_name = adata_batches[i].obs.SampleID.values[0]
    lost_cells_df.loc[lost_cells_df.SampleID==batch_name, 'High CC counts'] = lost_cells_df['initial_cells'][lost_cells_df.SampleID==batch_name] - adata_batches[i].n_obs



# filter by num genes detected 
# at least 300 genes with unique count

for i in range(len(adata_batches)):

    adata_batches[i] = adata_batches[i][(adata_batches[i].obs.n_genes_by_counts > 300), ]

    # view data after filtering
    batch_name = adata_batches[i].obs.SampleID.values[0]
    lost_cells_df.loc[lost_cells_df.SampleID==batch_name, 'Low num genes'] = lost_cells_df['initial_cells'][lost_cells_df.SampleID==batch_name] - adata_batches[i].n_obs



fig, ax = plt.subplots(1,1,figsize=(15,10))
for adata_batch in adata_batches:
    # plot distribution of number of unique genes
    sns.histplot(adata_batch.obs, x='n_genes_by_counts', kde=True, ax=ax)
ax.vlines(300, 0, ax.get_ylim()[1], color='r', linestyle='dashed')
plt.tight_layout()
plt.show()



# filter doublets 
for i in tqdm(range(len(adata_batches))):
    # Initialize Scrublet with the raw counts
    scrub = scr.Scrublet(adata_batches[i].X)

    # Predict doublets
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # Add doublet scores and predictions to the AnnData object
    adata_batches[i].obs['doublet_score'] = doublet_scores
    adata_batches[i].obs['predicted_doublet'] = predicted_doublets

    # Set potential doublet to observations with doublet score >=0.2
    adata_batches[i].obs['potential_doublet'] = adata_batches[i].obs['doublet_score'] >= 0.2



fig, ax = plt.subplots(1,1,figsize=(15,10))
for adata_batch in adata_batches:
    # plot distribution of number of unique genes
    sns.histplot(adata_batch.obs, x='doublet_score', kde=True, ax=ax)
ax.vlines(0.2, 0, ax.get_ylim()[1], color='r', linestyle='dashed')
plt.tight_layout()
plt.show()



lost_cells_df['Doublet'] = 0

# remove potential doublets

for i in range(len(adata_batches)):

    adata_batches[i] = adata_batches[i][~adata_batches[i].obs['potential_doublet'], :]

    # view data after filtering
    batch_name = adata_batches[i].obs.SampleID.values[0]
    lost_cells_df.loc[lost_cells_df.SampleID==batch_name, 'Doublet'] = lost_cells_df['initial_cells'][lost_cells_df.SampleID==batch_name] - adata_batches[i].n_obs

# looking at how many cells each sample lost 
plot_lost_cells(lost_cells_df)



#changing the gene ids to be the ensembl id instead of the gene name (this needs to be moved to the gene filtering section)
adata.var["gene_name"] = adata.var.index
adata.var.index = adata.var.gene_ids


# filtering out mitochondrial genes  
# view data before filtering
print("Genes before/after filtering:", adata_batches[0].n_vars, ' -> ', end=' ')
# remove all mitochondrial genes because this is data from the nucleus
for i in range(len(adata_batches)):
    # filter
    adata_batches[i] = adata_batches[i][:, adata_batches[i].var['mt'] != True]

# view data after filtering
print(adata_batches[0].n_vars)




#filter out genes that are very lowly/not expressed, we filter to only include genes that are detected in >1% of cells in any batch

adata_batches[0].var["gene_name"] = adata_batches[0].var.index
adata_batches[0].var.index = adata_batches[0].var.gene_ids

adata_batches[1].var["gene_name"] = adata_batches[1].var.index
adata_batches[1].var.index = adata_batches[1].var.gene_ids

adata_batches[2].var["gene_name"] = adata_batches[2].var.index
adata_batches[2].var.index = adata_batches[2].var.gene_ids

adata_batches[3].var["gene_name"] = adata_batches[3].var.index
adata_batches[3].var.index = adata_batches[3].var.gene_ids

genes_list = set([])
for i in range(len(adata_batches)):
    # we keep genes that we find in >1% of cells in any batch
    batch_genes_mask, _ = sc.pp.filter_genes(adata_batches[i], min_cells=0.01*len(adata_batches[i].obs), inplace=False)

    batch_genes_list = list(adata_batches[i].var[batch_genes_mask].index)
    batch_genes_list.sort()

    genes_list = genes_list | set(batch_genes_list)

# we include all IL genes and all CD genes (because we're interested in immunology/cell communication)
genes_list = genes_list | set([g for g in adata_batches[0].var_names if g.startswith('IL')])
genes_list = genes_list | set([g for g in adata_batches[0].var_names if g.startswith('CD')])

combined_genes_list = list(genes_list)

# sort list of genes
combined_genes_list.sort()

# view data before filtering
print("Genes before/after filtering:", adata_batches[0].n_vars, ' -> ', end=' ')
# filter genes to only include genes in combined genes list
for i in range(len(adata_batches)):

    # filter genes in adata to include only genes in combined_genes_list
    adata_batches[i] = adata_batches[i][:, combined_genes_list]

# view data after filtering genes (should be around 10k)
print(adata_batches[0].n_vars)



# Normalization, Clustering, and UMAP

for adata_batch in adata_batches:
    adata_batch.layers['X_raw'] = adata_batch.X
    sc.pp.normalize_total(adata_batch, target_sum=1e4)
    sc.pp.log1p(adata_batch)
    sc.pp.scale(adata_batch, max_value=10)


# combine the preprocessed adata objects
adata_preprocessed = ad.concat(adata_batches, merge='same')

# view the combined preprocessed data
adata_preprocessed



# delete the batches to save memory
del adata_batches
gc.collect()



# find best number of PCs
sc.tl.pca(adata_preprocessed, n_comps=50, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_preprocessed, 50, log=True)



meta_data = adata_preprocessed.obs
data_mat = adata_preprocessed.obsm['X_pca']

harmony_out = hm.run_harmony(data_mat, meta_data, ['SampleID'])

adata_preprocessed.obsm['X_pca_old'] = adata_preprocessed.obsm['X_pca']
adata_preprocessed.obsm['X_pca'] = np.array(harmony_out.Z_corr.T)

gc.collect()



# create neighbors map
sc.pp.neighbors(adata_preprocessed, n_pcs=50, use_rep='X_pca')

# create leiden and umap
sc.tl.leiden(adata_preprocessed, resolution=1)
sc.tl.umap(adata_preprocessed)


#check and save
adata_preprocessed
adata_preprocessed.write_h5ad('combined_samples_almost_processed_data.h5ad')



##calling cell type

ml GCCcore/.11.3.0 Python/3.10.4
source /scRNA_data_processing/bin/activate 
pip install anndata
pip install scanpy
pip install celltypist

python

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import celltypist


adata = ad.read_h5ad('combined_samples_almost_processed_data.h5ad')

adata.layers['X_processed'] = adata.X
adata.X = adata.layers['X_raw'].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#switch from ensembl id to gene name 
adata.var.index = adata.var.gene_name

#we have an issue w/ duplicate gene names (celltypist is throwing that error) and this was chatgpts sol 
adata.var_names = adata.var_names.astype(str)
adata.var_names_make_unique()

predictions_low = celltypist.annotate(adata, model = 'Immune_All_Low.pkl')
predictions_high = celltypist.annotate(adata, model = 'Immune_All_High.pkl')


adata.X = adata.layers['X_processed'].copy()
del adata.layers['X_processed']

adata.obs['cell_type_low'] = predictions_low.predicted_labels
adata.obs['cell_type_high'] = predictions_high.predicted_labels
adata.obs['cell_type_low_prob'] = np.max(predictions_low.probability_matrix, axis=1)
adata.obs['cell_type_high_prob'] = np.max(predictions_high.probability_matrix, axis=1)

#switch index back to ensembl id bc it has a problem w/ saving for some reason  
adata.var.index = adata.var.gene_ids

#save 
adata.write_h5ad('combined_samples_processed_w_cell_type.h5ad')

exit() 

deactivate







