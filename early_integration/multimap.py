# program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

#   # rna_unit contain mix, ref and ref_scRNA

#   return(rna_unit)
# }
import MultiMAP
import sys
import scanpy as sc
import pandas as pd
import numpy as np
from os import listdir
from rpy2 import robjects
from rpy2.robjects.packages import importr
base = importr("base")

input = NAME_list_mix_ref
output_folder = NAME_output

# load data
# WRITE IT
# D_rna is t(mix_rna) in pandas dataframe
# T_rna is t(ref_bulkRNA)
# D_met is t(mix_met)
# T_met is t(ref_met)

n_sample = D_rna.shape[0]
n_ref = T_rna.shape[0]

data_rna = pd.concat([D_rna, T_rna], ignore_index=True, sort=False)
data_met = pd.concat([D_met, T_met], ignore_index=True, sort=False)

print("data loaded")

# Transfo anndata
adata_rna = sc.AnnData(data_rna, dtype=np.asarray(data_rna).dtype)
adata_met = sc.AnnData(data_met, dtype=np.asarray(data_met).dtype)
adata_rna.obs['source'] = 'RNA'
adata_rna.obs['type'] = ['sample']*n_sample + ['ref']*n_ref
adata_met.obs['source'] = 'MET'
adata_met.obs['type'] = ['sample']*n_sample + ['ref']*n_ref

adata_rna_pca = adata_rna.copy()
sc.pp.scale(adata_rna_pca)
sc.pp.pca(adata_rna_pca)
adata_rna.obsm['X_pca'] = adata_rna_pca.obsm['X_pca'].copy()
sc.pp.highly_variable_genes(adata_met, n_top_genes=20000, subset=True)
adata_met_pca = adata_met.copy()
sc.pp.pca(adata_met_pca)
adata_met.obsm['X_pca'] = adata_met_pca.obsm['X_pca'].copy()
print("data transformed")

# Run MultiMAP
adata = MultiMAP.Integration([adata_rna, adata_met], ['X_pca', 'X_pca'], scale=True, n_components=10)
print("data mapped")
del adata.uns
del adata.obsp

# Save
adata.write_h5ad(output_folder+'adata.h5ad')
