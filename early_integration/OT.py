# program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

#   # rna_unit contain mix, ref and ref_scRNA

#   return(rna_unit)
# }
# https://github.com/caokai1073/uniPort
# https://uniport.readthedocs.io/en/latest/
# pip3 install uniport

import warnings
import uniport as up
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import numpy as np
import sys
from os import listdir
from rpy2 import robjects
from rpy2.robjects.packages import importr
base = importr("base")

input = NAME_list_mix_ref
output_folder = NAME_output

# load data
# WRITE IT
# D_rna is mix_rna in pandas dataframe
# T_rna is ref_bulkRNA
# D_met is mix_met
# T_met is ref_met

# load data
data_rna_DT = pd.concat([D_rna, T_rna], ignore_index=True, sort=False, axis=1)
data_rna_D = D_rna

data_met_DT = pd.concat([D_met, T_met], ignore_index=True, sort=False, axis=1)
data_met_D = D_met
print("data loaded")

# Make anndata
n_sample = data_rna_D.shape[1]
n_ref = data_rna_DT.shape[1] - n_sample
adata_rna_DT = sc.AnnData(data_rna_DT.T, dtype=np.float32)
adata_met_DT = sc.AnnData(data_met_DT.T, dtype=np.float32)
adata_rna_DT.obs['source'] = 'RNA'
adata_rna_DT.obs['type'] = ['sample']*n_sample + ['ref']*n_ref
adata_met_DT.obs['source'] = 'MET'
adata_met_DT.obs['type'] = ['sample']*n_sample + ['ref']*n_ref

# Select 2k HVG in each
sc.pp.normalize_total(adata_rna_DT)
sc.pp.log1p(adata_rna_DT)
sc.pp.highly_variable_genes(adata_rna_DT, n_top_genes=2000, inplace=False, subset=True)
sc.pp.highly_variable_genes(adata_met_DT, n_top_genes=2000, inplace=False, subset=True)
sc.pp.normalize_total(adata_met_DT)
sc.pp.log1p(adata_met_DT)
adata_rna_DT.X = adata_rna_DT.X/adata_rna_DT.X.max()
adata_met_DT.X = adata_met_DT.X/adata_met_DT.X.max()
adata_rna_DT.X = csr_matrix(adata_rna_DT.X, dtype=np.float32)
adata_met_DT.X = csr_matrix(adata_met_DT.X, dtype=np.float32)

warnings.filterwarnings("ignore", message="Cannot set number of intraop threads after parallel work has started")

# OT MET->RNA
adata_rna_DT_OT = adata_rna_DT
print("Analyzing dataset")
adata_rna_DT_OT = up.Run(adatas=[adata_rna_DT, adata_met_DT], mode='v', iteration=500, batch_size=50)

# OT RNA->MET
adata_met_DT_OT = adata_met_DT
adata_met_DT_OT = up.Run(adatas=[adata_met_DT, adata_rna_DT], mode='v', iteration=500, batch_size=50)

# Compute UMAP w/ or w/o OT
sc.pp.pca(adata_rna_DT_OT)
sc.pp.neighbors(adata_rna_DT_OT, use_rep='X_pca', key_added='withoutOT')
sc.tl.umap(adata_rna_DT_OT, min_dist=0.1, neighbors_key='withoutOT', n_components=10)
adata_rna_DT_OT.obsm['umap_withoutOT'] = adata_rna_DT_OT.obsm['X_umap']
sc.pp.neighbors(adata_rna_DT_OT, use_rep='latent', key_added='withOT')
sc.tl.umap(adata_rna_DT_OT, min_dist=0.1, neighbors_key='withOT', n_components=10)
adata_rna_DT_OT.obsm['umap_withOT'] = adata_rna_DT_OT.obsm['X_umap']

sc.pp.pca(adata_met_DT_OT)
sc.pp.neighbors(adata_met_DT_OT, use_rep='X_pca', key_added='withoutOT')
sc.tl.umap(adata_met_DT_OT, min_dist=0.1, neighbors_key='withoutOT', n_components=10)
adata_met_DT_OT.obsm['umap_withoutOT'] = adata_met_DT_OT.obsm['X_umap']
sc.pp.neighbors(adata_met_DT_OT, use_rep='latent', key_added='withOT')
sc.tl.umap(adata_met_DT_OT, min_dist=0.1, neighbors_key='withOT', n_components=10)
adata_met_DT_OT.obsm['umap_withOT'] = adata_met_DT_OT.obsm['X_umap']

# Save UMAP coordinates w/ or w/o OT
np.savetxt(output_folder+'metUMAPwithOT.csv', adata_met_DT_OT.obsm['umap_withOT'],
           delimiter=',', fmt='%.2f')
np.savetxt(output_folder+'rnaUMAPwithOT.csv', adata_rna_DT_OT.obsm['umap_withOT'],
           delimiter=',', fmt='%.2f')
    
