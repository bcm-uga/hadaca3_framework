# program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

#   # rna_unit contain mix, ref and ref_scRNA

#   return(rna_unit)
# }
# import sys
# from os import listdir
# from rpy2 import robjects
# from rpy2.robjects.packages import importr
# base = importr("base")

import MultiMAP
import scanpy as sc
import pandas as pd
import numpy as np


def program_block_EI(rna_unit,met_unit,path_og_dataset=""):

    # Align by shared genes
    common_genes    = rna_unit["mix"].index.intersection(rna_unit["ref"].index)
    rna_unit["mix"] = rna_unit["mix"].loc[common_genes]
    rna_unit["ref"] = rna_unit["ref"].loc[common_genes]


    common_genes_met = met_unit["mix"].index.intersection(met_unit["ref"].index)
    met_unit["mix"]  = met_unit["mix"].loc[common_genes_met]
    met_unit["ref"]  = met_unit["ref"].loc[common_genes_met]

    # D_rna is t(mix_rna) in pandas dataframe
    # T_rna is t(ref_bulkRNA)
    # D_met is t(mix_met)
    # T_met is t(ref_met)


    D_rna = rna_unit['mix'].T
    T_rna = rna_unit['ref'].T 

    D_met = met_unit['mix'].T
    T_met = met_unit['ref'].T 

    n_sample = D_rna.shape[0]
    n_ref = T_rna.shape[0]

    data_rna = pd.concat([D_rna, T_rna], ignore_index=True, sort=False)
    data_met = pd.concat([D_met, T_met], ignore_index=True, sort=False)


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


    adata = MultiMAP.Integration(
        [adata_rna, adata_met],
        use_reps=['X_pca', 'X_pca'],
        scale=True,
        n_components=10
    )

    def adata_to_df(adata):
        """Convert AnnData object to dense pandas DataFrame"""
        X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
        return pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)

    adata_mix = adata[adata.obs['type'] == 'sample']
    adata_ref = adata[adata.obs['type'] == 'ref']

    unit = {
        "mix": adata_to_df(adata_mix).T,
        "ref": adata_to_df(adata_ref).T,
        "ref_scRNA" : rna_unit["ref_scRNA"]
    }

    return unit
