program_block_EI <- function(rna_unit,met_unit,path_dataset) { 
  mix_rna = rna_unit$mix 
  ref_rna = rna_unit$ref
  mix_met = met_unit$mix 
  ref_met = met_unit$ref

  warning("not really an integration step, but not really a FS step either")
  if (!("mixOmics" %in% installed.packages())) {
    BiocManager::install("mixOmics")
  }
  
  n_feat = min(nrow(mix_rna), nrow(mix_met), 1e4)
  res.spls <- mixOmics::spls(t(mix_rna[sample(1:nrow(mix_rna), size = n_feat), ]),
                              t(mix_met[sample(1:nrow(mix_met), size = n_feat), ]),
                              keepX = c(100, 100), keepY = c(100, 100))
  mix_rna = mix_rna[rownames(res.spls$loadings$X),]
  mix_met <- mix_met[rownames(res.spls$loadings$Y),]

  common_genes <- intersect(rownames(mix_rna), rownames(ref_rna))

  mix_rna_aligned <- mix_rna[common_genes, ]
  ref_rna_aligned <- ref_rna[common_genes, ]
  
  common_meth_probes = intersect(rownames(mix_met), rownames(ref_met))
  
  mix_met_aligned <- mix_met[common_meth_probes, ]
  ref_met_aligned <- ref_met[common_meth_probes, ]



  # ref_rna <- ref_rna[rownames(mix_rna_aligned),]
  # ref_met <- ref_met[rownames(mix_met_aligned),]
  
  mix = list(rna=mix_rna_aligned, met=mix_met_aligned)
  ref = list(rna=ref_rna_aligned, met=ref_met_aligned)
  
  res_unit = list(mix=mix, ref = ref  )
  return(res_unit)
  # return(rna_unit)
}
