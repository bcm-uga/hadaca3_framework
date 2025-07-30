program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

  
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
  ref_scRNA <- lapply(ref_scRNA, function(x) list(counts=x$counts[rownames(mix_rna),], metadata=x$metadata))
  ref_met <- ref_met[rownames(mix_met),]
  
  mix = list(rna=mix_rna, met=mix_met)
  ref = list(rna=ref_scRNA, met=ref_met)

  return(rna_unit)
}
