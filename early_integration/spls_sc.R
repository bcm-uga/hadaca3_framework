program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

  mix_rna = rna_unit$mix 
  ref_rna = rna_unit$ref
  mix_met = met_unit$mix 
  ref_met = met_unit$ref

  ref_scRNA= rna_unit$ref_scRNA

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
  ref_scRNA <- lapply(ref_scRNA, function(x) {
   common_genes <- intersect(rownames(x$counts), rownames(mix_rna))
   list(counts=x$counts[rownames(common_genes),], metadata=x$metadata)
  }
   )
  ref_met <- ref_met[rownames(mix_met),]
  
  mix = list(rna=mix_rna, met=mix_met)
  ref = list(rna=ref_scRNA, met=ref_met)


  res_unit = list(mix=mix, ref = ref  )
  return(res_unit)
  # return(rna_unit)
}
