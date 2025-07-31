program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

  mix_rna = rna_unit$mix 
  ref_rna = rna_unit$ref
  mix_met = met_unit$mix 
  ref_met = met_unit$ref


  if (!("omicade4" %in% installed.packages())) {
    BiocManager::install("omicade4")
  }
  library(omicade4)
  
  # add samples' names
  if (is.null(colnames(mix_rna))) {
    colnames(mix_rna) = paste0("Sample",seq(ncol(mix_rna)))
  }
  if (is.null(colnames(mix_met))) {
    colnames(mix_met) = paste0("Sample",seq(ncol(mix_met)))
  }
  mix_rna = mix_rna[,colnames(mix_met)]

  common_genes <- intersect(rownames(mix_rna), rownames(ref_rna))

  mix_rna_aligned <- mix_rna[common_genes, ]
  ref_rna_aligned <- ref_rna[common_genes, ]
  
  common_meth_probes = intersect(rownames(mix_met), rownames(ref_met))
  
  mix_met_aligned <- mix_met[common_meth_probes, ]
  ref_met_aligned <- ref_met[common_meth_probes, ]

  
  # run omicade4
  data_combi <- list(RNA=as.matrix(cbind(mix_rna_aligned,ref_rna_aligned)),
                      DNAm=as.matrix(cbind(mix_met_aligned,ref_met_aligned)))
  combi <- mcia(data_combi, cia.nf=10)
  
  # Latent space is in combi$mcoa$SynVar
  projection = t(combi$mcoa$SynVar)
  projection = projection - min(projection)
  colnames(projection) = colnames(cbind(mix_rna_aligned,ref_rna_aligned))
  
  # Retrieve D and T
  mix = projection[,colnames(mix_rna_aligned)]
  ref = projection[,colnames(ref_rna_aligned)]

  res_unit = list(mix=mix, ref = ref  )
  return(res_unit)
}
