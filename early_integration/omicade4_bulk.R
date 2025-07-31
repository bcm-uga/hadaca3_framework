program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

  mix_rna = rna_unit$mix 
  ref_rna = rna_unit$ref
  mix_met = met_unit$mix 
  ref_met = met_unit$ref


  if (!("omicade4" %in% installed.packages())) {
    BiocManager::install("omicade4")
  }
  
  # add samples' names
  if (is.null(colnames(mix_rna))) {
    colnames(mix_rna) = paste0("Sample",seq(ncol(mix_rna)))
  }
  if (is.null(colnames(mix_met))) {
    colnames(mix_met) = paste0("Sample",seq(ncol(mix_met)))
  }
  mix_rna = mix_rna[,colnames(mix_met)]
  
  # run omicade4
  library(omicade4)
  data_combi <- list(RNA=as.matrix(cbind(mix_rna,ref_rna)),
                      DNAm=as.matrix(cbind(mix_met,ref_met)))
  combi <- mcia(data_combi, cia.nf=10)
  
  # Latent space is in combi$mcoa$SynVar
  projection = t(combi$mcoa$SynVar)
  projection = projection - min(projection)
  colnames(projection) = colnames(cbind(mix_rna,ref_rna))
  
  # Retrieve D and T
  mix = projection[,colnames(mix_rna)]
  ref = projection[,colnames(ref_rna)]

  res_unit = list(mix=mix, ref = ref  )
  return(res_unit)
}
