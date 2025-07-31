program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

  mix_rna = rna_unit$mix 
  ref_rna = rna_unit$ref
  mix_met = met_unit$mix 
  ref_met = met_unit$ref


  if (!("MOFA2" %in% installed.packages())) {
    BiocManager::install("MOFA2")
  }
  library(MOFA2)
  
  # add samples' names
  if (is.null(colnames(mix_rna))) {
    colnames(mix_rna) = paste0("Sample",seq(ncol(mix_rna)))
  }
  if (is.null(colnames(mix_met))) {
    colnames(mix_met) = paste0("Sample",seq(ncol(mix_met)))
  }
  mix_rna = mix_rna[,colnames(mix_met)]
  
  # MOFA
  MOFA <- create_mofa(list("RNA"=as.matrix(cbind(mix_rna,ref_rna)),
                            "DNAm"=as.matrix(cbind(mix_met,ref_met))))
  model_opts <- get_default_model_options(MOFA)
  model_opts$num_factors <- min(5,ncol(ref_rna)+1)
  train_opts <- get_default_training_options(MOFA)
  train_opts$seed <- 12
  MOFA <- prepare_mofa(MOFA, model_options = model_opts, training_options = train_opts)
  dir.create("tmp_mofa")
  MOFA <- run_mofa(MOFA, save_data=F, outfile="tmp_mofa/model.hdf5", use_basilisk=T)
  
  # Latent space is in MOFA@expectations$Z$group1
  projection = t(MOFA@expectations$Z$group1)
  projection = projection - min(projection)
  
  # Retrieve D and T
  mix = projection[,colnames(mix_rna)]
  ref = projection[,colnames(ref_rna)]

  res_unit = list(mix=mix, ref = ref  )
  return(res_unit)
}
