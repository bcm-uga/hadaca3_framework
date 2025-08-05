program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

  mix_rna = rna_unit$mix 
  ref_rna = rna_unit$ref
  mix_met = met_unit$mix 
  ref_met = met_unit$ref


  if (!("MOFA2" %in% installed.packages())) {
    BiocManager::install("MOFA2")
  }
  library(reticulate)
  use_condaenv("hadaca3framework_env", required = TRUE)
  library(MOFA2)
  
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

  
  # MOFA
  MOFA <- create_mofa(list("RNA"=as.matrix(cbind(mix_rna_aligned,ref_rna_aligned)),
                            "DNAm"=as.matrix(cbind(mix_met_aligned,ref_met_aligned))))
  model_opts <- get_default_model_options(MOFA)
  model_opts$num_factors <- min(5,ncol(ref_rna_aligned)+1)
  train_opts <- get_default_training_options(MOFA)
  train_opts$seed <- 12
  MOFA <- prepare_mofa(MOFA, model_options = model_opts, training_options = train_opts)
  dir.create("tmp_mofa")
  MOFA <- run_mofa(MOFA, save_data=F, outfile="tmp_mofa/model.hdf5",  use_basilisk=T) # ,  use_basilisk = FALSE)
  
  # Latent space is in MOFA@expectations$Z$group1
  projection = t(MOFA@expectations$Z$group1)
  projection = projection - min(projection)
  
  # Retrieve D and T
  mix = projection[,colnames(mix_rna_aligned)]
  ref = projection[,colnames(ref_rna_aligned)]

  res_unit = list(mix=mix, ref = ref  )
  return(res_unit)
}
