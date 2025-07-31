program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

  mix_rna = rna_unit$mix 
  ref_rna = rna_unit$ref
  mix_met = met_unit$mix 
  ref_met = met_unit$ref

  if (!("mixKernel" %in% installed.packages())) {
    BiocManager::install("mixKernel")
  }
  library(mixKernel)

  # add columns' names and order in the same way
  if (is.null(colnames(mix_rna))) {
    colnames(mix_rna) = paste0("Sample",seq(ncol(mix_rna)))
  }
  if (is.null(colnames(mix_met))) {
    colnames(mix_met) = paste0("Sample",seq(ncol(mix_met)))
  }
  mix_rna = mix_rna[,colnames(mix_met)]
  
  # compute kernels
  kernel_rna <- compute.kernel(t(as.matrix(cbind(mix_rna,ref_rna))), kernel.func = "abundance")
  kernel_met <- compute.kernel(t(as.matrix(cbind(mix_met,ref_met))), kernel.func = "abundance")
  kernel_all <- combine.kernels(kernel_rna = kernel_rna, kernel_met = kernel_met)
  kernel_pca <- kernel.pca(kernel_all, ncomp = 10)
  
  # Latent space is in kernel_pca$variates$X
  projection = t(kernel_pca$variates$X)
  projection = projection - min(projection)
  mix = projection[,colnames(mix_rna)]
  ref = projection[,colnames(ref_rna)]

  res_unit = list(mix=mix, ref = ref  )
  return(res_unit)
}
