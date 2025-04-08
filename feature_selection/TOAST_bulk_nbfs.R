
program_block_FS <- function(multi_data) {
  
  nb_fs_rna = 1e3

  hvg <- TOAST::findRefinx(multi_data$ref$ref_bulkRNA, nmarker = min(nrow(multi_data$ref$ref_bulkRNA),nb_fs_rna))
  multi_data$mix$mix_rna = multi_data$mix$mix_rna[hvg,]
  multi_data$ref$ref_bulkRNA = multi_data$ref$ref_bulkRNA[hvg,]
  multi_data$ref$ref_scRNA = lapply(multi_data$ref$ref_scRNA, function(x) list(counts = x$counts[hvg,], metadata = x$metadata))
  return(multi_data) 
}

