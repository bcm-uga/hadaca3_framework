
program_block_FS <- function(data,path_og_dataset='') {
  
  nb_fs_rna = 1e3
  og_ref_bulkRNA  =  read_all_ref_hdf5(path_og_dataset$ref,to_read = 'ref_bulkRNA')$ref_bulkRNA

  hvg <- TOAST::findRefinx(og_ref_bulkRNA, nmarker = min(nrow(og_ref_bulkRNA),nb_fs_rna))

  # multi_data$mix$mix_rna = multi_data$mix$mix_rna[hvg,]
  # multi_data$ref$ref_bulkRNA = multi_data$ref$ref_bulkRNA[hvg,]
  # multi_data$ref$ref_scRNA = lapply(multi_data$ref$ref_scRNA, function(x) list(counts = x$counts[hvg,], metadata = x$metadata))

  if(is.list(data)){
    data = lapply(data, function(x) list(counts = x$counts[hvg,], metadata = x$metadata))
  }else{
    data = data[hvg,]
  }
  return(data) 
}

