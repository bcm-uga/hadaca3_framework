
program_block_FS <- function(data,path_og_dataset) {
  
  nb_fs_rna = 1e3

  og_ref_bulkRNA  =  read_all_ref_hdf5(path_og_dataset$ref,to_read = 'ref_bulkRNA')$ref_bulkRNA


  determine_variable_genes <- function(mat, n_genes=nb_fs_rna) {
    vst <- log2(mat + median(mat[mat > 0]))
    var_sorted_genes <- TOAST::findRefinx(vst, nmarker = min(nrow(og_ref_bulkRNA),n_genes))
    return(var_sorted_genes)
  }
  top_genes <- determine_variable_genes(og_ref_bulkRNA, nb_fs_rna)
  
  # multi_data$mix$mix_rna = multi_data$mix$mix_rna[top_genes,]
  # multi_data$ref$ref_bulkRNA = multi_data$ref$ref_bulkRNA[top_genes,]

  if(is.list(data)){
    data = lapply(data, function(x) list(counts = x$counts[top_genes,], metadata = x$metadata))
  }else{
    data = data[top_genes,]
  }
  
  return(data) 
}

