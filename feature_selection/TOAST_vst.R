
program_block_FS <- function(multi_data) {
  
  nb_fs_rna = 1e3

  determine_variable_genes <- function(mat, n_genes=nb_fs_rna) {
    vst <- log2(mat + median(mat[mat > 0]))
    var_sorted_genes <- TOAST::findRefinx(vst, nmarker = min(nrow(multi_data$ref$ref_bulkRNA),n_genes))
    return(var_sorted_genes)
  }
  top_genes <- determine_variable_genes(multi_data$ref$ref_bulkRNA, nb_fs_rna)
  
  multi_data$mix$mix_rna = multi_data$mix$mix_rna[top_genes,]
  multi_data$ref$ref_bulkRNA = multi_data$ref$ref_bulkRNA[top_genes,]
  multi_data$ref$ref_scRNA = lapply(multi_data$ref$ref_scRNA, function(x) list(counts = x$counts[top_genes,], metadata = x$metadata))
  
  return(multi_data) 
}

