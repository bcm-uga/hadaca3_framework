program_block_FS <- function(data, path_og_dataset) {
  
  nb_fs_rna = 1e3

  og_ref_bulkRNA = read_all_ref_hdf5(path_og_dataset$ref, to_read = 'ref_bulkRNA')$ref_bulkRNA

  determine_variable_genes <- function(mat, n_genes=nb_fs_rna) {
    vst <- log2(mat + median(mat[mat > 0]))
    var_sorted_genes <- TOAST::findRefinx(vst, nmarker = min(nrow(og_ref_bulkRNA), n_genes))
    return(var_sorted_genes)
  }
  top_genes <- determine_variable_genes(og_ref_bulkRNA, nb_fs_rna)
  top_gene_names <- rownames(og_ref_bulkRNA)[top_genes]

  # multi_data$mix$mix_rna = multi_data$mix$mix_rna[top_genes,]
  # multi_data$ref$ref_bulkRNA = multi_data$ref$ref_bulkRNA[top_genes,]

  convert_and_filter_sparse_matrix <- function(x, top_gene_names){
    dense_counts <- as.matrix(x$counts)
    common_genes <- intersect(top_gene_names, rownames(dense_counts))

    subset_counts <- dense_counts[common_genes, , drop = FALSE]
    sparse_counts <- Matrix(subset_counts, sparse = TRUE)

    return(list(counts = sparse_counts, metadata = x$metadata) )
  }

  if (is.list(data)) {
    data = lapply(data, function(x) convert_and_filter_sparse_matrix(x, top_gene_names))
  } else {
    data = data[ intersect(top_gene_names,rownames(data)),]
  }
  
  return(data) 
}

