program_block_PP <- function(data,path_og_dataset='',omic='') {
    if(omic == 'ref_scRNA' ){ # is.list(data)
    if (!"edgeR" %in% rownames(installed.packages())) {
      library(BiocManager)
      BiocManager::install("edgeR")
    }
    
    # get binary markers for each cell type in each sample in each dataset
    pseudobulks <- lapply(data,\(x){
      counts <- x$counts
      unique_clones <- with(subset(x$metadata, !is.na(cell_type)), unique(paste0(sample, " - ", cell_type)))
      clones_matrix_id <- matrix(nrow = ncol(counts), ncol = length(unique_clones),
                                  dimnames = list(colnames(counts), unique_clones),
                                  data = 0)
      for (pt in unique(x$metadata$sample)) {
        for (ct in unique(x$metadata$cell_type)) {
          column <- paste0(pt, " - ", ct)
          if (column %in% colnames(clones_matrix_id)) {
            clones_matrix_id[rownames(subset(x$metadata, sample == pt & cell_type == ct)), column] <- 1
          }
        }
      }
      counts %*% clones_matrix_id})
      
      ## Generate normalized matrix and associate metadata
    shared_genes <- Reduce(intersect, lapply(pseudobulks, rownames))
    pseudobulks_matrix <- do.call(cbind, lapply(pseudobulks, function(x) {
      x[shared_genes, ]}))
    pseudobulks_normalized <- pseudobulks_matrix %*% Matrix::Diagonal(x = edgeR::calcNormFactors(pseudobulks_matrix))
    dimnames(pseudobulks_normalized) <- dimnames(pseudobulks_matrix)
    pseudobulks_metadata <- data.frame(
      sample = sapply(strsplit(colnames(pseudobulks_matrix), " - "), "[", 1),
      cell_type = sapply(strsplit(colnames(pseudobulks_matrix), " - "), "[", 2),
      dataset = unlist(lapply(seq_along(data), function(x)
        rep(names(data)[x], ncol(pseudobulks[[x]])))))
      
    


      data = list(ref_binarypseudobulk_log=list(
        counts=log2(1+pseudobulks_normalized),
        metadata=pseudobulks_metadata))
      
    }else{
     data = log2(1+data)
    # mix_rna = log2(1+mix_rna)
    # ref_bulkRNA = log2(1+ref_bulkRNA)

    }
  
  return(data) 
}