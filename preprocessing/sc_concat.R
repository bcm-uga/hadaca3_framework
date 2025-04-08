program_block_PP <- function(multi_data) {

  warning("should be sample-normalized")
  metadata <- do.call(rbind, lapply(multi_data$ref$ref_scRNA, function(x) x$metadata))
  ref_scRNA <- lapply(multi_data$ref$ref_scRNA, function(x) as.matrix(x$counts))
  shared_genes <- Reduce(intersect, lapply(ref_scRNA, rownames))
  ref_scRNA <- lapply(ref_scRNA, function(x) x[shared_genes,])
  ref_scRNA <- lapply(ref_scRNA, function(x) as(x,'dgCMatrix'))
  ref_scRNA_all <- do.call(cbind, ref_scRNA)
  
  metadata$dataset = sapply(rownames(metadata), function(x) strsplit(x, ".", fixed=T)[[1]][1])
  ref_scRNA = list("ref_concat"=list(counts = ref_scRNA_all,
                                       metadata = metadata))
  
  return(multi_data) 
}