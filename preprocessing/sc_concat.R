program_block_PP <- function(data, path_og_dataset='', omic='') {

  #only for scRNA
  stopifnot(is.list(data))
  
  warning("should be sample-normalized")
  
  metadata <- do.call(rbind, lapply(data, function(x) x$metadata))
  ref_scRNA <- lapply(data, function(x) as.matrix(x$counts))
  shared_genes <- Reduce(intersect, lapply(ref_scRNA, rownames))
  ref_scRNA <- lapply(ref_scRNA, function(x) x[shared_genes,])
  ref_scRNA <- lapply(ref_scRNA, function(x) as(x,'dgCMatrix'))
  ref_scRNA_all <- do.call(cbind, ref_scRNA)
  
  metadata$dataset = sapply(rownames(metadata), function(x) strsplit(x, ".", fixed=T)[[1]][1])
  data = list("ref_concat"=list(counts = ref_scRNA_all,
                                metadata = metadata))
  
  return(data) 
}