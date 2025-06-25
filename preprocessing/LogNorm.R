program_block_PP <- function(data, path_og_dataset='', omic='') {
  
  if (omic == 'ref_scRNA') { # is.list(data)
    data <- lapply(data, function(x) list(counts=Seurat::LogNormalize(x$counts), metadata=x$metadata))
  } else {
    data = exp(Seurat::LogNormalize(data))
  }
  
  return(data) 
} 