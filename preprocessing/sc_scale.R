program_block_PP <- function(data, path_og_dataset='', omic='') {
  #only for scRNA
  stopifnot(is.list(data))
        
  data <- lapply(data, function(x) {
    logical_matrix = x$counts > 1000
    x$counts[logical_matrix] = 1000
    list(counts = x$counts,
         metadata = x$metadata)})
  
  return(data) 
}