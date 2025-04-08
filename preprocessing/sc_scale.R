program_block_PP <- function(multi_data) {


        
  multi_data$ref$ref_scRNA <- lapply(multi_data$ref$ref_scRNA, function(x) {
      logical_matrix = x$counts > 1000 ; 
      x$counts[logical_matrix] = 1000;
      list(counts = x$counts,
           metadata = x$metadata)})
  
  return(multi_data) 
}