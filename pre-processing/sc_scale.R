program_block_PP <- function(mix = NULL, ref = NULL, 
                    ...) {
  
  preprocess_matrix <- function(matrix) {
    matrix <- as.matrix(matrix)
    # Remove rows  with NA or Inf
    matrix <- matrix[complete.cases(matrix), , drop = FALSE]
    matrix <- matrix[rowSums(is.infinite(matrix)) == 0, , drop = FALSE]
    # Remove rows with all zeros
    matrix <- matrix[rowSums(matrix) > 0, ]  # Remove rows with all zeros
    # Remove rows with zero variance
    matrix <- matrix[apply(matrix, 1, var) > 0, , drop = FALSE]
    return(matrix)
  }
  
  # Preprocess all matrices
  mix$mix_rna <- preprocess_matrix(mix$mix_rna)
  ref$ref_bulkRNA <- preprocess_matrix(ref$ref_bulkRNA)
  ref$ref_scRNA <- lapply(ref$ref_scRNA, function(x) {
    list(counts=preprocess_matrix(x$counts),
         metadata=x$metadata)})                    


  ref$ref_scRNA <- lapply(ref$ref_scRNA, function(x) {
      list(counts = x$counts[x$counts > 1000] = 1000,
           metadata = x$metadata)})


  multi_data = list(mix = mix,
             ref = ref)
  
  return(multi_data) 
}