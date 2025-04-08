program_block_PP <- function(multi_data) {

  scale_matrix <- function(mat) {
    mat = sweep(mat, 2, colSums(mat), "/")
    return(mat)
  }
  multi_data$ref$ref_met = scale_matrix(multi_data$ref$ref_met)
  multi_data$mix$mix_met = scale_matrix(multi_data$mix$mix_met)
  
  return(multi_data) 
}