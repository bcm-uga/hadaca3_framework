program_block_PP <- function(data,path_og_dataset='') {

  scale_matrix <- function(mat) {
    mat = sweep(mat, 2, colSums(mat), "/")
    return(mat)
  }
  data = scale_matrix(data)
  # multi_data$mix$mix_met = scale_matrix(multi_data$mix$mix_met)
  
  return(data) 
}