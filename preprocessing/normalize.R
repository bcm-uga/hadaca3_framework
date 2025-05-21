program_block_PP <- function(data,path_og_dataset='') {    
    # Normalize input matrices
    #' Normalize a Matrix by Column Sums
    #'
    #' @param mat A numeric matrix where rows represent features and columns represent 
    #'            samples.
    #'
    #' @return A numeric matrix with the same dimensions as the input `mat`, where each 
    #'         column has been normalized so that its elements sum to one.
    #'
    normalize_matrix <- function(mat) {
      mat = sweep(mat, 2, colSums(mat), "/")
      return(mat)
    }
    # multi_data$ref$ref_bulkRNA = normalize_matrix(multi_data$ref$ref_bulkRNA)
    # multi_data$ref$ref_met = normalize_matrix(multi_data$ref$ref_met)
    data = normalize_matrix(data)

  return(data) 
}