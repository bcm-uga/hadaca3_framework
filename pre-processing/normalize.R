program_block_PP <- function(mix = NULL, ref = NULL, 
                    ...) {
  # source("scr/teamF_Source_prior_known_features.R")                   
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
    ref$ref_bulkRNA = normalize_matrix(ref$ref_bulkRNA)
    ref$ref_met = normalize_matrix(ref$ref_met)


  
 multi_data = list(mix = mix,
             ref = ref)
  
  return(multi_data) 
}