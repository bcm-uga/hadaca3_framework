program_block_PP <- function(data,path_og_dataset='') {

  scale_matrix <- function(mat) {
    mat = sweep(mat, 2, colSums(mat), "/")
    return(mat)
  }

  if(is.list(data)){
    data = lapply(data, function(x) {
      list(counts=scale_matrix(x$counts),
           metadata=x$metadata)
    })
  }else{
    data = scale_matrix(data)
  }
  # multi_data$mix$mix_met = scale_matrix(multi_data$mix$mix_met)
  
  # stopifnot()


  return(data) 
}