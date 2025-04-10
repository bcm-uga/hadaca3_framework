# source("utils/data_processing.R") to if you want to use outside the wrapper "04_..."

program_block_li <- function(prop1,prop2,path_dataset) { 

  # prop1 = l_pred$prop1
  # prop2 = l_pred$prop2
  # path_dataset = l_pred$last_dataset
  mix = read_all_hdf5(path_dataset,('mix'))$mix$mix_rna

  real = mix / colSums(mix)
  reconstructed = lapply(list(prop1, prop2), function(x)
    (ref %*% x) / colSums(ref %*% x))
  rmse = sapply(reconstructed, function(x) 
    sqrt(mean(((as.matrix(x)) - as.matrix(real))^2)))
  
  # Normalize RMSEs
  rmse1_norm = rmse[[1]] / Reduce(`+`, rmse)
  rmse2_norm = rmse[[2]] / Reduce(`+`, rmse)
  
  # RMSE as weights  
  prop <- as.matrix(rmse1_norm * prop2 + rmse2_norm  * prop1)
  # multi_pred <- Reduce("+",l_pred) / length(l_pred)

  return(prop)
}


