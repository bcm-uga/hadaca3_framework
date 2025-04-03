program_block <- function(l_pred) { 

  prop1 = l_pred$prop1
  prop2 = l_pred$prop2
  mix = l_pred$mix

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
