program_block <- function(l_pred) { 


  
  multi_pred <- Reduce("+",l_pred) / length(l_pred)

  return(multi_pred)
}
