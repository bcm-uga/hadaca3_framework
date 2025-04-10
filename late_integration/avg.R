program_block_li <- function(prop1,prop2,path_dataset) { 
  # prop1 = l_pred$prop1
  # prop2 = l_pred$prop2

  prop = Reduce(`+`, list(prop1,prop2)) / length(list(prop1,prop2))

  return(prop)
}
