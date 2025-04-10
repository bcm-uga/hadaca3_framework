program_block_li <- function(prop1,prop2,path_dataset) { 

  # prop1 = l_pred$prop1
  # prop2 = l_pred$prop2

  prop <- rbind(
    prop1[c('fibro'), , drop = FALSE] * 1.8 ,
    prop1[c('endo'), , drop = FALSE] *2.5,
    prop1[c('immune'), , drop = FALSE] ,
    prop2[c('classic'), , drop = FALSE] *3,
    prop2[c('basal'), , drop = FALSE] *3
  )
  prop <- apply(prop, 2, function(x) x / sum(x))

  return(prop)
}
