

program_block_DE <- function(uni_data,path_og_dataset='') {


    prop <- t(epidish(uni_data$mix,uni_data$ref,method="RPC")$estF)

  
  return(prop) 
}

