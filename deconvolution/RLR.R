

program_block_DE <- function(uni_data,path_og_dataset='') {

  idx_feat = intersect(rownames(uni_data$mix), rownames(uni_data$ref))
  uni_data$mix = uni_data$mix[idx_feat,]
  uni_data$ref = uni_data$ref[idx_feat,]
  
  
  prop <- t(EpiDISH::epidish(uni_data$mix,uni_data$ref,method="RPC")$estF)

  return(prop) 
}

