

program_block_DE <- function(uni_data,path_og_dataset='') {
  
  idx_feat = intersect(rownames(uni_data$mix), rownames(uni_data$ref))
  uni_data$mix = uni_data$mix[idx_feat,]
  uni_data$ref = uni_data$ref[idx_feat,]
  
  
  prop_rlr = t(EpiDISH::epidish(uni_data$mix,uni_data$ref,method="RPC")$estF)  
  
  prop_nnls = apply(uni_data$mix, 2, function(b, A) {
    tmp_prop = nnls::nnls(b=b, A=A)$x
    tmp_prop = tmp_prop / sum(tmp_prop)
    return(tmp_prop)}, A=uni_data$ref)  
  rownames(prop_rlr) = colnames(uni_data$ref)
  rownames(prop_nnls) = colnames(uni_data$ref)
  
  reconstructed_rlr = (uni_data$ref %*% prop_rlr) / colSums(uni_data$ref %*% prop_rlr)
  reconstructed_nnls = (uni_data$ref %*% prop_nnls) / colSums(uni_data$ref %*% prop_nnls)
  real = uni_data$mix / colSums(uni_data$mix)
  rmse_rlr = sqrt(mean(((as.matrix(reconstructed_rlr)) - as.matrix(real))^2))
  rmse_nnls = sqrt(mean(((as.matrix(reconstructed_nnls)) - as.matrix(real))^2))
    
  # Choose the best method based on RMSE
  if (rmse_rlr < rmse_nnls) {
    prop = prop_rlr
  } else {
    prop = prop_nnls
  }

  
  return(prop) 
}

