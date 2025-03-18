# uni_data = SP_data$RNA
get_prop_nnls <- function(mix, ref) {
  prop = apply(mix, 2, function(b, A) {
    tmp_prop = nnls::nnls(b=b, A=A)$x
    tmp_prop = tmp_prop / sum(tmp_prop) # Sum To One
    return(tmp_prop)
  }, A=ref)
  rownames(prop) <- colnames(ref)
  
  return(prop)
}

program_block_DE <- function(uni_data) {




  idx_feat = intersect(rownames(uni_data$mix), rownames(uni_data$ref$bulk))
  uni_data$mix = uni_data$mix[idx_feat,]
  uni_data$ref$bulk = uni_data$ref$bulk[idx_feat,]
  uni_pred = get_prop_nnls(uni_data$mix,  uni_data$ref$bulk)
  
  return(uni_pred) 
}

