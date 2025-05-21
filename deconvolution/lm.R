

# pred_RNA = program_blockDE(uni_data = SP_data$RNA)
# pred_met = program_blockDE(uni_data = SP_data$met)

program_block_DE = function(uni_data,path_og_dataset='') {

  ##
  ## YOUR CODE BEGINS HERE
  ##
  # idx_feat corresponds to the intersection of features present in the references and in the mixtures.
  
  mix = uni_data$mix
  ref = uni_data$ref
  idx_feat = intersect(rownames(mix), rownames(ref))


  
  # Estimation of proportions
  prop = apply(mix[idx_feat,], 2, function(b, A) {
    tmp_prop = lm(b ~ A - 1)$coefficients  # Using `-1` to remove the intercept
    tmp_prop[tmp_prop < 0] = 0
    tmp_prop = tmp_prop / sum(tmp_prop)    # Sum To One
    return(tmp_prop)
  }, A=ref[idx_feat,])

  # Labeling of estimated proportions 
  rownames(prop) = colnames(ref)
  return(prop)
  
  ##
  ## YOUR CODE ENDS HERE
  ##
}