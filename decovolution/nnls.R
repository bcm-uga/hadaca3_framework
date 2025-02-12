# uni_data = SP_data$RNA

program_block <- function(uni_data) {

    ## RNA:
  #if ( !( is.null(x = mix_rna) ) ) {
    
  idx_feat = intersect(rownames(uni_data$mix), rownames(uni_data$ref$bulk))
  uni_data$mix = uni_data$mix[idx_feat,]
  uni_data$ref$bulk = uni_data$ref$bulk[idx_feat,]
  uni_pred = get_prop_nnls(uni_data$mix,  uni_data$ref$bulk)
  
  return(uni_pred) # MAG : je ne sais pas quoi sortir 
}


# pred_RNA = program_blockDE(uni_data = SP_data$RNA)
# pred_met = program_blockDE(uni_data = SP_data$met)