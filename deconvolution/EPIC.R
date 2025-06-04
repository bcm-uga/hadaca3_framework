

# pred_RNA = program_blockDE(uni_data = SP_data$RNA)
# pred_met = program_blockDE(uni_data = SP_data$met)

program_block_DE = function(uni_data,path_og_dataset='') {
  
  
  prop <- t(EPIC::EPIC(as.matrix(uni_data$mix),
        reference = list(refProfiles = uni_data$ref, sigGenes = rownames(uni_data$mix)),
        scaleExprs = F, withOtherCells = F)$mRNAProportions)

  return(prop)

}