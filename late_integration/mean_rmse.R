# source("utils/data_processing.R") to if you want to use outside the wrapper "04_..."

program_block_li <- function(prop1,prop2,path_dataset) { 

  # prop1 = l_pred$prop1
  # prop2 = l_pred$prop2
  # path_dataset = l_pred$last_dataset
  mix = read_all_hdf5(path_dataset$mix,('mix'))$mix
  mix_rna = mix$mix_rna
  mix_met = mix$mix_met
  ref = read_all_hdf5(path_dataset$ref,('ref'))$ref
  ref_bulkRNA = ref$ref_bulkRNA
  ref_met = ref$ref_met
  
  idx_feat_rna = intersect(rownames(mix_rna), rownames(ref_bulkRNA))
  mix_rna = mix_rna[idx_feat_rna,]
  ref_bulkRNA = ref_bulkRNA[idx_feat_rna,]
  idx_feat_met = intersect(rownames(mix_met), rownames(ref_met))
  mix_met = mix_met[idx_feat_met,]
  ref_met = ref_met[idx_feat_met,]
  

  real_rna = mix_rna / colSums(mix_rna)
  reconstructed_rna = (ref_bulkRNA %*% prop1) / colSums(ref_bulkRNA %*% prop1)
  rmse_rna = sqrt(mean(((as.matrix(reconstructed_rna)) - as.matrix(real_rna))^2))
    
  real_met = mix_met / colSums(mix_met)
  reconstructed_met = (ref_met %*% prop2) / colSums(ref_met %*% prop2)
  rmse_met = sqrt(mean(((as.matrix(reconstructed_met)) - as.matrix(real_met))^2))
  
  # Normalize RMSEs
  rmse_rna_norm = rmse_rna / (rmse_rna + rmse_met)
  rmse_met_norm = rmse_met / (rmse_rna + rmse_met)
  
  # RMSE as weights  
  prop <- as.matrix( rmse_met_norm * prop1 + rmse_rna_norm  * prop2 )
  # if prop1 is prop_rna and prop2 is prop_met

  return(prop)
}


