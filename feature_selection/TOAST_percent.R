program_block_FS <- function(data, path_og_dataset='') {
  
  TOAST_percent_met = 0.8

  og_ref_met = read_all_ref_hdf5(path_og_dataset$ref, to_read = 'ref_met')$ref_met

  nmarker_percent = round(TOAST_percent_met*nrow(og_ref_met))
  hvp <- TOAST::findRefinx(og_ref_met, nmarker = nmarker_percent)

  # multi_data$mix$mix_met <- multi_data$mix$mix_met[hvp,]
  # multi_data$ref$ref_met <- multi_data$ref$ref_met[hvp,]
  
  return(data[hvp,]) 
}

