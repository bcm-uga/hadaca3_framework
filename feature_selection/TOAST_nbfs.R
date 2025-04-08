
program_block_FS <- function(multi_data) {
  
  nb_fs_met=nrow(multi_data$ref$ref_bulkRNA)

  hvp <- TOAST::findRefinx(multi_data$ref$ref_met, nmarker = min(nrow(multi_data$ref$ref_met),nb_fs_met))

  multi_data$mix$mix_met <- multi_data$mix$mix_met[hvp,]
  multi_data$ref$ref_met <- multi_data$ref$ref_met[hvp,]
 
  return(multi_data) 
}

