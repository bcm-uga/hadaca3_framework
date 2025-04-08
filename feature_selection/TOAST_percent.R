
program_block_FS <- function(multi_data) {
  
  TOAST_percent_met = 0.8

  nmarker_percent = round(TOAST_percent_met*nrow(multi_data$ref$ref_met))
  hvp <- TOAST::findRefinx(multi_data$ref$ref_met, nmarker = nmarker_percent)

  multi_data$mix$mix_met <- multi_data$mix$mix_met[hvp,]
  multi_data$ref$ref_met <- multi_data$ref$ref_met[hvp,]
 
  return(multi_data) 
}

