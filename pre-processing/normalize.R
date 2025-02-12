program_block <- function(mix = NULL, ref = NULL, 
                    ...) {
  source("scr/teamF_Source_prior_known_features.R")                   
   # Normalize input matrices

    ref$ref_bulkRNA = normalize_matrix(ref$ref_bulkRNA)
    ref$ref_met = normalize_matrix(ref$ref_met)


  
 multi_data = list(mix = mix,
             ref = ref)
  
  return(multi_data) 
}