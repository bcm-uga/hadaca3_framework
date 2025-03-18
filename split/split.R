
program_blockSP <- function(multi_data) {
   
   uni_data_RNA = list(mix = multi_data$mix$mix_rna,
                       ref = list(bulk = multi_data$ref$ref_bulkRNA,
                                  scRNA =  multi_data$ref$ref_scRNA))
   
      
   uni_data_met = list(mix = multi_data$mix$mix_met,
                       ref = list(bulk = multi_data$ref$ref_met,
                                  scRNA =  multi_data$ref$ref_scRNA))
   
  uni_data = list(RNA =  uni_data_RNA,
                  met = uni_data_met)
  return(uni_data) 
}
