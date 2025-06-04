
program_block_FS <- function(data,path_og_dataset='') {
  
  og_ref_met  =  read_all_ref_hdf5(path_og_dataset$ref,to_read = 'ref_MET')$ref_MET
  
  
  # highest methylation based on biological knowledge (cancer has high proportion of methylated probes)
  
  quantiles = apply(og_ref_met, 2, quantile, 0.75) # select 25% of the probes
  methylated_probes = Reduce(  union,
    lapply(1:length(quantiles), function(x)  
      rownames(og_ref_met)[which(og_ref_met[,x] > quantiles[x])]))

  data = data[methylated_probes,]
  # ref_met = ref_met[methylated_probes,]
 
  return(data) 
}

