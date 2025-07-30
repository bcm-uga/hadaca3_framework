program_block_EI <- function(rna_unit,met_unit,path_dataset) { 


  mix_rna = rna_unit$mix 
  ref_rna = rna_unit$ref
  mix_met = met_unit$mix 
  ref_met = met_unit$ref

  mix = rbind(mix_rna, mix_met)
  ref = rbind(ref_rna, ref_met)


  res_unit = list(mix=mix, ref = ref  )
  return(res_unit)
}
