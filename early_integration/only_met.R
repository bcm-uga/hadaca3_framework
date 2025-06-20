program_block_EI <- function(rna_unit,met_unit,path_dataset) { 
  # prop1 = l_pred$prop1
  # prop2 = l_pred$prop2

  # prop = Reduce(`+`, list(prop1,prop2)) / length(list(prop1,prop2))
  met_unit[["ref_scRNA"]] = rna_unit$ref_scRNA
  return(met_unit)
}
