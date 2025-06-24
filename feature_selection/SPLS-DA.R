program_block_FS <- function(data,path_og_dataset='') {
  
  og_ref_met  =  read_all_ref_hdf5(path_og_dataset$ref,to_read = 'ref_met')$ref_met
  
  highly_var_idx <- TOAST::findRefinx(og_ref_met, nmarker=1e4)
  og_ref_met = og_ref_met[highly_var_idx,]
  sda <- mixOmics::splsda(t(log(og_ref_met/(1-og_ref_met))), Y=colnames(og_ref_met),
                    keepX = c(1000,1000))
  choose_markers_met <- unique(c(mixOmics::selectVar(sda, comp = 1)$name,
                                  mixOmics::selectVar(sda, comp = 2)$name))
  data = data[choose_markers_met,]
  
  return(data) 
}

