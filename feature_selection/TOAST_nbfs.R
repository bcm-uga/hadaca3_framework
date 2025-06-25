program_block_FS <- function(data, path_og_dataset='') {
  
  og_ref_met = read_all_ref_hdf5(path_og_dataset$ref, to_read = 'ref_met')$ref_met
  nb_fs_met = nrow(og_ref_met)

  hvp <- TOAST::findRefinx(og_ref_met, nmarker = min(nrow(og_ref_met), nb_fs_met))

  return(data[hvp,]) 
}

