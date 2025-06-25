program_block_FS <- function(data, path_og_dataset='') {
  
  nb_fs_rna = 1e3
  og_ref_bulkRNA = read_all_ref_hdf5(path_og_dataset$ref, to_read = 'ref_bulkRNA')$ref_bulkRNA
  
  if (is.list(data)) {
    og_ref_bulkRNA = lapply(data, function(x) og_ref_bulkRNA[intersect(rownames(og_ref_bulkRNA),
                                                                       rownames(x$counts)),])
    data = mapply(function(x,y) list(counts = x$counts[rownames(y),],
                                     metadata = x$metadata),
                  data, og_ref_bulkRNA, SIMPLIFY=F)
    hvg <- lapply(og_ref_bulkRNA, function(x)
      TOAST::findRefinx(x, nmarker = min(nrow(x), nb_fs_rna)))
  } else {
    og_ref_bulkRNA = og_ref_bulkRNA[intersect(rownames(og_ref_bulkRNA),
                                              rownames(data)),]
    data = data[rownames(og_ref_bulkRNA),]
    hvg <- TOAST::findRefinx(og_ref_bulkRNA, nmarker = min(nrow(og_ref_bulkRNA), nb_fs_rna))
  }

  if (is.list(data)) {
    data = mapply(function(x, y) list(counts = x$counts[y,], metadata = x$metadata),
                  data, hvg, SIMPLIFY=F)
  } else {data = data[hvg,]}
  
  return(data) 
}

