program_block_FS <- function(data,path_og_dataset='') {
  library(dplyr)
  
  og_ref_met  =  read_all_ref_hdf5(path_og_dataset$ref,to_read = 'ref_met')$ref_met
  
  ## 1. For each probe, compute the absolute difference between the probe's methylation level in a cell type versus the global average in the matrix
  ## 2. Ranking probes for each cell type by the score computed above
  ## 3. Extract the top n features per cell type
  n_discriminant <- function(mat, n) {
    mean_ct <- rowMeans(mat)
    df_mat <- mat %>% apply(2, function(x) round((abs(x - mean_ct)), 4))
    colnames(df_mat) %>%
      lapply(function(x) {
        df_mat[, x] %>%
          sort(decreasing = TRUE) %>%
          head(n) %>% names()}) %>% Reduce(c, .)
  }
  ## Return the maximal nb of probes (up to n_max) s.t. there is no overlap across cell types
  least_n_discriminant <- function(mat, n_max) {
    nb_features <- c(1:n_max) %>% parallel::mclapply(function(i) {
      top_n_each <- n_discriminant(mat, i)
      return(max(table(top_n_each)))
    }, mc.cores = 12)
    least_n <- min(which(unlist(nb_features) > 1)) - 1
    return(unique(n_discriminant(mat, least_n)))
  }
  
  features_met <- least_n_discriminant(og_ref_met, 100)
  
  data = data[features_met,]

    return(data) 
}

