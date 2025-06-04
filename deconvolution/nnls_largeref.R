

program_block_DE <- function(uni_data,path_og_dataset='') {


      # test in case we have an extra cell type in the reference
  seuil_res = 1
  for (col_to_delete in (1:length(colnames(uni_data$ref)))) {
    ref_met_missing <-uni_data$ref[ ,-col_to_delete]
    # nnls with truncated ref
    prop <- apply(uni_data$mix, 2, function(b, A) {
      tmp_prop <- nnls::nnls(b = b, A = A)$x
      tmp_prop[tmp_prop < 0] <- 0
      tmp_prop <- tmp_prop / sum(tmp_prop) # Sum To One
      return(tmp_prop)}, A = ref_met_missing)
    residu <- mean(prop)
    if (residu < seuil_res) {
      seuil_res <- residu
      n_mod <- col_to_delete
      best_mod <- prop
    }
  }
  rownames(best_mod) <- colnames(uni_data$ref)[-n_mod]
  best_mod <- rbind(best_mod, rep(0, ncol(best_mod)))
  rownames(best_mod)[ncol(uni_data$ref)] <- colnames(uni_data$ref)[n_mod]
  
  prop = apply(uni_data$mix, 2, function(b, A) {
    tmp_prop = nnls::nnls(b=b,A=A)$x  
    tmp_prop = tmp_prop / sum(tmp_prop)
    return(tmp_prop)}, A=uni_data$ref)
  rownames(prop) = colnames(uni_data$ref)
  
  rmse_trunc = sqrt(mean((uni_data$mix - uni_data$ref %*% best_mod)^2))
  rmse = sqrt(mean((uni_data$mix - uni_data$ref %*% prop)^2))
  
  if (rmse_trunc < rmse) {prop <- best_mod}

  
  return(prop) 
}

