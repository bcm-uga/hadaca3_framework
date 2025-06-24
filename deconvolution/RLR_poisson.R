

program_block_DE <- function(uni_data,path_og_dataset='') {

  idx_feat = intersect(rownames(uni_data$mix), rownames(uni_data$ref))
  uni_data$mix = uni_data$mix[idx_feat,]
  uni_data$ref = uni_data$ref[idx_feat,]

  get_weights_poisson <- function(ref) {
    weights = 1/apply(ref, 1, mean)
    return(weights/sum(weights))
  }
  compute_rlr_weighted <- function(beta.m, ref.m, weights) {
    est.m <- matrix(nrow = ncol(beta.m), ncol = ncol(ref.m))
    colnames(est.m) <- colnames(ref.m)
    rownames(est.m) <- colnames(beta.m)
    for (s in seq_len(ncol(beta.m))) {
      rlm.o <- MASS::rlm(beta.m[, s] ~ ref.m, maxit = 50, weights = weights)
      coef.v <- summary(rlm.o)$coef[2:(ncol(ref.m) + 1), 1]
      weight_rlm <- data.frame(GeneNames = rownames(beta.m), w = rlm.o$w)
      relevant_weights <- weight_rlm$w[weight_rlm$w ==1]
      kept_genes <- weight_rlm$GeneNames[weight_rlm$w ==1]
      # normalisation: remove negative coefficients and enforce the unit-simplex constraint
      coef.v[which(coef.v < 0)] <- 0
      total <- sum(coef.v)
      coef.v <- coef.v/total
      est.m[s, ] <- coef.v
    }
    return(t(est.m))
  }
  
  prop = compute_rlr_weighted(beta.m = uni_data$mix, ref.m = uni_data$ref, weights = get_weights_poisson(uni_data$ref))

  rownames(prop) <- colnames(uni_data$ref)


  return(prop) 
}

