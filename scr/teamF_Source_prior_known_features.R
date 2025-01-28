#' Compute Weights for RNA Data Based on Poisson Dispersion
#'
#' @param ref_bulkRNA A reference matrix for RNA data, where rows represent genes 
#'                    and columns represent samples or cell types.
#'
#' @return A numeric array of weights for genes. Each weight corresponds to a 
#'         gene in the `ref_bulkRNA` matrix and is normalized to sum to one.
#'
#' @details The weights are calculated as the reciprocal of the mean expression 
#'          levels for each gene (1/mean), based on the assumption of Poisson 
#'          dispersion. This approach gives more weight to genes with lower 
#'          mean expression, which are often considered to have higher 
#'          variability in Poisson models.
#'
get_weights_poisson <- function(ref_bulkRNA) {
  ## Compute weights 1/var, or 1/mean given poisson dispersion hypothesis
  weights = 1.0/apply(ref_bulkRNA, MARGIN = 1, FUN = mean)
  weights = weights/sum(weights)
  
  return(weights)
}

#' Estimate Proportion of Cell Types Using NNLS Method
#'
#' @param mix A matrix of mixed signals where rows represent features (e.g., genes) 
#'            and columns represent samples.
#' @param ref A reference matrix where rows represent the same features as in 
#'            the `mix` matrix, and columns represent pure cell types. Each column 
#'            is the expression profile of a specific cell type.
#'
#' @return A matrix where rows correspond to the cell types in the `ref` matrix, 
#'         and columns correspond to the samples in the `mix` matrix. Each entry 
#'         represents the estimated proportion of a cell type in a given sample.
#'
get_prop_nnls <- function(mix, ref) {
  prop = apply(mix, 2, function(b, A) {
    tmp_prop = nnls::nnls(b=b, A=A)$x
    tmp_prop = tmp_prop / sum(tmp_prop) # Sum To One
    return(tmp_prop)
  }, A=ref)
  rownames(prop) <- colnames(ref)
  
  return(prop)
}

#' Estimate Proportion of Cell Types Using RLR Method
#'
#' @param mix A matrix of mixed signals where rows represent features (e.g., genes) 
#'            and columns represent samples.
#' @param ref A reference matrix where rows represent the same features as in 
#'            the `mix` matrix, and columns represent pure cell types. Each column 
#'            is the expression profile of a specific cell type.
#' @param weights An optional vector of weights for the robust linear regression. 
#'                If `NULL`, equal weights are used. Default is `NULL`.
#'
#' @return A matrix where rows correspond to the samples in the `mix` matrix, and 
#'         columns correspond to the cell types in the `ref` matrix. Each entry 
#'         represents the estimated proportion of a cell type in a given sample.
#'         
get_prop_rlr <- function(mix, ref, weights = NULL) {
  rlr_epidish <- compute_epidish_rpc (beta.m = mix, ref.m = ref, maxit = 50, weights = weights)
  prop <- t(rlr_epidish$estimated_ratios)
  rownames(prop) <- colnames(ref)
  
  return(prop)
}

#' Deconvolove using a possibly weighted linear model
#'
#' with post hoc enforcement of frequency constraints
#' @param mix 
#' @param ref 
#' @param weights weights for the linear model fitting
#'
#' @return A matrix where rows correspond to the samples in the `mix` matrix, and 
#'         columns correspond to the cell types in the `ref` matrix. Each entry 
#'         represents the estimated proportion of a cell type in a given sample.
#' @export
#'
#' @examples
get_prop_lm <- function(mix, ref, weights = NULL){
  # Create unity weights by default
  if(is.null(weights)){
    weights<-seq(1.0,nrow(mix))
  }
  # Feed it to a weighted LM fit
  idx_feat = intersect(rownames(mix), rownames(ref))
  prop = apply(mix[idx_feat,], 2, function(b, A) {
    tmp_prop = lm.wfit(y=b, x=A, w = weights, offset = 0)$coefficients
    tmp_prop[tmp_prop < 0] = 0
    tmp_prop = tmp_prop / sum(tmp_prop) # Sum To One
    return(tmp_prop)
  }, A=ref[idx_feat,])  
  rownames(prop) = colnames(ref)
}

#' Normalize a Matrix by Column Sums
#'
#' @param mat A numeric matrix where rows represent features and columns represent 
#'            samples.
#'
#' @return A numeric matrix with the same dimensions as the input `mat`, where each 
#'         column has been normalized so that its elements sum to one.
#'
normalize_matrix <- function(mat) {
  mat = sweep(mat, 2, colSums(mat), "/")
  return(mat)
}


#' Calculate the Root Mean Squared Error (RMSE) Between Real and Estimated Matrices
#'
#' @param Y A matrix containing the real data. Rows represent features, and 
#'          columns represent samples.
#' @param Y_estimate A matrix containing the estimated data, with the same 
#'                   dimensions as `Y`.
#' @param weights An optional vector of weights to be applied to the squared errors. 
#'                It should have a length equal to the number of rows in `Y`.
#'                If not provided, equal weights are used.
#'
#' @return A vector representing the normalized RMSE score for every sample.
#'
rmse_mat <- function(Y=NULL, Y_estimate=NULL, weights=NULL) {
  n_sample = ncol(Y)
  rmse = rep(0, n_sample)
  
  if(is.null(weights)){
    # Use unity weights by default
    weights=rep(1.0, nrow(Y))
  }
  
  # compute rmse for each sample
  for (i in 1:n_sample)
  {
    rmse[i] = sqrt(weighted.mean(x = (Y[, i] - Y_estimate[, i])^2, w = weights))
  }
  
  # normalize the rmse values : BUGGED
  #range_rmse = max(rmse) - min(rmse)
  #if (range_rmse != 0)
  #{
  #  rmse = (rmse - min(rmse)) / (max(rmse) - min(rmse))
  #}
  
  return(rmse)
}


#' Copy-paste the EPIdish rpc score, relyong on robust linear regression
#' In EPIdish, they simply use the `MASS::rlm()` function with standard parameters, assuming an intercept and 
#' post-normalising the cellular ratios to enforce the unit-simplex constraint. 
#'
#' @param beta.m bulk matrix (transcriptomic or methylation)
#' @param ref.m reference matrix (transcriptomic or methylation. For single cell, expected pseudobulk)
#' @param weights weights for the different genes
#' @param maxit if convergence is not reached, maximal allowed number of iterations
#' @param verbose If TRUE, return the proportion of genes expressed. 
#' @param plot_weight_density Plot the weight density distribution?
#' @param ggplot_title Name of the dataset or experiment?
#' 
#' @return A list with `estimated_ratios` the cellular ratios, 
#' `kept_genes` a list of genes considered as expressed and 
#' `weight_ggplot` the ggplot2 distribution of weights after ILWS iteration. 

compute_epidish_rpc <- function(beta.m, ref.m, weights = NULL, maxit = 50, verbose = FALSE,
                                plot_weight_density = FALSE, ggplot_title ="Dataset name") {
  
  est.m <- matrix(nrow = ncol(beta.m), ncol = ncol(ref.m))
  colnames(est.m) <- colnames(ref.m)
  rownames(est.m) <- colnames(beta.m)
  for (s in seq_len(ncol(beta.m))) {
    rlm.o <- MASS::rlm(beta.m[, s] ~ ref.m, maxit = maxit, weights = weights)
    # retrieve coefficients, except the intercept
    coef.v <- summary(rlm.o)$coef[2:(ncol(ref.m) + 1), 1]
    
    weight_rlm <- data.frame(GeneNames = rownames(beta.m), 
                             w = rlm.o$w) # for ggplot2 tidy formats
    relevant_weights <- weight_rlm$w[weight_rlm$w ==1] # keep the genes observed
    kept_genes <- weight_rlm$GeneNames[weight_rlm$w ==1]
    
    if (verbose) {
      message(paste(format(100 * length(relevant_weights)/nrow(weight_rlm), digits = 3), 
                    "% of genes exhbit a weight of 1, and will be considered as expressed."))
    }
    
    # report to https://r-charts.com/distribution/histogram-binwidth-ggplot2/ for higher customisation
    # plot weight density distribution
    if (plot_weight_density) {
      weight_ggplot <- ggplot(weight_rlm, mapping = aes(x=w)) +         
        geom_histogram(aes(y = after_stat(density)), bins = 100) +
        geom_density(colour = "red", linewidth = 1) +
        ggtitle("Density plot and histogram of the posterior weights computed by RLR") +
        labs(x = "Weights", y = "Density distribution of weights") +
        annotate("text", x = 0.5, y = 40,
                 label = paste(format(100 * relevant_weights/nrow(weight_rlm), digits = 3), 
                               "% of genes exhbit a weight of 1."),
                 color = "blue", size = 4, fontface = "bold") +
        
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
        )
    }
    else {
      weight_ggplot <- NULL
    }
    
    # normalisation: remove negative coefficients and enforce the unit-simplex constraint
    coef.v[which(coef.v < 0)] <- 0
    total <- sum(coef.v)
    coef.v <- coef.v/total
    est.m[s, ] <- coef.v
  }
  return(list(estimated_ratios = est.m, 
              kept_genes = kept_genes, 
              weight_ggplot = weight_ggplot))
}


#' Enhanced RLM version
#'
#' Same parameters as \code{\link{compute_epidish_rpc}}. 
#' Recursive call to the function to remove absent cell types (ones associated with negative scores).
#' Stringent threshold to assign null weights to clearly outliers or irrelevant or noisy genes. 
#'
#' @inheritParams compute_epidish_rpc

compute_epidish_enhanced <- function(beta.m, ref.m, maxit = 50, verbose = FALSE) {
  
  est.m <- matrix(nrow = ncol(beta.m), ncol = ncol(ref.m))
  colnames(est.m) <- colnames(ref.m)
  rownames(est.m) <- colnames(beta.m)
  for (s in seq_len(ncol(beta.m))) {
    # print(paste("we are at", s, "iteration."))
    coef.v.restrained <- .compute_epidish_enhanced_internal (y = beta.m[, s],
                                                             X = ref.m,
                                                             maxit = maxit)
    
    # normalisation: enforce the unit-simplex constraint
    # ensure that removed cell populations are associated with a null value
    coef.v <- stats::setNames(rep(0, ncol(ref.m)), 
                              colnames(ref.m))
    coef.v[names(coef.v.restrained)] <- coef.v.restrained[names(coef.v.restrained)]
    #  to avoid null divisions, we can't remain on the simplex
    if (all(coef.v==0)) {
      coef.v <- stats::setNames(rep(1/ncol(ref.m), ncol(ref.m)),
                                colnames(ref.m))
    }
    
    total <- sum(coef.v)
    coef.v <- coef.v/total
    est.m[s, ] <- coef.v
  }
  return(list(estimated_ratios = est.m))
}


.compute_epidish_enhanced_internal <- function(y, X, maxit = maxit) {
  # new model assuming no intercept (alternatively, we could use anova or fisher-test)
  # and stronger and more stringent results
  tryCatch(
    {
      rlm.o <- suppressWarnings(MASS::rlm(y ~ X -1, 
                                          maxit = maxit, 
                                          method = "MM", 
                                          psi = psi.bisquare))
      # rlm.o <- suppressWarnings(MASS::rlm(y ~ x -1,
      #                                     maxit = maxit,
      #                                     psi = psi.hampel, 
      #                                     k = c(1, 2, 3)))  # Hampel tuning
      # retrieve coefficients
      coef.v <- stats::setNames(summary(rlm.o)$coef[, 1],
                                colnames(X))
    },
    error = function(e) {
      # Handle the error in case of bad convergence
      # Assume this case equal proportions -> better to use another RLM approach.
      message("An error occurred: ", e$message)
      coef.v <- stats::setNames(rep(1/ncol(X), ncol(X)),
                                colnames(X))
      # browser()  # Enter debugging mode for error investigation
    }
  )
  
  # if (all(coef.v==0)) {
  #   browser()
  # }
  
  if ((any(coef.v < 0)) & (length(coef.v) >=2)) {
    existing_cell_populations <- which(coef.v >0 )
    X <- X[,which(coef.v >0 )] # keep only cell populations present in the sample
    coef.v <- .compute_epidish_enhanced_internal (y = y,
                                                  X = X,
                                                  maxit = maxit)
  }
  # coef.v <- stats::setNames(coef.v, colnames(X))
  return(coef.v)
}

