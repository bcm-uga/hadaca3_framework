#' The function to estimate the A matrix
#' In the provided example, we use basic non-negative least squares (package "nnls"), which consists of minimizing the error term $||Mix - Ref \times Prop||^2$ with only positive entries in the prop matrix.
#' For methylation data, we source link_gene_CpG.R and probes.features.rds to use only CpG sites attached to a gene.
#'
#' @param mix_rna the list of bulk matrix associated to the transcriptome data set
#' @param mix_met the list of bulk matrix associated to the methylation data set
#' @param ref_bulkRNA the reference matrix associated to the transcriptome data set
#' @param ref_met the reference matrix associated to the methylation data set
#' @param ref_scRNA the reference list associated to the scRNA data set
#' 
#' @return the estimated A matrix
#' @examples
#' 
program <- function(mix_rna=NULL, ref_bulkRNA=NULL, 
                    mix_met=NULL, ref_met=NULL, ref_scRNA=NULL) {
  ##
  ## YOUR CODE BEGINS HERE
  ##
  
  
  source("baselines/attachement/teamF_Source_prior_known_features.R")
  
  ##############################################################################
  #                           Parameter tuning                                 #
  ##############################################################################
  
  rna_method = 'rlr' # or 'nnls'
  met_method = 'rlr' # or 'nnls'
  
  integration = 'rmse'
  
  weight = TRUE # weights for the rna data
  normalisation = TRUE 
  
  print(paste("rna method =", rna_method))
  print(paste("met_method =", met_method))
  print(paste("integration =", integration))
  print(paste("weighting =", weight))
  print(paste("normalisation =", normalisation))
  
  ##############################################################################
  #                           Preprocessing data                               #
  ##############################################################################
  
  # Normalize input matrices
  if (normalisation)
  {
    ref_bulkRNA = normalize_matrix(ref_bulkRNA)
    ref_met = normalize_matrix(ref_met)
  }
  
  ## Drop all 0 references
  drop_null_ref_cols <- function(ref_matrix){
    non_null_rows = apply(ref_matrix != 0,MARGIN = 1, FUN = any, simplify = TRUE)
    return(ref_matrix[non_null_rows,])
  }
  ref_bulkRNA <- drop_null_ref_cols(ref_bulkRNA)
  
  ##############################################################################
  #                 Proportion estimation for each omic                        #
  ##############################################################################
  
  ## RNA:
  if ( !( is.null(x = mix_rna) ) ) {
    
    idx_feat = intersect(rownames(mix_rna), rownames(ref_bulkRNA))
    mix_rna = mix_rna[idx_feat,]
    ref_bulkRNA = ref_bulkRNA[idx_feat,]
    
    # Compute gene weights for Poisson like regression
    bulkRNA_weights <- if (weight) get_weights_poisson(ref_bulkRNA) else NULL
    
    if (rna_method == 'nnls')
    {
      if (weight) warning("bulk RNA reweighting not availaible for nnls, use lmw or rlr. Weights will be used only for late RMSE fusion.")
      prop_rna = get_prop_nnls(mix_rna, ref_bulkRNA)
    }
    else if(rna_method == 'lmw'){
      prop_rna = get_prop_lm(mix_rna, ref_bulkRNA, bulkRNA_weights)  
    }
    else if (rna_method == 'rlr')
    {
      prop_rna = get_prop_rlr(mix_rna, ref_bulkRNA, bulkRNA_weights)
    }
    else
    {
      print("Parametrisation error : you have to indicate which method to use for bulk RNA.")
      prop_rna = NULL
    }
  }
  
  ## Methylation:
  if ( !( is.null(mix_met) ) ) {
    
    idx_feat = intersect(rownames(mix_met), rownames(ref_met))
    mix_met = mix_met[idx_feat,]
    ref_met = ref_met[idx_feat,]
    
    if (met_method == 'nnls')
    {
      prop_met = get_prop_nnls(mix_met, ref_met)
    }
    else if (met_method == 'rlr')
    {
      ## with robust linear regression approach ----
      prop_met = get_prop_rlr(mix_met, ref_met)
    }
    else
    {
      print("Parametrisation error : you have to indicate which method to use for methylation.")
      prop_met = NULL
    }
    
  }
  
  
  ##############################################################################
  #                           Late integration                                 #
  ##############################################################################
  
  ## we compute the mean with the RMSE weights of all the estimated A matrices
  ## as the final A matrix :
  
  if ( !is.null(x = mix_met) ) {
    if ( !is.null(x = mix_rna) ) {
      
      stopifnot( identical(x = dim(prop_rna), y = dim(prop_met)) ) # stop if not same number of cell type or samples
      
      if (integration == 'rmse')
      {
        ## compute RMSE for rna and methylation estimation
        rmse_rna = rmse_mat(Y = as.matrix(mix_rna),
                            Y_estimate = ref_bulkRNA %*% prop_rna,
                            weights = bulkRNA_weights)
        rmse_met = rmse_mat(as.matrix(mix_met), ref_met %*% prop_met)
        
        ## initialize with basic values we will modify
        prop <- (prop_rna + prop_met) / 2
        for (i in 1:length(rmse_met))
        {
          prop_val = (prop_rna[,i] / rmse_rna[i] + prop_met[,i] / rmse_met[i]) / (1 / rmse_rna[i] + 1 / rmse_met[i])
          prop[,i] = prop_val
        }
      }
      else if (integration == 'EnsDeconv') {
        source("baselines/attachement/teamF_CTS_EnsDeconv_LS.R")
        list_proportions <- list("RNA" = t(prop_rna),
                                 "Met" = t(prop_met))
        
        prop <- t(CTS_EnsDeconv_LS(What = list_proportions, maxit = 1000)$W)
      }
      
      else
      {
        prop <- (prop_rna + prop_met) / 2
      }
      
    } else {
      prop <- prop_met
    }
  } else {
    prop <- prop_rna
  }
  
  if (any(colSums(prop) != 1)) { # Sum To One 
    prop <- sapply(
      1:ncol(prop), 
      function(col_i) { prop[,col_i] / sum(prop[,col_i]) }
    )
  }
  
  
  ##############################################################################
  #                           Return result                                    #
  ##############################################################################
  
  return(prop)
  
  ##
  ## YOUR CODE ENDS HERE
  ##
}

