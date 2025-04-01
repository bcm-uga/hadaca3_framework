source("utils/data_processing.R")

# if (!exists("groundtruth_file")) {groundtruth_file = "01_groundtruth/groundtruth1_insilicodirichletEMFA_pdac.rds"} 
# if (!exists("prediction_file"))  {prediction_file = "02_prediction/prediction1_nnls_insilicodirichletEMFA_pdac.rds"}   
# if (!exists("score_file"))     {score_file = "03_prediction/score1_nnls_insilicodirichletEMFA_pdac.rds"}   

# Scoring functions - generate the worst RMSE/MAE possible to normalize metrics

weighgeomMean <- function(x,w) {
  prod(mapply(function(a,b) a^b, a=x, b=w), na.rm=T)^(1/sum(w))
}

#########################
# homogenization function to find the best match between real and predicted A matrices for unsupervised methods. -USELESS
homogeneized_cor_mat = function(A_real, A_pred) {
    cmat = cor(t(A_real),t(A_pred))
    pvec <- c(clue::solve_LSAP((1+cmat)^2,maximum = TRUE))
    return(A_pred[pvec,])
}

#########################
# Pre-treatment of predicted A
#not useful because reference based algorithm
prepare_A <- function(A_real, A_pred) {
  N <- ncol(A_real)
  K <- nrow(A_real)
    
  stopifnot(K > 1)
  stopifnot(ncol(A_pred) == N)
  stopifnot(!anyNA(A_pred))
    
  ### STEP 1 : matching the number of predicted components to the real number of cell types K
    
  ## if predicting too few cell types
  if (nrow(A_pred) < K) {
    #set random positive values closed to 0 for missing rows
    set.seed(1)
    random_data = abs(jitter(matrix(data = 0, nrow = K - nrow(A_pred), ncol = N), factor = 0.01))
    A_pred <- rbind(A_pred, random_data)
    print("Add rows of 0s to match K")
  }
           
  ## if predicting too many cell types
  
  ### strategy 1: keep the most abundant cell types
  if (nrow(A_pred) > K) {
    A_pred <- A_pred[order(rowSums(A_pred), decreasing = TRUE)[1:K],]
    print("Number of predicted cell types exceeds K, only the K most abundant cell types are kept for scoring.")
  }
  
  ### strategy 2: clustering of similar cell types to match K components -> SLIM
        
  ### STEP 2 : reordering predicted components to find the best match -> SLIM
  A_pred_best_cor = as.matrix(homogeneized_cor_mat(A_real, A_pred))
  
  return(A_pred_best_cor)
}

#########################
# Global Pearson/Spearman correlation coefficients
correlationP_tot = function(A_real, A_pred) {
  if (all(c(var(A_pred))==0)) {return(-1)} #worst case scenario where all cell types are predicted with the same proportions
  return(cor(c(A_real), c(A_pred), method = "pearson"))
}
correlationS_tot = function(A_real, A_pred) {
  if (all(c(var(A_pred))==0)) {return(-1)} #worst case scenario where all cell types are predicted with the same proportions
  return(cor(c(A_real), c(A_pred), method = "spearman"))
}

#########################
# Mean column/sample Pearson/Spearman correlation coefficients
correlationP_col = function(A_real, A_pred) {
  res = c()
  for (i in seq(ncol(A_real))) {
    if (sd(A_pred[, i]) > 0 & sd(A_real[, i]) > 0) {
      res[i] = cor(A_real[, i], A_pred[, i], method = "pearson")
    }
  }
  res = res[!is.na(res)]
  print(paste0(length(res), " columns are kept for correlation analysis"))
  if (length(res)==0) {return(-1)}
  return(mean(res))
}
correlationS_col = function(A_real, A_pred) {
  res = c()
  for (i in seq(ncol(A_real))) {
    if (sd(A_pred[, i]) > 0 & sd(A_real[, i]) > 0) {
      res[i] = cor(A_real[, i], A_pred[, i], method = "spearman")
    }
  }
  res = res[!is.na(res)]
  print(paste0(length(res), " columns are kept for correlation analysis"))
  if (length(res)==0) {return(-1)}
  return(mean(res))
}

#########################
# Mean row/cell type Pearson/Spearman correlation coefficients
correlationP_row = function (A_real, A_pred) {
  res = c()
  for (i in seq(nrow(A_real))) {
    if (sd(A_pred[i, ]) > 0 & sd(A_real[i, ]) > 0) {
      res[i] = cor(A_real[i, ], A_pred[i, ], method = "pearson")
    }
  }
  res = res[!is.na(res)]
  print(paste0(length(res), " rows are kept for correlation analysis"))
  if (length(res)==0) {return(-1)}
  return(mean(res))
}
correlationS_row = function (A_real, A_pred) {
  res = c()
  for (i in seq(nrow(A_real))) {
    if (sd(A_pred[i, ]) > 0 & sd(A_real[i, ]) > 0) {
      res[i] = cor(A_real[i, ], A_pred[i, ], method = "spearman")
    }
  }
  res = res[!is.na(res)]
  print(paste0(length(res), " rows are kept for correlation analysis"))
  if (length(res)==0) {return(-1)}
  return(mean(res))
}

#########################
# Aitchison distance
eval_Aitchison = function(A_real, A_pred, min = 1e-9) {
  # Aitchison dist doesn't like 0 
  A_real[A_real < min] = min
  A_pred[A_pred < min] = min
  
  # Aitchison apply a transformation to manage compositional data, so only sample by sample comparison have sens
  res = c()
  for (i in seq(ncol(A_real))) {
    res[i] = coda.base::dist(rbind(t(A_real[,i]), t(A_pred[,i])), method = "aitchison")[1]
  }
  res = res[!is.na(res)]
  return(mean(res))
}

#########################
# SDID with a twist : i don't change the sign as described in MPRA paper
# https://mpra.ub.uni-muenchen.de/84387/
eval_SDID = function(A_real, A_pred) {
  cos_angle = sapply(seq(ncol(A_real)), function(x)
    A_real[,x]%*%A_pred[,x]/norm(A_real[,x],type="2")/norm(A_pred[,x],type="2"))
  sin_angle = sapply(cos_angle, function(x) sqrt(1-min(1,x^2)))
  sdid = sqrt(sin_angle)
  return(mean(sdid))
}

#########################
# angular distance
eval_AID = function(A_real, A_pred) {
  angle_radian = sapply(seq(ncol(A_real)), function(x)
    acos(min(1,(A_real[,x]%*%A_pred[,x])/norm(A_real[,x],type="2")/norm(A_pred[,x],type="2"))))
  angle_degree = angle_radian * 90 / (pi/2)
  return(mean(angle_degree))
}

#########################
# Jensen Shannon divergence
eval_JSD = function(A_real, A_pred, min = 1e-9) {
  jsd = sapply(seq(ncol(A_real)), function(x) {
    pred = A_pred[,x]
    pred[pred < min] = min
    n <- 0.5 * (A_real[,x] + pred)
    JS <- 0.5 * (sum(A_real[,x] * log(A_real[,x] / n)) + sum(pred * log(pred / n)))})
  return(mean(jsd))
}

#########################
# RMSE
eval_RMSE = function(A_real, A_pred) {
  return(sqrt(mean((A_real - A_pred)^2)))
}

#########################
# MAE
eval_MAE = function (A_real, A_pred){
  return(mean(abs(A_real - A_pred)))
}

#########################
# Scoring function 
scoring_function <- function(A_real, A_pred) {
  # pre-treatment of predicted A
  #A_pred = prepare_A(A_real = A_real, A_pred = A_pred)
  
  # scoring with different metrics
  if (nrow(A_pred)==nrow(A_real)) {
    
    # Match cell types
    if (all(rownames(A_pred) %in% rownames(A_real))) {
      A_pred = A_pred[rownames(A_real),]
    }
    else {stop(paste0("Cell types names do not match:\n Predicted cell types are ",paste(rownames(A_pred), collapse=', '),"\n Real cell types are ",paste(rownames(A_real), collapse=', ')))}
    mae = eval_MAE(A_real, A_pred)
    rmse = eval_RMSE(A_real, A_pred)
    aitchison = eval_Aitchison(A_real, A_pred)
    sdid = eval_SDID(A_real, A_pred)
    aid = eval_AID(A_real, A_pred)
    jsd = eval_JSD(A_real, A_pred)
    pearson_tot = correlationP_tot(A_real, A_pred)
    pearson_col = correlationP_col(A_real, A_pred)
    pearson_row = correlationP_row(A_real, A_pred)
    spearman_tot = correlationS_tot(A_real, A_pred)
    spearman_col = correlationS_col(A_real, A_pred)
    spearman_row = correlationS_row(A_real, A_pred)
  }
  else if (nrow(A_pred)==(nrow(A_real)+1)) { # In case of an extra cell type in the reference matrix compared to the source/ground truth
    # Match cell types
    if (all(rownames(A_real) %in% rownames(A_pred))) {
      extra_type = matrix(0,ncol=ncol(A_real),nrow=1)
      rownames(extra_type) = setdiff(rownames(A_pred),rownames(A_real))
      A_real = rbind(A_real,extra_type)
      A_pred = A_pred[rownames(A_real),]
    }
    else {stop(paste0("Cell types names do not match:\n Predicted cell types are ",paste(rownames(A_pred), collapse=', '),"\n Real cell types are ",paste(rownames(A_real), collapse=', ')))}
    mae = eval_MAE(A_real, A_pred)
    rmse = eval_RMSE(A_real, A_pred)
    aitchison = eval_Aitchison(A_real, A_pred)
    sdid = eval_SDID(A_real, A_pred)
    aid = eval_AID(A_real, A_pred)
    jsd = eval_JSD(A_real, A_pred)
    pearson_tot = correlationP_tot(A_real, A_pred)
    pearson_col = correlationP_col(A_real, A_pred)
    pearson_row = correlationP_row(A_real, A_pred)
    spearman_tot = correlationS_tot(A_real, A_pred)
    spearman_col = correlationS_col(A_real, A_pred)
    spearman_row = correlationS_row(A_real, A_pred)
  }
  else if (nrow(A_pred)==(nrow(A_real)-1)) { # In case of a missing cell type in the reference matrix compared to the source/ground truth
    # Match cell types
    if (all(rownames(A_pred) %in% rownames(A_real))) {
      A_real = A_real[rownames(A_pred),]
    }
    else {stop(paste0("Cell types names do not match:\n Predicted cell types are ",paste(rownames(A_pred), collapse=', '),"\n Real cell types are ",paste(rownames(A_real), collapse=', ')))}
    mae = eval_MAE(A_real, A_pred)
    rmse = eval_RMSE(A_real, A_pred)
    aitchison = eval_Aitchison(A_real, A_pred)
    sdid = eval_SDID(A_real, A_pred)
    aid = eval_AID(A_real, A_pred)
    jsd = eval_JSD(A_real, A_pred)
    pearson_tot = correlationP_tot(A_real, A_pred)
    pearson_col = correlationP_col(A_real, A_pred)
    pearson_row = correlationP_row(A_real, A_pred)
    spearman_tot = correlationS_tot(A_real, A_pred)
    spearman_col = correlationS_col(A_real, A_pred)
    spearman_row = correlationS_row(A_real, A_pred)
  }
  else if (nrow(A_pred) > nrow(A_real) & setequal(rownames(A_real), c("basal",'classic'))) { # partial ground truth only for the in vivo dataset
    rmse = NA
    mae = NA
    aitchison = NA
    sdid = NA
    aid = NA
    jsd = NA
    pearson_tot = NA
    pearson_col = NA
    pearson_row = correlationP_row(A_real, A_pred[rownames(A_real),])
    spearman_tot = NA
    spearman_col = NA
    spearman_row = correlationS_row(A_real, A_pred[rownames(A_real),])
    weights_spec = c(1/2,1/2,rep(0,10))
  }
  all_judges = data.frame("pearson_row"=pearson_row,
                          "spearman_row"=spearman_row,
                          "pearson_tot"=pearson_tot,
                          "pearson_col"=pearson_col,
                          "spearman_tot"=spearman_tot,
                          "spearman_col"=spearman_col,
                          "rmse"=rmse,
                          "mae"=mae,
                          "aitchison"=aitchison,
                          "jsd"=jsd,
                          "sdid"=sdid,
                          "aid"=aid)
  
  # generate best/worst possible metrics
  fake_worst_pred = apply(A_real, 2, function(prop) { 
    tmp = rep(1e-9, length(prop))
    tmp[which.min(prop)] = 1
    return(tmp)})
  judge_candidate = rbind(all_judges,
                          data.frame("pearson_row"=1,
                                     "spearman_row"=1,
                                     "pearson_tot"=1,
                                     "pearson_col"=1,
                                     "spearman_tot"=1,
                                     "spearman_col"=1,
                                     "rmse"=0,
                                     "mae"=0,
                                     "aitchison"=0,
                                     "jsd"=0,
                                     "sdid"=0,
                                     "aid"=0),
                          data.frame("pearson_row"=-1,
                                     "spearman_row"=-1,
                                     "pearson_tot"=-1,
                                     "pearson_col"=-1,
                                     "spearman_tot"=-1,
                                     "spearman_col"=-1,
                                     "rmse"=eval_RMSE(A_real, fake_worst_pred),
                                     "mae"=eval_MAE(A_real, fake_worst_pred),
                                     "aitchison"=eval_Aitchison(A_real, fake_worst_pred), #verify it's the correct worst distance
                                     "jsd"=eval_JSD(A_real, fake_worst_pred), #verify it's the correct worst divergence
                                     "sdid"=1,
                                     "aid"=90 #verify it's the correct worst angle in multi-D
                                     ))
  
  # normalize scores with a linear shift between [0,1]
  # strategy when normalizing only one method
  CenterScaleNorm <- function(x) {
    tr1 = x - min(x, na.rm=T)
    tr2 = tr1/(max(tr1, na.rm=T))
    return(tr2)
  }
  judge_candidate_norm = apply(judge_candidate, 2, CenterScaleNorm)
  
  # transform scores s.t. 1 is the best score
  judge_candidate_norm = 1 - judge_candidate_norm
  judge_candidate_norm[,grep("pearson",colnames(judge_candidate_norm))] = 1 - judge_candidate_norm[,grep("pearson",colnames(judge_candidate_norm))]
  judge_candidate_norm[,grep("spearman",colnames(judge_candidate_norm))] = 1 - judge_candidate_norm[,grep("spearman",colnames(judge_candidate_norm))]

  # Average over judges with the geometric mean for the candidate of interest
  #score_aggreg = exp(mean(log(judge_candidate_norm[1,]),na.rm=T))
  weights = c(rep(1/3*1/2,2),
              rep(1/3*1/4,4),
              rep(1/3*1/4,4),
              rep(1/3*1/2,2))
  if (nrow(A_pred) > nrow(A_real) & setequal(rownames(A_real), c("basal",'classic'))) {weights = weights_spec}
  score_aggreg = weighgeomMean(judge_candidate_norm[1,],
                               weights)
  all_scores = c(score_aggreg,
                 pearson_row,spearman_row,
                 pearson_tot,pearson_col,spearman_tot,spearman_col,
                 rmse,mae,aitchison,jsd,sdid,aid,
                 judge_candidate_norm[1,])
  names(all_scores) = c("score_aggreg",
                       rep(c("pearson_row","spearman_row",
                 "pearson_tot","pearson_col","spearman_tot","spearman_col",
                 "rmse","mae","aitchison","jsd","sdid","aid"),times=2))
  names(all_scores)[(2+length(c("pearson_row","spearman_row",
                 "pearson_tot","pearson_col","spearman_tot","spearman_col",
                 "rmse","mae","aitchison","jsd","sdid","aid"))):length(all_scores)] = paste0(c("pearson_row","spearman_row",
                 "pearson_tot","pearson_col","spearman_tot","spearman_col",
                 "rmse","mae","aitchison","jsd","sdid","aid"),'_norm')
  return(all_scores)
}

# Run Scoring


prediction = read_hdf5(paste0(prediction_file))$pred
ground_truth = read_hdf5(paste0(groundtruth_file))$groundtruth

A_real = as.matrix(ground_truth)
A_pred = as.matrix(prediction)

scores_full = scoring_function(A_real=A_real,  A_pred=A_pred)
scores = as.numeric(scores_full)
names(scores) = names(scores_full)



write_global_hdf5(paste0(score_file),scores)

