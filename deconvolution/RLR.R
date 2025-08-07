

program_block_DE <- function(uni_data,path_og_dataset='') {

  library(caret)


  idx_feat = intersect(rownames(uni_data$mix), rownames(uni_data$ref))
  uni_data$mix = uni_data$mix[idx_feat,]
  uni_data$ref = uni_data$ref[idx_feat,]
  
  get_epidish_with_fallback <- function(beta.m, ref.m, cutoff = 0.99) {
    safe_epidish <- function(beta, ref) {
      tryCatch({
        return(t(EpiDISH::epidish(beta, ref, method = "RPC")$estF))
      }, error = function(e) {
        warning("EpiDISH failed, attempting to remove correlated columns...")

        cor_matrix <- cor(ref)
        to_remove <- caret::findCorrelation(cor_matrix, cutoff = cutoff)
        if (length(to_remove) == 0) stop("No columns removed, but EpiDISH still fails.")

        message("Removing correlated columns: ", paste(colnames(ref)[to_remove], collapse = ", "))
        ref_clean <- ref[, -to_remove, drop = FALSE]

        return(t(EpiDISH::epidish(beta, ref_clean, method = "RPC")$estF))
      })
    }

    return(safe_epidish(beta.m, ref.m))
  }

  prop <- get_epidish_with_fallback(uni_data$mix,uni_data$ref)

  # prop <- t(EpiDISH::epidish(uni_data$mix,uni_data$ref,method="RPC")$estF)

  return(prop) 
}

