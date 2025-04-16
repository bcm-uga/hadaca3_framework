
program_block_FS <- function(data,path_og_dataset='') {
   
   drop_null_ref_cols <- function(ref_matrix){
    non_null_rows = apply(ref_matrix != 0,MARGIN = 1, FUN = any, simplify = TRUE)
    return(ref_matrix[non_null_rows,])
   }

  # ref_matrix = multi_data$ref$ref_bulkRNA
  # multi_data$ref$ref_bulkRNA <- drop_null_ref_cols(ref_matrix)
  
  data <- drop_null_ref_cols(data)

  return(data) 
}

# FS_data = program_blockFS(multi_data = PP_data) 