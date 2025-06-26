program_block_PP <- function(data, path_og_dataset='', omic='') {
  
  warning("This method uses a single-cell signature from team H")
  avg_expression_df <- read.csv("teamH_scSignature.csv", row.names = 1)
  
  if (is.list(data)) {
    data <- list("ref_teamH"=list(counts=SeuratObject::as.sparse(avg_expression_df),
                                  metadata=data.frame(cell_type=colnames(avg_expression_df),
                                                      sample=NA)))
  } else {data = data[intersect(rownames(avg_expression_df), rownames(data)),]}
 
  return(data) 
}
