
program_block_FS <- function(data,path_og_dataset='') {


  og_ref_bulkRNA  =  read_all_ref_hdf5(path_og_dataset$ref,to_read = 'ref_bulkRNA')$ref_bulkRNA

  warning("This method uses a single-cell signature from team H")
  avg_expression_df <- read.csv("baselines/attachement/teamH_scSignature.csv", row.names = 1)
  common_genes <- intersect(rownames(avg_expression_df), rownames(mix_rna))



  # mix_rna <- mix_rna[common_genes,]
  # og_ref_bulkRNA <- og_ref_bulkRNA[common_genes,]
          

  if(is.list(data)){
  data <- list("ref_teamH"=list(counts=as.matrix(avg_expression_df[common_genes,colnames(og_ref_bulkRNA)]),
                    metadata=data.frame(cell_type=colnames(og_ref_bulkRNA),
    # data = lapply(data, function(x) list(counts = x$counts[top_genes,], metadata = x$metadata))
  }else{
    data = data[common_genes,]
  }
 
  return(data) 
}

