program_block_PP <- function(data,path_og_dataset='') {
  
    warning("maybe integrate the data first...")
    # calculate average profile for each cell type from the single cell data
    common_genes = Reduce(intersect, lapply(data, function(x) rownames(x$counts)))
    data = lapply(data, function(sc) {
      list(counts = sc$counts[common_genes,], metadata = sc$metadata)})
    ref_scRNA_centro = lapply(data, function(sc)
      t(pbapply::pbapply(sc$counts, 1, function(x)
        tapply(x, sc$metadata$cell_type, mean))))
    cell_types = unique(unlist(lapply(ref_scRNA_centro, colnames)))
    ref_scRNA_centro = sapply(cell_types, function(type) {
      list_sc = lapply(ref_scRNA_centro, function(sc) {
        if (type %in% colnames(sc)) {sc[,type]}})
      list_sc = list_sc[!sapply(list_sc,is.null)]
      Reduce("+", list_sc) / length(list_sc)})

      
    ref_scRNA_centro = ref_scRNA_centro[,colnames(ref_bulkRNA)]
    
    # scale to have the same library size as ref_bulkRNA
    ref_bulkRNA = ref_bulkRNA[intersect(rownames(ref_bulkRNA),rownames(ref_scRNA_centro)),]


    ref_scRNA_centro = ref_scRNA_centro[rownames(ref_bulkRNA),]
    median_ref = median(apply(ref_bulkRNA, 2, sum))
    dn = dimnames(ref_scRNA_centro)
    ref_scRNA_centro = ref_scRNA_centro %*% diag(median_ref / colSums(ref_scRNA_centro))
    dimnames(ref_scRNA_centro) = dn
    data = list(ref_pseudobulk = list(counts=ref_scRNA_centro,
                                           metadata=data.frame(cell_type=colnames(ref_scRNA_centro),
                                                               sample='average_pseudobulk')))
  return(data) 
}