program_block_PP <- function(data,path_og_dataset='',omic='') {
  
  
  # Sequence depth normalization for RNA
    seq_depth_normalization <- function(mat) {
      sweep(mat, 2, colSums(mat), "/") * 10^6
    }
    mix_rna_normed <- seq_depth_normalization(mix_rna)
    ref_bulkRNA_normed <- seq_depth_normalization(ref_bulkRNA)
    ref_scRNA_normed <- lapply(ref_scRNA, function(x) {
      list(counts=seq_depth_normalization(x$counts),
           metadata=x$metadata)})
    
    # Now scale the rows so ref and mix are on the same scale
    shared_genes <- intersect(Reduce(intersect, lapply(ref_scRNA_normed, function(x) rownames(x$counts))),
                              intersect(rownames(mix_rna_normed), rownames(ref_bulkRNA_normed)))
    mix_rna_normed = mix_rna_normed[shared_genes,]
    ref_bulkRNA_normed = ref_bulkRNA_normed[shared_genes,]
    ref_scRNA_normed <- lapply(ref_scRNA, function(x) {
      list(counts=x$counts[shared_genes,],
           metadata=x$metadata)})
    row_mean <- rowMeans(mix_rna_normed) + rowMeans(ref_bulkRNA_normed) + rowMeans(sapply(ref_scRNA_normed, function(x) rowMeans(x$counts))) + 1e-4
    mix_rna_normed <- sweep(mix_rna_normed, 1, row_mean, '/')
    ref_bulkRNA_normed <- sweep(ref_bulkRNA_normed, 1, row_mean, '/')
    ref_scRNA_normed <- lapply(ref_scRNA_normed, function(x)
      list(counts=sweep(x$counts, 1, row_mean, '/'),
           metadata=x$metadata))
    
    mix_rna=mix_rna_normed
    ref_bulkRNA=ref_bulkRNA_normed
    ref_scRNA=ref_scRNA_normed

  return(data) 
}