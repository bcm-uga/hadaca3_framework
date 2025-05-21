
# TODO
program_block_FS <- function(data,path_og_dataset='') {
    
    
    
    if (!("ref_cluster" %in% names(ref_scRNA))) {stop("This FS method requires to run the PP with option_sc set to sc_cluster")} 
    # Wilcoxon test one cluster vs other cluster 
    sc_markers = FindAllMarkers(ref_scRNA$ref_cluster$seurat_clustered, assay = NULL, features = NULL,
                          logfc.threshold = 0.1, test.use = "wilcox", slot = "data")
    sc_markers = sc_markers[which(sc_markers$p_val_adj < 0.05 & sc_markers$pct.1>0.6 & sc_markers$pct.2<0.3), ]$gene
    
    # select common genes with single cells
    # mix_rna = mix_rna[sc_markers,]
    # ref_bulkRNA = ref_bulkRNA[sc_markers,]

    if(is.list(data)){
      data = lapply(data, function(x) list(counts = x$counts[sc_markers,], metadata = x$metadata))
    }else{
      data = data[sc_markers,]
    }

 
  return(data) 
}

