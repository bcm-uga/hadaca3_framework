program_block_FS <- function(data,path_og_dataset='') {
    

    library("Seurat")

    if(is.list(data)){
      
      if (!("seurat_clustered" %in% names(data$ref_cluster))) {stop("This FS method requires the reference sc_cluster")} 
      
      sc_markers = FindAllMarkers(data$ref_cluster$seurat_clustered, assay = NULL, features = NULL,
                          logfc.threshold = 0.1, test.use = "wilcox", slot = "data")
      sc_markers = sc_markers[which(sc_markers$p_val_adj < 0.05 & sc_markers$pct.1>0.6 & sc_markers$pct.2<0.3), ]$gene
    

      data = lapply(data, function(x) list(counts = x$counts[sc_markers,], metadata = x$metadata))
    }else{

    og_ref = read_hdf5(path_og_dataset$ref)$ref_scRNA$ref_cluster

    if (!("seurat_clustered" %in% names(og_ref))) {stop("This FS method  requires the reference sc_cluster")} 
    
    # Wilcoxon test one cluster vs other cluster 
    sc_markers = FindAllMarkers(og_ref$seurat_clustered, assay = NULL, features = NULL,
                          logfc.threshold = 0.1, test.use = "wilcox", slot = "data")
    sc_markers = sc_markers[which(sc_markers$p_val_adj < 0.05 & sc_markers$pct.1>0.6 & sc_markers$pct.2<0.3), ]$gene
    
    data = data[sc_markers[sc_markers %in% rownames(data)],]
    
    }

 
  return(data) 
}