program_block_FS <- function(data,path_og_dataset='') {
    library("Seurat")

    if (is.list(data)) {
      sc = data$ref_cluster
    } else {sc = read_hdf5(path_og_dataset$ref)$ref_scRNA$ref_cluster}

    sc_markers = FindAllMarkers(sc$seurat_clustered, assay = NULL, features = NULL,
                              logfc.threshold = 0.1, test.use = "wilcox", slot = "data")
    sc_markers = sc_markers[which(sc_markers$p_val_adj < 0.05 & sc_markers$pct.1>0.6 & sc_markers$pct.2<0.3), ]$gene

    if (is.list(data)) {
      data = lapply(data, function(x) list(counts = x$counts[sc_markers,], metadata = x$metadata))
    } else {data = data[sc_markers,]}

  return(data) 
}
