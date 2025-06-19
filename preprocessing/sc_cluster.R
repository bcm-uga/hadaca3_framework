program_block_PP <- function(data,path_og_dataset='',omic='') {
  
    library(Seurat)
    #Error is coming from this line. 
    seurat_single = lapply(data, function(x) CreateSeuratObject(x$counts, meta.data=x$metadata))
    
    # Merge all
    seurat_merge = Reduce(merge, seurat_single)
    DefaultAssay(seurat_merge) <- "RNA"
    seurat_merge <- JoinLayers(seurat_merge, overwrite = T)
    seurat_merge@meta.data$dataset = unlist(lapply(seq_along(data), function(x)
      rep(names(data)[x], ncol(data[[x]]$counts))))
    
    options(future.globals.maxSize = 1024 * 1024 * 1024)  # 1 GB
    
    process_scrna = function(obj) {
      obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000, verbose=F) 
      obj <- ScaleData(object = obj, features = rownames(obj@assays$RNA)) 
      obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose=F)
      obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca", verbose=F)
      obj <- FindClusters(obj, resolution = 0.5, cluster.name = "clust_res.5")
      obj@meta.data[['clust_res.5']] = Idents(obj)
      Idents(obj) = obj@meta.data$cell_type[match(Idents(obj), obj@meta.data$clust_res.5)]  
      return(obj)
    }
    seurat_merge = process_scrna(seurat_merge)
    data = list("ref_cluster"=list(counts=seurat_merge@assays$RNA$counts,
                                        metadata=data.frame(cell_type=seurat_merge$cell_type,
                                                            sample=seurat_merge$sample,
                                                            dataset=seurat_merge$dataset)),
                  "seurat_clustered"=seurat_merge)
  return(data) 
}