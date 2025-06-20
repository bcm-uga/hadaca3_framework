program_block_PP <- function(data,path_og_dataset='',omic='') {
  
#  multi_data = list(mix = mix,
#              ref = ref)
  require(Seurat)
    cluster <- function(ds) {
      ds <- Seurat::NormalizeData(ds, verbose=F)
      ds <- Seurat::FindVariableFeatures(ds, verbose=F)
      ds <- Seurat::ScaleData(ds)
      ds <- Seurat::RunPCA(ds, verbose=F) 
      ds <- Seurat::FindNeighbors(ds, reduction = "pca", dims = 1:10, verbose=F) 
      ds <- Seurat::FindClusters(ds, resolution = 0.5)
      return(ds)
    }
    
    sc_list <- lapply(seq_along(data), function(x) {
      # the next line has issue : 
      sc_dataset <- CreateSeuratObject(counts = data[[x]]$counts, project = names(data)[x])
      sc_dataset$cell_type <- data[[x]]$metadata$cell_type
      sc_dataset$subjectname <- paste0(names(data)[x],data[[x]]$metadata$sample)
      sc_dataset <- cluster(sc_dataset)
    })
    sc_merge <- Reduce(merge, sc_list)
    sc_merge <- cluster(sc_merge)
    
    sc_integrated <- IntegrateLayers(object = sc_merge , method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = F)
    sc_integrated <- JoinLayers(sc_integrated)
    data = list(ref_integrated = list(counts=sc_integrated@assays$RNA$counts,
                                           metadata=data.frame(cell_type=sc_integrated$cell_type,
                                                               sample=sc_integrated$subjectname,
                                                               dataset=unlist(mapply(function(a,b) {rep(a,b)}, names(ref_scRNA), sapply(ref_scRNA, function(x) ncol(x$counts)))))))
  return (data)
} 