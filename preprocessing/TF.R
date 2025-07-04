program_block_PP <- function(data, path_og_dataset='', omic='') {
  
   if (omic == 'ref_bulkRNA') {
    data = readRDS("teamHtfrna_ref_modules.rds")
  }else{
    library(decoupleR)
  net = decoupleR::get_collectri(organism = "human", split_complexes = F)
  warning("This method uses a priori knowledge from team H, should be coupled with matching FS")

  compute.TFs.activity <- function(RNA.counts, universe, min_targets_size = 3) {
    tfs2viper_regulons <- function(df) {
      regulon_list <- split(df, df$source)
      regulons <- lapply(regulon_list, function(regulon) {
        tfmode <- stats::setNames(regulon$mor, regulon$target)
        list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
      })
      return(regulons)}
    net_regulons = tfs2viper_regulons(universe)
    sample_acts <- viper::viper(as.matrix(RNA.counts), net_regulons, minsize = min_targets_size, verbose=F, method = "scale")
    return(data.frame(t(sample_acts)))
  }
  create_tfs_modules = function(TF.matrix, network_tfs) {
    require(dplyr)
    tfs.modules = TF.matrix %>%
      t() %>%
      data.frame() %>%
      dplyr::mutate(Module = "na")
    for (i in 1:length(network_tfs[[3]])) {
      tfs.modules$Module[which(rownames(tfs.modules) %in% network_tfs[[3]][[i]])] = names(network_tfs[[3]])[i]
    }
    tfs_colors = tfs.modules %>%
      dplyr::pull(Module)
    MEList = WGCNA::moduleEigengenes(TF.matrix, colors = tfs_colors, scale = F) #Data already scale
    MEs = MEList$eigengenes
    MEs =  WGCNA::orderMEs(MEs)
    return(MEs)
  }
  minMax <- function(x) {
    #columns: features
    x = data.matrix(x)
    for (i in 1:ncol(x)) {
      x[,i] = (x[,i] - min(x[,i], na.rm = T)) / (max(x[,i], na.rm = T) - min(x[,i], na.rm = T))
    }
    return(x)
  }
  
  network = readRDS("teamHtfrna_network_modules.rds")

    if (is.list(data)) {
      metadata = lapply(data, function(x) x$metadata)
      data = lapply(data, function(x)
        counts = t(minMax(create_tfs_modules(compute.TFs.activity(ADImpute::NormalizeTPM(x$counts, log = T),
                                                                  universe = net),
                                            network))))
      data = lapply(seq_along(data), function(x)
        list(counts = SeuratObject::as.sparse(data[[x]]), metadata = metadata[[x]]))

    } else { 
      data = t(minMax(create_tfs_modules(compute.TFs.activity(ADImpute::NormalizeTPM(data, log = T),
                                                              universe = net),
                                          network)))
    }

  }
  
    
  return(data) 
}