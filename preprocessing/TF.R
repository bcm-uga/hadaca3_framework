program_block_PP <- function(data,path_og_dataset='',omic='') {
  
  # og_ref_scrna  =  read_all_ref_hdf5(path_og_dataset$ref,to_read = 'ref_scRNA')$ref_scRNA

    warning("This method uses a priori knowledge from team H")
    library(decoupleR)
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
      library(dplyr)
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
    




    if(omic == 'ref_scRNA'){
      net = decoupleR::get_collectri(organism = "human", split_complexes = F)
      
      ref_scRNA_metadata = lapply(data, function(x) x$metadata)
      data = lapply(data, function(x) ADImpute::NormalizeTPM(x$counts, log = T))
      data = lapply(data, function(x) compute.TFs.activity(x, universe = net))
      data = lapply(data, function(x) create_tfs_modules(x, network))
      data = lapply(data, function(x) t(minMax(x)))
      data = lapply(seq_along(data), function(x)
        list(counts = data[[x]], metadata = ref_scRNA_metadata[[x]]))
      # data = lapply(data, function(x) list(counts = x$counts[hvg,], metadata = x$metadata))
    
    }else if(omic == 'ref_bulkRNA'){
      ref_bulkRNA = readRDS("teamHtfrna_ref_modules.rds")
    }else{ 
      net = decoupleR::get_collectri(organism = "human", split_complexes = F)

      # network = readRDS("preprocessing/attachement/teamHtfrna_network_modules.rds")
      network = readRDS("teamHtfrna_network_modules.rds")
      data = ADImpute::NormalizeTPM(data, log = T)
      data = compute.TFs.activity(data, universe = net)
      data = create_tfs_modules(data, network)
      data = t(minMax(data))
    }
    
  return(data) 
}