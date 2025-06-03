program_block_PP <- function(data,path_og_dataset='') {
  

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
    
    network = readRDS("baselines/attachement/teamHtfrna_network_modules.rds")
    
    ref_scRNA_metadata = lapply(ref_scRNA, function(x) x$metadata)
    net = decoupleR::get_collectri(organism = "human", split_complexes = F)
    mix_rna = ADImpute::NormalizeTPM(mix_rna, log = T)
    ref_scRNA = lapply(ref_scRNA, function(x) ADImpute::NormalizeTPM(x$counts, log = T))
    mix_rna = compute.TFs.activity(mix_rna, universe = net)
    ref_scRNA = lapply(ref_scRNA, function(x) compute.TFs.activity(x, universe = net))
    mix_rna = create_tfs_modules(mix_rna, network)
    ref_scRNA = lapply(ref_scRNA, function(x) create_tfs_modules(x, network))
    mix_rna = t(minMax(mix_rna))
    ref_scRNA = lapply(ref_scRNA, function(x) t(minMax(x)))
    ref_scRNA = lapply(seq_along(ref_scRNA), function(x)
      list(counts = ref_scRNA[[x]], metadata = ref_scRNA_metadata[[x]]))
    
    ref_bulkRNA = readRDS("baselines/attachement/teamHtfrna_ref_modules.rds")
  return(data) 
}