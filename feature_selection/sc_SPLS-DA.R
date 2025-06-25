program_block_FS <- function(data, path_og_dataset='') {
    
    og_ref_rna  =  read_all_ref_hdf5(path_og_dataset$ref, to_read = 'ref_met')$ref_met
    
    if (!(any(c("ref_concat","ref_integrated","ref_cluster","ref_binarypseudobulk_log") %in% names(data)))) {stop("This FS method requires to run the PP set to concat, CCAintegration, cluster or binarypseudobulk_log")}
    ### SPLS-DA on sc
    sc_data = data[[1]]
    splsda.model <- mixOmics::mint.splsda(t(sc_data$counts), sc_data$metadata$cell_type, 
                                          study = sc_data$metadata$dataset, ncomp = 5,
                                          keepX = rep(400,5))
    choose_markers_scRNA <- unique(c(mixOmics::selectVar(splsda.model, comp = 1)$name,
                                     mixOmics::selectVar(splsda.model, comp = 2)$name,
                                     mixOmics::selectVar(splsda.model, comp = 3)$name,
                                     mixOmics::selectVar(splsda.model, comp = 4)$name,
                                     mixOmics::selectVar(splsda.model, comp = 5)$name))
    mix_rna = mix_rna[choose_markers_scRNA,]
    ref_bulkRNA = ref_bulkRNA[choose_markers_scRNA,]
    ref_scRNA = lapply(ref_scRNA, function(x) list(counts = x$counts[choose_markers_scRNA,], metadata = x$metadata))
    
    return(data)
}

