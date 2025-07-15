program_block_FS <- function(data, path_og_dataset='') {

    if (!is.list(data)) {
    sc = read_hdf5(path_og_dataset$ref)
    } else {sc = data}
    
    #if (!(any(c("ref_concat","ref_integrated","ref_cluster","ref_binarypseudobulk_log") %in% names(sc)))) {stop("This FS method requires to run the PP set to concat, CCAintegration, cluster or binarypseudobulk_log")}
    ### SPLS-DA on sc
    sc_data = sc[[1]]
    splsda.model <- mixOmics::mint.splsda(t(sc_data$counts), sc_data$metadata$cell_type, 
                                          study = sc_data$metadata$dataset, ncomp = 5,
                                          keepX = rep(400,5))
    choose_markers_scRNA <- unique(c(mixOmics::selectVar(splsda.model, comp = 1)$name,
                                     mixOmics::selectVar(splsda.model, comp = 2)$name,
                                     mixOmics::selectVar(splsda.model, comp = 3)$name,
                                     mixOmics::selectVar(splsda.model, comp = 4)$name,
                                     mixOmics::selectVar(splsda.model, comp = 5)$name))
    if (is.list(data)) {
        data = lapply(data, function(x) list(counts = x$counts[choose_markers_scRNA,], metadata = x$metadata))
    } else {data = data[choose_markers_scRNA,]}
    
    return(data)
}

