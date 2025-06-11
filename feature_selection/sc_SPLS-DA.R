
program_block_FS <- function(data,path_og_dataset='') {
   
    if (!(any(c("ref_concat","ref_integrated","ref_cluster","ref_binarypseudobulk_log") %in% names(ref_scRNA)))) {stop("This FS method requires to run the PP with option_sc set to concat, CCAintegration, cluster or binarypseudobulk_log")}
    ### SPLS-DA on sc
    sc_data = ref_scRNA[[1]]
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

