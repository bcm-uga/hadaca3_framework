program_block_FS <- function(data, path_og_dataset='') {
  
  if (is.list(data)) {
    sc = data
  } else {sc = read_hdf5(path_og_dataset$ref)$ref_scRNA}
  
  #if (!(any(c("ref_concat","ref_integrated","ref_cluster","ref_binarypseudobulk_log") %in% names(sc)))) {stop("This FS method requires to run the PP set to concat, CCAintegration, cluster or binarypseudobulk_log")}
  
  # select markers based on t statistic on sc data
  metadata = sc[[1]]$metadata
  counts = sc[[1]]$counts
  
  t_statistics <- lapply(combn(unique(metadata$cell_type), 2, simplify = F), function(cts) {
      lapply(unique(metadata$dataset), function(ds) {
        idx_col = metadata$dataset == ds & metadata$cell_type %in% cts
        # if (!any(idx_col) ) {
        #   # message(cts," ", ds ,"NULL")
        #   return(NULL)
        # }
        counts_trunc = counts[, which(idx_col)]

      if (is.null(dim(counts_trunc)) || any(dim(counts_trunc) == 0)) {
        return(NULL)
      }

        nb_celltype = length(unique(metadata[colnames(counts_trunc), "cell_type"]))

        apply(counts_trunc, 1, function(x) {

          if (nb_celltype > 1 & length(unique(x) > 0)) {
              res <- try(t.test(x ~ metadata[colnames(counts_trunc), "cell_type"])$statistic, silent = T)
            if (class(res) == "try-error") {return(0)}
            return(res)
          } else {return(NULL)}
        })
      })
    })


  n_top_genes = 20
  top_genes <- unique(do.call(c, lapply(t_statistics, function(x) {
    x <- x[!sapply(x, is.null)]
    Reduce(intersect, lapply(x, function(y) {
      y <- y[is.finite(y)]
      y <- y[order(y)]
      c(names(head(y, n_top_genes)), names(tail(y, n_top_genes)))}))
  })))
  
  res <- aggregate(as.data.frame(t(as.matrix(counts[top_genes, ]))),
                   list(metadata[, "cell_type"]),
                   mean)
  rownames(res) <- res[, 1]
  res <- t(res[-1])
  
  if (is.list(data)) {
    data = lapply(data, function(x) list(counts = x$counts[rownames(res),], metadata = x$metadata))
  } else {
    # which(rownames(res) %in% rownames(data))
    data = data[which(rownames(res) %in% rownames(data)),]

    # data = data[rownames(res),]
    }
  
  return(data) 
}

