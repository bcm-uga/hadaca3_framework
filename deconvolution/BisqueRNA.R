

program_block_DE = function(uni_data,path_og_dataset='') {

  if (!("BisqueRNA" %in% installed.packages())) {
    # install.packages(pkgs = "BisqueRNA", repos = "https://cloud.r-project.org")
    devtools::install_github("cozygene/bisque")
  }
  library(BisqueRNA)
  if (is.null(metadata)) {stop("BisqueRNA with EI requires to do EI with the spls_sc option")}
  # generate bulk matrix
  bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(uni_data$mix))
  # generate reference sc matrix
  sce <- as.matrix(uni_data$ref)
  phenoData <- data.frame(SubjectName = metadata$sample,
               cellType = metadata$cell_type)
               
  rownames(phenoData) = rownames(metadata)
  sc.eset <- Biobase::ExpressionSet(assayData = sce,
                          phenoData = Biobase::AnnotatedDataFrame(phenoData))
  
  # do deconvolution
  prop <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=F)$bulk.props

 
}