

source("utils/data_processing.R")
source(cleaner)

ref = read_all_ref_hdf5(reference_file)


ref$ref_bulkRNA <- preprocess_matrix(ref$ref_bulkRNA)
ref$ref_met <- preprocess_matrix(ref$ref_met)
ref$ref_scRNA <- lapply(ref$ref_scRNA, function(x) {
  list(counts=as(preprocess_matrix(x$counts),'dgCMatrix'),
        metadata=x$metadata)})          


write_all_ref_hdf5(output_file,ref)