


source("utils/data_processing.R")
source(script_file);


# pred_RNA = readRDS(input_file_rna)
# pred_met = readRDS(input_file_met)

pred_RNA = read_hdf5(input_file_rna)$pred
pred_met = read_hdf5(input_file_met)$pred





output= program_block( list(pred_RNA,pred_met)) 

# library(arrow)

write_global_hdf5(output_file, list(pred=output))