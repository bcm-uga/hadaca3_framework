

source(script_file);
pred_RNA = readRDS(input_file_rna)
pred_met = readRDS(input_file_met)
output= program_block( c(pred_RNA,pred_met)) 

# library(arrow)

saveRDS(output, output_file)