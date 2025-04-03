


source("utils/data_processing.R")
source(script_file);


# pred_RNA = readRDS(input_file_rna)
# pred_met = readRDS(input_file_met)

# print(input_need)

# prior_knowledge = list()

# if( input_needed != 'None'){
#     input_needed = strsplit(input_needed,' +')
#     prior_knowledge = read_all_hdf5(last_dataset,input_needed)

# }

pred_RNA = read_hdf5(input_file_rna)$pred
pred_met = read_hdf5(input_file_met)$pred

output= program_block_li( list(prop1= pred_RNA,prop2 = pred_met,last_dataset=last_dataset )  ) 

# library(arrow)

write_global_hdf5(output_file, list(pred=output))