

source("utils/data_processing.R")

mix = read_hdf5(input_file_mix)$mix_rna
ref = read_hdf5(input_file_rna)$ref_bulkRNA
scRNA= read_hdf5(input_file_sc)$ref_scRNA

path_og_dataset= list(mix =path_ogmix,ref = path_ogref )

rna_unit = list(mix= mix,ref =ref,ref_scRNA=  scRNA  )

source(script_de_rna)
pred = program_block_DE(rna_unit,path_og_dataset)
assert(length(pred) != 0, paste("prediction  in Pipeline A RNA unit is empty, script_file = ",script_de_rna ) )


write_global_hdf5(output_file,list(pred=pred))