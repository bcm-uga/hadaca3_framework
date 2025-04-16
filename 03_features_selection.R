source("utils/data_processing.R")



omic_name = omic2list_name[[get_omic(output_file)]]
path_og_dataset= list(mix =path_ogmix,ref = path_ogref )


data = read_hdf5(input_file)[[omic_name]]



source(script_file)
block = program_block_FS(data,path_og_dataset)

assert(length(block) != 0, paste("result data in Features selection is empty, script_file = ",script_file ) )


res =list()
res[[omic_name]] =  block
# write_global(output_file ,ppblock ) 
write_global_hdf5(output_file ,res ) 
# write_all_hdf5(output_file,block)