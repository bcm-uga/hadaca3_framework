source("utils/data_processing.R")



omic_name = omic2list_name[[get_omic(output_file)]]


data = read_hdf5(input_file)[[omic_name]]



source(script_file)
block = program_block_FS(data)



res =list()
res[[omic_name]] =  block
# write_global(output_file ,ppblock ) 
write_global_hdf5(output_file ,res ) 
# write_all_hdf5(output_file,block)