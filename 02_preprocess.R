source("utils/data_processing.R")
source(script_file)



path_og_dataset= list(mix =mixes_file,ref = reference_file )

omic_name = omic2list_name[[get_omic(output_file)]]

if (mixes_file!= ''){
    data = read_mix_hdf5(mixes_file)[[omic_name]]
}else{
    data = read_all_ref_hdf5(reference_file,to_read = omic_name)[[omic_name]]
}


print(mixes_file)
print(reference_file)

ppblock = program_block_PP(data,path_og_dataset)

assert(length(ppblock) != 0, paste("result data in Preprocess is empty, script_file = ",script_file ) )


res =list()
res[[omic_name]] =  ppblock
# write_global(output_file ,ppblock ) 
write_global_hdf5(output_file ,res ) 
# write_all_hdf5(output_file ,ppblock ) 
# write_global(output_file ,ppblock ) 