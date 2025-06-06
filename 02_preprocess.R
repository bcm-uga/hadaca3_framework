source(utils_script)
source(script_file)



path_og_dataset= list(mix =mixes_file,ref = reference_file )

omic_name = omic2list_name[[omic]]
# omic_name = omic

## Read either mix or ref. 
split_result <- strsplit(mixes_file, "/")[[1]]
if (any(split_result != 'none')){
    data = read_mix_hdf5(mixes_file)[[omic_name]]
}else{
    data = read_all_ref_hdf5(reference_file,to_read = omic_name)[[omic_name]]
}



ppblock = program_block_PP(data,path_og_dataset)

if(length(ppblock) == 0){
    stop(paste("result data in Preprocess is empty, script_file = ",script_file ))
}

res =list()
res[[omic_name]] =  ppblock
# write_global(output_file ,ppblock ) 
write_global_hdf5(output_file ,res ) 
# write_all_hdf5(output_file ,ppblock ) 
# write_global(output_file ,ppblock ) 