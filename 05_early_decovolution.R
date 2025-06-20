
source(utils_script)

# not implemented but should not be usefulll ? or ? 
# path_og_dataset= list(mix =path_ogmix,ref = path_ogref )
path_og_dataset = ""


source(script_de)

uni_data = read_hdf5(path_uni_data)

pred = program_block_DE(uni_data,path_og_dataset)

if(length(pred) == 0){
    stop(paste("prediction  in Pipeline A RNA unit is empty, script_file = ",script_de ) )
} 
# assert(length(pred) != 0, paste("prediction  in Pipeline A RNA unit is empty, script_file = ",script_de_rna ) )


write_global_hdf5(output_file,list(pred=pred))