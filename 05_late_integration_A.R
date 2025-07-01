
source(utils_script)

source(script_file);


path_og_dataset= list(mix =path_ogmix,ref = path_ogref )

if(input_file_met=='none_unit'){
    pred_met= NULL
}else{
    pred_met = read_hdf5(input_file_met)$pred
}


if(input_file_rna=='none_unit'){
    pred_RNA = NULL
}else {
    pred_RNA = read_hdf5(input_file_rna)$pred
}



pred= program_block_li(prop1 = pred_RNA,prop2 = pred_met,path_dataset = path_og_dataset )  

if(length(pred) == 0){ stop(paste("Late integration result is empty, script_file = ",script_file) ) }

write_global_hdf5(output_file, list(pred=pred))