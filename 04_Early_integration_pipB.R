

source(utils_script)

mix = read_hdf5(input_file_mix)$mix_met
ref = read_hdf5(input_file_met)$ref_met
# scRNA= read_hdf5(input_file_sc)$ref_scRNA
path_og_dataset= list(mix =path_ogmix,ref = path_ogref )

met_unit = list(mix= mix,ref =ref)#,ref_scRNA=  scRNA  )

source(script_de_met)
pred = program_block_DE(met_unit,path_og_dataset)

if(length(pred) == 0){
stop(paste("prediction  in Pipeline A MET unit is empty, script_file = ",script_de_met ) )
} 

write_global_hdf5(output_file,list(pred=pred))
