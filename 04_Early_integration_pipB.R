

source(utils_script)

mix_met = read_hdf5(input_file_mix_met)$mix_met
ref_met = read_hdf5(input_file_met)$ref_met

met_unit = list(mix= mix_met,ref =ref_met)#,ref_scRNA=  scRNA  )

mix_rna = read_hdf5(input_file_mix_rna)$mix_rna
ref_rna = read_hdf5(input_file_rna)$ref_bulkRNA
scRNA= read_hdf5(input_file_sc)$ref_scRNA

rna_unit = list(mix= mix_rna,ref =ref_rna,ref_scRNA=  scRNA  )


path_og_dataset= list(mix =path_ogmix,ref = path_ogref )


source(script_file)
integrated_unit = program_block_EI(rna_unit,met_unit,path_og_dataset)

if(length(integrated_unit) == 0){
    stop(paste("Early integration is empty, script_file = ",script_file ) )
} 

write_global_hdf5(output_file,integrated_unit)
