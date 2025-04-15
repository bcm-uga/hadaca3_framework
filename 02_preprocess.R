source("utils/data_processing.R")
source(script_file)

# ref = read_all_ref_hdf5(reference_file)
# mix = read_mix_hdf5(mixes_file)
# write_global_hdf5('output/tmp2.h5',ref)

# multi_data = list(mix = mix,
#              ref = ref)


omic_name = omic2list_name[[get_omic(output_file)]]

if (mixes_file!= ''){
    data = read_mix_hdf5(mixes_file)[[omic_name]]
}else{
    data = read_all_ref_hdf5(reference_file,to_read = omic_name)[[omic_name]]
}


print(mixes_file)
print(reference_file)

ppblock = program_block_PP(data)

res =list()
res[[omic_name]] =  ppblock
# write_global(output_file ,ppblock ) 
write_global_hdf5(output_file ,res ) 
# write_all_hdf5(output_file ,ppblock ) 
# write_global(output_file ,ppblock ) 