source("utils/data_processing.R")
source(script_file)

mix = read_mix_hdf5(mixes_file)
ref = read_all_ref_hdf5(reference_file)

multi_data = list(mix = mix,
             ref = ref)

ppblock = program_block_PP(multi_data)

write_all_hdf5(output_file ,ppblock ) 