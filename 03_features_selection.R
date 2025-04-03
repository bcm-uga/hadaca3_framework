source("utils/data_processing.R")
data = read_all_hdf5(input_file)

source(script_file)
block = program_block_FS(data)

write_all_hdf5(output_file,block)