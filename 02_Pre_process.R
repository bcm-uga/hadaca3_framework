# args <- commandArgs(trailingOnly = TRUE)

# print(length(args))
# print(args)


# print(reference_file)
# print(mixes_file)
# if (!exists("reference_file"))      {reference_file = "data/reference_pdac.rds"} 
# if (!exists("mixes_file"))          {mixes_file = "data/mixes_test_dataset.rds"} 

source("utils/data_processing.R")
source(script_file)


# script_file = script_file
# mix = readRDS(mixes_file)
# ref = readRDS(reference_file)

mix = read_mix_hdf5(mixes_file)
ref = read_all_ref_hdf5(reference_file)


ppblock = program_block(mix,ref)

# print(ppblock)

# library(arrow)

write_all_hdf5(output_file ,ppblock ) 