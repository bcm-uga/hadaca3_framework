# args <- commandArgs(trailingOnly = TRUE)

# print(length(args))
# print(args)


# print(reference_file)
# # print(mixes_file)
# if (!exists("reference_file"))      {reference_file = "data/reference_pdac.rds"} 
# if (!exists("mixes_file"))          {mixes_file = "data/mixes_test_dataset.rds"} 



# # script_file = script_file
# mix = readRDS(mixes_file)
# ref = readRDS(reference_file)

data = readRDS(input_file)

source(script_file)

block = program_block(data)

# print(ppblock)

# library(arrow)

saveRDS(block, output_file)