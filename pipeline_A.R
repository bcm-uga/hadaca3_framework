


multi_data = readRDS(input_file)

source(script_split)
source(script_de_met)
source(script_de_rna)



block = program_block(data)

# print(ppblock)

# library(arrow)

saveRDS(block, output_file)