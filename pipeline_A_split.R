

source("utils/data_processing.R")
source(script_split)

# multi_data contains  : $mix and $ref 
#  mix contains $mix_rna and $mix_met
# ref contains $ref_bulkRNA $ref_met $ref_scRNA

multi_data = read_all_hdf5(input_file)
split = program_blockSP(multi_data)

# split contains $RNA and $met
# split$RNA contains mix and ref 
#split$RNA$ref contains bulk and scRNA
# split$met contains

if (!exists("script_de_rna")) {
    #We decovolua the met
    source(script_de_met)
    pred = program_block_DE(split$met)

}else{

    source(script_de_rna)
    pred = program_block_DE(split$RNA)

}

# print(ppblock)

# library(arrow)

write_global_hdf5(output_file,list(pred=pred))