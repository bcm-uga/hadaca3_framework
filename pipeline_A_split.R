

source("utils/data_processing.R")
source(script_split)


multi_data = read_all_hdf5(input_file)
split = program_blockSP(multi_data)

if (!exists("script_de_rna")) {
    #We decovolua the met
    source(script_de_met)
    pred = program_block(split$met)

}else{

    source(script_de_rna)
    pred = program_block(split$RNA)

}

# print(ppblock)

# library(arrow)

write_global_hdf5(output_file,list(pred=pred))