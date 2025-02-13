


multi_data = readRDS(input_file)

source(script_split)
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

saveRDS(pred, output_file)