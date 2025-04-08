source("utils/data_processing.R")
source(cleaner)

mix = read_mix_hdf5(mixes_file)
  
mix$mix_met <- preprocess_matrix(mix$mix_met)
mix$mix_rna <- preprocess_matrix(mix$mix_rna)

write_mix_hdf5(output_file,mix)