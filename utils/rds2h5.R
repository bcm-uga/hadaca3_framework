source("utils/data_processing.R")

domain = 'pdac'

# l_datasets = c('invitro', 'invivo', 'insilicopseudobulk', 'insilicodirichletEMFA', 'insilicodirichletCopule')

l_datasets = c( 'insilicodirichletCopule')
# l_datasets = c( 'insilicodirichletEMFA', 'insilicodirichletCopule')


path_old_data = 'old_datasets/'
path_H5_data = 'data/' 

for (dataset in l_datasets){
    print(dataset)
    mix = paste0("mixes1_",dataset,"_",domain)
    groundtruth = paste0("groundtruth1_",dataset,"_",domain)


    r_mix = readRDS(paste0(path_old_data,"mixes/filtered",mix,'.rds'))
    r_groundtruth = readRDS(paste0(path_old_data,"groundtruth/",groundtruth,'.rds'))
    
    
    f_mix_h5 = paste0(path_H5_data,mix,'.h5') 
    write_mix_hdf5(f_mix_h5,r_mix)


    f_ground_h5 = paste0(path_H5_data,groundtruth,'.h5') 
    write_global_hdf5(f_ground_h5,list(groundtruth=r_groundtruth))
    

} 

