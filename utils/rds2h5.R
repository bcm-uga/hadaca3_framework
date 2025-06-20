source("utils/data_processing.R")

domain = 'pdac'

# l_datasets = c('invitro', 'invivo', 'insilicopseudobulk', 'insilicodirichletEMFA', 'insilicodirichletCopule')

# l_datasets = c("insilicodirichletNoDep", "insilicodirichletNoDep4CTsource", "insilicodirichletNoDep6CTsource",  "insilicodirichletEMFAImmuneLowProp" )

l_datasets = c('invitro', 'invivo', 'insilicopseudobulk', 'insilicodirichletEMFA', 'insilicodirichletCopule',
"insilicodirichletNoDep", "insilicodirichletNoDep4CTsource", "insilicodirichletNoDep6CTsource",  "insilicodirichletEMFAImmuneLowProp" )

# "insilicodirichletNoDep",  "SDN5"
# "insilicodirichletNoDep4CTsource",  "SDN4"
# "insilicodirichletNoDep6CTsource",  "SDN6"
# "insilicodirichletEMFAImmuneLowProp",  "SDEL"


# l_datasets = c( 'insilicodirichletCopule')
# l_datasets = c( 'insilicodirichletEMFA', 'insilicodirichletCopule')


path_old_data = 'data/old_datasets/'
path_H5_data = 'data/' 

for (dataset in l_datasets){
    print(dataset)
    # mix = paste0("mixes1_",dataset,"_",domain)
    # groundtruth = paste0("groundtruth1_",dataset,"_",domain)
    
    mix = paste0("mixes2_",dataset,"_",domain)
    groundtruth = paste0("groundtruth2_",dataset,"_",domain)


    # r_mix = readRDS(paste0(path_old_data,"mixes/filtered",mix,'.rds'))
    # r_groundtruth = readRDS(paste0(path_old_data,"groundtruth/",groundtruth,'.rds'))

    r_mix = readRDS(paste0(path_old_data,mix,'.rds'))
    r_groundtruth = readRDS(paste0(path_old_data,groundtruth,'.rds'))
    
    
    f_mix_h5 = paste0(path_H5_data,mix,'.h5') 
    write_mix_hdf5(f_mix_h5,r_mix)


    f_ground_h5 = paste0(path_H5_data,groundtruth,'.h5') 
    write_global_hdf5(f_ground_h5,list(groundtruth=r_groundtruth))
    

} 

