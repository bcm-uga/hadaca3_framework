
# library(assertr)
library(rhdf5)
library(Matrix) 



write_all_ref_hdf5 <- function(path, ref_all) {
  # assert( all(colnames(ref_peng$counts)== rownames(ref_peng$metadata)))
  h5createFile(path)

  h5createGroup(path, "ref_bulkRNA")
  h5write(ref_all$ref_bulkRNA, file=path , name=paste0( "ref_bulkRNA/data"  ))
  h5write(colnames(ref_all$ref_bulkRNA), file=path , name ="ref_bulkRNA/cell_types"   )
  h5write(rownames(ref_all$ref_bulkRNA), file=path , name ="ref_bulkRNA/genes"   )



  h5createGroup(path, "ref_met")
  h5write(ref_all$ref_met, file=path , name=paste0( "ref_met/data"  ))
  h5write(colnames(ref_all$ref_met), file=path , name ="ref_met/cell_types"   )
  h5write(rownames(ref_all$ref_met), file=path , name ="ref_met/CpG_sites"   )


  h5createGroup(path, 'ref_scRNA/')

  for (dataset in c( "ref_sc_peng"    ,  "ref_sc_baron"   , "ref_sc_raghavan")){
    counts = ref_all$ref_scRNA[[dataset]]$counts
    meta =   ref_all$ref_scRNA[[dataset]]$metadata
    
    
    group =  paste0( 'ref_scRNA/' , dataset  )

    h5createGroup(path, group)


    h5createDataset(path, dataset =  paste0(group,'/data'), dims = length(counts@x) , storage.mode = "integer")
    h5write(counts@x, file=path, name=paste0(group, "/data"))
    h5write(dim(counts), file=path, name=paste0(group, "/shape"))
    h5createDataset(path, dataset = paste0(group,'/indices'), dims = length(counts@i), storage.mode = "integer")
    h5write(counts@i, file=path, name=paste0(group, "/indices"))
    h5write(counts@p, file=path, name=paste0(group, "/indptr"))

    if (!is.null(rownames(counts))) {
      h5write(rownames(counts), file=path, name=paste0(group, "/genes"))
    }

    if (!is.null(colnames(counts))) {
      h5write(colnames(counts), file=path, name=paste0(group, "/cell"))
    }

    h5write(meta, file=path , name=paste0( group,"/meta"  ))
  }

  return(NULL)
}

# file = 'sparse_matrix_R2.h5'

write_mix_hdf5 <- function(path, mix) {
  # assert( all(colnames(ref_peng$counts)== rownames(ref_peng$metadata)))
  
  h5createFile(path)

  h5createGroup(path, "mix_rna")
  h5write(mix$mix_rna, file=path , name=paste0( "mix_rna/data"  ))

  if(length(colnames(mix$mix_rna) )){
    h5write(colnames(mix$mix_rna), file=path , name ="mix_rna/samples" )
  }
  h5write(rownames(mix$mix_rna), file=path , name ="mix_rna/genes" )

  h5createGroup(path, "mix_met")
  h5write(mix$mix_met, file=path , name=paste0( "mix_met/data"  ))
  if(length(colnames(mix$mix_met) ) ){
    h5write(colnames(mix$mix_met), file=path , name ="mix_met/samples"   )
  }
  h5write(rownames(mix$mix_met), file=path , name ="mix_met/CpG_sites"   )

  return(NULL)
}

# mix_rds = readRDS("data/mixe_test.rds")
# write_mix_hdf5("data/mix_test.h5",mix_rds)


read_mix_hdf5 <- function(path) {
  # Read ref_bulkRNA data
  mix_rna <- h5read(path, "mix_rna/data")
  file_structure <- h5dump(path,load=FALSE)
  if(!is.null(file_structure$mix_rna$samples)){
    colnames(mix_rna) <- h5read(path, "mix_rna/samples")
  }
  rownames(mix_rna) <- h5read(path, "mix_rna/genes")

  # Read ref_met data
  mix_met <- h5read(path, "mix_met/data")
  if(!is.null(file_structure$mix_met$samples)){
    colnames(mix_met) <- h5read(path, "mix_met/samples")
  }
  rownames(mix_met) <- h5read(path, "mix_met/CpG_sites")

  mix <- list(
    mix_rna = mix_rna,
    mix_met = mix_met
  )
  return(mix)
}

# mix_h5 = read_mix_hdf5("data/mix_test.h5")
read_all_ref_hdf5 <- function(path,to_read=c('ref_bulkRNA','ref_met','ref_scRNA')) {
  # Read ref_bulkRNA data
  ref_bulkRNA = list()
  if('ref_bulkRNA' %in% to_read){
    ref_bulkRNA <- h5read(path, "ref_bulkRNA/data")
    colnames(ref_bulkRNA) <- h5read(path, "ref_bulkRNA/cell_types")
    rownames(ref_bulkRNA) <- h5read(path, "ref_bulkRNA/genes")
  }
    

  # Read ref_met data
  ref_met = list()
  if('ref_met' %in% to_read){
    ref_met <- h5read(path, "ref_met/data")
    colnames(ref_met) <- h5read(path, "ref_met/cell_types")
    rownames(ref_met) <- h5read(path, "ref_met/CpG_sites")
  }

  # Read ref_scRNA data
  ref_scRNA <- list()
  if('ref_scRNA' %in% to_read){
    datasets <- c("ref_sc_peng", "ref_sc_baron", "ref_sc_raghavan")
    for (dataset in datasets) {
      group <- paste0('ref_scRNA/', dataset)

      counts_data <- as.numeric(h5read(path, paste0(group, "/data")))
      counts_shape <- as.integer(h5read(path, paste0(group, "/shape")))
      counts_indices <- as.integer(h5read(path, paste0(group, "/indices")))
      counts_indptr <- as.integer(h5read(path, paste0(group, "/indptr")))

      cells = h5read(path, paste0(group, "/cell"))
      counts <- new("dgCMatrix",
                    x = counts_data,
                    i = counts_indices,
                    p = counts_indptr,
                    Dim = counts_shape,
                    Dimnames = list(
                      h5read(path, paste0(group, "/genes")),
                      cells
                    ))

      meta <- h5read(path, paste0(group, "/meta"))
      rownames(meta) = cells
      ref_scRNA[[dataset]] <- list(
        counts = counts,
        metadata = meta
      )
    }
  }

  # Combine all data into a single list
  ref_all <- list(
    ref_bulkRNA = ref_bulkRNA,
    ref_met = ref_met,
    ref_scRNA = ref_scRNA
  )

  return(ref_all)
}

read_all_hdf5 <- function(path,to_read=c('mix','ref')) {
  
  mix =list()
  if('mix' %in% to_read ){
    mix = read_mix_hdf5(path)
  }
  ref= list()
  if('ref' %in% to_read ){
    ref = read_all_ref_hdf5(path)
  }
  return (list(
    mix = mix, 
    ref = ref
  ))

}

write_all_hdf5 <- function(path, combined_data){
  write_mix_hdf5(path, combined_data$mix)
  write_all_ref_hdf5(path, combined_data$ref)
}



write_global_hdf5 <- function(path, data_list) {
  # Create the HDF5 file
  h5createFile(path)

  # Iterate over each named element in the data_list
  for (name in names(data_list)) {
    # Create a group for each named element
    h5createGroup(path, name)

    # Write the data, column names, and row names to the HDF5 file
    h5write(data_list[[name]], file = path, name = paste0(name, "/data"))
    if(length(colnames(data_list[[name]]))){
      h5write(colnames(data_list[[name]]), file = path, name = paste0(name, "/samples"))
    }
    if(length(rownames(data_list[[name]]))){
      h5write(rownames(data_list[[name]]), file = path, name = paste0(name, "/genes"))
    }
  }
}

read_hdf5 <- function(path) {
  # Initialize a list to store the data
  # file_structure <- h5ls(path,recursive=FALSE)
  file_structure <- h5dump(path,load=FALSE)

  # Extract unique group names (assuming groups end with "/data")
  # group_names <- unique(sub("/data$", "", grep("/data$", file_structure$name, value = TRUE)))
  group_names= names(file_structure)

  # Initialize a list to store the data
  data_list <- list()

  # Iterate over each group name
  for (name in group_names) {
    # Read the data, column names, and row names from the HDF5 file
    data <- h5read(file = path, name = paste0("/",name, "/data"))
    
    if(!is.null(file_structure[[name]]$samples)){
      samples <- h5read(file = path, name = paste0("/",name, "/samples"))
      colnames(data) <- samples
    }
    if(!is.null(file_structure[[name]]$genes)){
      genes <- h5read(file = path, name = paste0("/",name, "/genes"))
      rownames(data) <- genes
    }

    data_list[[name]] <- data
    # Store the data in the list with the group name
  }

  return(data_list)
}

