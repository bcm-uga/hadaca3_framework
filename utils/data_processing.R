
# library(assertr)
library(rhdf5)
library(Matrix) 


get_omic <- function(path) {
  path_parts <- unlist(strsplit(path, "/"))
  return(path_parts[length(path_parts) - 1])
}

omic2list_name = list(mixRNA = "mix_rna" , 
mixMET =  "mix_met",
MET='ref_met',
RNA= 'ref_bulkRNA',
scRNA= 'ref_scRNA' ) 


write_data_frame <- function(name, path,data){
    h5createGroup(path, name)
    h5write(data, file = path, name = paste0(name, "/data"))
    if(length(colnames(data))){
      h5write(colnames(data), file = path, name = paste0(name, "/samples"))
    }
    if(length(rownames(data))){
      h5write(rownames(data), file = path, name = paste0(name, "/genes"))
    }
}



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
    
    write_sparse_matrix(paste0( 'ref_scRNA/' , dataset  ), path,counts,meta)

  }

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
  file_structure <- h5dump(path, load = FALSE)
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
      scRNA = read_sparse_matrix(paste0('ref_scRNA/',dataset), file_structure[[dataset]],path)

      ref_scRNA[[dataset]] <- list(
        counts = scRNA$counts,
        metadata = scRNA$meta
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


write_sparse_matrix <- function(group, path, counts, meta, data = NULL, scale_data = NULL) {
  h5createGroup(path, group)

  # Write counts
  h5createDataset(path, paste0(group, '/data'), dims = length(counts@x), storage.mode = "double")
  h5write(counts@x, path, paste0(group, "/data"))
  h5write(dim(counts), path, paste0(group, "/shape"))
  h5createDataset(path, paste0(group, '/indices'), dims = length(counts@i), storage.mode = "integer")
  h5write(counts@i, path, paste0(group, "/indices"))
  h5write(counts@p, path, paste0(group, "/indptr"))
  if (!is.null(rownames(counts))) {
    h5write(rownames(counts), path, paste0(group, "/genes"))
  }
  if (!is.null(colnames(counts))) {
    h5write(colnames(counts), path, paste0(group, "/cell"))
  }

  # Write metadata
  h5write(meta, path, paste0(group, "/meta"))

  # # Optionally write "data" and "scale.data" layers if provided
  if (!is.null(data)) {
    h5createDataset(path, paste0(group, '/normalized_data'), dims = length(data@x), storage.mode = "double")
    h5write(data@x, path, paste0(group, "/normalized_data"))
  }
  # if (!is.null(scale_data)) {
  #   if (inherits(scale_data, "dgCMatrix")) {
  #     h5createDataset(path, paste0(group, '/scale_data'), dims = length(scale_data@x), storage.mode = "double")
  #     h5write(scale_data@x, path, paste0(group, "/scale_data"))
  #   } else if (is.matrix(scale_data)) {
  #     h5write(as.numeric(scale_data), path, paste0(group, "/scale_data"))
  #     h5write(dim(scale_data), path, paste0(group, "/scale_data_shape"))
  #     h5write(rownames(scale_data), path, paste0(group, "/scale_data_genes"))
  #     h5write(colnames(scale_data), path, paste0(group, "/scale_data_cells"))
  #   }
  #   # h5write(scale_data@x, path, paste0(group, "/scale_data"))
  # }
}


write_global_hdf5 <- function(path, data_list) {
  # Create the HDF5 file
  h5createFile(path)

  for (name in names(data_list)) {

    if (name == "ref_scRNA") {
      h5createGroup(path, 'ref_scRNA/')

      for (dataset in names(data_list[[name]])) {
        group = paste0('ref_scRNA/', dataset)
        h5createGroup(path, group)

        # entry <- data_list[[name]][[dataset]]
        seurat_obj <- NULL
        seurat_field_name <- NULL
        group_seurat = group
        if(inherits(data_list[[name]][[dataset]], "Seurat")){
          seurat_obj <- data_list[[name]][[dataset]]
          seurat_field_name <- dataset
          
        }

        ##Check if there is a Seurat Object 
        for (field in names(data_list[[name]][[dataset]])) {
          if (inherits(data_list[[name]][[dataset]][[field]], "Seurat")) {
            seurat_obj <- data_list[[name]][[dataset]][[field]]
            seurat_field_name <- field
            group_seurat = paste0('ref_scRNA/', dataset,'/',field)
            # h5createGroup(path,paste0(path,group), field)
            h5createGroup(path,group_seurat)

            break
          }
        }

        if (!is.null(seurat_obj)) {
          library(Seurat)

          # counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
          # data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
          # scale_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "scale.data")
          # meta <- seurat_obj@meta.data

          h5write("seurat", path, paste0(group_seurat, "/object_type"))
          h5write(seurat_field_name, path, paste0(group_seurat, "/seurat_field_name"))

          write_sparse_matrix(group_seurat, path,
            counts = GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),
            meta = seurat_obj@meta.data , 
            data = GetAssayData(seurat_obj, assay = "RNA", layer = "data"), 
            scale_data= GetAssayData(seurat_obj, assay = "RNA", layer = "scale.data")
          )

        }
      
        if("counts" %in% names(data_list[[name]][[dataset]]) && "metadata" %in% names(data_list[[name]][[dataset]])) {
          write_sparse_matrix(group, path, data_list[[name]][[dataset]]$counts, data_list[[name]][[dataset]]$metadata)
        }
      }
    } else {
      write_data_frame(name, path, data_list[[name]])
    }
  }
}



read_data_frame <- function(path,name,file_structure){
    data <- h5read(file = path, name = paste0("/",name, "/data"))
    
    if(!is.null(file_structure[[name]]$samples)){
      samples <- h5read(file = path, name = paste0("/",name, "/samples"))
      colnames(data) <- samples
    }
    if(!is.null(file_structure[[name]]$genes)){
      genes <- h5read(file = path, name = paste0("/",name, "/genes"))
      rownames(data) <- genes
    }
    return(data)

}

read_sparse_matrix <- function(group,group_structure, path) {
  counts_data <- as.numeric(h5read(path, paste0(group, "/data")))
  counts_shape <- as.integer(h5read(path, paste0(group, "/shape")))
  counts_indices <- as.integer(h5read(path, paste0(group, "/indices")))
  counts_indptr <- as.integer(h5read(path, paste0(group, "/indptr")))
  gene_names <- as.character(h5read(path, paste0(group, "/genes")))
  cell_names <- as.character(h5read(path, paste0(group, "/cell")))

  counts <- new("dgCMatrix",
                x = counts_data,
                i = counts_indices,
                p = counts_indptr,
                Dim = counts_shape,
                Dimnames = list(gene_names, cell_names))

  meta <- h5read(path, paste0(group, "/meta"))
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  rownames(meta) <- cell_names

  # Detect if the original object was Seurat
  is_seurat <- FALSE
  # if (paste0(group, "/object_type") %in% h5ls(path, recursive = TRUE)$name) {
  if ( "object_type" %in% names(group_structure)) {
    object_type <- h5read(path, paste0(group, "/object_type"))
    is_seurat <- object_type == "seurat"
  }

  if (is_seurat) {
    library(Seurat)
    seurat_obj <- CreateSeuratObject(counts = counts, meta.data = meta)
    DefaultAssay(seurat_obj) <- "RNA"
    if ("normalized_data" %in% names(group_structure)) {
      # data_values <- as.numeric(h5read(path, paste0(group, "/normalized_data")))
      data <- counts
      data@x <- as.numeric(h5read(path, paste0(group, "/normalized_data")))
      seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], layer = "data", new.data = data)

      # seurat_obj[["RNA"]]@data <- data
    }

    if ("scale_data" %in% names(group_structure)) {
      scale_values <- as.numeric(h5read(path, paste0(group, "/scale_data")))
      scale_data <- counts
      scale_data@x <- scale_values
      # seurat_obj[["RNA"]]@scale.data <- scale_data
      seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], layer = "scale.data", new.data = scaled_data)

    }

    return(seurat_obj)
  } else {
    return(list(counts = counts, metadata = meta))
  }
}



read_hdf5 <- function(path) {
  file_structure <- h5dump(path, load = FALSE)
  group_names <- names(file_structure)
  data_list <- list()

  for (name in group_names) {
    if (name == "ref_scRNA") {
      ref_scRNA <- list()

      for (dataset in names(file_structure[[name]])) {
        seuratobj_or_counts_n_metadata <- read_sparse_matrix(paste0('ref_scRNA/', dataset), file_structure[[name]][[dataset]] , path)
        
        ## check if there is another data inside here !  TODO change this with a metadata list in the root of this file... 
        expected_names <- c("cell", "data", "genes", "indices", "indptr", "meta", "shape")
        actual_names <- names(file_structure[[name]][[dataset]])
        unexpected_names = actual_names[!actual_names %in% expected_names]

        for(field in unexpected_names ){
          seuratobj_or_counts_n_metadata[[field]] = read_sparse_matrix(paste0('ref_scRNA/', dataset,'/',field), file_structure[[name]][[dataset]][[field]] , path)
        }
        ref_scRNA[[dataset]] = seuratobj_or_counts_n_metadata
      }

      data_list[[name]] <- ref_scRNA
    } else {
      data_list[[name]] <- read_data_frame(path, name, file_structure)
    }
  }

  return(data_list)
}

# read_hdf5 <- function(path) {
#   file_structure <- h5dump(path, load = FALSE)
#   group_names <- names(file_structure)
#   data_list <- list()

#   for (name in group_names) {
#     if (name == "ref_scRNA") {
#       ref_scRNA <- list()

#       for (dataset in names(file_structure[[name]])) {
#         group <- paste0("ref_scRNA/", dataset)

#         object_type <- tryCatch(h5read(path, paste0(group, "/object_type")), error = function(e) NULL)
#         seurat_field_name <- tryCatch(h5read(path, paste0(group, "/seurat_field_name")), error = function(e) NULL)

#         result <- list()

#         # Always read raw counts + metadata
#         matrix_info <- read_sparse_matrix(group, path)
#         result$counts <- matrix_info$counts
#         result$metadata <- matrix_info$meta

#         # If it was originally a Seurat object
#         if (!is.null(object_type) && object_type == "seurat") {
#           library(Seurat)

#           # Optionally try to read data & scale_data
#           # data <- tryCatch(h5read(path, paste0(group, "/data")), error = function(e) NULL)
#           # scale_data <- tryCatch(h5read(path, paste0(group, "/scale_data")), error = function(e) NULL)

#           seurat_obj <- CreateSeuratObject(counts = matrix_info$counts, meta.data = matrix_info$meta)

#           # if (!is.null(data)) {
#           #   seurat_obj[["RNA"]]@data <- as(Matrix::Matrix(data, sparse = TRUE), "dgCMatrix")
#           # }

#           # if (!is.null(scale_data)) {
#           #   seurat_obj[["RNA"]]@scale.data <- scale_data
#           # }

#           # Save under the original field name or default to "seurat"
#           seurat_field <- if (!is.null(seurat_field_name)) seurat_field_name else "seurat"
#           result[[seurat_field]] <- seurat_obj
#         }

#         ref_scRNA[[dataset]] <- result
#       }

#       data_list[[name]] <- ref_scRNA
#     } else {
#       data_list[[name]] <- read_data_frame(path, name, file_structure)
#     }
#   }

#   return(data_list)
# }