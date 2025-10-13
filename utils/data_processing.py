import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
import anndata as ad


import numpy as np
import h5py
from scipy.sparse import csc_matrix, issparse
import pandas as pd

compression_lvl = 6
compression_type = "gzip"


# file_path = "sparse_matrix.h5"



# def save_sparse_hdf5(file_path, counts_adata, metadata):
#     """
#     Saves a sparse matrix (counts) and metadata (pandas DataFrame) to an HDF5 file.
#     """

#     counts =  csc_matrix(counts_adata.X.T)
#     with h5py.File(file_path, "w") as f:
#         # Save sparse matrix
#         group_counts = f.create_group("counts")
#         group_counts.create_dataset("data", data=counts.data, compression=compression_type, compression_opts=compression_lvl , chunks= (len(counts.data))  ,shuffle=True) #False)
#         group_counts.create_dataset("shape", data=counts.shape, compression=compression_type, compression_opts=compression_lvl , chunks= (len(counts.shape))  ,shuffle=True) #False)
#         group_counts.create_dataset("indices", data=counts.indices, compression=compression_type, compression_opts=compression_lvl , chunks= (len(counts.indices))  ,shuffle=True) #False)
#         group_counts.create_dataset("indptr", data=counts.indptr, compression=compression_type, compression_opts=compression_lvl , chunks= (len(counts.indptr))  ,shuffle=True) #False)

#         # Save row and column names
#         # if counts.shape[0] == len(counts_adata.var.index):  # Ensure rownames match
#         # group_counts.create_dataset("genes", data=np.array(counts_adata.var.index), compression=compression_type, compression_opts=compression_lvl, chunks=True)  # Convert to bytes
#         # if counts.shape[1] == len(counts_adata.var.index):  # Ensure colnames match
#         # group_counts.create_dataset("", data=np.array(counts_adata.obs.index), compression=compression_type, compression_opts=compression_lvl, chunks=True)  
        
#         group_counts.create_dataset("genes", data=np.array([x.encode('utf-8') for x in counts_adata.var.index]), compression=compression_type, compression_opts=compression_lvl,shuffle=True)
#         group_counts.create_dataset("cells", data=np.array([x.encode('utf-8') for x in counts_adata.obs.index]), compression=compression_type, compression_opts=compression_lvl,shuffle=True)   

#         # Save metadata (as a structured array)
#         # group_meta = f.create_group("meta")
#         # group_meta.create_dataset("/meta",data= metadata, compression=compression_type, compression_opts=compression_lvl  ,shuffle=True)
#         # metadata.to_hdf(file_path, key="meta/meta", mode="a")  # Append mode

# file_path_to_save = "sparse_matrix_py3.h5"


# # # Example Sparse Matrix
# # counts = csr_matrix(np.random.randint(0, 5, size=(1000, 500)))  # 1000 x 500 sparse matrix
# # metadata = pd.DataFrame({"cell_type": ["A"] * 500, "batch": ["1"] * 500}, index=[f"Cell{i}" for i in range(500)])

# # # Save to HDF5
# save_sparse_hdf5(file_path_to_save, adata, metadata)




# def load_sparse_hdf5(file_path):
#     # """
#     # Loads a sparse matrix (counts) and metadata from an HDF5 file.
#     # """
#     with h5py.File(file_path, "r") as f:
#         # Load sparse matrix
#         data = f["counts/data"][:]
#         indices = f["counts/indices"][:]
#         indptr = f["counts/indptr"][:]
#         shape = tuple(f["counts/shape"][:])
#         adata = ad.AnnData(csc_matrix((data, indices, indptr), shape=shape,dtype=np.int32).T)

#         # Load row and column names
#         rownames = [x.decode("utf-8") for x in f["counts/rownames"][:]] if "counts/rownames" in f else None
#         colnames = [x.decode("utf-8") for x in f["counts/colnames"][:]] if "counts/colnames" in f else None

#     if rownames is not None:
#         adata.var_names = rownames

#         # Set the column names (gene names) in the AnnData object
#     if colnames is not None:
#         adata.obs_names  = colnames

#         # Load metadata
#     metadata = pd.read_hdf(file_path, key="metameta")
#     metadata.index = colnames  # Restore rownames

#     return adata, metadata

# # Load Data
# adata, metadata = load_sparse_hdf5(file_path)

# print(sparse_restored)
# print(meta_restored.head())

# def encode_strings(data):
#     """Encode all string-like objects in a DataFrame or array to byte strings."""
#     if isinstance(data, pd.DataFrame):
#         return data.applymap(lambda x: x.encode('utf-8') if isinstance(x, str) else x)
#     elif isinstance(data, np.ndarray):
#         return np.vectorize(lambda x: x.encode('utf-8') if isinstance(x, str) else x)(data)
#     else:
#         return data

def read_all_ref_hdf5(path):
    with h5py.File(path, "r") as f:
        # Read ref_bulkRNA data
        ref_bulkRNA = pd.DataFrame(f["ref_bulkRNA/data"][:]).T
        ref_bulkRNA.columns = [x.decode('utf-8') for x in f["ref_bulkRNA/cell_types"][:]]
        ref_bulkRNA.index = [x.decode('utf-8') for x in f["ref_bulkRNA/genes"][:]]


        # Read ref_met data
        ref_met =  pd.DataFrame(f["ref_met/data"][:]).T
        ref_met.columns = [x.decode('utf-8') for x in f["ref_met/cell_types"][:]]
        ref_met.index = [x.decode('utf-8') for x in f["ref_met/CpG_sites"][:]]


        # Read ref_scRNA data
        ref_scRNA = {}
        datasets = ["ref_sc_peng", "ref_sc_baron", "ref_sc_raghavan"]

        for dataset in datasets:
            group = f"ref_scRNA/{dataset}"

            counts_data = f[f"{group}/data"][:]
            counts_shape = tuple(f[f"{group}/shape"][:])
            counts_indices = f[f"{group}/indices"][:]
            counts_indptr = f[f"{group}/indptr"][:]

            genes = [x.decode('utf-8') for x in f[f"{group}/genes"][:]] if f"{group}/genes" in f else None
            cells = [x.decode('utf-8') for x in f[f"{group}/cell"][:]] if f"{group}/cell" in f else None

            adata = ad.AnnData(csc_matrix((counts_data, counts_indices, counts_indptr), shape=counts_shape,dtype=np.int32).T)

            if genes is not None:
                adata.var_names = genes

            if cells is not None:
                adata.obs_names  = cells

            meta = pd.DataFrame(f[f"{group}/meta"][:])
            meta.index = cells
            for col in meta.columns:
                if meta[col].dtype == 'object':  # Check if the column contains byte strings
                    meta[col] = meta[col].apply(lambda x: x.decode('utf-8') if isinstance(x, bytes) else x)
            meta = meta.convert_dtypes()
            # meta = encode_strings(meta)

            ref_scRNA[dataset] = {
                "counts": adata,
                "metadata": meta
            }

        ref_all = {
            "ref_bulkRNA": ref_bulkRNA,
            "ref_met": ref_met,
            "ref_scRNA": ref_scRNA
        }

    return ref_all

# file = 'sparse_matrix_R2.h5'

# ref = read_all_ref_hdf5(file)



# def save_sparse_hdf5(file_path, counts_adata, metadata):
#     """
#     Saves a sparse matrix (counts) and metadata (pandas DataFrame) to an HDF5 file.
#     """

#     counts =  csc_matrix(counts_adata.X.T)
#     with h5py.File(file_path, "w") as f:
#         # Save sparse matrix
#         group_counts = f.create_group("counts")
#         group_counts.create_dataset("data", data=counts.data, compression=compression_type, compression_opts=compression_lvl , chunks= (len(counts.data))  ,shuffle=True) #False)
#         group_counts.create_dataset("shape", data=counts.shape, compression=compression_type, compression_opts=compression_lvl , chunks= (len(counts.shape))  ,shuffle=True) #False)
#         group_counts.create_dataset("indices", data=counts.indices, compression=compression_type, compression_opts=compression_lvl , chunks= (len(counts.indices))  ,shuffle=True) #False)
#         group_counts.create_dataset("indptr", data=counts.indptr, compression=compression_type, compression_opts=compression_lvl , chunks= (len(counts.indptr))  ,shuffle=True) #False)

#         # Save row and column names
#         # if counts.shape[0] == len(counts_adata.var.index):  # Ensure rownames match
#         # group_counts.create_dataset("genes", data=np.array(counts_adata.var.index), compression=compression_type, compression_opts=compression_lvl, chunks=True)  # Convert to bytes
#         # if counts.shape[1] == len(counts_adata.var.index):  # Ensure colnames match
#         # group_counts.create_dataset("", data=np.array(counts_adata.obs.index), compression=compression_type, compression_opts=compression_lvl, chunks=True)  
        
#         group_counts.create_dataset("genes", data=np.array([x.encode('utf-8') for x in counts_adata.var.index]), compression=compression_type, compression_opts=compression_lvl,shuffle=True)
#         group_counts.create_dataset("cells", data=np.array([x.encode('utf-8') for x in counts_adata.obs.index]), compression=compression_type, compression_opts=compression_lvl,shuffle=True)   

#         # Save metadata (as a structured array)
#         # group_meta = f.create_group("meta")
#         # group_meta.create_dataset("/meta",data= metadata, compression=compression_type, compression_opts=compression_lvl  ,shuffle=True)
#         # metadata.to_hdf(file_path, key="meta/meta", mode="a")  # Append mode

# file_path_to_save = "sparse_matrix_py3.h5"


# # Example Sparse Matrix
# counts = csr_matrix(np.random.randint(0, 5, size=(1000, 500)))  # 1000 x 500 sparse matrix
# metadata = pd.DataFrame({"cell_type": ["A"] * 500, "batch": ["1"] * 500}, index=[f"Cell{i}" for i in range(500)])

# # Save to HDF5
# save_sparse_hdf5(file_path_to_save, adata, metadata)



def write_all_ref_hdf5(path, ref_all):
    with h5py.File(path, "w") as f:
        #### Write ref_bulkRNA data
        group_bulk = f.create_group("ref_bulkRNA")
        group_bulk.create_dataset("data", data=ref_all['ref_bulkRNA'].T  , chunks = ref_all['ref_bulkRNA'].T.shape, 
                                  compression=compression_type, compression_opts=compression_lvl , shuffle=True  )
        group_bulk.create_dataset("cell_types", data=np.array(ref_all['ref_bulkRNA'].columns, dtype='S'),  chunks = len(ref_all['ref_bulkRNA'].columns),
                                  compression=compression_type, compression_opts=compression_lvl , shuffle=True    )
        group_bulk.create_dataset("genes", data=np.array(ref_all['ref_bulkRNA'].index, dtype='S'), chunks = len(ref_all['ref_bulkRNA'].index) ,
                                  compression=compression_type, compression_opts=compression_lvl , shuffle=True )#, chunks = ref_all['ref_bulkRNA'].T.shape) 

        # Write ref_met data
        group_met = f.create_group("ref_met")
        group_met.create_dataset("data", data=ref_all['ref_met'].T,
                                 compression=compression_type, compression_opts=compression_lvl , shuffle=True, chunks =ref_all['ref_met'].T.shape)
        group_met.create_dataset("cell_types", data=np.array(ref_all['ref_met'].columns, dtype='S' ), chunks = len(ref_all['ref_met'].columns) ,
                                 compression=compression_type, compression_opts=compression_lvl , shuffle=True   )
        group_met.create_dataset("CpG_sites", data=np.array(ref_all['ref_met'].index , dtype='S'), chunks = len(ref_all['ref_met'].index),
                                 compression=compression_type, compression_opts=compression_lvl , shuffle=True)

        ##### Write ref_scRNA data
        datasets = ["ref_sc_peng", "ref_sc_raghavan" , "ref_sc_baron"]
        for dataset in datasets:
            group_sc = f.create_group(f"ref_scRNA/{dataset}")
            counts =  csc_matrix(ref_all['ref_scRNA'][dataset]['counts'].X.T)   
            meta = ref_all['ref_scRNA'][dataset]['metadata']

            group_sc.create_dataset("data", data=counts.data,compression=compression_type ,chunks= (len(counts.data))  ,
                                    compression_opts=compression_lvl ,shuffle=True)
            group_sc.create_dataset("shape", data=counts.shape,compression=compression_type,
                                    compression_opts=compression_lvl   ,shuffle=True)
            group_sc.create_dataset("indices", data=counts.indices,compression=compression_type, chunks = len(counts.indices),
                                    compression_opts=compression_lvl   ,shuffle=True)
            group_sc.create_dataset("indptr", data=counts.indptr,compression=compression_type, chunks = len(counts.indptr),
                                    compression_opts=compression_lvl   ,shuffle=True)
            
            var_names = ref_all['ref_scRNA'][dataset]['counts'].var_names
            if counts.shape[0] == len(var_names):
                group_sc.create_dataset("genes", data= np.array(var_names, dtype='S'),chunks= len(var_names) ,
                                        compression=compression_type, compression_opts=compression_lvl,shuffle=True  )
           
            obs_name = ref_all['ref_scRNA'][dataset]['counts'].obs_names
            if counts.shape[1] == len(obs_name):
                group_sc.create_dataset("cell", data=np.array(obs_name, dtype='S'), chunks = len(obs_name) ,
                                        compression=compression_type, compression_opts=compression_lvl,shuffle=True) 


            dtype_list = []
            for col in meta.columns:
                if pd.api.types.is_integer_dtype(meta[col]):  # If the column is integer, store as int32
                    dtype_list.append((col, np.int32))
                else:  # Otherwise, store as a UTF-8 string
                    dtype_list.append((col, h5py.string_dtype('utf-8')))

            # Convert DataFrame to a structured NumPy array with the correct dtypes
            meta_np = np.array(list(meta.itertuples(index=False, name=None)), dtype=dtype_list)   

            group_sc.create_dataset("meta", data=meta_np,  compression=compression_type, 
                                    compression_opts=compression_lvl   ,shuffle=True , chunks= len (meta_np)) 



# file = 'ref_sc_peng.h5'
# write_all_ref_hdf5(file, ref)
# meta = ref['ref_scRNA']['ref_sc_baron']['metadata']

# for col in meta.columns:
#     if pd.api.types.is_integer_dtype(meta[col]):  # If the column is integer, store as int32
#         dtype_list.append((col, np.int32))
#     else:  # Otherwise, store as a UTF-8 string
#         dtype_list.append((col, h5py.string_dtype('utf-8')))





def write_sparse_matrix(group, fw, counts, meta, data=None, scale_data=None,
                        compression_type="gzip", compression_lvl=6):


    # group = fw.create_group(group)

    # Write sparse matrix components
    group.create_dataset("data", data=counts.data, compression=compression_type, compression_opts=compression_lvl,chunks= (len(counts.data)),shuffle=True)
    group.create_dataset("shape",   data=counts.shape, compression=compression_type, compression_opts=compression_lvl,chunks=   (len(counts.shape)),shuffle=True)
    group.create_dataset("indices", data=counts.indices, compression=compression_type, compression_opts=compression_lvl,chunks= (len(counts.indices)),shuffle=True)
    group.create_dataset("indptr",  data=counts.indptr, compression=compression_type, compression_opts=compression_lvl,chunks=  (len(counts.indptr)),shuffle=True)

    if hasattr(counts, 'var_names'):
        group.create_dataset("genes", data=np.array(counts.var_names, dtype='S'), compression=compression_type,chunks=  (len(counts.var_names)),shuffle=True)
    if hasattr(counts, 'obs_names'):
        group.create_dataset("cell", data=np.array(counts.obs_names, dtype='S'), compression=compression_type,chunks= (len(counts.obs_names)),shuffle=True)

    # Write metadata
    dtype_list = []
    for col in meta.columns:
        if pd.api.types.is_integer_dtype(meta[col]):
            dtype_list.append((col, np.int32))
        else:
            dtype_list.append((col, h5py.string_dtype('utf-8')))
    meta_np = np.array(list(meta.itertuples(index=False, name=None)), dtype=dtype_list)
    group.create_dataset("meta", data=meta_np, compression=compression_type, compression_opts=compression_lvl,chunks= (len(meta_np)),shuffle=True)

    # Optionally write normalized data
    if data is not None and issparse(data):
        group.create_dataset("normalized_data", data=data.data, compression=compression_type, compression_opts=compression_lvl,chunks= (len(data.data)),shuffle=True)

    # Optionally write scale.data
    if scale_data is not None:
        if issparse(scale_data):
            group.create_dataset("scale_data", data=scale_data.data, compression=compression_type, compression_opts=compression_lvl,chunks= (len(scale_data.data)),shuffle=True)
        elif isinstance(scale_data, np.ndarray) or isinstance(scale_data, pd.DataFrame):
            scale_arr = np.array(scale_data, dtype=np.float32)
            group.create_dataset("scale_data", data=scale_arr.flatten(), compression=compression_type, compression_opts=compression_lvl,chunks= (len(scale_arr)),shuffle=True)
            group.create_dataset("scale_data_shape", data=scale_arr.shape, compression=compression_type, compression_opts=compression_lvl,chunks= (len(scale_arr.shape)),shuffle=True)
            if isinstance(scale_data, pd.DataFrame):
                group.create_dataset("scale_data_genes", data=np.array(scale_data.index, dtype='S'), compression=compression_type, compression_opts=compression_lvl,chunks= (len(scale_data.index)),shuffle=True)
                group.create_dataset("scale_data_cells", data=np.array(scale_data.columns, dtype='S'), compression=compression_type, compression_opts=compression_lvl,chunks= (len(scale_data.columns)),shuffle=True)


def write_global_hdf5(path, data_list, compression_type="gzip", compression_lvl=4):
    with h5py.File(path, "w") as fw:
        for name in data_list:
            if name == "ref_scRNA":
                group_ref = fw.create_group("ref_scRNA")

                for dataset in data_list[name]:
                    group_path = f"ref_scRNA/{dataset}"
                    group = fw.create_group(group_path)

                    # seurat_obj = None
                    # seurat_field_name = None
                    # group_seurat_path = group_path

                    # Check if dataset is Seurat-like object (dict with laye        rs)
                    if "counts" in data_list[name][dataset] and "metadata" in data_list[name][dataset]:
                        # counts = csc_matrix(data_list[name][dataset]["counts"].X.T)
                        counts = data_list[name][dataset]["counts"]
                        meta = data_list[name][dataset]["metadata"]
                        data = None
                        scale_data = None

                        if "data" in data_list[name][dataset]:
                            # data = csc_matrix(data_list[name][dataset]["data"].X.T)
                            data = data_list[name][dataset]["data"]
                        if "scale.data" in data_list[name][dataset]:
                            scale_data = data_list[name][dataset]["scale.data"]

                        # # Write marker that this is a Seurat-like object
                        # group.create_dataset("object_type", data="seurat")
                        # group.create_dataset("seurat_field_name", data=str(dataset))

                        write_sparse_matrix(group, fw, counts, meta, data=data, scale_data=scale_data,
                                            compression_type=compression_type, compression_lvl=compression_lvl)
                    
                    else:
                        # If not a Seurat object but still has "counts" and "metadata"
                        if "counts" in data_list[name][dataset] and "metadata" in data_list[name][dataset]:
                            counts = csc_matrix(data_list[name][dataset]["counts"].X.T)
                            meta = data_list[name][dataset]["metadata"]
                            write_sparse_matrix(group_path, fw, counts, meta, compression_type=compression_type, compression_lvl=compression_lvl)

            else:
                # Fallback: save as flat table
                df = data_list[name]
                df_group = fw.create_group(name)
                df_group.create_dataset("data",data=df.T, compression=compression_type, compression_opts=compression_lvl)
                df_group.create_dataset("samples",data=list(map(str,df.columns)), compression=compression_type, compression_opts=compression_lvl)
                df_group.create_dataset("genes",data=list(map(str,df.index)), compression=compression_type, compression_opts=compression_lvl)
                # for col in df.columns:
                #     df_group.create_dataset(col, data=np.array(df[col]), compression=compression_type, compression_opts=compression_lvl)


def set_dataframe_index_and_columns(group, group_name, df):
    def get_string_array(field_name):
        if field_name in group_name :
            raw = group[field_name][()]
            if raw.dtype.kind in {'S', 'O'}:
                return [x.decode('utf-8') if isinstance(x, bytes) else x for x in raw]
            return list(raw)
        return None

    # Assign columns
    for field in ["samples", "cell_types"]:
        values = get_string_array(field)
        if values is not None:
            df.columns = values
            break  # Prioritize 'samples' over 'cell_types' like in R

    # Assign rows
    for field in ["genes", "CpG_sites"]:
        values = get_string_array(field)
        if values is not None:
            df.index = values
            break  # Prioritize 'genes' over 'CpG_sites'
    return df

def get_h5_structure(h5group):
    structure = {}
    for key in h5group:
        if isinstance(h5group[key], h5py.Group):
            structure[key] = get_h5_structure(h5group[key])
        else:
            structure[key] = "dataset"
    return structure

def read_data_frame(f, group_name, file_structure):
    # with h5py.File(path, "r") as f:
    group = f[group_name]
    df = pd.DataFrame(group["data"][()]).T

    df = set_dataframe_index_and_columns(group, file_structure[group_name],  df)

    return df
    


def read_sparse_matrix(group_path, group_structure, h5file):
    # Load sparse matrix components
    counts_data = np.array(h5file[f"{group_path}/data"])
    counts_shape = tuple(h5file[f"{group_path}/shape"])
    counts_indices = np.array(h5file[f"{group_path}/indices"])
    counts_indptr = np.array(h5file[f"{group_path}/indptr"])

    gene_names = [g.decode("utf-8") for g in h5file[f"{group_path}/genes"][:]]
    cell_names = [c.decode("utf-8") for c in h5file[f"{group_path}/cell"][:]]

    # Reconstruct sparse matrix
    counts = csc_matrix((counts_data, counts_indices, counts_indptr), shape=counts_shape)
    counts.var_names = gene_names
    counts.obs_names = cell_names

    # Load metadata
    meta_data_raw = h5file[f"{group_path}/meta"]
    meta = pd.DataFrame({key: meta_data_raw[key][:] for key in meta_data_raw.dtype.names})
    meta.columns = meta_data_raw.dtype.names
    meta.index = cell_names

    # Check if it was a Seurat object
    is_seurat = "object_type" in group_structure and h5file[f"{group_path}/object_type"][()].decode("utf-8") == "seurat"

    if is_seurat:
        # Add normalized data layer
        normalized_data = None
        if "normalized_data" in group_structure:
            normalized_data = csc_matrix((h5file[f"{group_path}/normalized_data"][:], counts_indices, counts_indptr), shape=counts_shape)

        # Add scaled data
        scale_data = None
        if "scale_data" in group_structure:
            scale_data = csc_matrix((h5file[f"{group_path}/scale_data"][:], counts_indices, counts_indptr), shape=counts_shape)

        return {
            "counts": counts,
            "metadata": meta,
            "data": normalized_data,
            "scale.data": scale_data,
            "is_seurat": True
        }

    else:
        return {
            "counts": counts,
            "metadata": meta
        }



def read_hdf5(path):
    data_list = {}

    with h5py.File(path, "r") as f:
        file_structure = get_h5_structure(f)

        for name in file_structure:
            if name == "ref_scRNA":
                ref_scRNA = {}
                for dataset in file_structure[name]:
                    group_path = f"{name}/{dataset}"
                    seurat_data = read_sparse_matrix(group_path, file_structure[name][dataset], f)

                    # Look for nested groups (fields like scaled_data, data, etc.)
                    expected = {"cell", "data", "genes", "indices", "indptr", "meta", "shape", "object_type", "seurat_field_name", "normalized_data", "scale_data"}
                    actual = set(file_structure[name][dataset].keys())
                    nested_fields = actual - expected

                    for field in nested_fields:
                        sub_group_path = f"{name}/{dataset}/{field}"
                        seurat_data[field] = read_sparse_matrix(sub_group_path, file_structure[name][dataset][field], f)

                    ref_scRNA[dataset] = seurat_data

                data_list[name] = ref_scRNA

            else:
                data_list[name] = read_data_frame(f, name, file_structure)


    return data_list




# file= "data/ref.h5"
# r_ref = read_hdf5(file)


# # file= "data/mixes1_insilicodirichletCopule_pdac.h5"
# file= "data/mixes1_invitro_pdac.h5"
# r =  read_hdf5(file)


# # path = file



# test_f = 'tmp.h5'

# write_global_hdf5(test_f,r_ref )
