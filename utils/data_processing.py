import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
import anndata as ad

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

file = 'sparse_matrix_R2.h5'

ref = read_all_ref_hdf5(file)



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



file = 'ref_sc_peng.h5'
# write_all_ref_hdf5(file, ref)
# meta = ref['ref_scRNA']['ref_sc_baron']['metadata']

# for col in meta.columns:
#     if pd.api.types.is_integer_dtype(meta[col]):  # If the column is integer, store as int32
#         dtype_list.append((col, np.int32))
#     else:  # Otherwise, store as a UTF-8 string
#         dtype_list.append((col, h5py.string_dtype('utf-8')))