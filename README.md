# Hadaca3 Framework

A framework to collectively develop multi-omic deconvolution methods.

## How to start?

```
cd ~/projects
git clone git@github.com:bcm-uga/hadaca3_framework.git
cd hadaca3_framework
```
## Conda environement

Set up your conda environement as follow:

```{Automatic CI-conda}
conda create -y -n hadaca3framework_env
conda activate hadaca3framework_env


mamba install -y  -c bioconda -c conda-forge -c r snakemake python r-base r-rmarkdown r-nnls r-seurat bioconductor-rhdf5 r-quadprog r-coda.base r-dt bioconductor-toast  psutil nextflow=24.10.5 r-lubridate

```
<!-- h5py   -->


<!-- r-clue r-coda.base r-ggpubr bioconductor-complexheatmap bioconductor-mofa2 r-viridis r-magrittr r-dplyr r-nnls graphviz r-tictoc  graphviz python-kaleido tenacity plotly r-bisquerna r-extraDistr r-MASS r-EPIC r-fmsb bioconductor-toast bioconductor-omicade4 r-mixomics r-mixkernel rpy2 scikit-learn keras tensorflow bioconductor-viper bioconductor-ADImpute r-WGCNA r-see r-ggfortify -->

## Getting otiginal data

The section describes which data are needed to execute the entire pipeline and provide the code to download it.

```{Automatic CI-data}
mkdir -p ~/projects/hadaca3_framework/data
cd ~/projects/hadaca3_framework/data

# from CIMENT/GRICAD cluster using rsync and cargo node
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicodirichletCopule_pdac.h5 .
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicodirichletEMFA_pdac.h5 .
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicopseudobulk_pdac.h5 . 
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_invitro_pdac.h5 .
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_invivo_pdac.h5 .
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicodirichletCopule_pdac.h5 .
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicodirichletEMFA_pdac.h5 .
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicopseudobulk_pdac.h5 .
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_invitro_pdac.h5 .
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_invivo_pdac.h5 .
rsync -auvP cargo:/bettik/hombergn/projects/hadaca3_framework/data/ref.h5 .

# from internet using wget
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_insilicodirichletCopule_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_insilicodirichletEMFA_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_insilicopseudobulk_pdac.h5 
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_invitro_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_invivo_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_insilicodirichletCopule_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_insilicodirichletEMFA_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_insilicopseudobulk_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_invitro_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_invivo_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/ref.h5
```
## Execute the pipeline: 

### snakemake

```
cd ~/projects/hadaca3_framework
snakemake --cores 1 -s 00_run_pipeline.smk -p clean  # keep it clean, keep it green!
snakemake --cores 4 -s 00_run_pipeline.smk -pn       # dry-run
```

This pipeline can be visualised by generating its DAG:
```
snakemake --forceall --dag -s 00_run_pipeline.smk | dot -Tpdf > dag.pdf
```

Run with another setup: 
```
snakemake -s 00_run_pipeline.smk  --cores 4  --config setup_folder='benchmark/setup/1/' 
```



### N E X T F L O W  

```
nextflow run 00_run_pipeline.nf
nextflow run 00_run_pipeline.nf -stub -resume  #continue and dry run
```

To create a a full report of the pipeline these options could be passed: 
```
nextflow run 00_run_pipeline.nf -with-dag -with-report -with-trace -with-timeline 
```

Run with another setup: 
```
nextflow run 00_run_pipeline.nf -resume --setup_folder benchmark/setup/1/
```



## Blocks description


This framework contains several blocks

- **preprocessing** :  This block is responsible for preparing the raw data for analysis. It may include tasks such as cleaning the data (handling missing values, removing duplicates), normalizing or scaling features, encoding categorical variables, and other transformations to make the data suitable for modeling. This block takes as input multi_data and return multi-data (see below for details).
  
- **feature_selection** : This block focuses on selecting the most relevant features (genes or Cpg sites) from the dataset to use in the model. It helps in reducing the dimensionality of the data, improving model performance, and reducing overfitting by eliminating irrelevant or redundant features. This block takes as input multi_data and return multi-data (see below for details).
  
- **deconvolution** : This block contain the algorithm that deconvoluate, such as lm, rlr, nnls...  This block takes as input uni-data and return a prediction (see below for details).
  
- **split** : This block only split multiomics (metylation and RNA) data to only metylation and only RNA.  This block takes as input multi_data and return sort of uni-data (see below for details).

- **early_int** : This block involves combining multiple omics data types (e.g., RNA, MET) into a unified dataset before applying the deconvolution method. It is part of pipeline B. This block takes as input multi_data and return uni-data (see below for details).
  
- **late_int** : This block focuses on integrating the results from multiple omics analyses into a single, cohesive prediction. It is part of pipeline A. This block takes as input a list of several(2) predictions and return one prediction (see below for details).
  
- **intermediate_int** :  This block combines both integration and deconvolution processes from multiple omics data types. It is part of pipeline C. This block takes as input multi_data and return a prediction (see below for details).


## Data types : 
The input and output of each block  
We have two types of data : 
 * Multi_data :  This format contains several mulli-omics such as metylation mixes and rna mix. This data organisation looks like this : 
```
multi_data
├── mix
│   ├── mix_rna
│   └── mix_met
└── ref
    ├── ref_bulkRNA
    ├── ref_met
    └── ref_scRNA
        ├── ref_sc_peng
        ├── ref_sc_baron
        └── ref_sc_raghavan
```
 * Uni-data : contains only one type of omics. 
```
uni_data
├── mix
└── ref
    ├── ref_bulkRNA
    ├── ref_met
    └── ref_scRNA
        ├── ref_sc_peng
        ├── ref_sc_baron
        └── ref_sc_raghavan
```
prediction : contains only the prediction table. 


## Snakemake shenanigan.

The snakeme make will create combinaison between compatible functions inside each block. 

To complete

## HDF5 format.

Hierarchical Data Format (HDF) is a set of file formats (HDF4, HDF5) designed to store and organise large amounts of data. 

It behaves like the OS file system with groups as folders and can handle symlink internaly. 
There are tools such as https://h5web.panosc.eu/ to visualise data. This tool also exists as a VS code extension by the name *H5Web*. 

There is a Linux program that read HDF5 data in a terminal : 
https://support.hdfgroup.org/documentation/hdf5/latest/_view_tools_view.html
You can install it with : 
`sudo apt-get install hdf5-tools`

### How to read and write H5 files ?  

In this hadaca3_framework project, Python and R libraries are provided to read and write data. 
They are named *data_processing* and are located in the *utils* folder.  

Useful functions: 
- *read_all_hdf5(path_of_file,..)* returns the a multi_data (cd data types). The second **optionnal** argument is *to_read=c('mix','ref')* which precise which data will be read. By default, all data of multi_data are read. 

- *write_all_hdf5(path,multi_data)* write to *path*, multi_data (cf data types). 

- *read_hdf5(path)* returns a data_list. This function browse the file from the file on path and browse read all subfolder inside this path. For instance, if the file "exemple.h5" contains /prop1 and /prop2 *read_hdf5(path)* return a list(prop1, prop2).

- *write_global_hdf5(path,data_list)* this function write all sub-data inside the data_list. 
For instance, data_list contain prop1 and prop2, *write_global_hdf5(file,data_list)* will wrote both prop into the file. 


All data should have HDF5 format with a compression level set to 6 and 'gzip' as the compression algorithm. Furthermore, to reduce storage footprints, the data are shuffled and written in one single chunk (chunk size = length(data)). *HDF5 shuffling does not impact order of the uncompressed file*


# Benchmark 

There is an attempt to perfom a benchmark of snakemake vs nextflow. 
The motivation behind the developpement of nextflow is the mandatory step of DAG creation in snakemake which was very time consuming. 

See the README.md inside benchmark folder. 

## TODO

* improve handling of hdf5 files to not rewrite unmodified data

* how to deal with function such as sc_cluster in Fs that requiere a specific pp. 


Nextflow : 
