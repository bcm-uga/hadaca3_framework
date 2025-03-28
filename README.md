# hadaca3_framework

## Purpose

A framework to collectively develop multi-omic deconvolution methods.

The framework contains several blocks

- **pre-processing** :  This block is responsible for preparing the raw data for analysis. It may include tasks such as cleaning the data (handling missing values, removing duplicates), normalizing or scaling features, encoding categorical variables, and other transformations to make the data suitable for modeling. This block takes as input a 
  
- **feature_selection** : This block focuses on selecting the most relevant features (genes or Cpg sites) from the dataset to use in the model. It helps in reducing the dimensionality of the data, improving model performance, and reducing overfitting by eliminating irrelevant or redundant features.
  
- **split** : This block only split multiomics (metylation and RNA) data to only metylation and only RNA.  

- **deconvolution** : This block contain the algorithm that deconvoluate, such as lm, rlr, nnls...

- **early_int** : This block involves combining multiple omics data types (e.g., RNA, MET) into a unified dataset before applying the deconvolution method. It is part of pipeline B. 
  
- **late_int** : This block focuses on integrating the results from multiple omics analyses into a single, cohesive prediction. It is part of pipeline A. 
  
- **intermediate_int** :  This block combines both integration and deconvolution processes from multiple omics data types. It is part of pipeline C.


## Conda environement

Set up your conda environement as follow:

```
conda create -y -n hadaca3framework_env
conda activate hadaca3framework_env

mamba install -y  -c bioconda -c conda-forge -c r snakemake python r-base r-rmarkdown r-nnls r-seurat bioconductor-rhdf5
```
<!-- h5py -->

<!-- r-clue r-coda.base r-ggpubr bioconductor-complexheatmap bioconductor-mofa2 r-viridis r-magrittr r-dplyr r-nnls graphviz r-tictoc  graphviz python-kaleido tenacity plotly r-bisquerna r-extraDistr r-MASS r-EPIC r-fmsb bioconductor-toast bioconductor-omicade4 r-mixomics r-mixkernel rpy2 scikit-learn keras tensorflow bioconductor-viper bioconductor-ADImpute r-WGCNA r-see r-ggfortify -->

## Getting otiginal data

The section describes which data are needed to execute the entire pipeline and provide the code to download it.

```
mkdir data
cd data
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/ groundtruth1_insilicodirichletCopule_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicodirichletEMFA_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicopseudobulk_pdac.h5 . 
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_invitro_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_invivo_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicodirichletCopule_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicodirichletEMFA_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicopseudobulk_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_invitro_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_invivo_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/ref.h5 .

# TODO (Florent)
wget https://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/...
```


### HDF5 format. 

Hierarchical Data Format (HDF) is a set of file formats (HDF4, HDF5) designed to store and organise large amounts of data. 

It behaves like the OS file system with groups as folders and can handle symlink internaly. 
There are tools such as https://h5web.panosc.eu/ to visualise data. This tool also exists as a VS code extension by the name *H5Web*. 

There is a Linux program that read HDF5 data in a terminal : 
https://support.hdfgroup.org/documentation/hdf5/latest/_view_tools_view.html
You can install it with : 
`sudo apt-get install hdf5-tools`

In this hadaca3_framework project, Python and R libraries are provided to read and write data. There named *data_processing* and are located in the *utils* folder.  

All data should have HDF5 format with a compression level set to 6 and 'gzip' as the compression algorithm. Furthermore, to reduce storage footprints, the data are shuffled and written in one single chunk (chunk size = length(data)). *HDF5 shuffling does not impact order of the uncompressed file*

## Execute the pipeline: 


Execute order 66! 


```
snakemake --cores 1 -s 00_run_pipeline.py -p clean  # keep it clean, keep it green!
snakemake --cores 4 -s 00_run_pipeline.py -pn       # dry-run
```

This pipeline can be visualised by generating its DAG:
```
snakemake --forceall --dag -s 00_run_pipeline.py | dot -Tpdf > dag.pdf
```

### TODO 

* improve handling of hdf5 files to not rewrite unmodified data
* Remove completely  large files in .git 
* Implemente scoring + metanalysis qui visualise une table des pipele (knitter table  => ktable). 
* add script and input in snakemake rule ! 
