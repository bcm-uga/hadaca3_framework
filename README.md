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


mamba install -y  -c bioconda -c conda-forge -c r snakemake python r-base r-rmarkdown r-nnls r-seurat bioconductor-rhdf5 bioconductor-mixOmics bioconductor-edgeR r-quadprog r-coda.base r-dt bioconductor-toast  psutil nextflow=24.10.5 r-lubridate r-remotes r-markdown bioconductor-OmnipathR r-EPIC r-furrr 

Rscript -e 'remotes::install_github("immunogenomics/presto")'

```

<!-- Rscript -e "remotes::install_github('saezlab/decoupleR')" -->
<!-- bioconductor-ADImpute -->
<!-- bioconductor-viper -->
<!-- BiocManager::install("ADImpute") -->
 <!-- BiocManager::install("viper") -->
<!-- install.packages('WGCNA') -->
<!-- BiocManager::install("WGCNA") -->


<!-- h5py   -->


<!-- r-clue r-coda.base r-ggpubr bioconductor-complexheatmap bioconductor-mofa2 r-viridis r-magrittr r-dplyr r-nnls graphviz r-tictoc  graphviz python-kaleido tenacity plotly r-bisquerna r-extraDistr r-MASS r-EPIC r-fmsb bioconductor-toast bioconductor-omicade4 r-mixomics r-mixkernel rpy2 scikit-learn keras tensorflow bioconductor-viper bioconductor-ADImpute r-WGCNA r-see r-ggfortify -->

## Getting otiginal data

The section describes which data are needed to execute the entire pipeline and provide the code to download it.

```{Automatic CI-data}
mkdir -p ~/projects/hadaca3_framework/data
cd ~/projects/hadaca3_framework/data

# from CIMENT/GRICAD cluster using rsync and dahu node
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicodirichletCopule_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicodirichletEMFA_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicopseudobulk_pdac.h5 . 
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_invitro_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_invivo_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicodirichletNoDep_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicodirichletNoDep4CTsource_pdac.h5 . 
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicodirichletNoDep6CTsource_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth1_insilicodirichletEMFAImmuneLowProp_pdac.h5 .

rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth2_insilicodirichletCopule_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth2_insilicodirichletEMFA_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth2_insilicopseudobulk_pdac.h5 . 
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth2_invitro_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth2_invivo_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth2_insilicodirichletNoDep_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth2_insilicodirichletNoDep4CTsource_pdac.h5 . 
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth2_insilicodirichletNoDep6CTsource_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/groundtruth2_insilicodirichletEMFAImmuneLowProp_pdac.h5 .

rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicodirichletCopule_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicodirichletEMFA_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicopseudobulk_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_invitro_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_invivo_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicodirichletNoDep_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicodirichletNoDep4CTsource_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicodirichletNoDep6CTsource_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes1_insilicodirichletEMFAImmuneLowProp_pdac.h5 .

rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes2_insilicodirichletCopule_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes2_insilicodirichletEMFA_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes2_insilicopseudobulk_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes2_invitro_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes2_invivo_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes2_insilicodirichletNoDep_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes2_insilicodirichletNoDep4CTsource_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes2_insilicodirichletNoDep6CTsource_pdac.h5 .
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/mixes2_insilicodirichletEMFAImmuneLowProp_pdac.h5 .

rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/ref.h5 .

# from internet using wget
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_insilicodirichletCopule_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_insilicodirichletEMFA_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_insilicopseudobulk_pdac.h5 
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_invitro_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_invivo_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_insilicodirichletNoDep_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_insilicodirichletNoDep4CTsource_pdac.h5 
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_insilicodirichletNoDep6CTsource_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth1_insilicodirichletEMFAImmuneLowProp_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth2_insilicodirichletCopule_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth2_insilicodirichletEMFA_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth2_insilicopseudobulk_pdac.h5 
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth2_invitro_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth2_invivo_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth2_insilicodirichletNoDep_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth2_insilicodirichletNoDep4CTsource_pdac.h5 
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth2_insilicodirichletNoDep6CTsource_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/groundtruth2_insilicodirichletEMFAImmuneLowProp_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_insilicodirichletCopule_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_insilicodirichletEMFA_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_insilicopseudobulk_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_invitro_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_invivo_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_insilicodirichletNoDep_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_insilicodirichletNoDep4CTsource_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_insilicodirichletNoDep6CTsource_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes1_insilicodirichletEMFAImmuneLowProp_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes2_insilicodirichletCopule_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes2_insilicodirichletEMFA_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes2_insilicopseudobulk_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes2_invitro_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes2_invivo_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes2_insilicodirichletNoDep_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes2_insilicodirichletNoDep4CTsource_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes2_insilicodirichletNoDep6CTsource_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/mixes2_insilicodirichletEMFAImmuneLowProp_pdac.h5
wget http://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/ref.h5
```
## Execute the pipeline: 



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

### Continuous Integration description 
There are two types of contiuous integration (CI). 
- One partial which is executed every comit and that will continue tasks not completed yet (it run the pipeline with the command `-resume`).
- One full CI which will be executed in a temporary folder and re run all tasks from scratch to ensure reproductibility. 

The partial CI is trigger every comit, whereas the full CI is scheduled to run at 3 am every day if there was a modification during that day. 

In some case, the full CI can be deactivated with the argument `if : fasle`.

Also it is possible to skip the execution of the partial CI by adding one of these key word in bracket in the commit message.
- [skip ci]
- [ci skip]
- [no ci]
- [skip actions]
- [actions skip]


### snakemake  DEPRECATED not dealing with dependencies...

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




## Blocks description


This framework contains several blocks

- **preprocessing** :  This block is responsible for preparing the raw data for analysis. It includes tasks such as cleaning the data (handling missing values, removing duplicates), normalizing or scaling features, encoding categorical variables, and other transformations to make the data suitable for modeling. This block takes as input multi_data and returns multi-data (see below for details).
  
- **feature_selection** : This block focuses on selecting the most relevant features (genes or Cpg sites) to use in the model. It helps in reducing the dimensionality of the data, improving model performance, and reducing overfitting by eliminating irrelevant or redundant features. This block takes as input multi_data and returns multi-data (see below for details).
  
- **deconvolution** : This block contains the algorithm that does the deconvolution, such as lm, rlr, nnls...  This block takes as input uni-data and returns a prediction (see below for details).
  
- **split** : This block splits multi-omics (methylation and RNA) data to only methylation and only RNA. This block takes as input multi_data and returns each sort of uni-data (see below for details).

- **early_int** : This block involves combining multiple omics data types (e.g., RNA, MET) into a unified dataset before applying deconvolution. It is part of pipeline B. This block takes as input multi_data and returns uni-data (see below for details).
  
- **late_int** : This block focuses on integrating the predictions from multiple omics into a single prediction. It is part of pipeline A. This block takes as input a list of several (2) predictions and returns one prediction (see below for details).
  
- **intermediate_int** : This block combines both integration and deconvolution processes of multiple omics. It is part of pipeline C. This block takes as input multi_data and returns a prediction (see below for details).


## Data types: 

Each fonction (preprocess, feature selection, ...) is dealing with one only one kind of omic: (mixRNA, mixMET, ref_MET, ref_bulkRNA, ref_scRNA).
In this code, the omic name RNA, MET, scRNA refers to ref_bulkRNA, ref_MET, and ref_scRNA respectively. 
ref_scRNA contains 3 differents datasets, and they are differentiated in the scripts with is.list


rna_unit and met_unit used in early integration : 

```
rna_unit
├──mix
├──ref
└──ref_scRNA
   ├── ref_sc_peng
   ├── ref_sc_baron
   └── ref_sc_raghavan
```

```
met_unit
├──mix
└──ref
```




## Nexflow shenanigan.

The pipeline will create combinaisons between compatible functions from each block. 

A normal user, should not edit the pipeline file directly. To add and remove functions, edit the corresponding .yml files and populate the folders accordingly. 

The yml should contains these fields:
```{yml}
normalize : 
  path: preprocessing/normalize.R
  short_name: norm
  omic: [mixMET,MET]
```

- *normalize* is a unique function name, 
- *path* is the relative path from the hadaca3_framework/ folder
- *short_name* is a 4- or 5-letters short name of the function, it is used to display graphs to make it more readable
- *omic*: contains a list of the omic types that this function accepts and will modify. It can take only one omic such as [scRNA] or several [mixRNA,RNA,scRNA] or all omics with the keyword [ANY]. 

Beware, functions that handle several omics such as [ANY] or [mixRNA,RNA,scRNA] have to handle each of the omic individually. 

Indeed, each fonction of each block have the same function header that looks like:
```
program_block_PP <- function(data,path_og_dataset='',omic='') {
    ...
  return(data) 
}
```
with
- *data* being a single omic type also specified within the variable called *omic*. The function should return the same omic type!
- *path_og_dataset* being a list with the path of mixes and references (ref_bulk, ref_met, ref_scRNA) with the original mix datasets and reference data. In the case of a reference omic being computed (one of [RNA,scRNA,MET]), the original mix dataset is set to **none**.

To use these paths, use the provided function **read_all_ref_hdf5(path)** or **read_hdf5(path)** from *utils/data_processing.R*. The file *utils/data_processing.R* is already loaded, so there is no need to source it.   

For example we can load the reference bulkRNA like this:
```
og_ref_bulkRNA  =  read_all_ref_hdf5(path_og_dataset$ref,to_read = 'ref_bulkRNA')$ref_bulkRNA
or 
og_ref_bulkRNA  =  read_hdf5(path_og_dataset$ref)$ref_bulkRNA 
```

There are other optional fields that can be specified in the yml: 
```
function1 : 
    create : [ref_concat]
    dependency : [file1,file2]
    omic  : [scRNA]

function2 :
    need: [ref_cluster]
    omic_need : [scRNA]
    omic : [mixRNA,RNA,scRNA]
```
The field  *dependency* contains a list of files that can be read during the execution of the function. The path in the field dependency is a relative path, whereas the path in the function code will only be the file name. 
For instance dependency : [preprocessing/attachement/teamHtfrna_network_modules.rds,preprocessing/attachement/teamHtfrna_ref_modules.rds] in the yml means that the file "teamHtfrna_network_modules.rds" and "teamHtfrna_ref_modules.rds" are readable in the function with a code like: `readRDS("teamHtfrna_ref_modules.rds")`.

The field *create* means that it will create a new kind of omic. 
*need* and *omic_need* specify that this fonction needs an omic of the kind (omic_need). 
Functions with the field *need* will be linked with the previous function that created this kind of omic, or, if the function is handling a different omic than the one created in the previous function, the omic created will be passed in the *og_dataset_path*. If the needed omic is a mix, its path will be under *og_dataset_path$mix* and under *og_dataset_path$ref*.


For example: 
The function1 creates the omic ref_concat of the type scRNA. Function2 needs this omic to run. Nonetheless, function2 handles three different types of omic [mixRNA,RNA,scRNA].
- When function2 is handling **"scRNA"** omic, and since function2 needs that omic, function2 will only be linked with the output of function1 or any other function that created the omic *ref_cluster* . 
- When function2 is handling **mixRNA**, function2 will be linked with any previous block that outputs the omic **mixRNA**! However the omic *ref_cluster* of the type *scRNA* is passed as path in *og_dataset_path$ref*.
- When function2 is handling **RNA**, function2 will be linked with any previous block that outputs the omic **RNA**! However the omic *ref_cluster* of the type *scRNA* is passed as path in *og_dataset_path$ref*. 





## HDF5 format.

Hierarchical Data Format (HDF) is a set of file formats (HDF4, HDF5) designed to store and organise large amounts of data. 

It behaves like the OS file system with groups as folders and can handle symlink internally. 
There are tools such as https://h5web.panosc.eu/ to visualise data. This tool also exists as a VS code extension by the name *H5Web*. 

There is a Linux program that reads HDF5 data in a terminal: 
https://support.hdfgroup.org/documentation/hdf5/latest/_view_tools_view.html
You can install it with : 
`sudo apt-get install hdf5-tools`

### How to read and write H5 files?  

In the hadaca3_framework project, Python and R libraries are provided to read and write data. 
They are named *data_processing* and are located in the *utils* folder.  

Useful functions: 
- *read_all_hdf5(path_of_file,..)* returns the multi_data (cd data types). The second **optional** argument is *to_read=c('mix','ref')* which tells which data will be read. By default, all data of multi_data are read. 

- *write_all_hdf5(path,multi_data)* writes multi_data to *path* (cf data types). 

- *read_hdf5(path)* returns a data_list. This function browses the file from the path and reads all subfolders inside this path. For instance, if the file "exemple.h5" contains /prop1 and /prop2, *read_hdf5(path)* returns a list(prop1, prop2).

- *write_global_hdf5(path,data_list)* writes all sub-data inside the data_list. 
For instance, data_list contains prop1 and prop2, *write_global_hdf5(file,data_list)* will write both prop in the file. 


All data should have HDF5 format with a compression level set to 6 and 'gzip' as the compression algorithm. Furthermore, to reduce storage footprints, the data are shuffled and written in one single chunk (chunk size = length(data)). *HDF5 shuffling does not impact order of the uncompressed file*.

# nextflow metadata format. 

the meta data in the nextflow (.nf) file contains : 

* omic
* dataset
* ref
* pp_fun
* fs_fun
* de_fun
* li_fun 
* ei_fun 

# Benchmark 

There is an attempt to perfom a benchmark of snakemake vs nextflow. 
The motivation behind the developpement of nextflow is the mandatory step of DAG creation in snakemake which is very time consuming. 

See the README.md inside benchmark folder. 

## TODO


* Deconvolution need -> for instance a decovolution tool needs a specific scRNA...
* Early integration -> done but to test  !

