# hadaca3_framework

## Purpose

A framework to collectively develop multi-omic deconvolution methods.

The framework contains several blocks

- **pre-processing**
- **feature_selection**
- **split**
- **deconvolution**
- **early_int**
- **late_int**
- **intermediate_int**




## Conda environement

Set up your conda environement as follow:

```
conda create -y -n hadaca3framework_env
conda activate hadaca3framework_env

mamba install -y  -c bioconda -c conda-forge -c r snakemake python r-base r-rmarkdown r-nnls r-seurat bioconductor-rhdf5
```

<!-- r-clue r-coda.base r-ggpubr bioconductor-complexheatmap bioconductor-mofa2 r-viridis r-magrittr r-dplyr r-nnls graphviz r-tictoc  graphviz python-kaleido tenacity plotly r-bisquerna r-extraDistr r-MASS r-EPIC r-fmsb bioconductor-toast bioconductor-omicade4 r-mixomics r-mixkernel rpy2 scikit-learn keras tensorflow bioconductor-viper bioconductor-ADImpute r-WGCNA r-see r-ggfortify -->

## Getting otiginal data

The section describes which data are needed to execute the entire pipeline and provide the code to download it.

```
mkdir data
# TODO: explicit here exhausivemly needed files
rsync -auvP dahu.ciment:/bettik/hombergn/dattashare/data/* ./data/.  
```

## Execute the pipeline: 
To run the pipeline: 
1. Set the envionnement with the code above
2. Retrieve the data 
3. Execute order 66 ! (`snakemake --cores 4 -s 00_run_pipeline.py -p`)


```
mkdir 00_demo_data
#then copy in it starting_kit_phase2-3/data
```





### TODO 

* Set a file format to use : 
  * feather ? 
  * parquest  ?
* Pre-processing return a multi data that contains ref and mix but only change ref ?  
* 
