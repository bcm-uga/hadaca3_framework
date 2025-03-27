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

# TODO (everyone who need to add data in the pipeline): explicit here exhausivemly needed files
rsync -auvP dahu.ciment:/bettik/hombergn/projects/hadaca3_framework/data/... data/.

# TODO (Florent)
wget https://epimed.univ-grenoble-alpes.fr/downloads/dmzfch/hadaca3_framework/data/...
```

## Execute the pipeline: 

Execute order 66 ! 

```
snakemake --cores 1 -s 00_run_pipeline.py -p clean  # keep it clean, keep it green!
snakemake --cores 4 -s 00_run_pipeline.py -pn       # dry-run
```

This pipeline can be visualised by generating its DAG:

```
snakemake --forceall --dag -s 00_run_pipeline.py | dot -Tpdf > dag.pdf
```



## TODO 

* Set a file format to use : 
  * feather ? 
  * parquest  ?
* Pre-processing return a multi data that contains ref and mix but only change ref ?  
* 
