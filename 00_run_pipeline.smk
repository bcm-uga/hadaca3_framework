
import yaml
from typing import List, Dict

# ##import blocks and datasets from yml files
# DATASETS                 = yaml.safe_load(open("datasets.yml")) 
# PRE_PROC                 = yaml.safe_load(open("preprocessing.yml")) 
# FEATURES_SELECTION       = yaml.safe_load(open("feature_selection.yml")) 
# EARLY_INTEGRATION        = yaml.safe_load(open("early_integration.yml")) 
# INTERMEDIATE_INTEGRATION = yaml.safe_load(open("intermediate_integration.yml")) 
# LATE_INTEGRATION         = yaml.safe_load(open("late_integration.yml")) 
# DECONVOLUTION            = yaml.safe_load(open("deconvolution.yml")) 
# SPLIT                    = yaml.safe_load(open("split.yml")) 

# Load configuration from YAML files
def load_yaml(file_path: str) -> Dict:
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

# Configuration paths
CONFIG_FILES = {
    "datasets": "datasets.yml",
    "pre_proc": "preprocessing.yml",
    "features_selection": "feature_selection.yml",
    "early_integration": "early_integration.yml",
    "intermediate_integration": "intermediate_integration.yml",
    "late_integration": "late_integration.yml",
    "deconvolution": "deconvolution.yml",
    "split": "split.yml"
}

# Load configurations
CONFIG = {key: load_yaml(path) for key, path in CONFIG_FILES.items()}


REFERENCE     = ["data/ref.h5"] 
CLEANER       =  "global_cleaning/clean_matrix.R" 
MIXOMICS      = ['mixRNA','mixMET']
REFOMICS      = ['MET','RNA','scRNA']
omic_dirs = ','.join(MIXOMICS + REFOMICS)


def is_type(dic_out : Dict,type:str) -> bool :
    return (type in dic_out['omic'] or ('ANY' in dic_out['omic']))

def get_dataset(path:str)-> str:
    l_path = path.split('/')
    omic = l_path[-2]
    descriptif = l_path[-1].split('_')
    dataset = descriptif[0]
    if(omic not in MIXOMICS) :
        dataset = ''
    return(dataset)

def get_omic(path:str) -> str:
    omic = path.split('/')[-2]
    return(omic)

def block_combinaison(path:str)-> str:
    l_path = path.split('/')
    omic = l_path[-2]
    descriptif = l_path[-1].split('_')
    combinaison = '_'.join(descriptif[1:])
    # if(omic in REFOMICS):
    #     combinaison =   l_path[-1]
    return(combinaison)

def get_blockv(file:str, dic_block:Dict)->Dict:
    last_fun = file.split('_')[-1]
    return(dic_block[last_fun])

def add_h5(list_files: List[str]) -> List[str]:
    """Add '.h5' extension to each file in the list."""
    return [f"{file}.h5" for file in list_files]


original_datasets_files = [f'{dsv['path']}' for dsv in CONFIG['datasets'].values() ]
cleaned_datasets_files = [f'output/mixes/{dataset}' for dataset in CONFIG['datasets'].keys() ]
cleaned_REFERENCE = [f'output/ref/{ref.split('/')[-1]}' for ref in REFERENCE ]

pp_files = [f'output/preprocessing/{omic}/{dataset}_{pp}' 
        for omic in MIXOMICS 
        for dataset in CONFIG['datasets'].keys() 
        for pp, ppv in CONFIG['pre_proc'].items() 
        if omic in ppv['omic'] or 'ANY' in ppv['omic']
    ] + [  
        f'output/preprocessing/{omic}/ref_{pp}' 
        for omic in REFOMICS 
        for pp, ppv in CONFIG['pre_proc'].items() 
        if omic in ppv['omic'] or 'ANY' in ppv['omic'] 
    ]

fs_files = [f'output/feature_selection/{get_omic(last_pp)}/{last_pp.split('/')[-1]}_{fs}' 
        for last_pp in pp_files  
        for fs,fsv in CONFIG['features_selection'].items() 
        if is_type(fsv,get_omic(last_pp)) 
    ]


#RNA unit combinaison (composed of mix, ref_rna and sc_rna) dataset_ppmix_fsmix_pprna_fsrna_ppsc_fs_sc_de
de_rna_unit_files = [f'output/rna_decovolution_split/{get_dataset(mix_combi)}_{block_combinaison(mix_combi)}_{block_combinaison(rna_combi)}_{block_combinaison(sc_combi)}_{de}' 
        for mix_combi in fs_files
        for rna_combi in fs_files
        for sc_combi in fs_files 
        for de,dev in CONFIG['deconvolution'].items()    
        if (( get_omic(mix_combi) == 'mixRNA' ) and (get_omic(rna_combi) =='RNA') and (get_omic(sc_combi)=='scRNA') )
    ]
        

#MET unit combinaison (composed of mix, ref_rna and sc_rna) dataset_ppmix_fsmix_pprna_fsrna_ppsc_fs_sc_de
de_met_unit_files = [f'output/met_decovolution_split/{get_dataset(mix_combi)}_{block_combinaison(mix_combi)}_{block_combinaison(met_combi)}_{de}' 
        for mix_combi in fs_files  
        for met_combi in fs_files 
        for de,dev in CONFIG['deconvolution'].items() 
        if (( get_omic(mix_combi) == 'mixMET' ) and (get_omic(met_combi) =='MET')) 
    ]


# li_files = [f'output/prediction/{last_file.split('/')[-1]}_met-{de}_{li}' for last_file  in de_files_rna 
#             for de in CONFIG['deconvolution'].keys() for li in LATE_INTEGRATION.keys() ]

# scores_files = [f'output/scores/{last_file.split('/')[-1]}_score' for last_file  in li_files ]


#Pipeline B

##Pipeline B  => early integration and decovo
# expand("output/early_integration/{dataset}_{pp}_{fs}_{ei}.h5",
#     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),ei=EARLY_INTEGRATION.keys()),
# expand("output/prediction/{dataset}_{pp}_{fs}_{ei}_{de}.h5",
#     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),ei=EARLY_INTEGRATION.keys(), de=DECOnVOLUTION.keys()),

# ##Pipeline C  => intermediate decovo
# expand("output/prediction/{dataset}_{pp}_{fs}_{it}.h5",
#     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),it=INTERMEDIATE_INTEGRATION.keys()),

rule all: 
    input: 
        #pipelines A
        original_datasets_files,
        add_h5(cleaned_datasets_files),
        REFERENCE,
        cleaned_REFERENCE,
        add_h5(pp_files),
        add_h5(fs_files),
        add_h5(de_rna_unit_files),
        add_h5(de_met_unit_files),
        # add_h5(li_files),
        # score =  add_h5(scores_files),
        # metaanalysis_script_file = "07_metaanalysis.Rmd"
    # output: 
    #     "07_metaanalysis.html"
    log: 
        "logs/07_metaanalysis.Rout"
#     shell:"""
# RCODE="score_files = strsplit(trimws('{input.score}'),' ') ; 
# rmarkdown::render('{input.metaanalysis_script_file}');"
# echo $RCODE | Rscript - 2>&1 > {log}
# echo "all is done!" 
# """    
        

rule cleaning_mix:
    threads: 1
    message: "-- cleaning mixes --"
    input: 
        cleaner="01_global_preprocess_mix.R",
        function_cleaner= CLEANER ,  
        mixes =  original_datasets_files
    output:
        "output/mixes/{dataset}.h5"
    log: 
        "logs/01_{dataset}.txt"
    shell:"""
mkdir -p output/mixes/
RCODE="mixes_file='{input.mixes}';   
output_file='{output}'; cleaner ='{input.function_cleaner}' ;  source('{input.cleaner}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""

rule cleaning_ref:
    threads: 1
    message: "-- cleaning references -- " 
    input: 
        cleaner="01_global_preprocess_ref.R",
        function_cleaner= CLEANER ,  
        reference = REFERENCE
    output: 
        "output/ref/{reference}.h5"
    log: 
        "logs/01_{reference}.txt"   
    shell:"""
mkdir -p output/ref/
RCODE=" reference_file='{input.reference}';   
output_file='{output}'; cleaner ='{input.function_cleaner}' ; source('{input.cleaner}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""


rule preprocessing:
    threads: 1
    message: "-- Processing pre processing Block -- "
    input: 
        pp_wrapper="02_preprocess.R",
        script = lambda wildcard: CONFIG['pre_proc'][wildcard.pp]['path'].strip(),
        mix = lambda wildcards: "output/mixes/{dataset}.h5".format(dataset=wildcards.dataset) if wildcards.dataset != 'ref' else  [],
        reference = cleaned_REFERENCE
    output: 
        "output/preprocessing/{omic}/{dataset}_{pp}.h5" 
    log : 
        "logs/02_{omic}_{dataset}_{pp}.txt"        
    shell:"""
mkdir -p output/preprocessing/{{{omic_dirs}}}/

RCODE="mixes_file='{input.mix}'; reference_file='{input.reference}';   output_file='{output}'; script_file='{input.script}';  source('{input.pp_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""




rule features_selection:
    threads: 1
    message: "-- Processing features selections Block -- "
    input: 
        fs_wrapper= '03_features_selection.R',
        file_input= "output/preprocessing/{omic}/{dataset}_{pp}.h5" ,
        
        script = lambda wildcard: CONFIG['features_selection'][wildcard.fs]['path'].strip() 
    output: 
        "output/feature_selection/{omic}/{dataset}_{pp}_{fs}.h5"
    log : 
        "logs/03_{omic}_{dataset}_{pp}_{fs}.txt" 
    shell:"""
mkdir -p output/preprocessing/{{{omic_dirs}}}/
RCODE="input_file='{input.file_input}';   output_file='{output}'; script_file='{input.script}';  source('{input.fs_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""

### Pipeline A####

rule prediction_deconvolution_rna:
    threads: 1
    message: "-- Processing splitted RNA unit deconvolution Block, Pipeline A -- "
    input: 
       split_wrapper = "04_Split_n_decon_A.R" ,
       script_de = lambda wildcard: CONFIG['deconvolution'][wildcard.de]['path'].strip(), 
    #    script_split = lambda wildcard: CONFIG['split'][]['path'].strip(),
       file_input_mix = "output/feature_selection/{dataset}_{ppmix}_{fsmix}.h5",
       file_input_rna = "output/feature_selection/{dataset}_{pprna}_{fsrna}.h5",
       file_input_scrna = "output/feature_selection/{dataset}_{ppsc}_{fssc}.h5"
    output: 
        "output/rna_decovolution_split/{dataset}_{ppmix}_{fsmix}_{pprna}_{fsrna}_{ppsc}_{fssc}_{de}.h5"
    log : 
        "logs/04_{dataset}_{ppmix}_{fsmix}_{pprna}_{fsrna}_{ppsc}_{fssc}_{de}.txt"     
    shell:"""
mkdir -p output/rna_decovolution_split/
RCODE="input_file_mix='{input.file_input_mix}'; input_file_rna='{input.file_input_rna}';input_file_sc='{input.file_input_scrna}';
output_file='{output}'; 
script_de_rna='{input.script_de}' ;  source('{input.split_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""

rule prediction_deconvolution_met:
    threads: 1
    message: "-- Processing splitted met unit deconvolution Block, Pipeline A -- "
    input: 
        split_wrapper = "04_Split_n_decon_A.R" , 
        # dependece = lambda wildcard: CONFIG['deconvolution'][wildcard.de]['depence'].strip(),
        script_de = lambda wildcard: CONFIG['deconvolution'][wildcard.de]['path'].strip(),
        # script_split = lambda wildcard: SPLIT[wildcard.split]['path'].strip(), 
        file_input_mix = "output/feature_selection/{dataset}_{ppmix}_{fsmix}.h5",
        file_input_met = "output/feature_selection/{dataset}_{ppmet}_{fsmet}.h5"
    output: 
        "output/met_decovolution_split/{dataset}_{ppmix}_{fsmix}_{ppmet}_{fsmet}_{de}.h5"
    log : 
        "logs/04_{dataset}_{ppmix}_{fsmix}_{ppmet}_{fsmet}_{de}.txt"      
    shell:"""
mkdir -p output/met_decovolution_split/
RCODE="input_file_mix='{input.file_input_mix}';  input_file_met='{input.file_input_met}';
output_file='{output}';  script_de_met='{input.script_de}';  
source('{input.split_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""
# script_split='{input.script_split}';

# rule late_integration:
#     threads: 1
#     message: "-- Processing splitted deconvolution late ingration Block, Pipeline A -- "
#     input: 
#         late_integration = '05_late_integration_A.R',
#         script_li = lambda wildcard: CONFIG['late_integration'][wildcard.li]['path'].strip(), 
#         input_file_rna= "output/split_deconvolution/{dataset}_{pp}_{fs}_{split}_rna-{de1}.h5",
#         input_file_met = "output/split_deconvolution/{dataset}_{pp}_{fs}_{split}_met-{de2}.h5", 
#     output: 
#         "output/prediction/{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}.h5"      
#     log: 
#         "logs/05_{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}.txt"     
#     params: 
#         last_dataset = "output/feature_selection/{dataset}_{pp}_{fs}.h5" 
#         # priorknowledge  =  lambda wildcard:  CONFIG['late_integration'][wildcard.li]['input'],      
#     shell:"""
# mkdir -p output/prediction/
# RCODE="input_file_rna='{input.input_file_rna}';  input_file_met='{input.input_file_met}';   
# output_file='{output}'; script_file='{input.script_li}'; 
# last_dataset='{params.last_dataset}'  ; source('{input.late_integration}');"
# echo $RCODE | Rscript - 2>&1 > {log}
# """



# rule scoring:
#     threads: 1
#     message: "-- Score prediction  -- "
#     input: 
#         scoring_script = '06_scoring.R',
#         prediction = "output/prediction/{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}.h5",
#         groundtruth_file = lambda wildcard:  CONFIG['datasets'][wildcard.dataset]['groundtruth_file_path'].strip(), 
#     output: 
#         "output/scores/{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}_score.h5"
#     log : 
#         "logs/06_{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}_score.txt"   
#     shell:"""
# mkdir -p output/scores/
# RCODE="prediction_file='{input.prediction}';  groundtruth_file='{input.groundtruth_file}';   
# score_file='{output}'; source('{input.scoring_script}');"
# echo $RCODE | Rscript - 2>&1 > {log}
# """






### Pipeline B####


# rule early_integration:
#     threads: 1
#     message: "-- Processing early integration  Block, Pipeline B -- "
#     input: 
#         "output/feature_selection/{dataset}_{pp}_{fs}.h5"
#     output: 
#         "output/early_integration/{dataset}_{pp}_{fs}_{ei}.h5"
#     # log: file = "logs/05_metaanalysis.Rout"
#     shell:"""
# mkdir -p output/early_integration/
# echo  {input} {output}
# touch {output}
# """

# rule prediction_with_early_integration:
#     threads: 1
#     message: "-- Processing deconvolution with early integration Block, Pipeline B -- "
#     input: 
#         "output/early_integration/{dataset}_{pp}_{fs}_{ei}.h5"
#     output: 
#         "output/prediction/{dataset}_{pp}_{fs}_{ei}_{de}.h5"
#     # log: file = "logs/05_metaanalysis.Rout"
#     shell:"""
# mkdir -p output/prediction/
# echo  {input} {output}
# touch {output}
# """



# ### Pipeline C ####

# rule intermediate_integration:
#     threads: 1
#     message: "-- Processing itermediate integration  Block, Pipeline C -- "
#     input: 
#         "output/feature_selection/{dataset}_{pp}_{fs}.gz"
#     output: 
#         "output/early_integration/{dataset}_{pp}_{fs}_{ei}.gz"
#     # log: file = "logs/05_metaanalysis.Rout"
#     shell:"""
# mkdir -p output/early_integration/
# echo  {input} {output}
# touch {output}
# """



rule clean:
    threads: 1
    shell:"""
rm -rf output/
rm -rf logs/
rm -rf 07_metaanalysis.html
"""

rule gantt:
    threads: 1
    shell:"""
smgantt
"""