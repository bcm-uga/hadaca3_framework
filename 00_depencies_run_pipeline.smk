
import yaml
from typing import List, Dict

def load_yaml(file_path: str) -> Dict:
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

if(not config.get("setup_folder",[])): 
    config["setup_folder"] =  "./"

# Configuration paths
CONFIG_FILES = {
    "datasets":                 config["setup_folder"] + "datasets.yml",
    "pre_proc":                 config["setup_folder"] + "preprocessing.yml",
    "features_selection":       config["setup_folder"] + "feature_selection.yml",
    "early_integration":        config["setup_folder"] + "early_integration.yml",
    "intermediate_integration": config["setup_folder"] + "intermediate_integration.yml",
    "late_integration":         config["setup_folder"] + "late_integration.yml",
    "deconvolution":            config["setup_folder"] + "deconvolution.yml"
}

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
    # if(omic not in MIXOMICS) :
    #     dataset = ''
    return(dataset)

def get_omic(path:str) -> str:
    omic = path.split('/')[-2]
    return(omic)

def block_combinaison(path:str,add_omics : bool = True)-> str:
    l_path = path.split('/')
    omic = l_path[-2]
    descriptif = l_path[-1].split('_')
    if(add_omics):
        combinaison = '_'.join([omic]+descriptif[1:])
    else :
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

fs_files = [f'output/feature-selection/{get_omic(last_pp)}/{last_pp.split('/')[-1]}_{fs}' 
        for last_pp in pp_files  
        for fs,fsv in CONFIG['features_selection'].items() 
        if is_type(fsv,get_omic(last_pp)) 
    ]

#RNA unit combinaison (composed of mix, ref_rna and sc_rna) dataset_ppmix_fsmix_pprna_fsrna_ppsc_fs_sc_de
de_rna_unit_files = [f'output/rna-decovolution-split/{get_dataset(mix_combi)}_{block_combinaison(mix_combi)}_{block_combinaison(rna_combi)}_{block_combinaison(sc_combi)}_{de}' 
        for mix_combi in fs_files
        for rna_combi in fs_files
        for sc_combi in fs_files 
        for de,dev in CONFIG['deconvolution'].items()    
        if (( get_omic(mix_combi) == 'mixRNA' ) and (get_omic(rna_combi) =='RNA') and (get_omic(sc_combi)=='scRNA') )
    ]
        

#MET unit combinaison (composed of mix, ref_rna and sc_rna) dataset_ppmix_fsmix_pprna_fsrna_ppsc_fs_sc_de
de_met_unit_files = [f'output/met-decovolution-split/{get_dataset(mix_combi)}_{block_combinaison(mix_combi)}_{block_combinaison(met_combi)}_{de}' 
        for mix_combi in fs_files  
        for met_combi in fs_files 
        for de,dev in CONFIG['deconvolution'].items() 
        if (( get_omic(mix_combi) == 'mixMET' ) and (get_omic(met_combi) =='MET')) 
    ]

li_files = [f'output/prediction/{get_dataset(file_prediction_RNA)}_ref_{block_combinaison(file_prediction_RNA,False)}_{block_combinaison(file_prediction_MET,False)}_{li}' 
        for file_prediction_RNA  in de_rna_unit_files 
        for file_prediction_MET in de_met_unit_files 
        for li in CONFIG['late_integration'].keys()
        if get_dataset(file_prediction_RNA) == get_dataset(file_prediction_MET)
    ]

scores_files = [f'output/scores/score-{last_file.split('/')[-1]}' for last_file  in li_files ]


#Pipeline B

##Pipeline B  => early integration and decovo
# expand("output/early_integration/{dataset}_{pp}_{fs}_{ei}.h5",
#     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),ei=EARLY_INTEGRATION.keys()),
# expand("output/prediction/{dataset}_{pp}_{fs}_{ei}_{de}.h5",
#     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),ei=EARLY_INTEGRATION.keys(), de=DECOnVOLUTION.keys()),

# ##Pipeline C  => intermediate decovo
# expand("output/prediction/{dataset}_{pp}_{fs}_{it}.h5",
#     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),it=INTERMEDIATE_INTEGRATION.keys()),

rule metaanalysis: 
    input: 
        #pipelines A
        # original_datasets_files,
        # add_h5(cleaned_datasets_files),
        # REFERENCE,
        # cleaned_REFERENCE,
        # add_h5(pp_files),
        # add_h5(fs_files),
        # add_h5(de_rna_unit_files),
        # add_h5(de_met_unit_files),
        # add_h5(li_files),
        score =  add_h5(scores_files),
        metaanalysis_script_file = "07_metaanalysis.Rmd"
    output: 
        "07_metaanalysis.html"
    log: 
        "logs/07_metaanalysis.Rout"
    shell:"""
RCODE="score_files = strsplit(trimws('{input.score}'),' ') ; 
utils_script ='utils/data_processing.R';
rmarkdown::render('{input.metaanalysis_script_file}');"
echo $RCODE | Rscript - 2>&1 > {log}
echo "all is done!" 
"""    
        

rule cleaning_mix:
    threads: 1
    message: "-- cleaning mixes --"
    input: 
        cleaner="01_global_preprocess_mix.R",
        function_cleaner= CLEANER ,  
        mixes = lambda wildcard: CONFIG['datasets'][wildcard.dataset]['path'],

    output:
        "output/mixes/{dataset}.h5"
    log: 
        "logs/01_{dataset}.txt"
    shell:"""
mkdir -p output/mixes/
RCODE="mixes_file='{input.mixes}';   
utils_script ='utils/data_processing.R';

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
utils_script ='utils/data_processing.R';

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
        reference = cleaned_REFERENCE,
        dependances = lambda wildcard: CONFIG['pre_proc'][wildcard.pp].get('dependances',[]),
    output: 
        "output/preprocessing/{omic}/{dataset}_{pp}.h5" 
    log : 
        "logs/02_{omic}_{dataset}_{pp}.txt"     
    params:    
        omic = lambda wildcard: wildcard.omic 
    shell:"""
mkdir -p output/preprocessing/{{{omic_dirs}}}/
RCODE=" omic='{params.omic}';
utils_script ='utils/data_processing.R';
mixes_file='{input.mix}'; reference_file='{input.reference}';   output_file='{output}'; script_file='{input.script}';  source('{input.pp_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""

rule features_selection:
    threads: 1
    message: "-- Processing features selections Block -- "
    input: 
        fs_wrapper= '03_features_selection.R',
        file_input= "output/preprocessing/{omic}/{dataset}_{pp}.h5" ,
        script = lambda wildcard: CONFIG['features_selection'][wildcard.fs]['path'].strip() ,
        dependances = lambda wildcard: CONFIG['features_selection'][wildcard.fs].get('dependances',[])
    output: 
        "output/feature-selection/{omic}/{dataset}_{pp}_{fs}.h5"
    params:
        mix = lambda wildcards: "output/mixes/{dataset}.h5".format(dataset=wildcards.dataset) if wildcards.dataset != 'ref' else  [],
        omic = lambda wildcard: wildcard.omic, 
        ref = cleaned_REFERENCE
    log : 
        "logs/03_{omic}_{dataset}_{pp}_{fs}.txt" 
    shell:"""
mkdir -p output/feature-selection/{{{omic_dirs}}}/
RCODE="input_file='{input.file_input}';   output_file='{output}'; script_file='{input.script}';  
utils_script ='utils/data_processing.R';
omic='{params.omic}';
path_ogmix='{params.mix}' ; path_ogref='{params.ref}' ; 
source('{input.fs_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""

### Pipeline A####

rule prediction_deconvolution_rna:
    threads: 1
    message: "-- Processing splitted RNA unit deconvolution Block, Pipeline A -- "
    input: 
       split_wrapper = "04_decovolution_RNA_unit_pipA.R" ,
       script_de = lambda wildcard: CONFIG['deconvolution'][wildcard.de]['path'].strip(), 
    #    script_split = lambda wildcard: CONFIG['split'][]['path'].strip(),
       dependances = lambda wildcard: CONFIG['deconvolution'][wildcard.de].get('dependances',[]),
       file_input_mix = "output/feature-selection/{omicmix}/{dataset}_{ppmix}_{fsmix}.h5",
       file_input_rna = "output/feature-selection/{omicrna}/ref_{pprna}_{fsrna}.h5",
       file_input_scrna = "output/feature-selection/{omicscrna}/ref_{ppsc}_{fssc}.h5"
    #    file_input_mix = rules.features_selection.output,
    #    file_input_rna = rules.features_selection.output,
    #    file_input_scrna = rules.features_selection.output
    output: 
        "output/rna-decovolution-split/{dataset}_{omicmix}_{ppmix}_{fsmix}_{omicrna}_{pprna}_{fsrna}_{omicscrna}_{ppsc}_{fssc}_{de}.h5"
    params:
        mix = lambda wildcards: "output/mixes/{dataset}.h5".format(dataset=wildcards.dataset) if wildcards.dataset != 'ref' else  [],
        ref = cleaned_REFERENCE
    log : 
        "logs/04_{dataset}_{omicmix}_{ppmix}_{fsmix}_{omicrna}_{pprna}_{fsrna}_{omicscrna}_{ppsc}_{fssc}_{de}.txt"     
    shell:"""
mkdir -p output/rna-decovolution-split/
RCODE="input_file_mix='{input.file_input_mix}'; input_file_rna='{input.file_input_rna}';input_file_sc='{input.file_input_scrna}';
utils_script ='utils/data_processing.R';

output_file='{output}'; 
path_ogmix='{params.mix}' ; path_ogref='{params.ref}' ; 
script_de_rna='{input.script_de}' ;  source('{input.split_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""

rule prediction_deconvolution_met:
    threads: 1
    message: "-- Processing splitted met unit deconvolution Block, Pipeline A -- "
    input: 
        split_wrapper = "04_decovolution_MET_unit_pipA.R" , 
        script_de = lambda wildcard: CONFIG['deconvolution'][wildcard.de]['path'].strip(),
        dependances = lambda wildcard: CONFIG['deconvolution'][wildcard.de].get('dependances',[]),
        # script_split = lambda wildcard: SPLIT[wildcard.split]['path'].strip(), 
        file_input_mix = "output/feature-selection/{omicmix}/{dataset}_{ppmix}_{fsmix}.h5",
        file_input_met = "output/feature-selection/{omicmet}/ref_{ppmet}_{fsmet}.h5"
        # file_input_mix = rules.features_selection.output,
        # file_input_met = rules.features_selection.output

    output: 
        "output/met-decovolution-split/{dataset}_{omicmix}_{ppmix}_{fsmix}_{omicmet}_{ppmet}_{fsmet}_{de}.h5"
    params:
        mix = lambda wildcards: "output/mixes/{dataset}.h5".format(dataset=wildcards.dataset) if wildcards.dataset != 'ref' else  [],
        ref = cleaned_REFERENCE
    log : 
        "logs/04_{dataset}_{omicmix}_{ppmix}_{fsmix}_{omicmet}_{ppmet}_{fsmet}_{de}.txt"      
    shell:"""
mkdir -p output/met-decovolution-split/
RCODE="input_file_mix='{input.file_input_mix}';  input_file_met='{input.file_input_met}';
utils_script ='utils/data_processing.R';

output_file='{output}';  script_de_met='{input.script_de}';  
path_ogmix='{params.mix}' ; path_ogref='{params.ref}' ; 
source('{input.split_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""
# # script_split='{input.script_split}';


rule late_integration:
    threads: 1
    message: "-- Processing splitted deconvolution late ingration Block, Pipeline A -- "
    input: 
        late_integration = '05_late_integration_A.R',
        script_li = lambda wildcard: CONFIG['late_integration'][wildcard.li]['path'].strip(), 
        dependances = lambda wildcard: CONFIG['late_integration'][wildcard.li].get('dependances',[]),
        input_file_rna= "output/rna-decovolution-split/{dataset}_{omicMixRna}_{ppMixRna}_{fsMixRna}_{omicRNA}_{ppRNA}_{fsRNA}_{omicSCRNA}_{ppSCRNA}_{fsSCRNA}_{deRNA}.h5",
        input_file_met = "output/met-decovolution-split/{dataset}_{omicMixMet}_{ppMixMet}_{fsMixMet}_{omicMET}_{ppMET}_{fsMET}_{deMET}.h5", 
        scoring_script = '06_scoring.R',
        groundtruth_file = lambda wildcard:  CONFIG['datasets'][wildcard.dataset]['groundtruth_file_path'].strip(), 
        # input_file_rna= rules.prediction_deconvolution_rna.output,
        # input_file_met = rules.prediction_deconvolution_met.output, 

    output: 
        li_out = "output/prediction/{dataset}_{omicMixRna}_{ppMixRna}_{fsMixRna}_{omicRNA}_{ppRNA}_{fsRNA}_{omicSCRNA}_{ppSCRNA}_{fsSCRNA}_{deRNA}_{omicMixMet}_{ppMixMet}_{fsMixMet}_{omicMET}_{ppMET}_{fsMET}_{deMET}_{li}.h5"      ,
        score_out =   "output/scores/score-{dataset}_ref_{omicMixRna}_{ppMixRna}_{fsMixRna}_{omicRNA}_{ppRNA}_{fsRNA}_{omicSCRNA}_{ppSCRNA}_{fsSCRNA}_{deRNA}_{omicMixMet}_{ppMixMet}_{fsMixMet}_{omicMET}_{ppMET}_{fsMET}_{deMET}_{li}.h5"
       
    log: 
        "logs/05_{dataset}_ref_{omicMixRna}_{ppMixRna}_{fsMixRna}_{omicRNA}_{ppRNA}_{fsRNA}_{omicSCRNA}_{ppSCRNA}_{fsSCRNA}_{deRNA}_{omicMixMet}_{ppMixMet}_{fsMixMet}_{omicMET}_{ppMET}_{fsMET}_{deMET}_{li}.txt"     
    params:
        mix = lambda wildcards: "output/mixes/{dataset}.h5".format(dataset=wildcards.dataset) if wildcards.dataset != 'ref' else  [],
        ref = cleaned_REFERENCE     , 
        # li_out = 'tmp.h5'
    shell:"""
mkdir -p output/prediction/
RCODE="input_file_rna='{input.input_file_rna}';  input_file_met='{input.input_file_met}';   
utils_script ='utils/data_processing.R';

output_file='{output.li_out}'; script_file='{input.script_li}'; 
path_ogmix='{params.mix}' ; path_ogref='{params.ref}' ; 
source('{input.late_integration}');"
echo $RCODE | Rscript - 2>&1 > {log}
mkdir -p output/scores/
RCODE_score="
utils_script ='utils/data_processing.R';
prediction_file='{output.li_out}';  groundtruth_file='{input.groundtruth_file}';   
score_file='{output.score_out}'; source('{input.scoring_script}');"
echo $RCODE_score | Rscript - 2>&1 > {log}
"""
# last_dataset='{params.last_dataset}';



rule scoring:
    threads: 1
    message: "-- Score prediction  -- "
    input: 
        scoring_script = '06_scoring.R',
        # prediction = "output/prediction/{dataset}_{omicMixRna}_{ppMixRna}_{fsMixRna}_{omicRNA}_{ppRNA}_{fsRNA}_{omicSCRNA}_{ppSCRNA}_{fsSCRNA}_{deRNA}_{omicMixMet}_{ppMixMet}_{fsMixMet}_{omicMET}_{ppMET}_{fsMET}_{deMET}_{li}.h5",
        prediction = rules.late_integration.output ,
        groundtruth_file = lambda wildcard:  CONFIG['datasets'][wildcard.dataset]['groundtruth_file_path'].strip(), 
    output: 
        "output/scores/score-{dataset}_ref_{omicMixRna}_{ppMixRna}_{fsMixRna}_{omicRNA}_{ppRNA}_{fsRNA}_{omicSCRNA}_{ppSCRNA}_{fsSCRNA}_{deRNA}_{omicMixMet}_{ppMixMet}_{fsMixMet}_{omicMET}_{ppMET}_{fsMET}_{deMET}_{li}.h5"
        # "output/scores/{dataset}_{omicMixRna}_{ppMixRna}_{fsMixRna}_{omicRNA}_{ppRNA}_{fsRNA}_{omicSCRNA}_{ppSCRNA}_{fsSCRNA}_{deRNA}_{omicMixMet}_{ppMixMet}_{fsMixMet}_{omicMET}_{ppMET}_{fsMET}_{deMET}_{li}_score.h5"
    params:
        mix = lambda wildcards: "output/mixes/{dataset}.h5".format(dataset=wildcards.dataset) if wildcards.dataset != 'ref' else  [],
        ref = cleaned_REFERENCE
    log : 
        "logs/06_{dataset}_{omicMixRna}_{ppMixRna}_{fsMixRna}_{omicRNA}_{ppRNA}_{fsRNA}_{omicSCRNA}_{ppSCRNA}_{fsSCRNA}_{deRNA}_{omicMixMet}_{ppMixMet}_{fsMixMet}_{omicMET}_{ppMET}_{fsMET}_{deMET}_{li}_score.txt"   
    shell:"""
mkdir -p output/scores/
RCODE="
utils_script ='utils/data_processing.R';
prediction_file='{input.prediction}';  groundtruth_file='{input.groundtruth_file}';   
score_file='{output}'; source('{input.scoring_script}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""






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