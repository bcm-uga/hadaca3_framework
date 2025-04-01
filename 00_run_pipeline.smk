
import yaml

##import blocks and datasets from yml files
DATASETS                 = yaml.safe_load(open("datasets.yml")) 
REFERENCE                = ["data/ref.h5"]
PRE_PROC                 = yaml.safe_load(open("preprocessing.yml")) 
FEATURES_SELECTION       = yaml.safe_load(open("feature_selection.yml")) 
EARLY_INTEGRATION        = yaml.safe_load(open("early_integration.yml")) 
INTERMEDIATE_INTEGRATION = yaml.safe_load(open("intermediate_integration.yml")) 
LATE_INTEGRATION         = yaml.safe_load(open("late_integration.yml")) 
DECONVOLUTION            = yaml.safe_load(open("deconvolution.yml")) 
SPLIT                    = yaml.safe_load(open("split.yml")) 


def compare_input_output(dic_out,dic_in):
    return (set(dic_out['output']) == set(dic_in['input']) or 'ANY' in dic_in['input'] or ('ANY' in dic_out['output']))


def get_blockv(file, dic_block):
    last_fun = file.split('_')[-1]
    return(dic_block[last_fun])

def add_h5(list_files):
    return (list(map( (lambda x : x+'.h5'), list_files)))


datasets_files = [f'{dsv['path']}' for dsv in DATASETS.values() ]
pp_files = [f'output/preprocessing/{dataset}_{pp}' for dataset in DATASETS.keys() for pp in PRE_PROC.keys()  ]
fs_files = [f'output/feature_selection/{last_file.split('/')[-1]}_{fs}' for  last_file in pp_files  
            for fs,fsv in FEATURES_SELECTION.items() if compare_input_output(get_blockv(last_file,PRE_PROC),fsv) ]

#Create combinaison for pipeline A. 
de_files_rna = [f'output/split_deconvolution/{last_file.split('/')[-1]}_{split}_rna-{de}' for last_file in fs_files 
               for split in SPLIT.keys() for de in DECONVOLUTION.keys() 
               if ('RNA' in get_blockv(last_file, FEATURES_SELECTION)["output"] or 'ANY' in get_blockv(last_file, FEATURES_SELECTION)["output"] ) ]
de_files_met = [f'output/split_deconvolution/{last_file.split('/')[-1]}_{split}_met-{de}' for last_file in fs_files 
               for split in SPLIT.keys() for de in DECONVOLUTION.keys() 
               if ('MET' in get_blockv(last_file, FEATURES_SELECTION)["output"] or 'ANY' in get_blockv(last_file, FEATURES_SELECTION)["output"])]
li_files = [f'output/prediction/{last_file.split('/')[-1]}_met-{de}_{li}' for last_file  in de_files_rna 
            for de in DECONVOLUTION.keys() for li in LATE_INTEGRATION.keys() ]
scores_files = [f'output/scores/{last_file.split('/')[-1]}_score' for last_file  in li_files ]

#Pipeline B




#Create the pipeline block combinaison.
rule all: 
    input: 
        #pipelines A
        datasets_files,
        add_h5(pp_files),
        add_h5(fs_files),
        add_h5(de_files_rna),
        add_h5(de_files_met),
        add_h5(li_files),
        score =  add_h5(scores_files),
        script_file ='06_metaanalysis.Rmd'
    output: 
        "06_metaanalysis.html"
    log : 
        "logs/06_metaanalysis.Rout"
    shell:"""
RCODE="score_files = strsplit(trimws('{input.score}'),' ') ; 
rmarkdown::render('{input.script_file}');"
echo $RCODE | Rscript - 2>&1 > {log}
echo "all is done!" 
"""    
        
        ##Pipeline B  => early integration and decovo
        # expand("output/early_integration/{dataset}_{pp}_{fs}_{ei}.h5",
        #     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),ei=EARLY_INTEGRATION.keys()),
        # expand("output/prediction/{dataset}_{pp}_{fs}_{ei}_{de}.h5",
        #     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),ei=EARLY_INTEGRATION.keys(), de=DECOnVOLUTION.keys()),
        
        # ##Pipeline C  => intermediate decovo
        # expand("output/prediction/{dataset}_{pp}_{fs}_{it}.h5",
        #     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),it=INTERMEDIATE_INTEGRATION.keys()),


rule preprocessing:
    threads: 1
    message: "-- Processing pre processing Block -- "
    input: 
        pp_wrapper="02_preprocess.R",
        script = lambda wildcard: PRE_PROC[wildcard.pp]['path'].strip(),
        mix = lambda wildcard: DATASETS[wildcard.dataset]['path'].strip(), # "data/mixes1_{dataset}_pdac.h5" ,
        reference = REFERENCE[0]
    output: 
        "output/preprocessing/{dataset}_{pp}.h5"
    log : 
        "logs/02_{dataset}_{pp}.h5"        
    shell:"""
mkdir -p output/preprocessing/
RCODE="mixes_file='{input.mix}'; reference_file='{input.reference}';   output_file='{output}'; script_file='{input.script}';  source('{input.pp_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""


rule features_selection:
    threads: 1
    message: "-- Processing features selections Block -- "
    input: 
        fs_wrapper= '03_features_selection.R',
        file_input= "output/preprocessing/{dataset}_{pp}.h5" ,
        script = lambda wildcard: FEATURES_SELECTION[wildcard.fs]['path'].strip() 
    output: 
        "output/feature_selection/{dataset}_{pp}_{fs}.h5"
    log : 
        "logs/03_{dataset}_{pp}_{fs}.h5" 
    shell:"""
mkdir -p output/feature_selection/
RCODE="input_file='{input.file_input}';   output_file='{output}'; script_file='{input.script}';  source('{input.fs_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""


### Pipeline A####

rule prediction_deconvolution_rna:
    threads: 1
    message: "-- Processing splitted rna deconvolution Block, Pipeline A -- "
    input: 
       split_wrapper = "04_pipeline_A_split.R" ,
       script_de = lambda wildcard: DECONVOLUTION[wildcard.de]['path'].strip(), 
       script_split = lambda wildcard: SPLIT[wildcard.split]['path'].strip(),
       file_input= "output/feature_selection/{dataset}_{pp}_{fs}.h5"
    output: 
        "output/split_deconvolution/{dataset}_{pp}_{fs}_{split}_rna-{de}.h5"
    log : 
        "logs/04_{dataset}_{pp}_{fs}_{split}_rna-{de}.h5"     
    shell:"""
mkdir -p output/split_deconvolution/
RCODE="input_file='{input.file_input}';   output_file='{output}'; script_split='{input.script_split}'; 
script_de_rna='{input.script_de}' ;  source('{input.split_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""

rule prediction_deconvolution_met:
    threads: 1
    message: "-- Processing splitted met deconvolution Block, Pipeline A -- "
    input: 
        split_wrapper = "04_pipeline_A_split.R" , 
        script_de = lambda wildcard: DECONVOLUTION[wildcard.de]['path'].strip(),
        script_split = lambda wildcard: SPLIT[wildcard.split]['path'].strip(), 
        file_input= "output/feature_selection/{dataset}_{pp}_{fs}.h5"
    output: 
        "output/split_deconvolution/{dataset}_{pp}_{fs}_{split}_met-{de}.h5"
    log : 
        "logs/04_{dataset}_{pp}_{fs}_{split}_met-{de}.h5"      
    shell:"""
mkdir -p output/split_deconvolution/
RCODE="input_file='{input.file_input}';   output_file='{output}'; 
script_split='{input.script_split}'; script_de_met='{input.script_de}';  
source('{input.split_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""


rule late_integration:
    threads: 1
    message: "-- Processing splitted deconvolution late ingration Block, Pipeline A -- "
    input: 
        merge_wrapper = '04_pipeline_A_merge.R',
        script_li = lambda wildcard: LATE_INTEGRATION[wildcard.li]['path'].strip(), 
        input_file_rna= "output/split_deconvolution/{dataset}_{pp}_{fs}_{split}_rna-{de1}.h5",
        input_file_met = "output/split_deconvolution/{dataset}_{pp}_{fs}_{split}_met-{de2}.h5"
    output: 
        "output/prediction/{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}.h5"      
    log : 
        "logs/04_{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}.h5"                
    shell:"""
mkdir -p output/prediction/
RCODE="input_file_rna='{input.input_file_rna}';  input_file_met='{input.input_file_met}';   
output_file='{output}'; script_file='{input.script_li}'; source('{input.merge_wrapper}');"
echo $RCODE | Rscript - 2>&1 > {log}
"""



rule scoring:
    threads: 1
    message: "-- Score prediction  -- "
    input: 
        scoring_script = '05_scoring.R',
        prediction = "output/prediction/{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}.h5",
        groundtruth_file = lambda wildcard: DATASETS[wildcard.dataset]['groundtruth_file_path'].strip(), 
    output: 
        "output/scores/{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}_score.h5"
    log : 
        "logs/05_{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}_score.h5"   
    shell:"""
mkdir -p output/scores/
RCODE="prediction_file='{input.prediction}';  groundtruth_file='{input.groundtruth_file}';   
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
rm 06_metaanalysis.html
"""

rule gantt:
    threads: 1
    shell:"""
smgantt
"""