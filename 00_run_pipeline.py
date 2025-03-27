
import yaml

##import blocks and datasets from yml files
DATASETS                 = yaml.safe_load(open("datasets.yml")) 
REFERENCE                = ["data/ref.h5"]
PRE_PROC                 = yaml.safe_load(open("preprocessing.yml")) 
FEATURES_SELECTION       = yaml.safe_load(open("feature_selection.yml")) 
EALRY_INTEGRATION        = yaml.safe_load(open("early_integration.yml")) 
INTERMEDIATE_INTEGRATION = yaml.safe_load(open("intermediate_integration.yml")) 
LATE_INTEGRATION         = yaml.safe_load(open("late_integration.yml")) 
DECOVOLUTION             = yaml.safe_load(open("decovolution.yml")) 
SPLIT                    = yaml.safe_load(open("split.yml")) 


def compare_input_output(dic_out,dic_in):
    return (set(dic_out['output']) == set(dic_in['input']) or 'ANY' in dic_in['input'] or ('ANY' in dic_out['output']))


def get_blockv(file, dic_block):
    last_fun = file.split('_')[-1]
    return(dic_block[last_fun])

def add_h5(list_files):
    return (list(map( (lambda x : x+'.h5'), list_files)))


#Create combinaison for pipeline A. 
datasets_files = [f'{dsv['path']}' for dsv in DATASETS.values() ]
pp_files = [f'output/preprocessing/{dataset}_{pp}' for dataset in DATASETS.keys() for pp in PRE_PROC.keys()  ]
fs_files = [f'output/feature_selection/{last_file.split('/')[-1]}_{fs}' for  last_file in pp_files  
            for fs,fsv in FEATURES_SELECTION.items() if compare_input_output(get_blockv(last_file,PRE_PROC),fsv) ]

de_files_rna = [f'output/split_decovolution/{last_file.split('/')[-1]}_{split}_rna-{de}' for last_file in fs_files 
               for split in SPLIT.keys() for de in DECOVOLUTION.keys() 
               if ('RNA' in get_blockv(last_file, FEATURES_SELECTION)["output"] or 'ANY' in get_blockv(last_file, FEATURES_SELECTION)["output"] ) ]
de_files_met = [f'output/split_decovolution/{last_file.split('/')[-1]}_{split}_met-{de}' for last_file in fs_files 
               for split in SPLIT.keys() for de in DECOVOLUTION.keys() 
               if ('MET' in get_blockv(last_file, FEATURES_SELECTION)["output"] or 'ANY' in get_blockv(last_file, FEATURES_SELECTION)["output"])]
# Decovolution tools accept both met and RNA so we do not need to 
li_files = [f'output/prediction/{last_file.split('/')[-1]}_met-{de}_{li}' for last_file  in de_files_rna 
            for de in DECOVOLUTION.keys() for li in LATE_INTEGRATION.keys() ]

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
        add_h5(li_files)
        

        # expand("data/{dataset}.h5",dataset= DATASETS.keys()),
        # expand("output/preprocessing/{dataset}_{pp}.h5",dataset= DATASETS.keys(),pp =PRE_PROC.keys()),
        # expand("output/feature_selection/{dataset}_{pp}_{fs}.h5",dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys()),
        
        # ##Pipeline A => split and decovo
        # expand("output/split_decovolution/{dataset}_{pp}_{fs}_{split}_rna-{de}.h5",
        #     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),split= SPLIT.keys(),  de=DECOVOLUTION.keys()),
        # expand("output/split_decovolution/{dataset}_{pp}_{fs}_{split}_met-{de}.h5",
        #     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),split= SPLIT.keys(),  de=DECOVOLUTION.keys()),
        # expand("output/prediction/{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}.h5",
        #     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),split= SPLIT.keys(),  de1=DECOVOLUTION.keys(),de2=DECOVOLUTION.keys(),li=LATE_INTEGRATION.keys())
        
        ##Pipeline B  => early integration and decovo
        # expand("output/early_integration/{dataset}_{pp}_{fs}_{ei}.h5",
        #     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),ei=EALRY_INTEGRATION.keys()),
        # expand("output/prediction/{dataset}_{pp}_{fs}_{ei}_{de}.h5",
        #     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),ei=EALRY_INTEGRATION.keys(), de=DECOVOLUTION.keys()),
        
        # ##Pipeline C  => intermediate decovo
        # expand("output/prediction/{dataset}_{pp}_{fs}_{it}.h5",
        #     dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),it=INTERMEDIATE_INTEGRATION.keys()),


# rule generate_data:
#     threads: 1
#     message: "-- generation of data -- "
#     input: 
#         "data/{dataset}.h5" 
#     output:
#         "data/{dataset}.h5" 
#     shell:"""
# python 01_data_TODO.py {input} {output}
# """


rule preprocessing:
    threads: 1
    message: "-- Processing pre processing Block -- "
    input: 
        mix = "data/mixes1_{dataset}_pdac.h5" ,
        reference = REFERENCE[0]
    output: 
        "output/preprocessing/{dataset}_{pp}.h5"
    params:
      script = lambda wildcard: PRE_PROC[wildcard.pp]['path'].strip()  # get_script
    # log: file = "logs/05_metaanalysis.Rout"
    shell:"""
mkdir -p output/preprocessing/
RCODE="mixes_file='{input.mix}'; reference_file='{input.reference}';   output_file='{output}'; script_file='{params.script}';  source('02_preprocess.R');"
echo $RCODE | Rscript -
"""


rule features_selection:
    threads: 1
    message: "-- Processing features selections Block -- "
    input: 
        "output/preprocessing/{dataset}_{pp}.h5"
    output: 
        "output/feature_selection/{dataset}_{pp}_{fs}.h5"
    params:
      script = lambda wildcard: FEATURES_SELECTION[wildcard.fs]['path'].strip()  # get_script

    # log: file = "logs/05_metaanalysis.Rout"
    shell:"""
mkdir -p output/feature_selection/
RCODE="input_file='{input}';   output_file='{output}'; script_file='{params.script}';  source('03_features_selection.R');"
echo $RCODE | Rscript -
"""


### Pipeline A####

rule prediction_decovolution_rna:
    threads: 1
    message: "-- Processing splitted rna decovolution Block, Pipeline A -- "
    input: 
        "output/feature_selection/{dataset}_{pp}_{fs}.h5"
    output: 
        "output/split_decovolution/{dataset}_{pp}_{fs}_{split}_rna-{de}.h5"
    # log: file = "logs/05_metaanalysis.Rout"
    params:
      script_de = lambda wildcard: DECOVOLUTION[wildcard.de]['path'].strip(), 
    #   script_de2 = lambda wildcard: DECOVOLUTION[wildcard.de2] ,
      script_split = lambda wildcard: SPLIT[wildcard.split]['path'].strip()
    shell:"""
mkdir -p output/split_decovolution/
RCODE="input_file='{input}';   output_file='{output}'; script_split='{params.script_split}'; 
script_de_rna='{params.script_de}' ;  source('pipeline_A_split.R');"
echo $RCODE | Rscript -
"""

rule prediction_decovolution_met:
    threads: 1
    message: "-- Processing splitted met decovolution Block, Pipeline A -- "
    input: 
        "output/feature_selection/{dataset}_{pp}_{fs}.h5"
    output: 
        "output/split_decovolution/{dataset}_{pp}_{fs}_{split}_met-{de}.h5"
    # log: file = "logs/05_metaanalysis.Rout"
    params:
      script_de = lambda wildcard: DECOVOLUTION[wildcard.de]['path'].strip(), 
    #   script_de2 = lambda wildcard: DECOVOLUTION[wildcard.de2] ,
      script_split = lambda wildcard: SPLIT[wildcard.split]['path'].strip() 
    shell:"""
mkdir -p output/split_decovolution/
RCODE="input_file='{input}';   output_file='{output}'; 
script_split='{params.script_split}'; script_de_met='{params.script_de}';  
source('pipeline_A_split.R');"
echo $RCODE | Rscript -
"""


rule late_integration:
    threads: 1
    message: "-- Processing splitted decovolution late ingration Block, Pipeline A -- "
    input: 
        input_file_rna= "output/split_decovolution/{dataset}_{pp}_{fs}_{split}_rna-{de1}.h5",
        input_file_met = "output/split_decovolution/{dataset}_{pp}_{fs}_{split}_met-{de2}.h5"
    output: 
        "output/prediction/{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}.h5"
    # log: file = "logs/05_metaanalysis.Rout"
    params:
      script_li = lambda wildcard: LATE_INTEGRATION[wildcard.li]['path'].strip(), 
    shell:"""
mkdir -p output/prediction/
RCODE="input_file_rna='{input.input_file_rna}';  input_file_met='{input.input_file_met}';   output_file='{output}'; script_file='{params.script_li}'; source('pipeline_A_merge.R');

 "
echo $RCODE | Rscript -
"""





### Pipeline B####


rule early_integration:
    threads: 1
    message: "-- Processing early integration  Block, Pipeline B -- "
    input: 
        "output/feature_selection/{dataset}_{pp}_{fs}.h5"
    output: 
        "output/early_integration/{dataset}_{pp}_{fs}_{ei}.h5"
    # log: file = "logs/05_metaanalysis.Rout"
    shell:"""
mkdir -p output/early_integration/
echo  {input} {output}
touch {output}
"""

rule prediction_with_early_integration:
    threads: 1
    message: "-- Processing decovolution with early integration Block, Pipeline B -- "
    input: 
        "output/early_integration/{dataset}_{pp}_{fs}_{ei}.h5"
    output: 
        "output/prediction/{dataset}_{pp}_{fs}_{ei}_{de}.h5"
    # log: file = "logs/05_metaanalysis.Rout"
    shell:"""
mkdir -p output/prediction/
echo  {input} {output}
touch {output}
"""



### Pipeline C ####

rule intermediate_integration:
    threads: 1
    message: "-- Processing itermediate integration  Block, Pipeline C -- "
    input: 
        "output/feature_selection/{dataset}_{pp}_{fs}.gz"
    output: 
        "output/early_integration/{dataset}_{pp}_{fs}_{ei}.gz"
    # log: file = "logs/05_metaanalysis.Rout"
    shell:"""
mkdir -p output/early_integration/
echo  {input} {output}
touch {output}
"""



rule clean:
    threads: 1
    shell:"""
rm -rf output/
"""
# rm -rf tmp_all/
# rm -rf 01_*
# rm -rf 02_prediction/
# rm -rf 03_scores/
# rm -rf 04_visu/
# """

rule gantt:
    threads: 1
    shell:"""
smgantt
"""