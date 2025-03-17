
import yaml

# ## read csv parameters files to import blocks to compute
# DATASETS                 = {row[0].strip(): row[1].strip() for row in csv.reader(open("datasets.csv")) if not row[0].startswith("#")}
# REFERENCE                = ["data/reference_pdac.h5"]
# PRE_PROC                 = {row[0].strip(): row[1:] for row in csv.reader(open("pre-processing.csv")) if not row[0].startswith("#")} 
# FEATURES_SELECTION       = {row[0].strip(): row[1:] for row in csv.reader(open("feature_selection.csv")) if not row[0].startswith("#")} 
# EALRY_INTEGRATION        = {row[0].strip(): row[1:] for row in csv.reader(open("early_integration.csv")) if not row[0].startswith("#")} 
# INTERMEDIATE_INTEGRATION = {row[0].strip(): row[1:] for row in csv.reader(open("intermediate_integration.csv")) if not row[0].startswith("#")} 
# LATE_INTEGRATION         = {row[0].strip(): row[1:] for row in csv.reader(open("late_integration.csv")) if not row[0].startswith("#")} 
# DECOVOLUTION             = {row[0].strip(): row[1:] for row in csv.reader(open("decovolution.csv")) if not row[0].startswith("#")} 
# SPLIT                    = {row[0].strip(): row[1:] for row in csv.reader(open("split.csv")) if not row[0].startswith("#")} 


##import blocks and datasets from yml files
DATASETS                 = yaml.safe_load(open("datasets.yml")) 
REFERENCE                = ["data/ref.h5"]
PRE_PROC                 = yaml.safe_load(open("pre-processing.yml")) 
FEATURES_SELECTION       = yaml.safe_load(open("feature_selection.yml")) 
EALRY_INTEGRATION        = yaml.safe_load(open("early_integration.yml")) 
INTERMEDIATE_INTEGRATION = yaml.safe_load(open("intermediate_integration.yml")) 
LATE_INTEGRATION         = yaml.safe_load(open("late_integration.yml")) 
DECOVOLUTION             = yaml.safe_load(open("decovolution.yml")) 
SPLIT                    = yaml.safe_load(open("split.yml")) 







#Create the pipeline block combinaison.
rule all: 
    input: 
        expand("data/{dataset}.h5",dataset= DATASETS.keys()),
        expand("output/pre-processing/{dataset}_{pp}.h5",dataset= DATASETS.keys(),pp =PRE_PROC.keys()),
        expand("output/feature_selection/{dataset}_{pp}_{fs}.h5",dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys()),
        
        ##Pipeline A => split and decovo
        expand("output/split_decovolution/{dataset}_{pp}_{fs}_{split}_rna-{de}.h5",
            dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),split= SPLIT.keys(),  de=DECOVOLUTION.keys()),
        expand("output/split_decovolution/{dataset}_{pp}_{fs}_{split}_met-{de}.h5",
            dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),split= SPLIT.keys(),  de=DECOVOLUTION.keys()),
        expand("output/prediction/{dataset}_{pp}_{fs}_{split}_rna-{de1}_met-{de2}_{li}.h5",
            dataset= DATASETS.keys(),pp =PRE_PROC.keys(),fs = FEATURES_SELECTION.keys(),split= SPLIT.keys(),  de1=DECOVOLUTION.keys(),de2=DECOVOLUTION.keys(),li=LATE_INTEGRATION.keys())
        
        ##Pipeline B  => early integrategration and decovo
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


rule pre_processing:
    threads: 1
    message: "-- Processing pre processing Block -- "
    input: 
        mix = "data/{dataset}.h5" ,
        reference = REFERENCE[0]
    output: 
        "output/pre-processing/{dataset}_{pp}.h5"
    params:
      script = lambda wildcah5: PRE_PROC[wildcah5.pp]['path'].strip()  # get_script
    # log: file = "logs/05_metaanalysis.Rout"
    shell:"""
mkdir -p output/pre-processing/
RCODE="mixes_file='{input.mix}'; reference_file='{input.reference}';   output_file='{output}'; script_file='{params.script}';  source('02_Pre_process.R');"
echo $RCODE | Rscript -
"""


rule features_selection:
    threads: 1
    message: "-- Processing features selections Block -- "
    input: 
        "output/pre-processing/{dataset}_{pp}.h5"
    output: 
        "output/feature_selection/{dataset}_{pp}_{fs}.h5"
    params:
      script = lambda wildcah5: FEATURES_SELECTION[wildcah5.fs]['path'].strip()  # get_script

    # log: file = "logs/05_metaanalysis.Rout"
    shell:"""
mkdir -p output/feature_selection/
RCODE="input_file='{input}';   output_file='{output}'; script_file='{params.script}';  source('middle_man.R');"
echo $RCODE | Rscript -
"""


### Pipeline A####

rule prediction_split_decovolution_rna:
    threads: 1
    message: "-- Processing splitted rna decovolution Block, Pipeline A -- "
    input: 
        "output/feature_selection/{dataset}_{pp}_{fs}.h5"
    output: 
        "output/split_decovolution/{dataset}_{pp}_{fs}_{split}_rna-{de}.h5"
    # log: file = "logs/05_metaanalysis.Rout"
    params:
      script_de = lambda wildcah5: DECOVOLUTION[wildcah5.de]['path'].strip(), 
    #   script_de2 = lambda wildcah5: DECOVOLUTION[wildcah5.de2] ,
      script_split = lambda wildcah5: SPLIT[wildcah5.split]['path'].strip()
    shell:"""
mkdir -p output/split_decovolution/
RCODE="input_file='{input}';   output_file='{output}'; script_split='{params.script_split}'; 
script_de_rna='{params.script_de}' ;  source('pipeline_A_split.R');"
echo $RCODE | Rscript -
"""

rule prediction_split_decovolution_met:
    threads: 1
    message: "-- Processing splitted met decovolution Block, Pipeline A -- "
    input: 
        "output/feature_selection/{dataset}_{pp}_{fs}.h5"
    output: 
        "output/split_decovolution/{dataset}_{pp}_{fs}_{split}_met-{de}.h5"
    # log: file = "logs/05_metaanalysis.Rout"
    params:
      script_de = lambda wildcah5: DECOVOLUTION[wildcah5.de]['path'].strip(), 
    #   script_de2 = lambda wildcah5: DECOVOLUTION[wildcah5.de2] ,
      script_split = lambda wildcah5: SPLIT[wildcah5.split]['path'].strip() 
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
      script_li = lambda wildcah5: LATE_INTEGRATION[wildcah5.li]['path'].strip(), 
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