#!/usr/bin/env nextflow

import groovy.yaml.YamlSlurper

nextflow.enable.dsl=2

workDir='~/project/hadaca3_framework'

params.config_files = [
    datasets: 'datasets.yml',
    pre_proc: 'preprocessing.yml',
    features_selection: 'feature_selection.yml',
    early_integration: 'early_integration.yml',
    intermediate_integration: 'intermediate_integration.yml',
    late_integration: 'late_integration.yml',
    deconvolution: 'deconvolution.yml'
]

params.reference = ['data/ref.h5']
params.cleaner = 'global_cleaning/clean_matrix.R'
params.mixomics = ['mixRNA', 'mixMET']
params.refomics = ['MET', 'RNA', 'scRNA']
params.omic_dirs = params.mixomics + params.refomics


params.wrapper = [
    script_01: '01_global_preprocess_mix.R',
    script_02 : '01_global_preprocess_ref.R'
]



def get_omic(path) {
    return path.split('/')[-2]
}

def get_dataset(path) {
    return path.split('/')[-1].split('_')[0]
}

def block_combinaison(path, add_omics = true) {
    def l_path = path.split('/')
    def omic = l_path[-2]
    def descriptif = l_path[-1].split('_')
    return add_omics ? [omic] + descriptif[1..-1] : descriptif[1..-1]
}

workflow { 
     def CONFIG = [:]

    params.config_files.each { key, filePath ->
        def parsedContent = new YamlSlurper().parse(filePath as File)
        CONFIG[key] = parsedContent
    }


    original_datasets_files = CONFIG['datasets'].collect { k, v -> v['path'] }
    cleaned_datasets_files = CONFIG['datasets'].keySet().collect { "output/mixes/${it}" }
    cleaned_REFERENCE = params.reference.collect { "output/ref/${it.split('/')[-1]}" }

    pp_files = []
    CONFIG['pre_proc'].each { pp, ppv ->
        params.mixomics.each { omic ->
            if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')) {
                CONFIG['datasets'].keySet().each { dataset ->
                    pp_files.add("output/preprocessing/${omic}/${dataset}_${pp}")
                }
            }
        }
        params.refomics.each { omic ->
            if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')) {
                pp_files.add("output/preprocessing/${omic}/ref_${pp}")
            }
        }
    }

    fs_files = []
    pp_files.each { last_pp ->
        omic = last_pp.split('/')[-2]
        CONFIG['features_selection'].each { fs, fsv ->
            if (fsv['omic'].contains(omic) || fsv['omic'].contains('ANY')) {
                fs_files.add("output/feature-selection/${omic}/${last_pp.split('/')[-1]}_${fs}")
            }
        }
    }



    def de_rna_unit_files = []
    fs_files.each { mix_combi ->
        fs_files.each { rna_combi ->
            fs_files.each { sc_combi ->
                if (get_omic(mix_combi) == 'mixRNA' && get_omic(rna_combi) == 'RNA' && get_omic(sc_combi) == 'scRNA') {
                    CONFIG['deconvolution'].each { de, dev ->
                        de_rna_unit_files.add("output/rna-decovolution-split/${get_dataset(mix_combi)}_${block_combinaison(mix_combi).join('_')}_${block_combinaison(rna_combi).join('_')}_${block_combinaison(sc_combi).join('_')}_${de}")
                    }
                }
            }
        }
    }

    // Generate de_met_unit_files
    def de_met_unit_files = []
    fs_files.each { mix_combi ->
        fs_files.each { met_combi ->
            if (get_omic(mix_combi) == 'mixMET' && get_omic(met_combi) == 'MET') {
                CONFIG['deconvolution'].each { de, dev ->
                    de_met_unit_files.add("output/met-decovolution-split/${get_dataset(mix_combi)}_${block_combinaison(mix_combi).join('_')}_${block_combinaison(met_combi).join('_')}_${de}")
                }
            }
        }
    }

    
    li_files = []
    de_rna_unit_files.each { file_prediction_RNA ->
        de_met_unit_files.each { file_prediction_MET ->
            if (get_dataset(file_prediction_RNA) == get_dataset(file_prediction_MET)) {
                CONFIG['late_integration'].keySet().each { li ->
                    li_files.add("output/prediction/${get_dataset(file_prediction_RNA)}_${block_combinaison(file_prediction_RNA, false).join('_')}_${block_combinaison(file_prediction_MET, false).join('_')}_${li}")
                }
            }
        }
    }

    scores_files = li_files.collect { "output/scores/${it.split('/')[-1]}_score" }


    // println "scores: ${li_files.sort()}"
    // println "scores: ${scores_files.size()}"


     cleaning_mix( Channel.fromPath(original_datasets_files),Channel.fromPath(params.cleaner) ,
     Channel.fromPath(params.wrapper['script_01'])  )
    // cleaning_ref(params.reference)
    // preprocessing(pp_files)
//     features_selection(fs_files)
//     prediction_deconvolution_rna(de_rna_unit_files)
//     prediction_deconvolution_met(de_met_unit_files)
//     late_integration(li_files)
//     scoring(scores_files)
}

process cleaning_mix {
    input:
    path mixes
    path cleaner
    path wrapper01
     

    output:
    path "output/mixes/${mixes.baseName}.h5"

    script:
    """
    mkdir -p output/mixes/
    RCODE="mixes_file='${mixes}'; output_file='output/mixes/${mixes.baseName}.h5'; 
    cleaner='${cleaner}'; source('${wrapper01}');"
    echo \$RCODE | Rscript -
    """

    workingDir : '.'
}

// process cleaning_mix {
//     input:
//     path mixes

//     output:
//     path "output/mixes/${mixes.baseName}.h5"

//     script:
//     """
//     mkdir -p output/mixes/
//     RCODE="mixes_file='${mixes}'; output_file='${output}'; cleaner='${params.cleaner}'; source('01_global_preprocess_mix.R');"
//     echo \$RCODE | Rscript - 2>&1 > logs/01_${mixes.baseName}.txt
//     """
// }

process cleaning_ref{
    input:
    path reference

    output:
    path "output/ref/${reference.baseName}.h5"

    script:
    """
    mkdir -p output/ref/
    RCODE="reference_file='${reference}'; output_file='${output}'; cleaner='${params.cleaner}'; source('01_global_preprocess_ref.R');"
    echo \$RCODE | Rscript - 2>&1 > logs/01_${reference.baseName}.txt
    """
}

process preprocessing {
    input:
    path mix, reference
    path script

    output:
    path "output/preprocessing/${omic}/${dataset}_${pp}.h5"

    script:
    """
    mkdir -p output/preprocessing/${params.omic_dirs.join(',')}/
    RCODE="mixes_file='${mix}'; reference_file='${reference}'; output_file='${output}'; script_file='${script}'; source('02_preprocess.R');"
    echo \$RCODE | Rscript - 2>&1 > logs/02_${omic}_${dataset}_${pp}.txt
    """
}

process features_selection {
    input:
    path file_input
    path script

    output:
    path "output/feature-selection/${omic}/${dataset}_${pp}_${fs}.h5"

    script:
    """
    mkdir -p output/feature-selection/${params.omic_dirs.join(',')}/
    RCODE="input_file='${file_input}'; output_file='${output}'; script_file='${script}'; source('03_features_selection.R');"
    echo \$RCODE | Rscript - 2>&1 > logs/03_${omic}_${dataset}_${pp}_${fs}.txt
    """
}

process prediction_deconvolution_rna {
    input:
    path file_input_mix, file_input_rna, file_input_scrna
    path script_de

    output:
    path "output/rna-decovolution-split/${dataset}_${omicmix}_${ppmix}_${fsmix}_${omicrna}_${pprna}_${fsrna}_${omicscrna}_${ppsc}_${fssc}_${de}.h5"

    script:
    """
    mkdir -p output/rna-decovolution-split/
    RCODE="input_file_mix='${file_input_mix}'; input_file_rna='${file_input_rna}'; input_file_sc='${file_input_scrna}'; output_file='${output}'; script_de_rna='${script_de}'; source('04_decovolution_RNA_unit_pipA.R');"
    echo \$RCODE | Rscript - 2>&1 > logs/04_${dataset}_${omicmix}_${ppmix}_${fsmix}_${omicrna}_${pprna}_${fsrna}_${omicscrna}_${ppsc}_${fssc}_${de}.txt
    """
}

process prediction_deconvolution_met {
    input:
    path file_input_mix, file_input_met
    path script_de

    output:
    path "output/met-decovolution-split/${dataset}_${omicmix}_${ppmix}_${fsmix}_${omicmet}_${ppmet}_${fsmet}_${de}.h5"

    script:
    """
    mkdir -p output/met-decovolution-split/
    RCODE="input_file_mix='${file_input_mix}'; input_file_met='${file_input_met}'; output_file='${output}'; script_de_met='${script_de}'; source('04_decovolution_MET_unit_pipA.R');"
    echo \$RCODE | Rscript - 2>&1 > logs/04_${dataset}_${omicmix}_${ppmix}_${fsmix}_${omicmet}_${ppmet}_${fsmet}_${de}.txt
    """
}

process late_integration {
    input:
    path input_file_rna, input_file_met
    path script_li

    output:
    path "output/prediction/${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}.h5"

    script:
    """
    mkdir -p output/prediction/
    RCODE="input_file_rna='${input_file_rna}'; input_file_met='${input_file_met}'; output_file='${output}'; script_file='${script_li}'; source('05_late_integration_A.R');"
    echo \$RCODE | Rscript - 2>&1 > logs/05_${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}.txt
    """
}

process scoring {
    input:
    path prediction
    path groundtruth_file
    path scoring_script

    output:
    path "output/scores/${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}_score.h5"

    script:
    """
    mkdir -p output/scores/
    RCODE="prediction_file='${prediction}'; groundtruth_file='${groundtruth_file}'; score_file='${output}'; source('${scoring_script}');"
    echo \$RCODE | Rscript - 2>&1 > logs/06_${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}_score.txt
    """
}



// / #!/usr/bin/env nextflow

// nextflow.enable.dsl=2

// params.config_files = [
//     datasets: 'datasets.yml',
//     pre_proc: 'preprocessing.yml',
//     features_selection: 'feature_selection.yml',
//     early_integration: 'early_integration.yml',
//     intermediate_integration: 'intermediate_integration.yml',
//     late_integration: 'late_integration.yml',
//     deconvolution: 'deconvolution.yml'
// ]

// params.reference = ['data/ref.h5']
// params.cleaner = 'global_cleaning/clean_matrix.R'
// params.mixomics = ['mixRNA', 'mixMET']
// params.refomics = ['MET', 'RNA', 'scRNA']
// params.omic_dirs = params.mixomics + params.refomics

// process loadConfig {
//     input:
//     path config_file

//     output:
//     path "config.yml"

//     script:
//     """
//     cp ${config_file} config.yml
//     """
// }

// workflow {
//     Channel.fromPath(params.config_files.values()) | loadConfig

//     CONFIG = Channel.fromPath('config.yml').map { yaml ->
//         return groovy.yaml.YamlSlurper().parseText(yaml.text)
//     }

//     original_datasets_files = CONFIG['datasets'].collect { k, v -> v['path'] }
//     cleaned_datasets_files = CONFIG['datasets'].keySet().collect { "output/mixes/${it}" }
//     cleaned_REFERENCE = params.reference.collect { "output/ref/${it.split('/')[-1]}" }

//     pp_files = []
//     CONFIG['pre_proc'].each { pp, ppv ->
//         params.mixomics.each { omic ->
//             if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')) {
//                 CONFIG['datasets'].keySet().each { dataset ->
//                     pp_files.add("output/preprocessing/${omic}/${dataset}_${pp}")
//                 }
//             }
//         }
//         params.refomics.each { omic ->
//             if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')) {
//                 pp_files.add("output/preprocessing/${omic}/ref_${pp}")
//             }
//         }
//     }

//     fs_files = []
//     pp_files.each { last_pp ->
//         omic = last_pp.split('/')[-2]
//         CONFIG['features_selection'].each { fs, fsv ->
//             if (fsv['omic'].contains(omic) || fsv['omic'].contains('ANY')) {
//                 fs_files.add("output/feature-selection/${omic}/${last_pp.split('/')[-1]}_${fs}")
//             }
//         }
//     }

//     de_rna_unit_files = []
//     fs_files.combinations(3).each { mix_combi, rna_combi, sc_combi ->
//         if (get_omic(mix_combi) == 'mixRNA' && get_omic(rna_combi) == 'RNA' && get_omic(sc_combi) == 'scRNA') {
//             CONFIG['deconvolution'].each { de, dev ->
//                 de_rna_unit_files.add("output/rna-decovolution-split/${get_dataset(mix_combi)}_${block_combinaison(mix_combi)}_${block_combinaison(rna_combi)}_${block_combinaison(sc_combi)}_${de}")
//             }
//         }
//     }

//     de_met_unit_files = []
//     fs_files.combinations(2).each { mix_combi, met_combi ->
//         if (get_omic(mix_combi) == 'mixMET' && get_omic(met_combi) == 'MET') {
//             CONFIG['deconvolution'].each { de, dev ->
//                 de_met_unit_files.add("output/met-decovolution-split/${get_dataset(mix_combi)}_${block_combinaison(mix_combi)}_${block_combinaison(met_combi)}_${de}")
//             }
//         }
//     }

//     li_files = []
//     de_rna_unit_files.each { file_prediction_RNA ->
//         de_met_unit_files.each { file_prediction_MET ->
//             if (get_dataset(file_prediction_RNA) == get_dataset(file_prediction_MET)) {
//                 CONFIG['late_integration'].keySet().each { li ->
//                     li_files.add("output/prediction/${get_dataset(file_prediction_RNA)}_${block_combinaison(file_prediction_RNA, false)}_${block_combinaison(file_prediction_MET, false)}_${li}")
//                 }
//             }
//         }
//     }

//     scores_files = li_files.collect { "output/scores/${it.split('/')[-1]}_score" }

//     cleaning_mix(original_datasets_files)
//     cleaning_ref(params.reference)
//     preprocessing(pp_files)
//     features_selection(fs_files)
//     prediction_deconvolution_rna(de_rna_unit_files)
//     prediction_deconvolution_met(de_met_unit_files)
//     late_integration(li_files)
//     scoring(scores_files)
// }

// def get_omic(path) {
//     return path.split('/')[-2]
// }

// def get_dataset(path) {
//     return path.split('/')[-1].split('_')[0]
// }

// def block_combinaison(path, add_omics = true) {
//     def descriptif = path.split('/')[-1].split('_')
//     def omic = path.split('/')[-2]
//     return add_omics ? [omic] + descriptif[1..-1] : descriptif[1..-1]
// }

// process cleaning_mix {
//     input:
//     path mixes

//     output:
//     path "output/mixes/${mixes.baseName}.h5"

//     script:
//     """
//     mkdir -p output/mixes/
//     RCODE="mixes_file='${mixes}'; output_file='${output}'; cleaner='${params.cleaner}'; source('01_global_preprocess_mix.R');"
//     echo \$RCODE | Rscript - 2>&1 > logs/01_${mixes.baseName}.txt
//     """
// }00_run_pipeline.nf

//     output:
//     path "output/ref/${reference.baseName}.h5"

//     script:
//     """
//     mkdir -p output/ref/
//     RCODE="reference_file='${reference}'; output_file='${output}'; cleaner='${params.cleaner}'; source('01_global_preprocess_ref.R');"
//     echo \$RCODE | Rscript - 2>&1 > logs/01_${reference.baseName}.txt
//     """
// }

// process preprocessing {
//     input:
//     path mix, reference
//     path script

//     output:
//     path "output/preprocessing/${omic}/${dataset}_${pp}.h5"

//     script:
//     """
//     mkdir -p output/preprocessing/${params.omic_dirs.join(',')}/
//     RCODE="mixes_file='${mix}'; reference_file='${reference}'; output_file='${output}'; script_file='${script}'; source('02_preprocess.R');"
//     echo \$RCODE | Rscript - 2>&1 > logs/02_${omic}_${dataset}_${pp}.txt
//     """
// }

// process features_selection {
//     input:
//     path file_input
//     path script

//     output:
//     path "output/feature-selection/${omic}/${dataset}_${pp}_${fs}.h5"

//     script:
//     """
//     mkdir -p output/feature-selection/${params.omic_dirs.join(',')}/
//     RCODE="input_file='${file_input}'; output_file='${output}'; script_file='${script}'; source('03_features_selection.R');"
//     echo \$RCODE | Rscript - 2>&1 > logs/03_${omic}_${dataset}_${pp}_${fs}.txt
//     """
// }

// process prediction_deconvolution_rna {
//     input:
//     path file_input_mix, file_input_rna, file_input_scrna
//     path script_de

//     output:
//     path "output/rna-decovolution-split/${dataset}_${omicmix}_${ppmix}_${fsmix}_${omicrna}_${pprna}_${fsrna}_${omicscrna}_${ppsc}_${fssc}_${de}.h5"

//     script:
//     """
//     mkdir -p output/rna-decovolution-split/
//     RCODE="input_file_mix='${file_input_mix}'; input_file_rna='${file_input_rna}'; input_file_sc='${file_input_scrna}'; output_file='${output}'; script_de_rna='${script_de}'; source('04_decovolution_RNA_unit_pipA.R');"
//     echo \$RCODE | Rscript - 2>&1 > logs/04_${dataset}_${omicmix}_${ppmix}_${fsmix}_${omicrna}_${pprna}_${fsrna}_${omicscrna}_${ppsc}_${fssc}_${de}.txt
//     """
// }

// process prediction_deconvolution_met {
//     input:
//     path file_input_mix, file_input_met
//     path script_de

//     output:
//     path "output/met-decovolution-split/${dataset}_${omicmix}_${ppmix}_${fsmix}_${omicmet}_${ppmet}_${fsmet}_${de}.h5"

//     script:
//     """
//     mkdir -p output/met-decovolution-split/
//     RCODE="input_file_mix='${file_input_mix}'; input_file_met='${file_input_met}'; output_file='${output}'; script_de_met='${script_de}'; source('04_decovolution_MET_unit_pipA.R');"
//     echo \$RCODE | Rscript - 2>&1 > logs/04_${dataset}_${omicmix}_${ppmix}_${fsmix}_${omicmet}_${ppmet}_${fsmet}_${de}.txt
//     """
// }

// process late_integration {
//     input:
//     path input_file_rna, input_file_met
//     path script_li

//     output:
//     path "output/prediction/${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}.h5"

//     script:
//     """
//     mkdir -p output/prediction/
//     RCODE="input_file_rna='${input_file_rna}'; input_file_met='${input_file_met}'; output_file='${output}'; script_file='${script_li}'; source('05_late_integration_A.R');"
//     echo \$RCODE | Rscript - 2>&1 > logs/05_${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}.txt
//     """
// }

// process scoring {
//     input:
//     path prediction
//     path groundtruth_file
//     path scoring_script

//     output:
//     path "output/scores/${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}_score.h5"

//     script:
//     """
//     mkdir -p output/scores/
//     RCODE="prediction_file='${prediction}'; groundtruth_file='${groundtruth_file}'; score_file='${output}'; source('${scoring_script}');"
//     echo \$RCODE | Rscript - 2>&1 > logs/06_${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}_score.txt
//     """
// }
