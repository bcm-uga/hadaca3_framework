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
    script_01_mix: '01_global_preprocess_mix.R',
    script_01_ref : '01_global_preprocess_ref.R', 
    script_02 : "02_preprocess.R", 
    script_03 : '03_features_selection.R',
    script_04_rna : '04_decovolution_RNA_unit_pipA.R',
    script_04_met : '04_decovolution_MET_unit_pipA.R',
    script_04_met : '04_decovolution_MET_unit_pipA.R',
    script_05 : '05_late_integration_A.R',
    script_06_scoring : '06_scoring.R',
    script_07_anal : '07_metaanalysis.Rmd'
]

params.utils = "utils/data_processing.R"

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


    // ################## Reading yml file and populating CONFIG file. 
     def CONFIG = [:]

    params.config_files.each { key, filePath ->
        def parsedContent = new YamlSlurper().parse(filePath as File)
        CONFIG[key] = parsedContent
    }
    original_datasets_files = CONFIG['datasets'].collect { k, v -> v['path'] }.flatten()
    cleaned_datasets_files = CONFIG['datasets'].keySet().collect { "output/mixes/${it}" }
    cleaned_REFERENCE = params.reference.collect { "output/ref/${it.split('/')[-1]}" }


    // // Generate de_met_unit_files
    // def de_met_unit_files = []
    // fs_files.each { mix_combi ->
    //     fs_files.each { met_combi ->
    //         if (get_omic(mix_combi) == 'mixMET' && get_omic(met_combi) == 'MET') {
    //             CONFIG['deconvolution'].each { de, dev ->
    //                 de_met_unit_files.add("output/met-decovolution-split/${get_dataset(mix_combi)}_${block_combinaison(mix_combi).join('_')}_${block_combinaison(met_combi).join('_')}_${de}")
    //             }
    //         }
    //     }
    // }

    
    // li_files = []
    // de_rna_unit_files.each { file_prediction_RNA ->
    //     de_met_unit_files.each { file_prediction_MET ->
    //         if (get_dataset(file_prediction_RNA) == get_dataset(file_prediction_MET)) {
    //             CONFIG['late_integration'].keySet().each { li ->
    //                 li_files.add("output/prediction/${get_dataset(file_prediction_RNA)}_${block_combinaison(file_prediction_RNA, false).join('_')}_${block_combinaison(file_prediction_MET, false).join('_')}_${li}")
    //             }
    //         }
    //     }
    // }

    // scores_files = li_files.collect { "output/scores/${it.split('/')[-1]}_score" }


    def utils_channel =  Channel.fromPath(params.utils)

    // computing ref 
    out_cleaned_ref =     cleaning_ref( 
        Channel.fromPath(params.reference),
        Channel.fromPath(params.cleaner) ,
        Channel.fromPath(params.wrapper['script_01_ref']), 
        Channel.fromPath(params.utils))
    
    out_ref_keyed = out_cleaned_ref.map { file ->
        def key = file.baseName 
        tuple(key, file)
    }

    // ################## Computing cleaned mix and ref cleand and generating a tuple with key as the dataset or ref name and output file. 

    // computing mixes
    Channel.fromPath(original_datasets_files)
    .combine (Channel.fromPath(params.cleaner))
    .combine (Channel.fromPath(params.wrapper['script_01_mix']))
    .combine(utils_channel)
    .set {  combined_inputs_cleaning_mix}
    
    out_mix = combined_inputs_cleaning_mix |cleaning_mix

    out_mix_keyed = out_mix.map { file ->
        def key = file.baseName 
        tuple(key, file)
    }.concat( Channel.of(tuple('none',file('none')) ) )  


    // ################## Generate combinaison and prediction for the preprocess 

    
    pp_mix_path = []
    CONFIG['pre_proc'].each { pp, ppv ->
        params.mixomics.each { omic ->
            if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')){
                pp_mix_path.add(tuple(omic,file(ppv['path'])))
            }
        }
    }
    pp_ref_path = []
    CONFIG['pre_proc'].each { pp, ppv ->
        params.refomics.each { omic ->
            if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')){
                pp_ref_path.add(tuple(omic,file(ppv['path']),file('none')))
            }
        }
    }


    Channel.fromList( pp_mix_path)
    .combine(out_mix)
    .combine(out_cleaned_ref)
    .combine(Channel.fromPath(params.wrapper.script_02))
    .combine(utils_channel)
    .set{pp_mix}

    Channel.fromList( pp_ref_path)
    .combine(out_cleaned_ref)
    .combine(Channel.fromPath(params.wrapper.script_02))
    .combine(utils_channel)
    .set{pp_ref}
    
    out_pp = pp_ref.concat(pp_mix) | preprocessing

    // ################## Generate combinaison and prediction for  features selection
    fs_files = out_pp.map { last_pp ->
        def omic = last_pp[0]
        def dataset = last_pp[1]
        def ref = last_pp[2]

        def results = []
        CONFIG['features_selection'].each { fs, fsv ->
            if (fsv['omic'].contains(omic) || fsv['omic'].contains('ANY')) {
                results.add( tuple(
                                dataset,
                                ref, 
                                omic,
                                file(fsv['path']),
                                file(last_pp[-1]),
                                file(params.wrapper.script_03),
                                file(params.utils)
                ))
            }   
         }
        return results
    }.flatten().collate(7)

    // Add the correct file from out_mix to each fs_files tuple
    complete_fs_files = fs_files
    .combine(out_mix_keyed)
    .filter { fs_dataset_key, fs_ref_key, omic, script, input_file, fs_script,utils, dataset_key, dataset_file ->
        fs_dataset_key == dataset_key 
    }
    .combine(out_ref_keyed)
    .filter {fs_dataset_key, fs_ref_key, omic, script, input_file, fs_script,utils, dataset_key, dataset_file,ref_key,ref_file ->
        fs_ref_key == ref_key 
    }
    .map { fs_dataset_key, fs_ref_key, omic, script, input_file, fs_script,utils, dataset_key, dataset_file,ref_key,ref_file->
        tuple(omic, script, input_file, dataset_file,ref_file,fs_script,utils)
    }

    out_fs = complete_fs_files | features_selection

    // ################## Generate combinaison for the RNA unit 

        // def de_rna_unit_files = []
    // fs_files.each { mix_combi ->
    //     fs_files.each { rna_combi ->
    //         fs_files.each { sc_combi ->
    //             if (get_omic(mix_combi) == 'mixRNA' && get_omic(rna_combi) == 'RNA' && get_omic(sc_combi) == 'scRNA') {
    //                 CONFIG['deconvolution'].each { de, dev ->
    //                     de_rna_unit_files.add("output/rna-decovolution-split/${get_dataset(mix_combi)}_${block_combinaison(mix_combi).join('_')}_${block_combinaison(rna_combi).join('_')}_${block_combinaison(sc_combi).join('_')}_${de}")
    //                 }
    //             }
    //         }
    //     }
    // }


    de_rna_unit_files = out_fs.map { last_fs ->
        def omic = last_fs[0]
        def dataset = last_fs[1]
        def ref = last_fs[2]
        def output_file = last_fs[-1]
        def results = []
        CONFIG['features_selection'].each { fs, fsv ->
            if (fsv['omic'].contains(omic) || fsv['omic'].contains('ANY')) {
                results.add( tuple(
                                dataset,
                                ref, 
                                omic,
                                file(fsv['path']),
                                file(last_pp[-1]),
                                file(params.wrapper.script_03),
                                file(params.utils)
                ))
            }   
         }
        return results
    }.flatten().collate(7)

    // fs_files = out_pp.map { last_pp ->
    //     def omic = last_pp[0]
    //     println omic
    //     def results = []
    //     params.CONFIG['features_selection'].each { fs, fsv ->
    //             results.add(file(fsv['path']))
    //         }
    //     }
    //     return results
    // }.flatten()

    // fs_files.view()

    // pp_mix.view()
    // pp_ref.concat(pp_mix).count().view()

    // pp_mix.combine(pp_ref.collect()).view()
    // out_pp_mix = pp_mix | preprocessing 
    // Channel.fromList([1, 2, 3]).set { A }
    // Channel.fromList(['a', 'b', 'c']).set { B }

    // A.combine((B.collect()))
    //     .view()



//     prediction_deconvolution_rna(de_rna_unit_files)
//     prediction_deconvolution_met(de_met_unit_files)
//     late_integration(li_files)
//     scoring(scores_files)


}

process cleaning_mix {
    input:
    tuple path(mixes), path(cleaner), path(wrapper01), path(utils)

    output:
    path "output/mixes/${mixes.baseName}.h5"

    script:
    """
    mkdir -p output/mixes/
    RCODE="mixes_file='${mixes}'; output_file='output/mixes/${mixes.baseName}.h5'; 
    utils_script='${utils}'; cleaner='${cleaner}';
     source('${wrapper01}');"
    echo \$RCODE | Rscript -
    """
    
    stub:
    """
    mkdir -p output/mixes/
    RCODE="mixes_file='${mixes}'; output_file='output/mixes/${mixes.baseName}.h5'; 
    utils_script='${utils}'; cleaner='${cleaner}';
     source('${wrapper01}');"
    echo \$RCODE
    touch output/mixes/${mixes.baseName}.h5
    """
    // workingDir : '.'
}


process cleaning_ref{
    input:
    path reference
    path cleaner
    path wrapper01
    path utils

    output:
    path  "output/ref/${reference.baseName}.h5"

    script:
    """
    mkdir -p output/ref/
    RCODE="reference_file='${reference}'; output_file='output/ref/${reference.baseName}.h5'; 
    utils_script='${utils}';cleaner='${cleaner}'; 
    source('${wrapper01}');"
    echo \$RCODE | Rscript -
    """

    stub : 
        """
    mkdir -p output/ref/
    RCODE="reference_file='${reference}'; output_file='output/ref/${reference.baseName}.h5'; 
    utils_script='${utils}';cleaner='${cleaner}'; 
    source('${wrapper01}');"
    echo \$RCODE 
    touch output/ref/${reference.baseName}.h5
    """
}
    // echo \$RCODE | Rscript - 2>&1 > logs/01_${reference.baseName}.txt

process preprocessing {
    cpus 1

    input:
    tuple val(omic),
        path(pp_script), 
        path(mix), 
        path(reference), 
        path(wrapper02),
        path(utils)
    output:
    tuple val(omic),val(mix.baseName),val(reference.baseName), path("output/preprocessing/${omic}/${reference.baseName}_${mix.baseName}_${pp_script.baseName}.h5")

    script:
    """
    mkdir -p output/preprocessing/${omic}/
    RCODE=" omic='${omic}'; 
    mixes_file='${mix}'; reference_file='${reference}'; 
    output_file='output/preprocessing/${omic}/${reference.baseName}_${mix.baseName}_${pp_script.baseName}.h5'; 
    utils_script='${utils}'; 
    script_file='${pp_script}'; 
    source('${wrapper02}');"
    echo \$RCODE | Rscript -
    """
    // 

    stub : 
        """
    mkdir -p output/preprocessing/${omic}/
    RCODE=" omic='${omic}';
    mixes_file='${mix}'; reference_file='${reference}'; 
    output_file='output/preprocessing/${omic}/${reference.baseName}_${mix.baseName}_${pp_script.baseName}.h5'; 
    utils_script='${utils}'; 
    script_file='${pp_script}'; 
    source('${wrapper02}');"
    echo \$RCODE 
    touch output/preprocessing/${omic}/${reference.baseName}_${mix.baseName}_${pp_script.baseName}.h5
    """
}
    // echo \$RCODE | Rscript - 2>&1 > logs/02_${omic}_${dataset}_${pp}.txt

process features_selection {
    cpus 1
    
    input:
        tuple val(omic),
        path(fs_script), 
        path(file_input),
        path(mix), 
        path(reference), 
        path(wrapper03),
        path(utils)


    output:
    tuple(val(omic),val(mix.baseName),val(reference.baseName),path("output/feature-selection/${omic}/${file_input.baseName}_${fs_script.baseName}.h5"))

    script:
    """
        mkdir -p output/feature-selection/${omic}/
        RCODE="omic='${omic}'; 
        input_file='${file_input}'; output_file='output/feature-selection/${omic}/${file_input.baseName}_${fs_script.baseName}.h5';
        path_ogmix='${mix}' ; path_ogref='${reference}' ; 
        script_file='${fs_script}'; 
        utils_script='${utils}'; 
        source('${wrapper03}');"
        echo \$RCODE | Rscript - 
    """

    stub:
    """
        mkdir -p output/feature-selection/${omic}/
        RCODE="omic='${omic}'; 
        input_file='${file_input}'; output_file='output/feature-selection/${omic}/${file_input.baseName}_${fs_script.baseName}.h5'; 
        path_ogmix='${mix}' ; path_ogref='${reference}' ; 
        script_file='${fs_script}'; 
        utils_script='${utils}'; 
        source('${wrapper03}');"
        echo \$RCODE 
        touch output/feature-selection/${omic}/${file_input.baseName}_${fs_script.baseName}.h5
    """
}

process prediction_deconvolution_rna {
    cpus 1
    
     input:
        tuple val(omic),
        path(fs_script), 
        path(file_input_mix),
        path(file_input_rna),
        path(file_input_scrna),
        path(mix), 
        path(reference), 
        path(wrapper03),
        path(utils)

    path file_input_mix, file_input_rna, file_input_scrna
    path script_de

    output:
    path "output/rna-decovolution-split/${dataset}_${omicmix}_${ppmix}_${fsmix}_${omicrna}_${pprna}_${fsrna}_${omicscrna}_${ppsc}_${fssc}_${de}.h5"

    output:
    tuple(val(omic),val(mix.baseName),val(reference.baseName),path("output/feature-selection/${omic}/${file_input.baseName}_${fs_script.baseName}.h5"))
    
    script:
    """
    mkdir -p output/rna-decovolution-split/
    RCODE="input_file_mix='${file_input_mix}'; input_file_rna='${file_input_rna}'; input_file_sc='${file_input_scrna}'; output_file='${output}'; script_de_rna='${script_de}'; source('04_decovolution_RNA_unit_pipA.R');"
    echo \$RCODE | Rscript - 2>&1 > logs/04_${dataset}_${omicmix}_${ppmix}_${fsmix}_${omicrna}_${pprna}_${fsrna}_${omicscrna}_${ppsc}_${fssc}_${de}.txt
    """


    script:
    """
        mkdir -p output/feature-selection/${omic}/
        RCODE="input_file='${file_input}'; output_file='output/feature-selection/${omic}/${file_input.baseName}_${fs_script.baseName}.h5';
        mixes_file='${mix}'; reference_file='${reference}'; 
        script_file='${fs_script}'; 
        utils_script='${utils}'; 
        source('${wrapper03}');"
        echo \$RCODE | Rscript - 
    """

    stub:
    """
        mkdir -p output/feature-selection/${omic}/
        RCODE="input_file='${file_input}'; output_file='output/feature-selection/${omic}/${file_input.baseName}_${fs_script.baseName}.h5'; 
        script_file='${fs_script}'; 
        utils_script='${utils}'; 
        source('${wrapper03}');"
        echo \$RCODE 
        touch output/feature-selection/${omic}/${file_input.baseName}_${fs_script.baseName}.h5
    """
}

process prediction_deconvolution_met {
    cpus 1
    
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
    cpus 1
    
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
    cpus 1
    
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
