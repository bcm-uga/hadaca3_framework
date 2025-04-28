#!/usr/bin/env nextflow

import groovy.yaml.YamlSlurper

nextflow.enable.dsl=2

workDir='~/project/hadaca3_framework'


params.EXECUTE = true

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
    script_05 : '05_late_integration_A.R',
    script_06 : '06_scoring.R',
    script_07 : '07_metaanalysis.Rmd'
]

params.utils = "utils/data_processing.R"

// def get_omic(path) {
//     return path.split('/')[-2]
// }

// def get_dataset(path) {
//     return path.split('/')[-1].split('_')[0]
// }

// def block_combinaison(path, add_omics = true) {
//     def l_path = path.split('/')
//     def omic = l_path[-2]
//     def descriptif = l_path[-1].split('_')
//     return add_omics ? [omic] + descriptif[1..-1] : descriptif[1..-1]
// }


workflow { 


    // ################## Reading yml file and populating CONFIG file. 
     def CONFIG = [:]

    params.config_files.each { key, filePath ->
        def parsedContent = new YamlSlurper().parse(filePath as File)
        CONFIG[key] = parsedContent
    }
    original_datasets_files = CONFIG['datasets'].collect { k, v -> v['path'] }.flatten()
    cleaned_datasets_files = CONFIG['datasets'].keySet().collect { "output/mixes/${it}" }


    // ################## Computing cleaned mix and ref cleand and generating a tuple with key as the dataset or ref name and output file. 

    def utils_channel =  Channel.fromPath(params.utils)

    // computing ref 


    ref_input = Channel.of (
        tuple(     
        [ id: "ref",
        // ref: "ref", 
        output : "cleaning-ref-ref.h5"
        ],
        file(params.reference[0]), 
        file(params.cleaner),  
        file(params.wrapper.script_01_ref),
        file(params.utils)
        )
    )

    
    out_cleaned_ref = ref_input |cleaning_ref 

    
    // out_ref_keyed = out_cleaned_ref.map { file ->
    //     def key = file.baseName 
    //     tuple(key, file)
    // }
    // computing mixes


    dataset_tuple = []
    CONFIG.datasets.each{dt,dtv-> 
        dataset_tuple.add( tuple( 
                [ id: dt,
                output : "cleaning-mix-"+dt+".h5"
                ],
                file(dtv.path),
                file(params.cleaner),
                file(params.wrapper.script_01_mix),
                file(params.utils)                
            ))
        }

    out_mix =  Channel.fromList(dataset_tuple) |cleaning_mix
    
    // out_mix.view()
        // out_mix = combined_inputs_cleaning_mix |cleaning_mix

//     out_mix_keyed = out_mix.map { file ->
//         def key = file.baseName 
//         tuple(key, file)
//     }
// .concat( Channel.of(tuple('none',file('none')) ) )  




//     // ################## Generate combinaison and prediction for the preprocess 


    // out_mix.map{  mix_meta,mix_file -> 
    // def pp_mix_path = []
    //     CONFIG['pre_proc'].each { pp, ppv ->
    //         params.mixomics.each { omic ->
    //             if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')){
    //                 pp_mix_path.add(tuple(
    //                     [ id: dt,
    //                     omic, 
    //                     output : "cleaned-"+dt+".h5"
    //                     ],
    //                     file(ppv['path'])))
    //             }
    //         }
    //     }


    // }

    pp_mix_path = []
    CONFIG['pre_proc'].each { pp, ppv ->
        params.mixomics.each { omic ->
            if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')){
                pp_mix_path.add(tuple(
                    [ pp_fun: pp,
                    omic: omic, 
                    ],
                    file(ppv['path'])))
            }
        }
    }

    pp_mix = Channel.fromList( pp_mix_path)
    .combine(out_mix)
    .combine(out_cleaned_ref)
    .map{pp_meta,pp_file,mix_meta,mix_file,ref_met,ref_file ->
        def dup_pp_meta = pp_meta.clone()
        dup_pp_meta['dataset'] = mix_meta.id
        dup_pp_meta['ref']=ref_met.id
        dup_pp_meta['output']= "out-prepross-"+[dup_pp_meta.omic,mix_meta.id, ref_met.id,dup_pp_meta.pp_fun ].join('_')+'.h5'
        tuple(dup_pp_meta,pp_file,mix_file,ref_file,file(params.wrapper.script_02),file(params.utils))
    }
    pp_ref_path = []
    CONFIG['pre_proc'].each { pp, ppv ->
        params.refomics.each { omic ->
            if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')){
                pp_ref_path.add(tuple(
                    [ pp_fun: pp,
                    omic: omic, 
                    ],
                file(ppv['path']),
                file('none')))
            }
        }
    }
    pp_ref =  Channel.fromList( pp_ref_path)
    .combine(out_cleaned_ref)
    .map{pp_meta,pp_file,mix_file,ref_met,ref_file ->
        def dup_pp_meta = pp_meta.clone() 
        dup_pp_meta['dataset'] = 'none'
        dup_pp_meta['ref']=ref_met.id
        dup_pp_meta['output']= "out-prepross-"+[dup_pp_meta.omic, dup_pp_meta.dataset, ref_met.id,   dup_pp_meta.pp_fun ].join('_')+'.h5'
        tuple(dup_pp_meta,pp_file,mix_file,ref_file,file(params.wrapper.script_02),file(params.utils))
    }
    
    out_pp = pp_ref.concat(pp_mix) | preprocessing


//     // ################## Generate combinaison and prediction for  features selection


    fs_files = out_pp.map { meta , last_pp_file ->
        def results = []
        CONFIG['features_selection'].each { fs, fsv ->
            def dup_meta = meta.clone() 

            if (fsv['omic'].contains(dup_meta.omic) || fsv['omic'].contains('ANY')) {
                dup_meta['fs_fun'] = fs
                results.add( 
                                [
                                    dup_meta,
                                file(fsv['path']),
                                file(last_pp_file),
                                file(params.wrapper.script_03),
                                file(params.utils)
                                ]
                )
            }   
         }
        return results
    }.flatMap()
    

    out_mix_with_none = out_mix.concat(Channel.of(tuple( [id:'none'],file('none'))))
    

    complete_fs_files = fs_files
    .combine(out_mix_with_none)
    .filter { fs_meta, a,b,c,d, dataset_meta, dataset_file ->
        fs_meta.dataset == dataset_meta.id 
    }
    .combine(out_cleaned_ref)
    .filter {fs_meta, a,b,c,d, dataset_meta, dataset_file, ref_meta, ref_file ->
        fs_meta.ref == ref_meta.id 
    }
    .map {fs_meta, a,b,c,d, dataset_meta, dataset_file, ref_meta, ref_file ->
        def dup_fs_meta = fs_meta.clone()
        dup_fs_meta.output ="out-fs-"+ [dup_fs_meta.omic,dup_fs_meta.dataset, dup_fs_meta.ref,dup_fs_meta.pp_fun, dup_fs_meta.fs_fun ].join('_')+'.h5'
        tuple(dup_fs_meta, a,b,c,d,dataset_file,ref_file )
    }

    out_fs = complete_fs_files | features_selection

    // out_fs.view()

    // complete_fs_files.view{ v -> v[0].dataset}


// ################## Generate combinaison for the RNA unit 

    deco_path =  []
    CONFIG.deconvolution.each {de,dev -> deco_path.add( [ [de_fun : de]  , file(dev.path)])}

    de_channel = Channel.fromList(deco_path)

    fs_mixRNA = out_fs.filter{meta,_ -> meta.omic=='mixRNA'  }
    fs_RNA =out_fs.filter{meta,_ -> meta.omic=='RNA'}
    fs_scRNA = out_fs.filter{meta,_ -> meta.omic=='scRNA'}

    de_rna_unit = 
    de_channel.combine(fs_mixRNA).combine(fs_RNA).combine(fs_scRNA).combine(out_mix).combine(out_cleaned_ref)
    .filter{ meta_de,  de_script, meta_mix,file_input_mix, meta_RNA,file_input_RNA,meta_scRNA,file_input_scRNA,   dataset_meta, dataset_file, ref_meta, ref_file ->
        meta_scRNA.ref == ref_meta.id && meta_RNA.ref == ref_meta.id && meta_mix.ref == ref_meta.id && meta_mix.dataset == dataset_meta.id
    }
    .map{ meta_de, de_script, meta_mix,file_input_mix, meta_RNA,file_input_RNA,meta_scRNA,file_input_scRNA,   dataset_meta, dataset_file, ref_meta, ref_file->
        // def meta_unit_rna = 
        def dup_meta_de =meta_de.clone()
        dup_meta_de['mixRNA'] = meta_mix
        dup_meta_de['RNA'] = meta_RNA
        dup_meta_de['scRNA'] = meta_scRNA
        dup_meta_de['dataset'] = dup_meta_de.mixRNA.dataset
        dup_meta_de['ref'] = dup_meta_de.mixRNA.ref
        def output_name = "out-derna-" + [dup_meta_de.dataset,dup_meta_de.ref].join('_')  +
                [dup_meta_de.mixRNA.omic, dup_meta_de.mixRNA.pp_fun, dup_meta_de.mixRNA.fs_fun ].join('_')   +
                [dup_meta_de.RNA.omic, dup_meta_de.RNA.pp_fun, dup_meta_de.RNA.fs_fun ].join('_')   +
                [dup_meta_de.scRNA.omic, dup_meta_de.scRNA.pp_fun, dup_meta_de.scRNA.fs_fun ].join('_')   +'.h5'
                
        dup_meta_de.output =output_name
        tuple(dup_meta_de,de_script,file_input_mix, file_input_RNA,file_input_scRNA, dataset_file,ref_file)
    }
    .combine(Channel.of(tuple(file(params.wrapper.script_04_rna),file(params.utils))))



    out_de_rna_unit = de_rna_unit | prediction_deconvolution_rna 

    

//     // ################## Generate combinaison for the MET unit 

    fs_mixMET = out_fs.filter{meta,_ -> meta.omic=='mixMET'  }
    fs_MET =out_fs.filter{meta,_ -> meta.omic=='MET'  }

    deco_path_met =  []
    CONFIG.deconvolution.each {de,dev -> deco_path_met.add( [ [de_fun : de]  , file(dev.path)])}

    de_channel_met = Channel.fromList(deco_path_met)


    de_met_unit = 
    de_channel_met.combine(fs_mixMET).combine(fs_MET).combine(out_mix).combine(out_cleaned_ref)
    .filter{ meta_de,  de_script, meta_mix, file_input_mix, meta_MET,file_input_MET,dataset_meta, dataset_file, ref_meta, ref_file ->
        meta_MET.ref == ref_meta.id && meta_mix.ref == ref_meta.id && meta_mix.dataset == dataset_meta.id
    }
    .map{meta_de, de_script, meta_mix,file_input_mix, meta_MET,file_input_MET,dataset_meta, dataset_file, ref_meta, ref_file ->
        def dup_meta_de = meta_de.clone()
        dup_meta_de['mixMET'] = meta_mix
        dup_meta_de['MET'] = meta_MET
        dup_meta_de['dataset'] = meta_mix.dataset
        dup_meta_de['ref'] = meta_mix.ref
        def output_name = "out-demet-" + [dup_meta_de.dataset,dup_meta_de.ref].join('_')  +
                [dup_meta_de.mixMET.omic, dup_meta_de.mixMET.pp_fun, dup_meta_de.mixMET.fs_fun ].join('_')   +
                [dup_meta_de.MET.omic, dup_meta_de.MET.pp_fun, dup_meta_de.MET.fs_fun ].join('_')   +'.h5'
                
        dup_meta_de.output =output_name
        tuple(dup_meta_de, de_script,file_input_mix, file_input_MET, dataset_file,ref_file)
    }
    .combine(Channel.of(tuple(file(params.wrapper.script_04_met),file(params.utils))))


    out_de_met_unit= de_met_unit | prediction_deconvolution_met




    // out_de_met_unit.count().view()
    // out_de_met_unit.last().view()

    // de_met_unit.filter{v ->
    //     v[0].mixMET == 'fsID'
    // }.count().view()
    
    // de_met_unit

// ################## Generate combinaison for late integration

    li_path =  []
    CONFIG.late_integration.each {li,liv -> li_path.add([ [li:li ] , file(liv.path)])}

    li_channel = Channel.fromList(li_path)

    li_combinaison = 
    li_channel.combine(out_de_rna_unit).combine(out_de_met_unit)
    .combine(out_mix).combine(out_cleaned_ref)
    // .filter{m, li_dic, meta_rna, file_rna,    meta_met, file_met         -> 
    //     meta_met.mixMET == 'fsID'
    // }
    // .filter{m, li_dic, meta_rna, file_rna,    meta_met, file_met         -> 
    //     meta_rna.dataset == meta_met.dataset
    // }
    

    .filter{m, li_dic, meta_rna, file_rna,    meta_met, file_met,           dataset_meta, dataset_file,   ref_meta, ref_file-> 
        meta_rna.dataset == dataset_meta.id &&         meta_met.dataset == dataset_meta.id && 
        meta_rna.ref == ref_meta.id         &&         meta_met.ref ==ref_meta.id
    }
    .map{
        li_meta,li_path,meta_rna, file_rna , meta_met, file_met ,dataset_meta, dataset_file, ref_meta, ref_file ->
        def dup_li_meta = li_meta.clone()
        dup_li_meta['dataset'] = meta_rna.dataset
        dup_li_meta['ref'] = meta_rna.ref
        dup_li_meta['rna_unit'] = meta_rna
        dup_li_meta['met_unit'] = meta_met

        def output_name = "out-li-" + [dup_li_meta.dataset,dup_li_meta.ref].join('_')  +
            [dup_li_meta.rna_unit.mixRNA.pp_fun, dup_li_meta.rna_unit.mixRNA.fs_fun ].join('_')   +
            [dup_li_meta.rna_unit.RNA.pp_fun, dup_li_meta.rna_unit.RNA.fs_fun ].join('_')           +
            [dup_li_meta.rna_unit.scRNA.pp_fun, dup_li_meta.rna_unit.scRNA.fs_fun ].join('_')      +
            [dup_li_meta.met_unit.mixMET.pp_fun, dup_li_meta.met_unit.mixMET.fs_fun ].join('_')   +
            [dup_li_meta.met_unit.MET.pp_fun, dup_li_meta.met_unit.MET.fs_fun ].join('_')           +'.h5'
        dup_li_meta["output"] = output_name
        tuple( dup_li_meta,li_path , file_rna , file_met, dataset_file,ref_file )
    }.combine(Channel.of(tuple(file(params.wrapper.script_05),file(params.utils))))


    // li_combinaison.view{ v->  """${v[2].dataset} ${v[4].dataset}
    // ${v[2].mixRNA}
    // ${v[2].RNA}
    // ${v[2].scRNA}
    
    // ${v[4].mixMET}
    // ${v[4].MET}    
    // \n\n"""  }
    

    
    li_combinaison.count().view()



    out_li = li_combinaison | late_integration 

    // out_li.last().view()
    // out_li.count().view()

//     // ################## Generate score 

    score_input = out_li.map{   meta,file_path  -> 

        def dup_meta = meta.clone()
        // dup_meta.output = 
        def output_name = "score-" + [dup_meta.dataset,dup_meta.ref].join('_')  +
            [dup_meta.rna_unit.mixRNA.pp_fun, dup_meta.rna_unit.mixRNA.fs_fun ].join('_')   +
            [dup_meta.rna_unit.RNA.pp_fun, dup_meta.rna_unit.RNA.fs_fun ].join('_')           +
            [dup_meta.rna_unit.scRNA.pp_fun, dup_meta.rna_unit.scRNA.fs_fun ].join('_')      +
            [dup_meta.met_unit.mixMET.pp_fun, dup_meta.met_unit.mixMET.fs_fun ].join('_')   +
            [dup_meta.met_unit.MET.pp_fun, dup_meta.met_unit.MET.fs_fun ].join('_')           +'.h5'
        dup_meta["output"] = output_name
        tuple(dup_meta, file_path, file(CONFIG.datasets[dup_meta.dataset].groundtruth_file_path),file(params.wrapper.script_06),file(params.utils))
     }   

    // score_input.last().view()
    // score_input.count().view()


    score_out = score_input | scoring

    // score_out.last().view()
    // score_out.count().view()

}

process cleaning_mix {
    input:
    tuple val(meta),  
    path(mixes), 
    path(cleaner), 
    path(wrapper01),
    path(utils)

    output:
    // path "${mixes.baseName}.h5"
    // tuple val(meta), path("cleaning-mix-*.h5")
    tuple val(meta), path("${meta.output}")
    

    script:
    """
    RCODE="mixes_file='${mixes}'; output_file='${meta.output}'; 
    utils_script='${utils}'; cleaner='${cleaner}';
     source('${wrapper01}');"
    echo \$RCODE | Rscript -
    """
    
    stub:
    """
    RCODE="mixes_file='${mixes}'; output_file='${meta.output}'; 
    utils_script='${utils}'; cleaner='${cleaner}';
     source('${wrapper01}');"
    echo \$RCODE
    touch ${meta.output}
    """
    // workingDir : '.'
}


process cleaning_ref{
    input:
    tuple val(meta), 
    path(reference),
    path(cleaner) ,
    path(wrapper01),
    path(utils)
    

    output:
    // tuple val(meta), path("cleaning-ref*.h5")
    tuple val(meta), path("${meta.output}")



    script:
    """
    RCODE="reference_file='${reference}'; output_file='${meta.output}'; 
    utils_script='${utils}';cleaner='${cleaner}'; 
    source('${wrapper01}');"
    echo \$RCODE | Rscript -
    """

    stub : 
        """
    RCODE="reference_file='${reference}'; output_file='${meta.output}'; 
    utils_script='${utils}';cleaner='${cleaner}'; 
    source('${wrapper01}');"
    echo \$RCODE 
    touch ${meta.output}
    """
}
    // echo \$RCODE | Rscript - 2>&1 > logs/01_${reference.baseName}.txt

process preprocessing {
    cpus 1

    input:
    tuple val(meta),
        path(pp_script), 
        path(mix), 
        path(reference), 
        path(wrapper02),
        path(utils)
    output:
    tuple val(meta), path("${meta.output}")

    script:
    """
    RCODE=" omic='${meta.omic}'; 
    mixes_file='${mix}'; reference_file='${reference}'; 
    output_file='${meta.output}'; 
    utils_script='${utils}'; 
    script_file='${pp_script}'; 
    source('${wrapper02}');"
    echo \$RCODE | Rscript -
    """

    stub : 
    """
    RCODE=" omic='${meta.omic}';
    mixes_file='${mix}'; reference_file='${reference}'; 
    output_file='${meta.output}'; 
    utils_script='${utils}'; 
    script_file='${pp_script}'; 
    source('${wrapper02}');"
    echo \$RCODE 
    touch ${meta.output}
    """
}
    // echo \$RCODE | Rscript - 2>&1 > logs/02_${omic}_${dataset}_${pp}.txt

process features_selection {
    cpus 1
    
    input:
        tuple( val(meta),
            path(fs_script), 
            path(file_input),
            path(wrapper03),
            path(utils),
            path(mix), 
            path(reference),
            )


    output:
    tuple val(meta), path("${meta.output}")

    // tuple(val(meta),path("out-fs-*.h5"))

    script:
    """
        RCODE="omic='${meta.omic}'; 
        input_file='${file_input}'; output_file='${meta.output}';
        path_ogmix='${mix}' ; path_ogref='${reference}' ; 
        script_file='${fs_script}'; 
        utils_script='${utils}'; 
        source('${wrapper03}');"
        echo \$RCODE | Rscript - 
    """

    stub:
    """
        RCODE="omic='${meta.omic}'; 
        input_file='${file_input}'; output_file='${meta.output}'; 
        path_ogmix='${mix}' ; path_ogref='${reference}' ; 
        script_file='${fs_script}'; 
        utils_script='${utils}'; 
        source('${wrapper03}');"
        echo \$RCODE 
        touch ${meta.output}
    """
}

process prediction_deconvolution_rna {
    cpus 1
    
     input:
        tuple val(meta),
        path(de_script), 
        path(file_input_mix),
        path(file_input_rna),
        path(file_input_scrna),
        path(mix), 
        path(reference), 
        path(wrapper04),
        path(utils)

    output:
    tuple val(meta), path("${meta.output}")

    // tuple(val(meta),path("out-derna-*.h5"))s

    // tuple(val(mix.baseName),val(reference.baseName),path("${meta.output}"))

    script:
    """
        RCODE="
        input_file_mix='${file_input_mix}'; input_file_rna='${file_input_rna}'; input_file_sc='${file_input_scrna}';
        output_file='${meta.output}';
        path_ogmix='${mix}' ; path_ogref='${reference}' ; 
        script_de_rna='${de_script}'; 
        utils_script='${utils}'; 
        source('${wrapper04}');"
        echo \$RCODE | Rscript - 
    """

    stub:
    """
        RCODE="
        input_file_mix='${file_input_mix}'; input_file_rna='${file_input_rna}'; input_file_sc='${file_input_scrna}';
        output_file='${meta.output}'; 
        path_ogmix='${mix}' ; path_ogref='${reference}' ; 
        script_de_rna='${de_script}'; 
        utils_script='${utils}'; 
        source('${wrapper04}');"
        echo \$RCODE 
        touch ${meta.output}
    """
}

process prediction_deconvolution_met {
    cpus 1
     input:
        tuple val(meta),
        path(de_script), 
        path(file_input_mix),
        path(file_input_met),
        path(mix), 
        path(reference), 
        path(wrapper04),
        path(utils)

    output:
    tuple val(meta), path("${meta.output}")

    // tuple(val(meta),path("out-demet*.h5"))

    script:
    """
        RCODE="
        input_file_mix='${file_input_mix}'; input_file_met='${file_input_met}';
        output_file='${meta.output}';
        path_ogmix='${mix}' ; path_ogref='${reference}' ; 
        script_de_met='${de_script}'; 
        utils_script='${utils}'; 
        source('${wrapper04}');"
        echo \$RCODE | Rscript - 
    """

    stub:
    """
        RCODE="
        input_file_mix='${file_input_mix}'; input_file_met='${file_input_met}';
        output_file='${meta.output}'; 
        path_ogmix='${mix}' ; path_ogref='${reference}' ; 
        script_de_met='${de_script}'; 
        utils_script='${utils}'; 
        source('${wrapper04}');"
        echo \$RCODE 
        touch ${meta.output}
    """
}

process late_integration {
    cpus 1
    
    input:
    tuple val(meta),
    path(script_li),
    path(input_file_rna), 
    path(input_file_met),
    path(mix), 
    path(reference), 
    path(wrapper04),
    path(utils)

    output:
    tuple val(meta), path("${meta.output}")

    // tuple(val(meta),path("out-li*.h5"))

    // path "${mix.baseName}_${reference.baseName}_${script_li.baseName}.h5"

    script:
    """
    RCODE="input_file_rna='${input_file_rna}'; input_file_met='${input_file_met}'; 
    output_file='${meta.output}'; 
    script_file='${script_li}'; 
    utils_script='${utils}'; 
    source('${wrapper04}');"
    echo \$RCODE | Rscript - 
    """

    stub:
    """
    RCODE="input_file_rna='${input_file_rna}'; input_file_met='${input_file_met}'; 
    output_file='${meta.output}'; 
    script_file='${script_li}'; 
    utils_script='${utils}'; 
    source('${wrapper04}');"
    echo \$RCODE
    touch ${meta.output}
    """
}
    // path "output/prediction/${mix.baseName}_${reference.baseName}_${input_file_rna.baseName}_${input_file_met.baseName}_${script_li.baseName}.h5"
    // 2>&1 > logs/05_${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}.txt

process scoring {
    cpus 1
    
    input:   
    tuple  val(meta) ,
    path(prediction),
    path(groundtruth_file),
    path(scoring_script),
    path(utils)

    output:
    tuple val(meta), path("${meta.output}")

    // path "score.h5"

    script:
    """
    RCODE="prediction_file='${prediction}'; groundtruth_file='${groundtruth_file}'; 
    score_file='${meta.output}'; 
    utils_script='${utils}'; 
    source('${scoring_script}');"
    echo \$RCODE | Rscript - 
    """

    stub:
    """
    RCODE="prediction_file='${prediction}'; groundtruth_file='${groundtruth_file}'; 
    score_file='${meta.output}'; 
    utils_script='${utils}'; 
    source('${scoring_script}');"
    echo \$RCODE 
    touch ${meta.output}
    """
}

    // 2>&1 > logs/06_${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}_score.txt
    // path "output/scores/${dataset}_${omicMixRna}_${ppMixRna}_${fsMixRna}_${omicRNA}_${ppRNA}_${fsRNA}_${omicSCRNA}_${ppSCRNA}_${fsSCRNA}_${deRNA}_${omicMixMet}_${ppMixMet}_${fsMixMet}_${omicMET}_${ppMET}_${fsMET}_${deMET}_${li}_score.h5"
