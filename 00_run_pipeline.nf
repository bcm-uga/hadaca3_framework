#!/usr/bin/env nextflow

import groovy.yaml.YamlSlurper

nextflow.enable.dsl=2

workDir='~/project/hadaca3_framework/'

// This should be overwritten 
params.setup_folder = './'


params.config_files = [
    datasets:                   params.setup_folder + 'datasets.yml',
    pre_proc:                   params.setup_folder + 'preprocessing.yml',
    features_selection:         params.setup_folder + 'feature_selection.yml',
    early_integration:          params.setup_folder + 'early_integration.yml',
    intermediate_integration:   params.setup_folder + 'intermediate_integration.yml',
    late_integration:           params.setup_folder + 'late_integration.yml',
    deconvolution:              params.setup_folder + 'deconvolution.yml'
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

    
    out_cleaned_ref = ref_input | Cleaning_ref 

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

    out_mix =  Channel.fromList(dataset_tuple) | Cleaning_mix
    



//     // ################## Generate combinaison and prediction for the preprocess 


    pp_mix_path = []
    CONFIG['pre_proc'].each { pp, ppv ->
        // println(ppv.getOrDefault('dependency','none'))
        params.mixomics.each { omic ->
            if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')){
                pp_mix_path.add(tuple(
                    [ pp_fun: pp,
                    omic: omic, 
                    pp_create : ppv.getOrDefault('create','none'),
                    pp_omics : ppv.omic             
                    ],
                    file(ppv['path']),
                    tuple(ppv.getOrDefault('dependency',['none_dep']).collect{f -> file(f)} )
                    // file(ppv.getOrDefault('dependency','none'))
                    // file("preprocessing/attachement/*")
                    ))
            }
        }
    }

    pp_mix = Channel.fromList( pp_mix_path)
    .combine(out_mix)
    .combine(out_cleaned_ref)
    .map{pp_meta,pp_file,file_dep,mix_meta,mix_file,ref_met,ref_file ->
        def dup_pp_meta = pp_meta.clone()
        dup_pp_meta['dataset'] = mix_meta.id
        dup_pp_meta['ref']=ref_met.id
        dup_pp_meta['output']= "out-prepross-"+[dup_pp_meta.omic,mix_meta.id, ref_met.id,dup_pp_meta.pp_fun ].join('_')+'.h5'
        tuple(dup_pp_meta,pp_file,file_dep,mix_file,ref_file,file(params.wrapper.script_02),file(params.utils))
    }
    pp_ref_path = []
    CONFIG['pre_proc'].each { pp, ppv ->
        params.refomics.each { omic ->
            if (ppv['omic'].contains(omic) || ppv['omic'].contains('ANY')){
                pp_ref_path.add(tuple(
                    [ pp_fun: pp,
                    omic: omic, 
                    pp_create : ppv.getOrDefault('create','none'),
                    pp_omics : ppv.omic
                    ],
                file(ppv['path']),
                tuple(ppv.getOrDefault('dependency', ['none_dep']).collect {f -> file(f)} ),

                file('none')
                ))
            }
        }
    }

    pp_ref =  Channel.fromList( pp_ref_path)
    .combine(out_cleaned_ref)
    .map{pp_meta,pp_file,file_dep,mix_file,ref_met,ref_file ->
        def dup_pp_meta = pp_meta.clone() 
        dup_pp_meta['dataset'] = 'none'
        dup_pp_meta['ref']=ref_met.id
        dup_pp_meta['output']= "out-prepross-"+[dup_pp_meta.omic, dup_pp_meta.dataset, ref_met.id,   dup_pp_meta.pp_fun ].join('_')+'.h5'
        tuple(dup_pp_meta,pp_file,file_dep,mix_file,ref_file,file(params.wrapper.script_02),file(params.utils))
    }
    
    // pp_ref.view()

    // out_pp = pp_ref.concat(pp_mix) | Preprocessing
    out_pp = pp_ref.mix(pp_mix) | Preprocessing


//     // ################## Generate combinaison and prediction for  features selection


    fs_files = out_pp.map { meta , last_pp_file ->
        def results = []
        def pp_create = meta.pp_create ;
        if (pp_create instanceof List) {
            pp_create = meta.pp_create[0] ; 
        }

        def pp_omics = meta.pp_omics ;  
        // println(pp_omics)
        CONFIG['features_selection'].each { fs, fsv ->
            def dup_meta = meta.clone() 
            // def pp_create = dup_meta.pp_create  ; 
            // def pp_omics = dup_meta.pp_omics ;  
            dup_meta["fs_need"] = fsv.getOrDefault('need','none')
            dup_meta["fs_omic_need"] = fsv.getOrDefault('omic_need',['none'])
            def fs_need = fsv.getOrDefault('need',['none'])
            def fs_omic_need = fsv.getOrDefault('omic_need',['none'])
            if (fsv['omic'].contains(dup_meta.omic) || fsv['omic'].contains('ANY')    ) {
                // println( dup_meta.omic + ' ' + dup_meta.pp_fun + ' ' + pp_create + ' ; ' + fs + ' ' + fs_need + ' ' +   fs_need.contains(pp_create) )
                // fs need contains pp_create (works also for none)
                // OR the omic being computed is not never created but fs need another kind of omic therfor not in pp_omicS and pp_create is not special
                if( 
                 (fs_need.contains(pp_create)  && (fs_omic_need.contains(dup_meta.omic) || fs_omic_need.contains('ANY') || fs_omic_need.contains('none') )) || 
                 (pp_create == 'none' && !fs_omic_need.contains(dup_meta.omic))
                 ){
                    dup_meta['fs_fun'] = fs
                    results.add( 
                        [
                            dup_meta,
                            file(fsv['path']),
                            file(last_pp_file),
                            file(params.wrapper.script_03),
                            file(params.utils),
                            tuple(fsv.getOrDefault('dependency', ['none_dep']).collect {f -> file(f)} ),
                        ]
                    )
                }
            }   
         }
        return results
    }.flatMap()


    

    out_mix_with_none = out_mix.concat(Channel.of(tuple( [id:'none'],file('none'))))
    complete_fs_files = fs_files
    .combine(out_mix_with_none)
    .filter { fs_meta, a,b,c,d,e, dataset_meta, dataset_file ->
        fs_meta.dataset == dataset_meta.id 
    }
    .combine(out_cleaned_ref)
    .filter {fs_meta, a,b,c,d,e, dataset_meta, dataset_file, ref_meta, ref_file ->
        fs_meta.ref == ref_meta.id 
    }
    .map {fs_meta, a,b,c,d,e, dataset_meta, dataset_file, ref_meta, ref_file ->
        def dup_fs_meta = fs_meta.clone()
        dup_fs_meta.output ="out-fs-"+ [dup_fs_meta.omic,dup_fs_meta.dataset, dup_fs_meta.ref,dup_fs_meta.pp_fun, dup_fs_meta.fs_fun ].join('_')+'.h5'
        tuple(dup_fs_meta, a,b,c,d,e,dataset_file,ref_file )
    }

    out_pp_create_filtered = out_pp.filter{pp_meta, out_file -> 
        pp_meta.pp_create !='none'
    }
    
    // out_pp_create_filtered.view()

    complete_fs_files.branch{fs_meta, a,b,c,d,e,dataset_file,ref_file  -> 
        simple_fs : (fs_meta.fs_need =='none'  ||   fs_meta.omic in fs_meta.fs_omic_need )
        fs_dependency_MET : 'mixMET' in  fs_meta.fs_need || 'MET' in  fs_meta.fs_need
        fs_dependency_RNA : true 
    }.set { fs_branch }

    // fs_branch.simple_fs.view{v-> v[0].fs_omic_need}
    // fs_branch.fs_dependency_MET.count()

    // out_pp_create_filtered.view()

    
    fs_RNA = fs_branch.fs_dependency_RNA
    .combine(out_pp_create_filtered)
    .filter{ fs_meta, a,b,c,d,e,dataset_file,ref_file, pp_meta,pp_file ->
        pp_meta.create in fs_meta.need   &&    fs_meta.omic !in fs_meta.fs_omic_need   
    }
    .map{fs_meta, a,b,c,d,e,dataset_file,ref_file, pp_meta,pp_file  ->
        if( "mixRNA" in  fs_meta.omic_need   ){
            tuple(fs_meta, a,b,c,d,e,pp_file,ref_file )
        }else { //place the ref in ref. 
            tuple(fs_meta, a,b,c,d,e,dataset_file,pp_file )
        }
    }

    fs_MET = fs_branch.fs_dependency_MET
    .combine(out_pp_create_filtered)
    .filter{ fs_meta, a,b,c,d,e,dataset_file,ref_file, pp_meta,pp_file ->
        pp_meta.create in fs_meta.need   &&    fs_meta.omic !in fs_meta.fs_omic_need   
    }
    .map{fs_meta, a,b,c,d,e,dataset_file,ref_file, pp_meta,pp_file  ->
        if( "mixMET" in  fs_meta.omic_need   ){
            tuple(fs_meta, a,b,c,d,e,pp_file,ref_file )
        }else { //place the pp_file in ref. 
            tuple(fs_meta, a,b,c,d,e,dataset_file,pp_file )
        }
    }


    // fs_RNA.view()
    
    // fs_RNA.count()



    // fs_branch.fs_dependency_RNA.view{fs_meta, a,b,c,d,e,dataset_file,ref_file, pp_meta,pp_file -> 
    // println(fs_meta.omic)
    // println(fs_meta.fs_omic_need)
    // println(fs_meta.omic !in fs_meta.fs_omic_need   )
    
    // }
    // .filter { fs_meta, a,b,c,d,e, dataset_meta, dataset_file ->
    //     fs_meta.dataset == dataset_meta.id 
    // }
    // .combine(out_cleaned_ref)
    // .filter {fs_meta, a,b,c,d,e, dataset_meta, dataset_file, ref_meta, ref_file ->
    //     fs_meta.ref == ref_meta.id 
    // }
    // .map {fs_meta, a,b,c,d,e, dataset_meta, dataset_file, ref_meta, ref_file ->
    //     def dup_fs_meta = fs_meta.clone()
    //     dup_fs_meta.output ="out-fs-"+ [dup_fs_meta.omic,dup_fs_meta.dataset, dup_fs_meta.ref,dup_fs_meta.pp_fun, dup_fs_meta.fs_fun ].join('_')+'.h5'
    //     tuple(dup_fs_meta, a,b,c,d,e,dataset_file,ref_file )
    // }


    // complete_fs_files.view{ v-> v[0].output } 


    out_fs =  fs_branch.simple_fs.mix(fs_MET , fs_RNA) | Features_selection
    // out_fs = complete_fs_files | Features_selection



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


    out_de_rna_unit = de_rna_unit | Prediction_deconvolution_rna 

    
// ################## Generate combinaison for the MET unit 

    fs_mixMET = out_fs.filter{meta,a -> meta.omic=='mixMET'  }
    fs_MET =out_fs.filter{meta,a -> meta.omic=='MET'  }

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


    out_de_met_unit= de_met_unit | Prediction_deconvolution_met

// ################## Generate combinaison for late integration

    li_path =  []
    CONFIG.late_integration.each {li,liv -> li_path.add([ [li_fun:li ] , file(liv.path)])}

    li_channel = Channel.fromList(li_path)

    li_combinaison = 
    li_channel.combine(out_de_rna_unit).combine(out_de_met_unit)
    .combine(out_mix).combine(out_cleaned_ref)
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
            [dup_li_meta.rna_unit.scRNA.pp_fun, dup_li_meta.rna_unit.scRNA.fs_fun,dup_li_meta.rna_unit.de_fun ].join('_')      +
            [dup_li_meta.met_unit.mixMET.pp_fun, dup_li_meta.met_unit.mixMET.fs_fun ].join('_')   +
            [dup_li_meta.met_unit.MET.pp_fun, dup_li_meta.met_unit.MET.fs_fun ,dup_li_meta.met_unit.de_fun].join('_')           +'.h5'
        dup_li_meta["output"] = output_name
        tuple( dup_li_meta,li_path , file_rna , file_met, dataset_file,ref_file )
    }.combine(Channel.of(tuple(file(params.wrapper.script_05),file(params.utils))))


    out_li = li_combinaison | Late_integration 


// ################## Generate score 

    score_input = out_li.map{   meta,file_path  -> 

        def dup_meta = meta.clone()
        def output_name = "score-" + [dup_meta.dataset,dup_meta.ref].join('_')  + '_' +
            ['mixRNA',  dup_meta.rna_unit.mixRNA.pp_fun, dup_meta.rna_unit.mixRNA.fs_fun ].join('_')   +   '_' +  //dup_meta.rna_unit.mixRNA.omic,
            ['RNA', dup_meta.rna_unit.RNA.pp_fun, dup_meta.rna_unit.RNA.fs_fun ].join('_')           +   '_' +  //dup_meta.rna_unit.RNA.omic,
            ['scRNA', dup_meta.rna_unit.scRNA.pp_fun, dup_meta.rna_unit.scRNA.fs_fun,dup_meta.rna_unit.de_fun ].join('_')   + '_'    +    //dup_meta.rna_unit.scRNA.omic,
            ['mixMET',  dup_meta.met_unit.mixMET.pp_fun, dup_meta.met_unit.mixMET.fs_fun ].join('_') + '_'  +    //dup_meta.rna_unit.mixMET.omic,
            ['MET', dup_meta.met_unit.MET.pp_fun, dup_meta.met_unit.MET.fs_fun ,dup_meta.met_unit.de_fun].join('_')  + '_'  +       //dup_meta.rna_unit.MET.omic,
            dup_meta.li_fun    +'.h5'
        dup_meta["output"] = output_name
        tuple(dup_meta, file_path, file(CONFIG.datasets[dup_meta.dataset].groundtruth_file_path),file(params.wrapper.script_06),file(params.utils))//,file(params.config_files.datasets))
     }   

    score_out = score_input | Scoring


// ################## Generate Metaanalysis

    // input_meta = score_out.collect(flat: false)

    score_out.map{ v-> 
    v[0] }.set{l_meta}

    score_out.map{ v-> 
    v[1] }.set{l_path}

    input_meta = l_meta.collect(flat: false).concat(l_path.collect(flat: false)).collect(flat: false)
    .combine(Channel.of(tuple(file(params.wrapper.script_07),file(params.utils))))


    // input_meta.count().view()
    // input_meta.view()

    // input.count().view()
    // input_meta.view { v-> v.size()  }
    // input_meta.view { v-> v[0].size() +' ' +v[1].size()   }
    // input_meta.view { v-> v[0] +' \n\n\n' +v[1][0]   }

    // input_meta.view { v-> v[0][1] +' ' +v[1][1]   }

    input_meta | Metaanalysis 

}

process Cleaning_mix {
    input:
    tuple val(meta),  
    path(mixes), 
    path(cleaner), 
    path(wrapper01),
    path(utils)

    output:
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
}


process Cleaning_ref{
    input:
    tuple val(meta), 
    path(reference),
    path(cleaner) ,
    path(wrapper01),
    path(utils)
    

    output:
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

process Preprocessing {
    cpus 1

    input:
    tuple val(meta),
        path(pp_script), 
        path(dependency),
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

process Features_selection {
    cpus 1
    
    input:
        tuple( val(meta),
            path(fs_script), 
            path(file_input),
            path(wrapper03),
            path(utils),
            path(files_dep),
            path(mix), 
            path(reference),
            )


    output:
    tuple val(meta), path("${meta.output}")

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

process Prediction_deconvolution_rna {
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

process Prediction_deconvolution_met {
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

process Late_integration {
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


    script:
    """
    RCODE="input_file_rna='${input_file_rna}'; input_file_met='${input_file_met}'; 
    path_ogmix='${mix}' ; path_ogref='${reference}' ; 
    output_file='${meta.output}'; 
    script_file='${script_li}'; 
    utils_script='${utils}'; 
    source('${wrapper04}');"
    echo \$RCODE | Rscript - 
    """

    stub:
    """
    RCODE="input_file_rna='${input_file_rna}'; input_file_met='${input_file_met}'; 
    path_ogmix='${mix}' ; path_ogref='${reference}' ; 
    output_file='${meta.output}'; 
    script_file='${script_li}'; 
    utils_script='${utils}'; 
    source('${wrapper04}');"
    echo \$RCODE
    touch ${meta.output}
    """
}

process Scoring {
    cpus 1
    
    input:   
    tuple  val(meta) ,
        path(prediction),
        path(groundtruth_file),
        path(scoring_script),
        path(utils)

    output:
    tuple val(meta), path("${meta.output}")


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


process Metaanalysis {
    publishDir '.' 
    cpus 1
    
    input:   
    tuple( 
        val(meta), 
        path(input_score), 
        path(meta_script), 
        path(utils), 
        // path(file_dataset)
    )
    // tuple( tuple(val(meta) ,path(input_score)),
    // path(meta_script),
    // path(utils))

    output:
    path("07_metaanalysis.html")


    script:
    """
    RCODE="
    score_files = strsplit(trimws('${input_score}'),' ') ; 
    utils_script ='${utils}';
    rmarkdown::render('${meta_script}');"
    echo \$RCODE | Rscript -
    """  
    // file_dataset = '${file_dataset}';

    stub:
    """
    RCODE="
    score_files = strsplit(trimws('${input_score}'),' ') ; 
    utils_script ='${utils}';
    rmarkdown::render('${meta_script}');"
    echo \$RCODE 
    touch 07_metaanalysis.html
    """

}
