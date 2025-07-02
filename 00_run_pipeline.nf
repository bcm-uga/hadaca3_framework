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
    script_04_early : '04_Early_integration_pipB.R',
    script_05_early_deco : '05_early_decovolution.R',
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


    // ################## Computing cleaned mix and ref clean and generating a tuple with key as the dataset or ref name and output file. 

    def utils_channel =  Channel.fromPath(params.utils)

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
    

    // Separate into branch if the og_path should be changed to a pp_output. 
    complete_fs_files.branch{fs_meta, a,b,c,d,e,dataset_file,ref_file  -> 
        simple_fs : (fs_meta.fs_need =='none'  ||   fs_meta.omic in fs_meta.fs_omic_need )
        fs_dependency_MET : 'mixMET' in  fs_meta.fs_need || 'MET' in  fs_meta.fs_need
        fs_dependency_RNA : true 
    }.set { fs_branch }

    
    fs_RNA = fs_branch.fs_dependency_RNA
    .combine(out_pp_create_filtered)
    .filter{ fs_meta, a,b,c,d,e,dataset_file,ref_file, pp_meta,pp_file ->
        pp_meta.create in fs_meta.need   &&    fs_meta.omic !in fs_meta.fs_omic_need   
    }
    .map{fs_meta, a,b,c,d,e,dataset_file,ref_file, pp_meta,pp_file  ->
        if( "mixRNA" in  fs_meta.omic_need   ){
            tuple(fs_meta, a,b,c,d,e,pp_file,ref_file )
        }else { //place the ref in ref of og_path.... 
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

    fs_mixRNA = out_fs.filter{meta, file -> meta.omic=='mixRNA'  }
    fs_RNA =out_fs.filter{meta, file -> meta.omic=='RNA'}
    fs_scRNA = out_fs.filter{meta, file -> meta.omic=='scRNA'}

    rna_unit = 
    fs_mixRNA.combine(fs_RNA).combine(fs_scRNA).combine(out_mix).combine(out_cleaned_ref)
    // filter on the dataset and ref (same dataset/ref for pp, fs and all omics)
    .filter{ meta_mix,file_input_mix, meta_RNA,file_input_RNA,meta_scRNA,file_input_scRNA,   dataset_meta, dataset_file, ref_meta, ref_file ->
        meta_scRNA.ref == ref_meta.id && meta_RNA.ref == ref_meta.id && meta_mix.ref == ref_meta.id && meta_mix.dataset == dataset_meta.id
    }
    //filter on same function for mix and bulk 
    .filter{ meta_mix,file_input_mix, meta_RNA,file_input_RNA,meta_scRNA,file_input_scRNA,   dataset_meta, dataset_file, ref_meta, ref_file ->
        // println(meta_mix.pp_fun + meta_mix.fs_fun)
        meta_mix.pp_fun == meta_RNA.pp_fun && meta_mix.fs_fun == meta_RNA.fs_fun
    }


    de_rna_unit = 
    de_channel.combine(rna_unit)
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
    


    out_de_rna_unit = 
    de_rna_unit.combine(Channel.of(tuple(file(params.wrapper.script_04_rna),file(params.utils)))) 
    | Prediction_deconvolution_rna 

    
// ################## Generate combinaison for the MET unit 

    fs_mixMET = out_fs.filter{meta,a -> meta.omic=='mixMET'  }
    fs_MET =out_fs.filter{meta,a -> meta.omic=='MET'  }

    deco_path_met =  []
    CONFIG.deconvolution.each {de,dev -> deco_path_met.add( [ [de_fun : de]  , file(dev.path)])}

    de_channel_met = Channel.fromList(deco_path_met)

    met_unit=  
    fs_mixMET.combine(fs_MET).combine(out_mix).combine(out_cleaned_ref)    
    .filter{meta_mix, file_input_mix, meta_MET,file_input_MET,dataset_meta, dataset_file, ref_meta, ref_file ->
        meta_MET.ref == ref_meta.id && meta_mix.ref == ref_meta.id && meta_mix.dataset == dataset_meta.id
    } 
    // filter mix and bulk to have the same functions : 
    .filter{meta_mix, file_input_mix, meta_MET,file_input_MET,dataset_meta, dataset_file, ref_meta, ref_file ->
        meta_mix.pp_fun == meta_MET.pp_fun && meta_mix.fs_fun == meta_MET.fs_fun
    }


    de_met_unit = 
    de_channel_met.combine(met_unit)
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
    


    out_de_met_unit= 
    de_met_unit.combine(Channel.of(tuple(file(params.wrapper.script_04_met),file(params.utils))))
     | Prediction_deconvolution_met

// ################## Generate combinaison for late integration

    //add none prediction for rna/met only. 



    li_path =  []
    CONFIG.late_integration.each {li,liv -> li_path.add([ [li_fun:li ] , file(liv.path)])}

    li_channel = Channel.fromList(li_path)

    none_unit = Channel.of( tuple( [dataset: "nodt", ref:'nore',
                            mixRNA:[pp_fun:"nopp",fs_fun:"nofs"],
                            RNA:[pp_fun:"nopp",fs_fun:"nofs"],
                            scRNA:[pp_fun:"nopp",fs_fun:"nofs"],
                            mixMET:[pp_fun:"nopp",fs_fun:"nofs"],
                            MET:[pp_fun:"nopp",fs_fun:"nofs"],
                            de_fun:"node"], 
                            file("none_unit") 
                            ))
    


    // out_de_met_unit.count().view()
    // out_de_rna_unit.count().view()
    
    li_combinaison = li_channel
    // .combine(out_de_rna_unit)  // inject none_unit into RNA
    // .combine(out_de_met_unit) 
    .combine(out_de_rna_unit.mix(none_unit))  // inject none_unit into RNA
    .combine(out_de_met_unit.mix(none_unit))  // inject none_unit into MET
    .combine(out_mix)
    .combine(out_cleaned_ref)
    // .filter{m, li_file, meta_rna, file_rna,    meta_met, file_met,           dataset_meta, dataset_file,   ref_meta, ref_file-> 
    //     meta_rna.dataset == dataset_meta.id &&         meta_met.dataset == dataset_meta.id && 
    //     meta_rna.ref == ref_meta.id         &&         meta_met.ref ==ref_meta.id
    // }
    .filter { m, li_file, meta_rna, file_rna, meta_met, file_met, dataset_meta, dataset_file, ref_meta, ref_file ->

        // Keep all valid matches OR explicitly allow 'none_unit'
        def isValidMatch = (
            meta_rna.dataset == dataset_meta.id &&
            meta_met.dataset == dataset_meta.id &&
            meta_rna.ref == ref_meta.id &&
            meta_met.ref == ref_meta.id && 
            m.li_fun !='OnlyRna' && 
            m.li_fun !='OnlyMet'
        )

        def isOnlyRna = (
            m.li_fun == 'OnlyRna' && 
            meta_rna.dataset == dataset_meta.id &&
            meta_rna.ref == ref_meta.id &&
            meta_rna.dataset != "nodt"  &&
            meta_met.dataset == "nodt"
        ) 

        def isOnlyMet = (
            m.li_fun =='OnlyMet' && 
            meta_met.dataset == dataset_meta.id &&
            meta_met.ref == ref_meta.id && 
            meta_met.dataset != "nodt" &&
            meta_rna.dataset == "nodt"
        )
        return (isValidMatch || isOnlyRna  || isOnlyMet)
    }
    .map{
        li_meta,li_file,meta_rna, file_rna , meta_met, file_met ,dataset_meta, dataset_file, ref_meta, ref_file ->
        def dup_li_meta = li_meta.clone()
        if(meta_rna.dataset !="nodt"){
            dup_li_meta['dataset'] = meta_rna.dataset
            dup_li_meta['ref'] = meta_rna.ref
        }else{
            dup_li_meta['dataset'] = meta_met.dataset
            dup_li_meta['ref'] = meta_met.ref
        }
        dup_li_meta['rna_unit'] = meta_rna
        dup_li_meta['met_unit'] = meta_met

        def output_name = "pred-li-" + [dup_li_meta.dataset,dup_li_meta.ref].join('_')                                      + '_' +
            [dup_li_meta.rna_unit.mixRNA.pp_fun, dup_li_meta.rna_unit.mixRNA.fs_fun ].join('_')                            + '_' +
            [dup_li_meta.rna_unit.RNA.pp_fun, dup_li_meta.rna_unit.RNA.fs_fun ].join('_')                                  + '_' +
            [dup_li_meta.rna_unit.scRNA.pp_fun, dup_li_meta.rna_unit.scRNA.fs_fun,dup_li_meta.rna_unit.de_fun ].join('_')  + '_' +
            [dup_li_meta.met_unit.mixMET.pp_fun, dup_li_meta.met_unit.mixMET.fs_fun ].join('_')                            + '_' +
            [dup_li_meta.met_unit.MET.pp_fun, dup_li_meta.met_unit.MET.fs_fun ,dup_li_meta.met_unit.de_fun,dup_li_meta.li_fun].join('_')          +'.h5'
        dup_li_meta["output"] = output_name
        tuple( dup_li_meta,li_file , file_rna , file_met, dataset_file,ref_file )
    }.combine(Channel.of(tuple(file(params.wrapper.script_05),file(params.utils))))

    // li_combinaison.count().view()
    // li_combinaison.count()


    out_li = li_combinaison | Late_integration 
    // out_li = Channel.of()

    // out_li.view()
// ################################  Early integration and deconvolution path


    ei_path =  []
    CONFIG.early_integration.each {ei,eiv -> ei_path.add([ [ei_fun:ei ] , file(eiv.path)])}
    

    ei_channel = Channel.fromList(ei_path)

    ei_combinaison = 
    ei_channel.combine(rna_unit).combine(met_unit)
    .filter{meta_ei, ei_file,     
            meta_mixRNA, file_input_mixRNA,    meta_RNA,file_input_RNA,    meta_scRNA,file_input_scRNA,    dataset_meta_rna, dataset_file_rna,    ref_meta_rna, ref_file_rna,   
            meta_mixMET, file_input_mixMET, meta_MET,file_input_MET,  dataset_meta_MET, dataset_file_MET, ref_meta_MET, ref_file_MET -> 
        
        dataset_meta_rna.id == dataset_meta_MET.id && ref_meta_rna.id == ref_meta_MET.id
    }.map{ meta_ei, ei_file,     
           meta_mixRNA, file_input_mixRNA,    meta_RNA,file_input_RNA,    meta_scRNA,file_input_scRNA,    dataset_meta_rna, dataset_file_rna,    ref_meta_rna, ref_file_rna,   
           meta_mixMET, file_input_mixMET, meta_MET,file_input_MET,  dataset_meta_MET, dataset_file_MET, ref_meta_MET, ref_file_MET -> 
        def dup_meta_ei = meta_ei.clone() 
        if(dataset_meta_rna.id!= "nodt"){
            dup_meta_ei['dataset'] = dataset_meta_rna.id
            dup_meta_ei['ref'] = ref_meta_rna.id
        }else{
            dup_meta_ei['dataset'] = dataset_meta_MET.id
            dup_meta_ei['ref'] = ref_meta_MET.id
        }

        dup_meta_ei['mixRNA'] = meta_mixRNA
        dup_meta_ei['RNA'] = meta_RNA
        dup_meta_ei['scRNA'] = meta_scRNA
        dup_meta_ei['mixMET'] = meta_mixMET
        dup_meta_ei['MET'] = meta_MET
        // dup_meta_ei['rna_unit'] = meta_rna
        // dup_meta_ei['met_unit'] = meta_met

        def output_name = "out-li-" + [dup_meta_ei.dataset,dup_meta_ei.ref].join('_')   + '_' + 
            [dup_meta_ei.mixRNA.pp_fun, dup_meta_ei.mixRNA.fs_fun ].join('_')           + '_' +
            [dup_meta_ei.RNA.pp_fun, dup_meta_ei.RNA.fs_fun ].join('_')                 + '_' +
            [dup_meta_ei.scRNA.pp_fun, dup_meta_ei.scRNA.fs_fun].join('_')              + '_' +
            [dup_meta_ei.mixMET.pp_fun, dup_meta_ei.mixMET.fs_fun ].join('_')           + '_' +
            [dup_meta_ei.MET.pp_fun, dup_meta_ei.MET.fs_fun , dup_meta_ei.ei_fun ].join('_')           +'.h5'
        dup_meta_ei["output"] = output_name

        tuple(
           dup_meta_ei, ei_file,     
        file_input_mixRNA,    file_input_RNA,    file_input_scRNA,      
           file_input_mixMET, file_input_MET,  dataset_file_MET, ref_file_MET
        )
    }
    .combine(Channel.of(tuple(file(params.wrapper.script_04_early),file(params.utils))))

    ei_out = ei_combinaison | early_integration

// ################## Deconvolution of the early integration. 


    ei_decon_path =  []
    CONFIG.deconvolution.each {de,dev -> ei_decon_path.add( [ [de_fun : de]  , file(dev.path)])}
    
    de_channel_EI = 
    Channel.fromList(ei_decon_path)
    .combine(ei_out)
    .map{de_meta , de_file, ei_meta, ei_file -> 
        def dup_de_meta = ei_meta.clone()
        dup_de_meta["de_fun"] = de_meta.de_fun
        def output_name = "out-li-" + [dup_de_meta.dataset,dup_de_meta.ref].join('_')   + '_' + 
            [dup_de_meta.mixRNA.pp_fun, dup_de_meta.mixRNA.fs_fun ].join('_')           + '_' +
            [dup_de_meta.RNA.pp_fun, dup_de_meta.RNA.fs_fun ].join('_')                 + '_' +
            [dup_de_meta.scRNA.pp_fun, dup_de_meta.scRNA.fs_fun].join('_')              + '_' +
            [dup_de_meta.mixMET.pp_fun, dup_de_meta.mixMET.fs_fun ].join('_')           + '_' +
            [dup_de_meta.MET.pp_fun, dup_de_meta.MET.fs_fun , dup_de_meta.ei_fun, dup_de_meta.de_fun ].join('_') +'.h5'

        dup_de_meta["output"] = output_name
        tuple(dup_de_meta,de_file,ei_file)
    }
    .combine(Channel.of(tuple(file(params.wrapper.script_05_early_deco),file(params.utils))))

    out_ei_de = de_channel_EI | decovolution_EI

    // ei_de_out.view()

// ################## Generate score 

    score_input_li = out_li.map{   meta,file_path  -> 
        def dup_meta = meta.clone()
        def output_name = "score-li-" + [dup_meta.dataset,dup_meta.ref].join('_')  + '_' +
            ['mixRNA',  dup_meta.rna_unit.mixRNA.pp_fun, dup_meta.rna_unit.mixRNA.fs_fun ].join('_')   +   '_' +  //dup_meta.rna_unit.mixRNA.omic,
            ['RNA', dup_meta.rna_unit.RNA.pp_fun, dup_meta.rna_unit.RNA.fs_fun ].join('_')           +   '_' +  //dup_meta.rna_unit.RNA.omic,
            ['scRNA', dup_meta.rna_unit.scRNA.pp_fun, dup_meta.rna_unit.scRNA.fs_fun,dup_meta.rna_unit.de_fun ].join('_')   + '_'    +    //dup_meta.rna_unit.scRNA.omic,
            ['mixMET',  dup_meta.met_unit.mixMET.pp_fun, dup_meta.met_unit.mixMET.fs_fun ].join('_') + '_'  +    //dup_meta.rna_unit.mixMET.omic,
            ['MET', dup_meta.met_unit.MET.pp_fun, dup_meta.met_unit.MET.fs_fun ,dup_meta.met_unit.de_fun].join('_')  + '_'  +       //dup_meta.rna_unit.MET.omic,
            dup_meta.li_fun    +'.h5'
        dup_meta["output"] = output_name
        tuple(dup_meta, file_path, file(CONFIG.datasets[dup_meta.dataset].groundtruth_file_path),file(params.wrapper.script_06),file(params.utils))//,file(params.config_files.datasets))
     }   

    score_input_ei = out_ei_de.map{ ei_meta, ei_file ->
        def dup_ei_de_meta = ei_meta.clone()
        def output_name = "score-ei-" + [dup_ei_de_meta.dataset,dup_ei_de_meta.ref].join('_') + '_' +
            ['mixRNA',dup_ei_de_meta.mixRNA.pp_fun, dup_ei_de_meta.mixRNA.fs_fun ].join('_')  + '_' +
            ['RNA',dup_ei_de_meta.RNA.pp_fun, dup_ei_de_meta.RNA.fs_fun ].join('_')           + '_' +
            ['scRNA',dup_ei_de_meta.scRNA.pp_fun, dup_ei_de_meta.scRNA.fs_fun].join('_')      + '_' +
            ['mixMET',dup_ei_de_meta.mixMET.pp_fun, dup_ei_de_meta.mixMET.fs_fun ].join('_')  + '_' +
            ['MET',dup_ei_de_meta.MET.pp_fun, dup_ei_de_meta.MET.fs_fun , dup_ei_de_meta.ei_fun, dup_ei_de_meta.de_fun ].join('_') +'.h5'

        dup_ei_de_meta["output"] = output_name
        tuple(dup_ei_de_meta,ei_file,file(CONFIG.datasets[dup_ei_de_meta.dataset].groundtruth_file_path),file(params.wrapper.script_06),file(params.utils))
    }

    score_out = score_input_li.mix(score_input_ei) | Scoring


// ################## Generate Metaanalysis

    // input_meta = score_out.collect(flat: false)

    score_out.map{ v-> 
    v[0] }.set{l_meta}

    score_out.map{ v-> 
    v[1] }.set{l_path}

    score_input_li.map{ v ->
    v[1]}.set{ pred_files }

    input_meta = l_meta.collect(flat: false).concat(l_path.collect(flat: false).concat( pred_files.collect(flat:false) )).collect(flat: false)
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

process early_integration {
    cpus 1
    
    input:
    tuple val(meta),
    path(script_ei),
    path(file_input_mix_rna),
    path(file_input_rna),
    path(file_input_scrna),
    path(file_input_mix_met),
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
    input_file_mix_rna='${file_input_mix_rna}'; input_file_rna='${file_input_rna}'; input_file_sc='${file_input_scrna}';
    input_file_mix_met='${file_input_mix_met}'; input_file_met='${file_input_met}';
    path_ogmix='${mix}' ; path_ogref='${reference}' ;
    output_file='${meta.output}'; 
    script_file='${script_ei}'; 
    utils_script='${utils}'; 
    source('${wrapper04}');"
    echo \$RCODE | Rscript - 
    """

    stub:
    """
    RCODE="
    input_file_mix_rna='${file_input_mix_rna}'; input_file_rna='${file_input_rna}'; input_file_sc='${file_input_scrna}';
    input_file_mix_met='${file_input_mix_met}'; input_file_met='${file_input_met}';
    path_ogmix='${mix}' ; path_ogref='${reference}' ; 
    output_file='${meta.output}'; 
    script_file='${script_ei}'; 
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
    path(wrapper05),
    path(utils)

    output:
    tuple val(meta), path("${meta.output}")


    script:
    """
    RCODE="
    
    input_file_rna='${input_file_rna}'; input_file_met='${input_file_met}'; 
    path_ogmix='${mix}' ; path_ogref='${reference}' ; 
    output_file='${meta.output}'; 
    script_file='${script_li}'; 
    utils_script='${utils}'; 
    source('${wrapper05}');"
    echo \$RCODE | Rscript - 
    """

    stub:
    """
    RCODE="input_file_rna='${input_file_rna}'; input_file_met='${input_file_met}'; 
    path_ogmix='${mix}' ; path_ogref='${reference}' ; 
    output_file='${meta.output}'; 
    script_file='${script_li}'; 
    utils_script='${utils}'; 
    source('${wrapper05}');"
    echo \$RCODE
    touch ${meta.output}
    """
}

process decovolution_EI {
    cpus 1
     input:
        tuple val(meta),
        path(de_script), 
        path(uni_data),
        path(wrapper05_deco),
        path(utils)

    output:
    tuple val(meta), path("${meta.output}")

    script:
    """
        RCODE="
        path_uni_data='${uni_data}';
        output_file='${meta.output}';
        script_de='${de_script}'; 
        utils_script='${utils}'; 
        source('${wrapper05_deco}');"
        echo \$RCODE | Rscript - 
    """

    stub:
    """
        RCODE="
        path_uni_data='${uni_data}';
        output_file='${meta.output}';
        script_de='${de_script}'; 
        utils_script='${utils}'; 
        source('${wrapper05_deco}');"
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
        path(pred_files),
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
    pred_files = strsplit(trimws('${pred_files}'),' ');
    utils_script ='${utils}';
    rmarkdown::render('${meta_script}');"
    echo \$RCODE | Rscript -
    """  
    // file_dataset = '${file_dataset}';

    stub:
    """
    RCODE="
    score_files = strsplit(trimws('${input_score}'),' ') ;
    pred_files = strsplit(trimws('${pred_files}'),' ');
    utils_script ='${utils}';
    rmarkdown::render('${meta_script}');"
    echo \$RCODE 
    touch 07_metaanalysis.html
    """

}
