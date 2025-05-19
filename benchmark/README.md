# Benchmaaaark 


## how to run this benchmark 

```
cd benchmark
conda activate hadaca3framework_env
python benchmark.py
```

## Run the analysis: 

```
cd benchmark
conda activate hadaca3framework_env
Rscript -e "rmarkdown::render('analysis.Rmd')"
```

## setup details : 

### setup1 
```
job                             count
----------------------------  -------
cleaning_mix                        1
cleaning_ref                        1
features_selection                  5
late_integration                    4
metaanalysis                        1
prediction_deconvolution_met        2
prediction_deconvolution_rna        2
preprocessing                       5
scoring                             4
total                              25
```


### setup2 
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                  7
late_integration                    8
metaanalysis                        1
prediction_deconvolution_met        4
prediction_deconvolution_rna        4
preprocessing                       7
scoring                             8
total                              42
```

### setup3
```
job                             count
----------------------------  -------
cleaning_mix                        4
cleaning_ref                        1
features_selection                 11
late_integration                   16
metaanalysis                        1
prediction_deconvolution_met        8
prediction_deconvolution_rna        8
preprocessing                      11
scoring                            16
total                              76

```

### setup4 
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 16
late_integration                  128
metaanalysis                        1
prediction_deconvolution_met       64
prediction_deconvolution_rna        4
preprocessing                      10
scoring                           128
total                             354

```

### setup5 

```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 16
late_integration                  384
metaanalysis                        1
prediction_deconvolution_met       64
prediction_deconvolution_rna        4
preprocessing                      10
scoring                           384
total                             866

```

### setup6
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 16
late_integration                  512
metaanalysis                        1
prediction_deconvolution_met       64
prediction_deconvolution_rna        4
preprocessing                      10
scoring                           512
total                            1122
```
### setup7
(4min45 dag creation dry run )
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 26
late_integration                 4608
metaanalysis                        1
prediction_deconvolution_met       96
prediction_deconvolution_rna       96
preprocessing                      13
scoring                          4608
total                            9451
```


### setup8
(13min dag creation dry run)
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 33
late_integration                10368
metaanalysis                        1
prediction_deconvolution_met      216
prediction_deconvolution_rna       96
preprocessing                      13
scoring                         10368
total                           21098
```

### setup9
(43min dag creation)
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 33
late_integration                31104
metaanalysis                        1
prediction_deconvolution_met      216
prediction_deconvolution_rna       96
preprocessing                      13
scoring                         31104
total                           62570
```enchmark/
(85min dag creation)
```
job                             count
----------------------------  -------
cleaning_mix                        2
cleaning_ref                        1
features_selection                 33
late_integration                62208
metaanalysis                        1
prediction_deconvolution_met      216
prediction_deconvolution_rna       96
preprocessing                      13
scoring                         62208
total                          124778
```
