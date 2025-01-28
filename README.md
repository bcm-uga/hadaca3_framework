# hadaca3_framework

A framework to collectively develop multi-omic deconvolution methods.

The framework contains several blocks

- **pre-processing** : takes as input a multimodal .txt.gz, outputs a multimodal txt.gz


- **feature_selection** : takes as input a multimodal .txt.gz, outputs a multimodal txt.gz

- **deconvolution** : takes as input a multimodal mix.txt.gz and a ref.txt.gz, output a prediction.txt.gz

- **early_int** : takes as input a multimodal mix.txt.gz and a ref.txt.gz, outputs an integrated mix.txt.gz and a ref.txt.gz

- **late_int** : takes as input a a multimodal prediction.txt.gz, outputs an integrated prediction

- **intermediate_int** : takes as input a multimodal mix.txt.gz and a ref.txt.gz, outputs an integrated prediction

To run the demo, execute 
```
mkdir 00_demo_data
#then copy in it starting_kit_phase2-3/data
```
