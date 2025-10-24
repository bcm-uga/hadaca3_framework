
# import utils.data_processing as dp 
# import data_processing as dp 
import importlib.util
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description="Early_python_Python version")
parser.add_argument("--input_mix_rna", required=True)
parser.add_argument("--input_rna", required=True)
parser.add_argument("--input_scrna", required=True)
parser.add_argument("--input_mix_met", required=True)
parser.add_argument("--input_met", required=True)
parser.add_argument("--path_ogmix", required=True)
parser.add_argument("--path_ogref", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--script_file", required=True)
parser.add_argument("--utils", required=True)
args = parser.parse_args()



# print("=== Running Prediction Deconvolution (Python) ===")
# print(f"Input Mix RNA: {args.input_mix_rna}")
# print(f"Input RNA: {args.input_rna}")
# print(f"Input scRNA: {args.input_scrna}")
# print(f"Input Mix Met: {args.input_mix_met}")
# print(f"Input Met: {args.input_met}")
# print(f"Mix Path: {args.path_ogmix}")
# print(f"Reference Path: {args.path_ogref}")
# print(f"Output File: {args.output}")
# print(f"script_file: {args.script_file}")
# print(f"Utils: {args.utils}")


path_og_dataset = ""



def load_module(module_symlink):
    path = Path(module_symlink).resolve()
    module_name =  path.stem 
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

dp =  load_module(args.utils)


mix_met = dp.read_hdf5(args.input_mix_met)['mix_met']
ref_met = dp.read_hdf5(args.input_met)['ref_met']

met_unit = {"mix": mix_met,"ref" :ref_met}#,ref_scRNA=  scRNA  )

mix_rna = dp.read_hdf5(args.input_mix_rna)['mix_rna']
ref_rna = dp.read_hdf5(args.input_rna)['ref_bulkRNA']
scRNA= dp.read_hdf5(args.input_scrna)['ref_scRNA']

rna_unit = {"mix" :mix_rna,"ref" :ref_rna,"ref_scRNA" :  scRNA }


path_og_dataset= {"mix" :args.path_ogmix,"ref" : args.path_ogref }



program_script = load_module(args.script_file)




# print(rna_unit)
# print(met_unit)
integrated_unit = program_script.program_block_EI(rna_unit,met_unit,path_og_dataset)

# print(integrated_unit)

if(len(integrated_unit) == 0):
    raise Exception("Early integration is empty, script_file = ",args.script_file )
#     (paste("Early integration is empty, script_file = ",args.script_file ) )


dp.write_global_hdf5(args.output,integrated_unit)
