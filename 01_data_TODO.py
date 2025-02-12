# import rpy2.robjects
import gzip
import sys

# readRDS = rpy2.robjects.r['readRDS']


def convert_to_gz(data_rds,output_gz):
    # lines = readRDS(data_rds)
    lines = "Placeholders"
    with gzip.open(output_gz, 'wb') as f:
        f.writelines(map(str.encode,lines))


dataset_name_rds = sys.argv[1]
dataset_name_gz = sys.argv[2]


convert_to_gz(dataset_name_rds ,dataset_name_gz)




