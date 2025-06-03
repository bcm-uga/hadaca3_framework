
import subprocess
import os
import psutil
import time
from datetime import datetime
import shutil

# nb_cores doit être modifié dans nextflow.config. 

# nextflow_cmd_stub = "nextflow run 00_run_pipeline.nf -stub -with-report --setup_folder "
# nextflow_cmd      = "nextflow run 00_run_pipeline.nf -with-report --setup_folder "

nextflow_cmd_stub = "nextflow run 00_run_pipeline.nf -stub --setup_folder "
nextflow_cmd      = "nextflow run 00_run_pipeline.nf --setup_folder "


smk_cmd_dry     = "snakemake --cores 4 --forceall -s 00_run_pipeline.smk -n --config setup_folder="
smk_dep_cmd_dry     = "snakemake --cores 4 --forceall -s 00_run_pipeline.smk -n --config setup_folder="
smk_cmd         = "snakemake --cores 4 --forceall -s 00_run_pipeline.smk --config setup_folder="
smk_cmd_clean   = "snakemake --cores 4 -s 00_run_pipeline.smk -p clean"

d_cmd = {"nextflow_stub":nextflow_cmd_stub,"snakemake_dry":smk_cmd_dry ,"nextflow":nextflow_cmd, "snakemake":smk_cmd,"snakemake_dep_dry":smk_dep_cmd_dry}
# d_cmd = { "snakemake":smk_cmd }
# d_cmd = {"nextflow_stub":nextflow_cmd_stub,"nextflow":nextflow_cmd,}
# d_cmd = {"snakemake_dry":smk_cmd_dry , "snakemake":smk_cmd}


path_setup = 'benchmark/setup/'
conda_env = "hadaca3framework_env"
meta_file_path = "07_metaanalysis.html"


conda_activate = "conda run -n "+ conda_env 


# setup_nb= range(1,3)
# setup_nb= range(2,3)
# setup_nb= range(5,11)
setup_nb= range(1,11)


nb_replica = 3


bench_path = "benchmark/results/"


##### change directory here
os.chdir("..")
print(os.getcwd())

if os.path.exists(meta_file_path):
    os.remove(meta_file_path)


def run_process(command, bench_path, process_name):
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"Current time: {current_time}")

    start_time = time.time()
    process = subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)
    # process = subprocess.Popen(command)
    parent = psutil.Process(process.pid)

    memory_usage = 0
    # with open(file_path_stdout, "w",buffering=1) as f_stdout, open(file_path_err, "w",buffering=1) as f_err:
    while process.poll() is None:
        try:
            # Get all children + the parent itself
            children = parent.children(recursive=True)
            all_processes = [parent] + children
            total_mem = sum(p.memory_info().rss for p in all_processes if p.is_running())
            memory_usage = max(memory_usage, total_mem / (1024 * 1024))  # in MB


        except psutil.NoSuchProcess:    
            break

    
    stdout, stderr = process.communicate()


    end_time = time.time()
    time.sleep(1)
    return (start_time, end_time, memory_usage)

file_path_res = bench_path + "data.txt"


f_res =  open(file_path_res, "w")
d_result = {}


def delete_work_folder():
    work_folder_path ="work"
    if os.path.exists(work_folder_path):
        shutil.rmtree(work_folder_path)
    nextflow_folder_path =".nextflow/"
    if os.path.exists(nextflow_folder_path):
        shutil.rmtree(nextflow_folder_path)


def check_and_delete_metaanalysis():
    if os.path.exists(meta_file_path):
        # Delete the file
        os.remove(meta_file_path)
        print(f"Metaanalysis has been deleted successfully, therefor the pipeline was a success.")
    else:
        raise ValueError(f"Metaanalysis does not exist,therefor the pipeline did not completed without erros.")
        


# for work_flow in [nextflow_cmd_stub,smk_cmd_dry]:
for i in setup_nb :
    for w_name,work_flow in    d_cmd.items():
        for rep in range(nb_replica):
            process_name = w_name+str(i)
            cmd = conda_activate +  ' '+ work_flow+path_setup+str(i)+'/'
            # print(cmd)
            print("replica :" + str(rep+1) +' '+ w_name+" setup : "+str(i))

            # print(cmd)
            start_time, end_time, memory_usage = run_process(cmd.split(' '), bench_path, process_name)
            
            if ('dry' not in w_name ):
                check_and_delete_metaanalysis()
            
            if("nextflow" in w_name):
                delete_work_folder()

            d_result[process_name] = (end_time - start_time,memory_usage)
            f_res.write( f"{process_name} : ({end_time - start_time},{memory_usage})\n")
            
            
            print(f"Time elapsed: {end_time - start_time:.4f} seconds")
            print(f"Memory usage: {memory_usage:.4f} MB\n")
            print('\n')

        f_res.write( f"\n")
        print('\n')

