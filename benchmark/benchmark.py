
import subprocess
import os
import psutil
import time


# nb_cores doit être modifié dans nextflow.config. 

# nextflow_cmd_stub = "nextflow run 00_run_pipeline.nf -stub -with-report --setup_folder "
# nextflow_cmd      = "nextflow run 00_run_pipeline.nf -with-report --setup_folder "

nextflow_cmd_stub = "nextflow run 00_run_pipeline.nf -stub --setup_folder "
nextflow_cmd      = "nextflow run 00_run_pipeline.nf --setup_folder "


smk_cmd_dry     = "snakemake --cores 4 --forceall -s 00_run_pipeline.smk -pn --config setup_folder="
smk_cmd         = "snakemake --cores 4 --forceall -s 00_run_pipeline.smk -p --config setup_folder="
smk_cmd_clean   = "snakemake --cores 4 -s 00_run_pipeline.smk -p clean"

d_cmd = {"nextflow_stub":nextflow_cmd_stub,"snakemake_dry":smk_cmd_dry ,"nextflow":nextflow_cmd, "snakemake":smk_cmd}
# d_cmd = { "snakemake":smk_cmd }
# d_cmd = {"nextflow_stub":nextflow_cmd_stub,"nextflow":nextflow_cmd,}


path_setup = 'benchmark/setup'
conda_env = "hadaca3framework_env"

conda_activate = "conda run -n "+ conda_env 

setup_nb= range(1,3)
# setup_nb= range(2,3)



bench_path = "benchmark/results/"

file_path_res = bench_path+"data.txt"

# os.chdir("~/projects/hadaca3_framework/")
os.chdir("..")
print(os.getcwd())

def run_process(command):
    start_time = time.time()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    parent = psutil.Process(process.pid)

    memory_usage = 0
    while process.poll() is None:
        try:
            # Get all children + the parent itself
            children = parent.children(recursive=True)
            all_processes = [parent] + children
            total_mem = sum(p.memory_info().rss for p in all_processes if p.is_running())
            memory_usage = max(memory_usage, total_mem / (1024 * 1024))  # in MB
        except psutil.NoSuchProcess:
            break
        time.sleep(0.1)

    # # Poll process for memory usage while it's running
    # while process.poll() is None:
    #     try:
    #         memory_info = psutil.Process(process.pid).memory_info()
    #         memory_usage = max(memory_usage, memory_info.rss / (1024 * 1024))
    #     except psutil.NoSuchProcess:
    #         break
    #     time.sleep(0.1)  # Sleep briefly to avoid high CPU usage

    stdout, stderr = process.communicate()
    end_time = time.time()

    return (start_time, end_time, memory_usage, stdout, stderr)


# compute mean conda activate time 
l_time_act =[]
l_mem_act=[]
for i in range(10) :
    cmd = conda_activate + 'ls'
    start_time, end_time, memory_usage, stdout, stderr = run_process(cmd.split(' '))
    l_time_act.append(end_time - start_time)
    l_mem_act.append(memory_usage)


f_res =  open(file_path_res, "w")
f_res.write( f"conda : ({str(l_time_act)},{str(l_mem_act)})\n")
#     f.write(str(d_result))
d_result = {"conda":(l_time_act,l_mem_act)}

# for work_flow in [nextflow_cmd_stub,smk_cmd_dry]:
for w_name,work_flow in    d_cmd.items():
    for i in setup_nb :
        cmd = conda_activate +  ' '+ work_flow+path_setup+str(i)+'/'
        print(cmd)

        start_time, end_time, memory_usage, stdout, stderr = run_process(cmd.split(' '))
        
        process_name = w_name+str(i)
        file_path_stdout = bench_path + "stdout"+process_name+".txt"
        with open(file_path_stdout, "w") as f_stdout:
            f_stdout.write(stdout)

        file_path_err = bench_path+"stderr"+process_name+".txt"
        with open(file_path_err, "w") as f_err:
            f_err.write(stderr)
        
        d_result[process_name] = (end_time - start_time,memory_usage)
        f_res.write( f"{process_name} : ({end_time - start_time},{memory_usage})\n")
        
        
        print(f"Time elapsed: {end_time - start_time:.4f} seconds")
        print(f"Memory usage: {memory_usage:.4f} MB")

    f_res.write( f"{process_name} : ({end_time - start_time},{memory_usage})\n")




# with open(file_path_res, "w") as f:
#     f.write(str(d_result))