
import subprocess
import os
import psutil
import time
from datetime import datetime

# nb_cores doit être modifié dans nextflow.config. 

# nextflow_cmd_stub = "nextflow run 00_run_pipeline.nf -stub -with-report --setup_folder "
# nextflow_cmd      = "nextflow run 00_run_pipeline.nf -with-report --setup_folder "

nextflow_cmd_stub = "nextflow run 00_run_pipeline.nf -stub --setup_folder "
nextflow_cmd      = "nextflow run 00_run_pipeline.nf --setup_folder "


smk_cmd_dry     = "snakemake --cores 4 --forceall -s 00_run_pipeline.smk -n --config setup_folder="
smk_cmd         = "snakemake --cores 4 --forceall -s 00_run_pipeline.smk --config setup_folder="
smk_cmd_clean   = "snakemake --cores 4 -s 00_run_pipeline.smk -p clean"

d_cmd = {"nextflow_stub":nextflow_cmd_stub,"snakemake_dry":smk_cmd_dry ,"nextflow":nextflow_cmd, "snakemake":smk_cmd}
# d_cmd = { "snakemake":smk_cmd }
# d_cmd = {"nextflow_stub":nextflow_cmd_stub,"nextflow":nextflow_cmd,}
# d_cmd = {"snakemake_dry":smk_cmd_dry , "snakemake":smk_cmd}


path_setup = 'benchmark/setup/'
conda_env = "hadaca3framework_env"

conda_activate = "conda run -n "+ conda_env 

# setup_nb= range(4,7)
# setup_nb= range(1,3)
# setup_nb= range(2,3)
setup_nb= range(5,11)



bench_path = "benchmark/results/"

file_path_res = bench_path+"data.txt"

# os.chdir("~/projects/hadaca3_framework/")
os.chdir("..")
print(os.getcwd())

def run_process(command, bench_path, process_name):
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"Current time: {current_time}")

    start_time = time.time()
    # process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)
    process = subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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

# def run_process(command, bench_path, process_name):
#     current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#     print(f"Current time: {current_time}")

#     start_time = time.time()
#     process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)
#     parent = psutil.Process(process.pid)

#     # file_path_stdout = bench_path  + process_name + "-stdout.txt"
#     # file_path_err = bench_path  + process_name + "-stderr.txt"

#     memory_usage = 0

#     # with open(file_path_stdout, "w",buffering=1) as f_stdout, open(file_path_err, "w",buffering=1) as f_err:
#     while process.poll() is None:
#         try:
#             # Get all children + the parent itself
#             children = parent.children(recursive=True)
#             all_processes = [parent] + children
#             total_mem = sum(p.memory_info().rss for p in all_processes if p.is_running())
#             memory_usage = max(memory_usage, total_mem / (1024 * 1024))  # in MB

#             # # Read and write stdout and stderr incrementally
#             # stdout_line = process.stdout.readline()
#             # stderr_line = process.stderr.readline()

#             # if stdout_line:
#             #     f_stdout.write(stdout_line)
#             # if stderr_line:
#             #     f_err.write(stderr_line)
#             # f_stdout.flush()
#             # f_err.flush()
#         except psutil.NoSuchProcess:    
#             break

#         # Ensure any remaining output is captured
#     stdout, stderr = process.communicate()
#         # if stdout:
#         #     f_stdout.write(stdout)
#         # if stderr:
#         #     f_err.write(stderr)

#     end_time = time.time()
#     time.sleep(1)
#     return (start_time, end_time, memory_usage)

f_res =  open(file_path_res, "w")
d_result = {}

# compute mean conda activate time 
# l_time_act =[]
# l_mem_act=[]
# for i in range(10) :
#     cmd = conda_activate + ' ls'
#     start_time, end_time, memory_usage, stdout, stderr = run_process(cmd.split(' '))
#     l_time_act.append(end_time - start_time)
#     l_mem_act.append(memory_usage)

# f_res.write( f"conda : ({str(l_time_act)},{str(l_mem_act)})\n\n")
# d_result["conda"]=(l_time_act,l_mem_act)






# for work_flow in [nextflow_cmd_stub,smk_cmd_dry]:
for i in setup_nb :
    for w_name,work_flow in    d_cmd.items():
        process_name = w_name+str(i)
        cmd = conda_activate +  ' '+ work_flow+path_setup+str(i)+'/'
        # print(cmd)
        print( w_name+" setup : "+str(i))

        start_time, end_time, memory_usage = run_process(cmd.split(' '), bench_path, process_name)
        

        # if ('snakemake_dry' in w_name ):
        #     file_path_stdout = bench_path + "stdout"+process_name+".txt"
        #     with open(file_path_stdout, "w") as f_stdout:
        #         f_stdout.write(stdout)

        #     file_path_err = bench_path+"stderr"+process_name+".txt"
        #     with open(file_path_err, "w") as f_err:
        #         f_err.write(stderr)
        
        d_result[process_name] = (end_time - start_time,memory_usage)
        f_res.write( f"{process_name} : ({end_time - start_time},{memory_usage})\n")
        
        
        print(f"Time elapsed: {end_time - start_time:.4f} seconds")
        print(f"Memory usage: {memory_usage:.4f} MB\n")

    f_res.write( f"\n")




# with open(file_path_res, "w") as f:
#     f.write(str(d_result))