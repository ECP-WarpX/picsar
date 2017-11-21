import os, sys, shutil
import argparse, re, time

# This script runs automated performance tests for WarpX.
# It runs tests in list test_list defined below, and write
# results in file performance_log.txt in warpx/performance_tests/

# ---- User's manual ----
# Before running performance tests, make sure you have the latest version 
# of performance_log.txt
# A typical execution reads:
# > python run_alltests_1node.py --no-recompile --compiler=gnu --machine=cori1 --mode=run --input_file=input_file.pixr --n_node=1 --log_file='my_performance_log.txt'
# These are default values, and will give the same result as 
# > python run_alltests.py
# To add a new test item, extent the test_list with a line like
# test_list.extend([[input_file, nprocx, nprocy, nprocz, nx, ny, nz, ntilex, ntiley, ntilez, n_omp]]*n_repeat)
# - my_input_file must be in picsar/performance_tests
# - the tests must take <10 min on average or they will timeout

# ---- Developer's manual ----
# This script can run in two modes:
# - 'run' mode: for each test item, a batch job is executed.
#     create folder '$SCRATCH/performance_warpx/'
#     recompile the code if option --recompile is used
#     loop over test_list and submit one batch script per item
#     Submit a batch job that executes the script in read mode
#     This last job runs once all others are completed
# - 'read' mode: Get performance data from all test items
#     create performance log file if does not exist
#     loop over test_file 
#         read initialization time and step time
#         write data into the performance log file
#         push file performance_log.txt on the repo

# Read command-line arguments
# ---------------------------

# Create parser and read arguments
parser = argparse.ArgumentParser(
    description='Run performance tests and write results in files')
parser.add_argument('--recompile', dest='recompile', action='store_true', default=False)
parser.add_argument('--no-recompile', dest='recompile', action='store_false', default=False)
parser.add_argument('--commit', dest='commit', action='store_true', default=False)
parser.add_argument( '--compiler', choices=['gnu', 'intel'], default='gnu',
    help='which compiler to use')
parser.add_argument( '--machine', choices=['cori1', 'cori2'], default='cori1',
    help='which NERSC machine')
parser.add_argument( '--mode', choices=['run', 'read'], default='run',
    help='whether to run perftests or read their perf output. run calls read')
parser.add_argument( '--log_file', dest = 'log_file', default='my_performance_log.txt',
    help='name of log file where data will be written. ignored if option --commit is used')
parser.add_argument('--n_node', dest='n_node', default=1, help='nomber of nodes for the runs')
parser.add_argument('--input_file', dest='input_file', default='input_file.pixr', 
    type=str, help='input file to run')

args = parser.parse_args()

log_file = args.log_file
if args.commit == True:
    log_file = 'performance_log.txt'

# Dictionaries
# compiler names. Used for WarpX executable name
compiler_name = {'intel': 'intel', 'gnu': 'gcc'}
# architecture. Used for WarpX executable name
module_name = {'cori1': 'haswell', 'cori2': 'mic-knl'}
# architecture. Used in batch scripts
module_Cname = {'cori1': 'haswell', 'cori2': 'knl,quad,cache'}
# Define environment variables
cwd = os.getcwd() + '/'
res_dir_base = os.environ['SCRATCH'] + '/performance_picsar/'
bin_dir = cwd + 'Bin/'
bin_name = 'picsar_' + args.machine
log_dir  = cwd

day = time.strftime('%d')
month = time.strftime('%m')
year = time.strftime('%Y')
run_name = args.input_file
# Initialize tests
# ----------------
if args.mode == 'run':
# Set default options for compilation and execution
    config_command = ''
    config_command += 'module unload darshan;' 
    config_command += 'module load craype-hugepages4M;'
    if args.machine == 'cori2':
        if args.compiler == 'intel':
            config_command += 'module unload PrgEnv-gnu;'
            config_command += 'module load PrgEnv-intel;'
        elif args.compiler == 'gnu':
            config_command += 'module unload PrgEnv-intel;'
            config_command += 'module load PrgEnv-gnu;'
        config_command += 'module unload craype-haswell;'
        config_command += 'module load craype-mic-knl;'
    elif args.machine == 'cori1':
        if args.compiler == 'intel':
            config_command += 'module unload PrgEnv-gnu;'
            config_command += 'module load PrgEnv-intel;'
        elif args.compiler == 'gnu':
            config_command += 'module unload PrgEnv-intel;'
            config_command += 'module load PrgEnv-gnu;'
        config_command += 'module unload craype-mic-knl;'
        config_command += 'module load craype-haswell;'
    # Create main result directory if does not exist
    if not os.path.exists(res_dir_base):
        os.mkdir(res_dir_base)    

# Recompile if requested
if args.recompile == True:
    with open(cwd + 'Makefile_perftest') as makefile_handler:
        makefile_text = makefile_handler.read()
    with open(cwd + 'Makefile_perftest', 'w') as makefile_handler:
        makefile_handler.write( makefile_text )
    os.system(config_command + " make -f Makefile_perftest clean ; " + " rm *.mod; make -f Makefile_perftest SYS=" + args.machine + " MODE=prod")

# Define functions to run a test and analyse results
# --------------------------------------------------
def run_batch_nnode(test_list, res_dir, n_node=1):
    # Clean res_dir
    if os.path.exists(res_dir):
        shutil.rmtree(res_dir)
    os.makedirs(res_dir)
    # Copy files to res_dir
    shutil.copyfile(bin_dir + bin_name, res_dir + bin_name)
    shutil.copyfile(cwd + args.input_file, res_dir + 'input_file.pixr')
    os.chdir(res_dir)
    # Calculate simulation time. Take 10 min + 10 min / simulation
    job_time_min = 10. + len(test_list)*10.
    job_time_str = str(int(job_time_min/60)) + ':' + str(int(job_time_min%60)) + ':00'
    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#SBATCH --job-name=' + res_dir + '\n'
    batch_string += '#SBATCH --time=' + job_time_str + '\n'
    batch_string += '#SBATCH -C ' + module_Cname[args.machine] + '\n'
    batch_string += '#SBATCH -N ' + str(n_node) + '\n'
    batch_string += '#SBATCH --partition=regular\n'
    batch_string += '#SBATCH --qos=normal\n'
    batch_string += '#SBATCH -e error.txt\n'
    batch_string += '#SBATCH --account=m2852\n'
    for count, test_item in enumerate(test_list):
        # test_item reads
        # [runname, nprocx (total), nprocy (total), nprocz (total),                                                                                                                                 
        #  nx, ny, nz, ntilex, ntiley, ntilez, int n_omp]                                                                                                                                                       
        runname = test_item[0];
        nprocx, nprocy, nprocz = test_item[1:4]
        n_mpi = (nprocx * nprocy * nprocz) / n_node
        nx, ny, nz = test_item[4:7]
        ntilex, ntiley, ntilez = test_item[7:10]
        n_omp = test_item[10]
        srun_string = ''
        srun_string += 'export OMP_NUM_THREADS=' + str(n_omp) + '\n'
        # number of logical cores per MPI process
        if args.machine == 'cori1':
            cflag_value = max(1, int(32/n_mpi) * 2) # Follow NERSC directives
        elif args.machine == 'cori2':
            cflag_value = max(1, int(64/n_mpi) * 4) # Follow NERSC directives
        output_filename = 'out_' + '_'.join([runname, str(n_node), str(nprocx), str(nprocy), \
                          str(nprocz), str(n_omp), str(nx), str(ny), str(nz), str(ntilex), \
                          str(ntiley), str(ntilez)]) + '.txt'
        srun_string += 'srun --cpu_bind=cores '+ \
                       ' -n ' + str(n_node*n_mpi) + \
                       ' -c ' + str(cflag_value)   + \
                       ' ./'  + bin_name + \
                       ' --nx ' + str(nx) + ' --ny ' + str(ny) + ' --nz ' + str(nz) + \
                       ' --nprocx ' + str(nprocx) + ' --nprocy ' + str(nprocy) + ' --nprocz ' + str(nprocz) + \
                       ' --ntilex ' + str(ntilex) + ' --ntiley ' + str(ntiley) + ' --ntilez ' + str(ntilez) + \
                       ' 2> ' + output_filename + '\n'
        batch_string += srun_string
    batch_file = 'slurm'
    f_exe = open(batch_file,'w')
    f_exe.write(batch_string)
    f_exe.close()
    os.system('chmod 700 ' + bin_name)
    os.system(config_command + 'sbatch ' + batch_file + ' >> ' + cwd + 'log_jobids_tmp.txt')
    return 0

# Read output file and return init time and 1-step time
def read_run_perf(filename):
    timing_list = []
    # Search inclusive time to get simulation step time
    partition_limit = 'Step part  min (s)  ave (s)  max (s)  per (%) /it (ms)'
    with open(filename) as file_handler:
        output_text = file_handler.read()
    search_area = output_text.partition(partition_limit)[2]
    # Get total simulation time
    line_match_looptime = re.search(' Total time:.*', search_area)
    time_wo_initialization = float(line_match_looptime.group(0).split()[3])
    timing_list += [str(time_wo_initialization/n_steps)]
    file_handler.close()
    with open(filename) as file_handler:
        output_text = file_handler.read()
    pattern_list = ['Particle pusher \+ field g.*',\
                    'Particle MPI bound\. cond.*',\
                    'Current deposition.*',\
                    'Push bfield fdtd.*',\
                    'Push efield fdtd.*',\
                    'Total time Maxwell solver.*',\
                    'Total time bound. cond.*',\
                    'Particle OpenMP bound.*']
    position_list = [6,5,3,4,4,5,5,5]
    for count, pattern in enumerate(pattern_list):
        timing = '0'
        line_match = re.search(pattern, search_area)
        if line_match is not None:
            timing = [str(float(line_match.group(0).split()[position_list[count]])/n_steps)]
        timing_list += timing
    print(timing_list)
    return timing_list

# Write time into logfile
def write_perf_logfile(log_file):
    log_line = ' '.join([year, month, day, runname, args.compiler,\
                         args.machine, str(n_node), str(nprocx), str(nprocy), str(nprocz), str(n_omp), \
                         str(nx), str(ny), str(nz), str(ntilex), str(ntiley), str(ntilez), \
                         str(n_steps)] +  timing_list + ['\n'])
    f_log = open(log_file, 'a')
    f_log.write(log_line)
    f_log.close()
    return 0

def get_nsteps(runname):
    with open(runname) as file_handler:
        runname_text = file_handler.read()
    line_match_nsteps = re.search('nstep.*', runname_text)
    nsteps = float(line_match_nsteps.group(0).split()[2])
    return nsteps

def process_analysis():
    dependencies = ''
    f_log = open(cwd + 'log_jobids_tmp.txt','r')
    line = f_log.readline()
    dependencies += line.split()[3] + ':'
    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#SBATCH --job-name=perftests_read\n'
    batch_string += '#SBATCH --time=00:05:00\n'
    batch_string += '#SBATCH -C ' + module_Cname[args.machine] + '\n'
    batch_string += '#SBATCH -N 1\n'
    batch_string += '#SBATCH -S 4\n'
    batch_string += '#SBATCH --partition=regular\n'
    batch_string += '#SBATCH --qos=normal\n'
    batch_string += '#SBATCH -e read_error.txt\n'
    batch_string += '#SBATCH -o read_output.txt\n'
    batch_string += '#SBATCH --mail-type=end\n'
    batch_string += '#SBATCH --account=m2852\n'
    batch_string += 'python ' + __file__ + ' --no-recompile --compiler=' + args.compiler + \
                    ' --machine=' + args.machine + ' --mode=read' + ' --log_file=' + log_file + \
                    ' --input_file=' + args.input_file
    if args.commit == True:
        batch_string += ' --commit'
    batch_string += '\n'
    batch_file = 'slurm_perfread'
    f_exe = open(batch_file,'w')
    f_exe.write(batch_string)
    f_exe.close()
    os.system('chmod 700 ' + batch_file)
    os.system('sbatch  --dependency afterok:' + dependencies[0:-1] + ' ' + batch_file)
    return 0
 
# Loop over the tests and return run time + details
# -------------------------------------------------

test_list = []
# TEMPORARY VERSION. USE N_REPEAT=1 AT THE MOMENT. OTHERWISE, SIMILAR RUNS WILL
# WRITE OUTPUT IN THE SAME FILE.
n_repeat = 1
filename1 = args.input_file
# each element of test_list contains
# [str runname, int nprocx (total), nprocy (total), nprocz (total),
#  nx, ny, nz, ntilex, ntiley, ntilez, int n_omp]
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 1, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 2, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 3, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 4, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 5, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 6, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 8, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 10, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 12, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 14, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 16, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 18, 4, 4, 2]]*n_repeat)
test_list.extend([[filename1, 4, 4, 4, 128, 128, 128, 20, 4, 4, 2]]*n_repeat)

n_tests   = len(test_list)

if args.mode == 'run':
    # Remove file log_jobids_tmp.txt if exists.
    # This file contains the jobid of every perf test
    # It is used to manage the analysis script dependencies
    if os.path.isfile(cwd + 'log_jobids_tmp.txt'):
        os.remove(cwd + 'log_jobids_tmp.txt')
    n_steps  = get_nsteps(cwd + run_name)
    res_dir = res_dir_base
    res_dir += '_'.join([year, month, day, args.compiler,\
                         args.machine, args.input_file]) + '/'
    run_batch_nnode(test_list, res_dir, n_node=int(args.n_node))
    os.chdir(cwd)
    process_analysis()

if args.mode == 'read':
    # Create log_file for performance tests if does not exist
    if not os.path.isfile(log_dir + log_file):
        log_line = '## year, month, day, runname, compiler, machine, n_node, nprocx, nprocy, nprocz, ' +\
                   'n_omp, nx, ny, nz, ntilex, ntiley, ntilez, n_steps, step, Particle pusher + field g, ' +\
                   'Particle MPI bound. cond, Current deposition, Push bfield fdtd, Push efield fdtd, '+\
                   'Total time Maxwell solver, Total time bound. cond, Particle OpenMP bound\n'
        f_log = open(log_dir + log_file, 'a')
        f_log.write(log_line)
        f_log.close()
    res_dir = res_dir_base
    res_dir += '_'.join([year, month, day, args.compiler,\
                         args.machine, args.input_file]) + '/'
    n_node   = args.n_node
    for count, test_item in enumerate(test_list):
        # Results folder
        print('read ' + str(test_item))
        runname = test_item[0];
        nprocx, nprocy, nprocz = test_item[1:4]
        n_mpi = (nprocx * nprocy * nprocz) / int(n_node)
        nx, ny, nz = test_item[4:7]
        ntilex, ntiley, ntilez = test_item[7:10]
        n_omp = test_item[10]
        output_filename = 'out_' + '_'.join([str(runname), str(n_node), str(nprocx), str(nprocy), \
                          str(nprocz), str(n_omp), str(nx), str(ny), str(nz), str(ntilex), \
                          str(ntiley), str(ntilez)]) + '.txt'
        n_steps  = get_nsteps(cwd  + args.input_file)
        print('n_steps = ' + str(n_steps))
        # Read performance data from the output file
        timing_list = read_run_perf(res_dir + output_filename)
        # Write performance data to the performance log file
        write_perf_logfile(log_dir + log_file)

    # Store test parameters fot record
    dir_record_base = './perf_warpx_record/'
    if not os.path.exists(dir_record_base):
        os.mkdir(dir_record_base)
    count = 0
    dir_record = dir_record_base + '_'.join([year, month, day]) + '_0'
    while os.path.exists(dir_record):
        count += 1
        dir_record = dir_record[:-1] + str(count)
    os.mkdir(dir_record)
    shutil.copy(__file__, dir_record)
    shutil.copy(log_dir + log_file, dir_record)
    for count, current_run in enumerate(test_list):
        shutil.copy(current_run[0], dir_record)

    # Commit results to the Repo
    if args.commit == True:
        os.system('git add ' + log_dir + log_file + ';'\
                  'git commit -m "performance tests";'\
                  'git push -u origin master')
        
