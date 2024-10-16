#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation

# Set project storage
#SBATCH -A naiss2024-1-13

#SBATCH -p shared

# The name of the script: sr_pf_muc_type
#SBATCH -J 1_0.55_0.0_SS

# 2 days wall-clock time will be given to this job
#SBATCH --time=23:59:00 #### h:m:s

# Number of nodes
#SBATCH -N 1

### Number of MPI processes.
###SBATCH -n 1

### Number of cores hosting OpenMP threads
#SBATCH -c 64

##SBATCH --dependency=afterok:403029 runscript.sh

#SBATCH -e error_file.e
#SBATCH -o output_file.o

# Set OMP_NUM_THREADS to the same value as -c
# with a fallback in case it isn't set.
# SLURM_CPUS_PER_TASK is set to the value of -c, but only if -c is explicitly set
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi

export OMP_NUM_THREADS=$omp_threads

echo $OMP_NUM_THREADS

# Clear the environment from any previously loaded modules
#module purge > /dev/null 2>&1

# Load the module environment suitable for the job

# Run the executable named myexe 
# and write the output into my_output_file

srun -n 1 ./HLGD1 > terminal_print.txt 2>&1
