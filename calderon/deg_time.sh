#!/bin/bash

# To do ----
# Adjust cluster resources
# Set directory
# Set script
# Submit job with the following code in the directory with this .sh file
# > sbatch [filename].sh

# Request cluster resources ----------------------------------------------------
# Adjust the below slurm commands SBATCH

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1

# Request processor cores
#SBATCH --ntasks=1

# Request CPUs per task
#SBATCH --cpus-per-task=4

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem-per-cpu=16G

# Reserve walltime -- hh:mm:ss -- 30 hrs max
#SBATCH --time=30:00:00

# Name this job
#SBATCH --job-name=deg_time

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=%x_%j.out 

# Email reports -- at beginning, end, and failure to my email
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tsoleary12@gmail.com

# ------------------------------------------------------------------------------
# Move to the directory where the script is located
cd ${HOME}/projects/heater/calderon

# Echo some useful and interesting information 
echo "Starting sbatch script at:`date`"
echo "  running host:    ${SLURMD_NODENAME}"
echo "  assigned nodes:  ${SLURM_JOB_NODELIST}"
echo "  partition used:  ${SLURM_JOB_PARTITION}"
echo "  jobid:           ${SLURM_JOBID}"

# Run the script on the galaxy2 cluster
Rscript deg_time.R