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

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=1G

# Reserve walltime -- hh:mm:ss -- 30 hrs below
#SBATCH --time=00:01:00

# Name this job
#SBATCH --job-name=singler_annot

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=%x_%j.out 

# ------------------------------------------------------------------------------
# Move to the directory where the script is located
cd ${HOME}/projects/heater/calderon/slurm_test

# Echo some useful and interesting information 
echo "Starting sbatch script at:`date`"
echo "  running host:    ${SLURMD_NODENAME}"
echo "  assigned nodes:  ${SLURM_JOB_NODELIST}"
echo "  partition used:  ${SLURM_JOB_PARTITION}"
echo "  jobid:           ${SLURM_JOBID}"

# Run the script on the galaxy2 cluster
Rscript myscript.R