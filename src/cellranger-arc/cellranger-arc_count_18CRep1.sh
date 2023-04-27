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
#SBATCH --nodes=8

# Request processor cores
#SBATCH --ntasks=8

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem-per-cpu=10G

# Reserve walltime -- hh:mm:ss -- 30 hrs below
#SBATCH --time=10:00:00

# Name this job
#SBATCH --job-name=cellrangerarc_18C_1

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=%x_%j.out 

# ------------------------------------------------------------------------------

# Echo some useful and interesting information 
echo "Starting sbatch script at:`date`"
echo "  running host:    ${SLURMD_NODENAME}"
echo "  assigned nodes:  ${SLURM_JOB_NODELIST}"
echo "  partition used:  ${SLURM_JOB_PARTITION}"
echo "  jobid:           ${SLURM_JOBID}"

# Run the script on the galaxy2 cluster
export PATH=/netfiles02/lockwood_lab/cellranger-arc/cellranger-arc-2.0.2:$PATH
cellranger-arc count --id=18C_Rep1 \
                       --reference=/netfiles02/lockwood_lab/heater/data/raw/seq/ref/BDGP6_32_filtered \
                       --libraries=/netfiles02/lockwood_lab/heater/data/raw/seq/libraries.csv \
                       --localcores=16 \
                       --localmem=64