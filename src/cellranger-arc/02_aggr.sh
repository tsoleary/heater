#!/bin/bash

# Aggregate together all samples -----------------------------------------------

# Directions for running this script -------------------------------------------
# 0. Before running run cellranger-arc count for all samples and create 
#      aggr_libraries.csv file pointing to atac_fragments.tsv.gz and 
#      gex_molecule_info.h5 files per 10X Genomics instructions.
# 1. In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=31000

# Reserve walltime -- hh:mm:ss -- 30 hrs max below
#SBATCH --time=30:00:00

# Name this job
#SBATCH --job-name=aggr_all

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=logs/%x_%j.out 

# Echo some useful and interesting information 
echo "Starting sbatch script at: `date`"
echo "  running host:    ${SLURMD_NODENAME}"
echo "  assigned nodes:  ${SLURM_JOB_NODELIST}"
echo "  partition used:  ${SLURM_JOB_PARTITION}"
echo "  jobid:           ${SLURM_JOBID}"

# Set up directories and variables ---------------------------------------------

# Move to the directory where files will be aggregated
cd /netfiles02/lockwood_lab/heater/data/processed/seq/

# Point towards installed cellranger-arc program
export PATH=/netfiles02/lockwood_lab/cellranger-arc/cellranger-arc-2.0.2:$PATH

# Set variables for the script
ref_path="/netfiles02/lockwood_lab/heater/data/raw/seq/ref/BDGP6_32_filtered"
libraries_path="/netfiles02/lockwood_lab/heater/data/processed/seq/samples/aggr_libraries.csv"

# Run cellranger-arc aggr to aggregate all samples together --------------------
cellranger-arc aggr --id=all \
  --csv=${libraries_path} \
  --reference=${ref_path} \
  --normalize=depth

