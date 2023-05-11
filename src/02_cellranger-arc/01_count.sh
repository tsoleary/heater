#!/bin/bash

# Directions for running this script -------------------------------------------
# 0. Be sure the following things have happened before beginning
#   - Demultiplex all samples
#   - Create reference package
#   - Create a separate libraries csv file pointing to demuxed fastq files
#   - Set up slurm.template file according to 10X Genomics instructions in 
#      cellranger-arc-2.0.2/external/martian/jobmanagers
# 1. In the command line, run the following command: sbatch path/to/this/file.sh

# Set up directories and variables ---------------------------------------------

# Move to the directory where you want the output files to be saved
cd /netfiles02/lockwood_lab/heater/data/processed/seq/samples/

# Point towards installed cellranger-arc program
export PATH=/netfiles02/lockwood_lab/cellranger-arc/cellranger-arc-2.0.2:$PATH

# Set variables for the script
ref_path="/netfiles02/lockwood_lab/heater/data/raw/seq/ref/BDGP6_32_filtered"
libraries_path="/netfiles02/lockwood_lab/heater/data/raw/seq/libraries/"
n_cores=32
mem=64

# Run cellranger-arc count for each sample -------------------------------------

# Sample 18°C Rep 1
cellranger-arc count --id=18C_Rep1 \
  --reference=${ref_path} \
  --libraries=${libraries_path}/18C_Rep1.csv \
  --localcores=${n_cores} \
  --localmem=${mem} \
  --jobmode=slurm
  
# Sample 18°C Rep 2
cellranger-arc count --id=18C_Rep2 \
  --reference=${ref_path} \
  --libraries=${libraries_path}/18C_Rep2.csv \
  --localcores=${n_cores} \
  --localmem=${mem} \
  --jobmode=slurm
 
# Sample 25°C Rep 1
cellranger-arc count --id=25C_Rep1 \
  --reference=${ref_path} \
  --libraries=${libraries_path}/25C_Rep1.csv \
  --localcores=${n_cores} \
  --localmem=${mem} \
  --jobmode=slurm
  
# Sample 25°C Rep 2
cellranger-arc count --id=25C_Rep2 \
  --reference=${ref_path} \
  --libraries=${libraries_path}/25C_Rep2.csv \
  --localcores=${n_cores} \
  --localmem=${mem} \
  --jobmode=slurm