#!/bin/bash

# Make a reference package for D. melanogaster ---------------------------------

# Directions for running this script -------------------------------------------
# 0. Be sure the following things have happened before beginning
#   - Download FASTA and GTF files for D. melanogaster from Ensembl
#   - Create a contig file according to 10X Genomics instructions
# 1. In the command line, run the following command: bash path/to/this/file.sh

# Set up directories and variables ---------------------------------------------

# Move to the directory where you want the output files to be saved
cd /netfiles02/lockwood_lab/heater/data/processed/seq/ref/

# Point towards installed cellranger-arc program
export PATH=/netfiles02/lockwood_lab/cellranger-arc/cellranger-arc-2.0.2:$PATH

# # Set variables for the script
# ref_path="/netfiles02/lockwood_lab/heater/data/processed/seq/ref/"
# libraries_path="/netfiles02/lockwood_lab/heater/data/raw/seq/libraries/"

# Create a raw unfiltered reference package ------------------------------------
cellranger-arc mkref --config=BDGP6.32.config

# Filter the GTF to restrict it to protein-coding, lncRNA, antisense, and 
#   immune-related genes as 10X Genomics has done for their human and mouse refs
cellranger-arc mkgtf Drosophila_melanogaster.BDGP6.32.109.gtf \
  Drosophila_melanogaster.BDGP6.32.109_filtered.gtf \
  --attribute=gene_biotype:protein_coding \
  --attribute=gene_biotype:lncRNA \
  --attribute=gene_biotype:antisense \
  --attribute=gene_biotype:IG_LV_gene \
  --attribute=gene_biotype:IG_V_gene \
  --attribute=gene_biotype:IG_V_pseudogene \
  --attribute=gene_biotype:IG_D_gene \
  --attribute=gene_biotype:IG_J_gene \
  --attribute=gene_biotype:IG_J_pseudogene \
  --attribute=gene_biotype:IG_C_gene \
  --attribute=gene_biotype:IG_C_pseudogene \
  --attribute=gene_biotype:TR_V_gene \
  --attribute=gene_biotype:TR_V_pseudogene \
  --attribute=gene_biotype:TR_D_gene \
  --attribute=gene_biotype:TR_J_gene \
  --attribute=gene_biotype:TR_J_pseudogene \
  --attribute=gene_biotype:TR_C_gene

# Create a filtered reference package ------------------------------------------
cellranger-arc mkref --config=BDGP6.32.filtered.config 