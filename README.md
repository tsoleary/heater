# Single-nuclei multiome ATAC and RNA sequencing reveals the molecular basis of thermal plasticity in Drosophila melanogaster embryos

### Project Goal

All code used for analysis for the manuscript available <INSERT LINK TO BioRxiv Preprint>.

### File structure

The files in this project are organized into the following directory structure:

- `data/` 
    - `raw/` - all raw data
    - `processed/` - all processed data
    
- `src/` - all source files
    - `00_pheno/` - acclimation phenotype
    - `01_nuclei/` - nuclei isolation
    - `02_cellranger-arc` - mapping sn-multiome reads
    - `03_seurat` - all analysis in Seurat
    - `04_scenic` - all analysis for the SCENIC+ pipeline
    - `05_plots` - scripts used to generate plots
    - `scenic_plus.ipynb` - python notebook for SCENIC+ analysis
    
- `output/` - all output files
    - `figs/`
    - `dars/`
    - `degs/`
    - `markers/`
    - `pheno/`
    - `tf_motif/`
    - `wgcna/`
    
- `docs/` - summary and public facing documents
    - `analysis.html` - summary of all analysis
    - `shiny/` - Shiny App for exploring the data
    
- `scratch/` - auxiliary data analysis that is not a part of the final product
  - `calderon/` information and analysis from Calderon _et al._ 2022
    - `nuclei/` info and analysis related to nuclei troubleshooting
        - `flow/` flow cytometry data and analysis
    - `pbmc/` practice pbmc data set
    - `writing` scraps methods sections scratch writing 

