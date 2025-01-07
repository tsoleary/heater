# Single-nuclei multiome ATAC and RNA sequencing reveals the molecular basis of thermal plasticity in Drosophila melanogaster embryos

All code used for analysis for the manuscript available **INSERT LINK TO BioRxiv Preprint**.

### File structure

The files in this project are organized into the following directory structure:
    
- `src/` - all source files
    - `00_pheno/` - acclimation phenotype
    - `01_nuclei/` - nuclei isolation
    - `02_cellranger-arc` - mapping sn-multiome reads
    - `03_seurat` - all analysis in Seurat
    - `04_scenic` - all analysis for the SCENIC+ pipeline
    - `05_plots` - scripts used to generate plots

**Untracked directories with large file sizes**

- `data/` 
    - `raw/` - raw data
    - `processed/` - processed data
    
- `output/` - all output files
    - `figs/`
    - `dars/`
    - `degs/`
    - `markers/`
    - `pheno/`
    - `tf_motif/`
    - `wgcna/`
    
- `docs/`
    - `analysis.html`
    - `shiny/` 
    
- `scratch/` - misc files

