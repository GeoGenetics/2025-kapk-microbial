# Ancient Metagenomics Analysis for the Kap København microbial manuscript

This repository contains the complete analysis pipeline and downstream analysis code for processing and analyzing ancient DNA (aDNA) metagenomic data. The project is organized into two main components: workflow pipelines (`wf/`) and analysis code (`analysis/`).

## Data Availability

The required data and results folders are available for download from:
- Data: http://files.metagenomics.eu/2025-kapk-microbial/data.tar.gz
- Results: http://files.metagenomics.eu/2025-kapk-microbial/results.tar.gz

Download and extract these archives into the `analysis/` directory:

```bash
cd analysis
wget http://files.metagenomics.eu/2025-kapk-microbial/data.tar.gz
wget http://files.metagenomics.eu/2025-kapk-microbial/results.tar.gz
tar xzf data.tar.gz
tar xzf results.tar.gz
```

The taxonomic database used for the taxonomic profiling can be downloaded [here](https://doi.org/n6ks)

## Project Structure

```
.
├── analysis/                    # Analysis code and data
│   ├── data/                   # Input data directory
│   │   ├── biomarkers/        # Biomarker analysis data
│   │   ├── briggs/            # Briggs analysis data
│   │   ├── cdata/             # Sample metadata
│   │   ├── function/          # Functional analysis data
│   │   ├── mag-distribution/  # MAG distribution data
│   │   ├── re-extractions/    # Re-extraction analysis data
│   │   ├── sourcetracker/     # SourceTracker data
│   │   ├── stats/             # Pipeline statistics
│   │   └── taxonomy/          # Taxonomic data
│   ├── libs/                  # R helper libraries
│   ├── results/               # Analysis results
│   └── scripts/               # R analysis scripts
└── wf/                        # Workflow pipelines
    └── analysis/              # Analysis workflows
        ├── function/          # Functional analysis pipeline
        │   ├── config/        # Pipeline configuration
        │   ├── envs/          # Conda environments
        │   └── rules/         # Snakemake rules
        └── taxonomy/          # Taxonomic analysis pipeline
            ├── config/        # Pipeline configuration
            ├── envs/          # Conda environments
            └── rules/         # Snakemake rules
```

## Components

### Workflow Pipelines (`wf/`)

The project includes two main Snakemake pipelines:

1. **Taxonomic Analysis Pipeline** (`wf/analysis/taxonomy/`)
   - Read preprocessing and quality control
   - Taxonomic profiling
   - Damage pattern analysis
   - See [Taxonomy Pipeline README](wf/analysis/taxonomy/README.md) for details

2. **Functional Analysis Pipeline** (`wf/analysis/function/`)
   - KEGG pathway analysis
   - CAZyme profiling
   - Integration with damage patterns
   - See [Function Pipeline README](wf/analysis/function/README.md) for details

### Analysis Code (`analysis/`)

The analysis directory contains R scripts and supporting files for downstream analysis:

1. **Scripts** (`analysis/scripts/`)
   - Sequential analysis scripts (numbered 01-11)
   - Control analysis through re-extraction analysis
   - See [Analysis README](analysis/README.md) for details

2. **Helper Libraries** (`analysis/libs/`)
   - Common functions and utilities
   - Virus-specific helper functions

### Data Organization

The `analysis/data/` directory contains all input data organized by analysis type:
- Sample metadata (`cdata/`)
- Taxonomic data (`taxonomy/`)
- Functional analysis data (`function/`)
- Various specialized analyses (biomarkers, MAG distribution, etc.)

## Getting Started

1. Install the required dependencies for both pipelines:
   - Follow the installation instructions in each pipeline's README
   - Set up the required conda environments

2. Download and prepare the data:
   - Download and extract the data and results archives
   - Verify the directory structure

3. Run the pipelines:
   - Execute the taxonomic and functional pipelines
   - Verify the output files

4. Run the analysis scripts:
   - Follow the sequence of R scripts
   - Check the results in the output directories

## Requirements

### Pipeline Dependencies
- Snakemake (≥6.0.0)
- Python (≥3.9)
- Conda/Mamba
- Various bioinformatics tools (installed via conda)

### Analysis Dependencies
- R (≥4.0.0)
- Required R packages (tidyverse, phyloseq, etc.)
- Sufficient computational resources for large datasets

## Citation

If you use this workflow, please cite:
https://doi.org/10.1101/2023.06.10.544454

## Contact

For questions or issues, please open an issue in the repository.