# R Analysis Code

This repository contains the R analysis code for processing and analyzing ancient DNA (aDNA) sequence data from Fernandez-Guerra et. al, 2025, with a focus on taxonomic and functional profiling, damage pattern analysis, and various ecological analyses.

## Overview

The analysis workflow consists of multiple R scripts that process the output from the Snakemake pipeline. Key features include:

- Taxonomic data processing and filtering
- Damage pattern analysis and threshold selection
- Source tracking analysis
- MAG distribution analysis
- Functional profiling
- Biomarker analysis
- Viral genome analysis

## Data Availability

The required data and results folders are available for download from:
- Data: http://files.metagenomics.eu/2025-kapk-microbial/data.tar.gz
- Results: http://files.metagenomics.eu/2025-kapk-microbial/results.tar.gz

Download and extract these archives to set up your working directory:

```bash
# Download data and results
wget http://files.metagenomics.eu/2025-kapk-microbial/data.tar.gz
wget http://files.metagenomics.eu/2025-kapk-microbial/results.tar.gz

# Extract archives
tar xzf data.tar.gz
tar xzf results.tar.gz
```

## Requirements

### R Dependencies

The analysis requires R (≥4.0.0) and the following main packages:

```r
# Core packages
library(tidyverse)
library(janitor)
library(ggthemr)
library(ggpubr)
library(showtext)

# Specialized packages
library(phyloseq)
library(biomformat)
library(lvplot)
library(ggdensity)
library(ggpol)
library(phytools)
library(microViz)
library(readxl)
library(unikn)
```

## Script Organization

The analysis consists of multiple R scripts that should be run in sequence:

1. `01--control-analysis.R`: Control sample analysis
2. `02--prepare-taxonomic-data.R`: Initial taxonomic data processing
3. `03--dmg-threshold-selection.R`: Damage pattern threshold selection
4. `04--taxonomic-profiling.R`: Main taxonomic profiling analysis
5. `05--sourcetracker.R`: Source tracking analysis
6. `06--taxonomic-mag-distribution.R`: MAG distribution analysis
7. `07--functional-profiling.R`: Functional profiling
8. `08--virome.R`: Viral genome analysis
9. `09--biomarkers.R`: Biomarker analysis
10. `10--briggs.R`: Briggs pattern analysis
11. `11--reextractions.R`: Re-extraction analysis

## Data Requirements

The analysis expects the following data structure:

```
.
├── data/                          # Input data directory
│   ├── biomarkers/               # Biomarker analysis data
│   ├── briggs/                   # Briggs analysis data
│   ├── cdata/                    # Sample metadata
│   ├── function/                 # Functional analysis data
│   │   ├── dbcan/               # CAZyme annotations
│   │   │   ├── all/             # All samples
│   │   │   └── damaged/         # Damaged samples only
│   │   └── kegg/                # KEGG annotations
│   │       ├── all/             # All samples
│   │       └── damaged/         # Damaged samples only
│   ├── mag-distribution/        # MAG distribution data
│   ├── re-extractions/          # Re-extraction analysis data
│   │   ├── stats/              # Statistics
│   │   └── taxonomy/           # Taxonomic data
│   ├── sourcetracker/          # SourceTracker data
│   │   └── cdata/             # SourceTracker metadata
│   ├── stats/                  # Pipeline statistics
│   └── taxonomy/               # Taxonomic data
│       └── viruses-aa/        # Viral protein analysis
├── libs/                       # Helper libraries
├── results/                    # Output directory
│   ├── controls/              # Control analysis results
│   ├── damage/                # Damage analysis results
│   ├── sourcetracker/         # SourceTracker results
│   └── taxonomy/              # Processed taxonomic results
│       └── anvio/            # Anvio visualization data
└── scripts/                   # Analysis scripts
```

## Usage

### Basic Workflow

1. Set up your working directory with the required data structure
2. Run the scripts in sequence:

```r
# Run each script in order
source("01--control-analysis.R")
source("02--prepare-taxonomic-data.R")
# ... continue with remaining scripts
```

### Key Functions

Each script contains several key functions for data processing and visualization. Common patterns include:

```r
# Data loading functions
load_metadata()
load_taxonomy()
load_abundance()

# Data processing functions
process_taxonomic_data()
calculate_sequence_stats()
analyze_damage_patterns()

# Visualization functions
plot_sequence_distribution()
plot_damage_patterns()
plot_taxonomic_proportions()
```

## Output

The analysis generates various outputs including:

- Filtered and processed taxonomic data
- Damage pattern thresholds and classifications
- Source tracking results
- MAG distribution analyses
- Functional profiles
- Publication-ready figures

## Additional Resources

### Helper Libraries

The analysis uses custom helper functions located in:
- `libs/lib.R`: Core helper functions
- `libs/lib-vir.R`: Virus-specific helper functions

### Constants and Parameters

Key parameters and thresholds are defined in the scripts:

```r
ABUNDANCE_THRESHOLD <- 0.01
SIGNIFICANCE_THRESHOLD <- 2
MIN_SEQUENCES <- 10e6
```

## Citation

If you use this analysis code, please cite:
https://doi.org/10.1101/2023.06.10.544454

## Contact

For questions or issues, please open an issue in the repository.