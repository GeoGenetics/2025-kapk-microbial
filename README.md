# Kap København 2-million-year-old microbial communities

Analysis code and workflows for Fernandez-Guerra et al. 2025 — ancient metagenomics of the Kap København Formation, Greenland.

## Repository structure

```
.
├── analysis/
│   ├── scripts/          # R scripts for figures and supplementary tables (01–15)
│   │   └── extra/        # Auxiliary scripts (simulations, benchmarks, sourcetracker curation)
│   ├── libs/             # Shared R helper functions
│   ├── beast2/           # BEAST2 XMLs and MCC trees (Methanoflorens stasis)
│   └── .scripts/         # Packaging utilities (stage_erda.sh)
└── wf/
    ├── analysis/
    │   ├── function/     # DART/AGP functional profiling pipeline (KEGG + CAZy)
    │   └── taxonomy/     # Taxonomic profiling + metaDMG authentication pipeline
    ├── binning/          # MAG binning quality assessment (CheckM2, GUNC, GTDB-Tk)
    └── stasis/           # Methanoflorens evolutionary stasis pipeline (phylogenomics + BEAST2)
```

## Data

Input data and pre-computed results are distributed via ERDA (not included in this repository).

Download and extract into the `analysis/` working directory:

```bash
cd analysis
wget https://erda.ku.dk/TODO/data.tar.gz
wget https://erda.ku.dk/TODO/results.tar.gz
wget https://erda.ku.dk/TODO/manuscript.tar.gz
tar xzf data.tar.gz
tar xzf results.tar.gz
tar xzf manuscript.tar.gz
```

The taxonomic database used for profiling is available at: https://doi.org/n6ks

## R analysis scripts

Run from the `analysis/` directory in order:

| Script | Description |
|--------|-------------|
| `01--control-analysis.R` | Control sample QC |
| `02--prepare-taxonomic-data.R` | Taxonomic data processing |
| `03--dmg-threshold-selection.R` | Damage threshold selection |
| `04--taxonomic-profiling.R` | Main taxonomic profiling |
| `05--sourcetracker.R` | Source tracking |
| `06--taxonomic-mag-distribution.R` | MAG distribution analysis |
| `07a--agp-functional.R` | DART/AGP functional profiling (KEGG + CAZy) |
| `08d--virome-figures.R` | Viral community analysis |
| `09--biomarkers.R` | Lipid biomarker analysis |
| `10--briggs.R` | Briggs damage pattern analysis |
| `11-reextractions.R` | Re-extraction comparison |
| `14-summary-tables.R` | Supplementary tables generation |
| `15--methanoflorens-stasis.R` | Methanoflorens evolutionary stasis figure |
| `generate_dart_tables.R` | DART supplementary tables |
| `generate_mag_table.R` | MAG quality table |
| `rerender_heatmap.R` | Functional heatmap figure |

## Workflows

All pipelines use Snakemake (≥7.0) with conda environments.

### Taxonomic profiling (`wf/analysis/taxonomy/`)
Read preprocessing, mapping to the hires taxonomic database, filterBAM filtering, and metaDMG damage authentication.

### Functional profiling — DART/AGP (`wf/analysis/function/`)
DART gene prediction on ancient reads, MMseqs2 search against KEGG/CAZy databases, damage annotation via AGP, and per-sample aggregation.

### MAG quality assessment (`wf/binning/`)
CheckM2 completeness/contamination, GUNC chimerism detection, and GTDB-Tk phylogenetic classification of recovered MAGs.

### Evolutionary stasis (`wf/stasis/`)
Phylogenomics pipeline for Methanoflorens: marker gene extraction, alignment, IQ-TREE ML trees, BEAST2 tip-dating, SNP matrix construction, and pN/pS calculation.

## Requirements

- R ≥ 4.0 with tidyverse, phyloseq, ape, ggtree, openxlsx, and related packages
- Snakemake ≥ 7.0, conda/mamba
- DART/AGP (see [DART repository](https://github.com/genomewalker/damage-aware-read-tagger))

## Citation

Fernandez-Guerra et al. 2025. https://doi.org/10.1101/2023.06.10.544454

## Contact

For questions or issues, please open an issue in this repository.
