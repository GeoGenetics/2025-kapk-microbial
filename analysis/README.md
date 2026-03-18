# R Analysis Scripts

R scripts for Fernandez-Guerra et al. 2025. Run from the `analysis/` directory.

## Data

Input data and pre-computed results are distributed via ERDA. Download and extract into `analysis/`:

```bash
wget https://erda.ku.dk/TODO/data.tar.gz
wget https://erda.ku.dk/TODO/results.tar.gz
wget https://erda.ku.dk/TODO/manuscript.tar.gz
tar xzf data.tar.gz
tar xzf results.tar.gz
tar xzf manuscript.tar.gz
```

## Scripts

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
| `generate_mag_table.R` | MAG quality table (CheckM2 + GUNC + GTDB-Tk) |
| `rerender_heatmap.R` | Functional heatmap figure |

Auxiliary scripts (simulations, benchmarks, sourcetracker curation) are in `scripts/extra/`.

## MIMAG quality thresholds

Applied in `generate_mag_table.R`:
- **HQ**: completeness ≥ 90%, contamination < 5%
- **MQ**: completeness ≥ 50%, contamination < 10%
- **LQ**: all others

Strict MIMAG also requires rRNA genes + ≥18 tRNAs for HQ; not applied here due to ancient DNA fragmentation. GUNC pass/fail is reported separately.

## Helper libraries

- `libs/lib.R`: shared helper functions
- `libs/lib-vir.R`: virus-specific helper functions

## Requirements

R ≥ 4.0 with: tidyverse, phyloseq, ape, ggtree, openxlsx, janitor, ggpubr, phytools, microViz, biomformat, and related packages.
