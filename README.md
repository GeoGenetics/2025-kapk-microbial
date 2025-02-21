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


## Miscellaneous
### Identifying and authenticating taxa in samples

To identify and authenticate the taxa present in environmental or ancient DNA samples, we use **filterBAM** to filter out low-confidence references and **metaDMG** to authenticate ancient DNA damage patterns. These steps ensure accurate taxonomic identification while distinguishing between ancient and modern microbial communities.

#### **Filtering References with filterBAM**

The **filterBAM filter** module is used to remove low-confidence references based on several metrics:

- **Gini Coefficient Filtering (`-g auto`)**: Measures the uniformity of read coverage across a reference genome, removing references with highly uneven distributions.
- **Entropy Filtering (`-e auto`)**: Measures the distribution of read coverage, removing references with highly uneven distributions.
- **Minimum Read Count (`-n 100`)**: Discards references with fewer than 100 mapped reads.
- **Breadth Expectation Ratio (`-b 0.75`)**: Compares observed breadth to expected breadth, filtering out references with significantly lower coverage.
- **Minimum Breadth (`--min-breadth 0.01`)**: Ensures that at least 1% of a reference genome is covered before inclusion.
- **ANI Filtering (`-A 94 -a 90`)**: Reads are required to have at least 94% ANI for inclusion, while reference-wide ANI must be above 90%.
- **Read Hit Frequency (`--read-hits-count`)**: Tracks how frequently reads map to each reference to refine filtering criteria.

One of the main applications of **bam-filter** is to reliably identify which potential organisms are present in a metagenomic ancient sample, and get relatively accurate taxonomic abundances, even when they are present in very low abundances. The resulting BAM file then can be used as input for [metaDMG](https://github.com/metaDMG-dev/metaDMG-cpp). We rely on several measures to discriminate between noise and a potential signal, analyzing the mapping results at two different levels:

- Is the observed breadth aligned with the expected one?
- Are the reads spread evenly across the reference, or are they clumped in a few regions?

To assess the first question, we use the concepts defined [here](https://doi.org/10.1016/0888-7543(88)90007-9). We estimate the ratio between the observed and expected breadth as a function of the coverage. If we get a **breadth_exp_ratio** close to 1, it means that the coverage we observed is close to the one we expect based on the calculated coverage. While this measure is already a strong indicator of a potential signal, we complement it with the metrics that measure the **normalized positional entropy** and the **normalized distribution inequality** (Gini coefficient) of the positions in the coverage. For details on how these are calculated, check [here](https://www.frontiersin.org/articles/10.3389/fmicb.2022.918015/full). These two metrics help identify cases where we get a high **breadth_exp_ratio**, but the coverage is not evenly distributed across the reference, instead clumping in a few regions.

To calculate these metrics, we bin the reference sequences using [numpy.histogram](https://numpy.org/doc/stable/reference/generated/numpy.histogram.html) to determine the [number of bins](https://numpy.org/doc/stable/reference/generated/numpy.histogram_bin_edges.html#numpy.histogram_bin_edges), either using the Sturges or the Freedman-Diaconis rule. Finally, we apply the [knee point detection algorithm](https://github.com/arvkevi/kneed) to identify optimal thresholds for filtering the Gini coefficient as a function of positional entropy.

Example command used in the manuscript:
```bash
filterBAM filter --bam sample.bam -N -g auto -e auto -n 100 -b 0.75 -c 0.4 -A 94 -a 90 --read-length-freqs --read-hits-count --include-low-detection --min-breadth 0.01
```

#### **Authenticating Damaged References with metaDMG**

To verify whether references contain authentic ancient DNA rather than modern contaminants, **metaDMG** is used in local mode to estimate DNA damage:

- **Post-Mortem Damage Analysis**: Evaluates characteristic damage patterns (e.g., cytosine deamination leading to C-to-T transitions at read termini).
- **Bayesian Estimation**: Uses probabilistic models to determine whether a reference exhibits significant ancient DNA damage.
- **Local Mode**: Assesses damage on an individual reference basis rather than relying on lowest common ancestor (LCA) methods.

By applying **metaDMG**, references that show statistically significant post-mortem damage patterns are retained as authentic ancient DNA, while non-damaged references are flagged as potential modern contaminants.

These two steps—**filterBAM filtering** and **metaDMG authentication**—ensure that only reliable, ancient microbial references are retained for downstream analysis.

## Citation

If you use this workflow, please cite:
https://doi.org/10.1101/2023.06.10.544454

## Contact

For questions or issues, please open an issue in the repository.