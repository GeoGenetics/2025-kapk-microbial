# Ancient Metagenomics Analysis Workflow: function

A Snakemake workflow for analyzing functional potential in ancient DNA (aDNA) sequences, with a focus on KEGG pathway analysis and CAZyme profiling.

## Overview

This pipeline provides a comprehensive workflow for processing and analyzing ancient DNA sequences with a focus on functional analysis. Key features include:

- Read preprocessing (quality control, extension, dereplication)
- KEGG pathway analysis
- CAZyme (dbCAN) profiling
- Damage pattern analysis
- Support for both ancient and modern DNA samples
- Extensive quality control and filtering options

## Requirements

### Software Dependencies

The pipeline requires the following main tools:

- Snakemake (≥6.0.0)
- Python (3.9+)
- MMseqs2
- SeqKit
- VSEARCH
- BBMap (Tadpole)
- Various other tools installed via conda environments

### Reference Databases

The pipeline requires several reference databases:
- KEGG database (specified in config as `kegg_db`)
- dbCAN database for CAZyme analysis
- Other specialized databases as needed

## Installation

1. Clone the repository:
```bash
git clone [repository-url]
```

2. Create and activate the required conda environments:
```bash
# Create environments from provided YAML files
conda env create -f envs/qc.yaml
conda env create -f envs/kegg-xfilter.yaml
conda env create -f envs/kegg-xfilter-anvio.yaml
conda env create -f envs/get-dmg-reads.yaml
```

3. Configure the pipeline by modifying `config.yaml` with your specific paths and parameters.

## Usage

### Basic Usage

1. Prepare your input data:
   - Place sequence files in your source directory
   - Create a sample sheet (TSV format) with required columns: `label` and `file`

2. Configure your run:
   - Modify `config.yaml` to specify input/output directories and parameters
   - Set appropriate computational resources

3. Run the pipeline:
```bash
snakemake -s Snakefile --use-conda -j [cores]
```

### Configuration

Key configuration options in `config.yaml`:

```yaml
# Directories
wdir: "/path/to/working/dir"
sdir: "/path/to/source/dir"
rdir: "/path/to/results/dir"

# Input parameters
sample_file_read: "path/to/samples.tsv"
read_minlen: 30

# Processing parameters
rename_reads: false
derep_reads: true
extend_reads: true
get_dmg_reads: true

# Tool-specific parameters
seqkit_threads: 16
mmseqs2_threads: 32
anvio_threads: 4
```

### Sample Sheet Format

The sample sheet (`samples.tsv`) should contain the following columns:
- `label`: Sample identifier
- `file`: Input file name
- `study_type`: "ancient" or "modern" (optional, defaults to "ancient")
- `libprep`: Library preparation method (optional, defaults to "double")

## Pipeline Steps

1. **Read Preprocessing**
   - Initial statistics collection (`stats-initial.smk`)
   - Optional read renaming (`read-rename.smk`)
   - Read extension for ancient samples (`read-extension.smk`)
   - Read dereplication (`read-derep.smk`)

2. **KEGG Analysis**
   - MMseqs2 profile search against KEGG database
   - Filtering and abundance calculations
   - KEGG module analysis using Anvi'o
   - Summarization of KEGG hits and modules

3. **CAZyme Analysis**
   - MMseqs2 profile search against dbCAN database
   - Filtering and abundance calculations
   - Summary statistics generation

4. **Damage Pattern Analysis**
   - Identification of damaged reads
   - Separate functional analysis of damaged reads
   - Integration with KEGG and CAZyme results

## Output Structure

The pipeline generates results in the following directory structure:

```
results/
├── stats/                    # Quality control statistics
├── read-renamed/            # Renamed read files
├── read-extension/          # Extended reads
├── read-derep/              # Dereplicated reads
├── kegg-reads/              # KEGG analysis results
├── kegg-reads-anvio/        # Anvi'o KEGG module analysis
├── dbcan-reads/             # CAZyme analysis results
├── read-dmg/                # Damage pattern analysis
└── *-summary/               # Summary statistics and reports
```

## Advanced Usage

### KEGG Analysis Parameters

Key parameters for KEGG analysis:

```yaml
mmseqs_kegg_parms: "--comp-bias-corr 0 --mask 0 -e 1e-5 --exact-kmer-matching 1 --sub-mat PAM30.out -s 3 -k 6 --spaced-kmer-mode 1 --spaced-kmer-pattern 11011101 --min-length 15 --format-mode 2 -c 0.8 --cov-mode 2 --min-seq-id 0.6"
xfilter_parms: "-n 25 -b 20 -e 1e-5 --breadth-expected-ratio 0 -f depth_evenness --depth-evenness 1"
```

### Read Extension Parameters

Parameters for ancient DNA read extension:

```yaml
tadpole_mode: "extend"
tadpole_length: 100
tadpole_k: 17
trimends: 9
ibb: false
ecc: true
ecco: false
conservative: true
```

## Troubleshooting

Common issues and solutions:

1. **Memory Issues**
   - Adjust memory parameters in config.yaml
   - Monitor resource usage with the benchmarking files

2. **Missing Files**
   - Verify all reference database paths
   - Check input file permissions and formats

3. **Pipeline Failures**
   - Check log files in the logs/ directory
   - Verify conda environment activation

## Citation

If you use this workflow, please cite:
https://doi.org/10.1101/2023.06.10.544454

## Contact

For questions or issues, please open an issue in the repository.