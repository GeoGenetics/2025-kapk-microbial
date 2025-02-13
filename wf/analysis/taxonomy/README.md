# Ancient Metagenomics Analysis Workflow: taxonomy

A Snakemake workflow for analyzing ancient DNA (aDNA) sequences with a focus on taxonomic profiling and damage pattern analysis.

## Overview

This pipeline provides a comprehensive workflow for processing and analyzing ancient DNA sequences. Key features include:

- Read preprocessing (renaming, extension, dereplication)
- Taxonomic profiling using multiple reference databases
- Damage pattern analysis with metaDMG
- Support for viral protein analysis
- Flexible input handling for both raw and preprocessed reads
- Extensive quality control and filtering options

## Requirements

### Software Dependencies

The pipeline requires the following main tools:

- Snakemake (≥6.0.0)
- Python (3.9)
- Bowtie2 
- Samtools
- Picard
- SeqKit
- metaDMG
- Various other tools installed via conda environments

### Reference Databases

The pipeline requires several reference databases:
- Taxonomic database (specified in config as `bowtie2_tp_db`)
- MAG databases for distribution analysis
- Viral protein databases
- Other specialized databases (archaeal viruses, borgs, etc.)

## Installation

1. Clone the repository:
```bash
git clone [repository-url]
```

2. Create and activate the required conda environments:
```bash
# Create environments from provided YAML files
conda create -n snakemake
```

3. Configure the pipeline by modifying `config.yaml` with your specific paths and parameters.

## Usage

### Basic Usage

1. Prepare your input data:
   - Place raw sequence files in your source directory
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
wdir: "path/to/working/dir"
sdir: "path/to/source/dir"
rdir: "path/to/results/dir"

# Input parameters
sample_file_read: "path/to/samples.tsv"
read_minlen: 30
extend_reads: true

# Processing parameters
derep_reads: true
use_vsearch: true
tp_bowtie2_k: [1000]
mdmg_weight: [1]
```

### Sample Sheet Format

The sample sheet (`samples.tsv`) should contain the following columns:
- `label`: Sample identifier
- `file`: Input file name
- `study_type`: "ancient" or "modern" (optional, defaults to "ancient")
- `libprep`: Library preparation method (optional, defaults to "double")

## Pipeline Steps

1. **Read Preprocessing**
   - Initial statistics collection
   - Read renaming (optional)
   - Read extension for ancient samples
   - Read dereplication

2. **Taxonomic Profiling**
   - Mapping against reference databases
   - Filtering of mapped reads
   - Coverage and abundance calculations

3. **Damage Pattern Analysis**
   - metaDMG analysis on filtered alignments
   - Damage pattern detection and quantification
   - Read classification based on damage patterns

4. **Specialized Analyses**
   - Viral protein mapping and analysis
   - MAG distribution analysis
   - Additional database searches

## Output Structure

The pipeline generates results in the following directory structure:

```
results/
├── stats/                    # Quality control statistics
├── read-renamed/            # Renamed read files
├── read-extension/          # Extended reads
├── read-derep/              # Dereplicated reads
├── taxonomic-profiling*/    # Taxonomic profiling results
├── read-dmg*/               # Damage pattern analysis
└── *-summary/               # Summary statistics and reports
```

## Advanced Usage

### Customizing Reference Databases

The pipeline supports multiple reference databases defined in `config.yaml`:

```yaml
db_mag_distribution:
  - database_name:
      bowtie2_db: "path/to/bowtie2/index"
      db_mag_sizes: "path/to/genome/sizes"

db_other:
  - custom_db:
      bowtie2_db: "path/to/custom/db"
      db_genome_sizes: "path/to/sizes"
```

### Performance Tuning

Key parameters for performance optimization:

```yaml
# Threading
seqkit_threads: 16
bowtie2_tp_threads: 24
mdmg_threads: 8

# Memory usage
tp_filter_sort_mem: 8G
picard_java_opts: "-Xms2g -Xmx100g"
```

## Troubleshooting

Common issues and solutions:

1. **Memory Issues**
   - Adjust memory parameters in config.yaml
   - Use appropriate thread counts for your system

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