# Batch Processing Script for Paired-End FASTQ Files

## Overview

This is a batch processing script for paired-end FASTQ files, leveraging `fastp` and `umierrorcorrect` tools for quality filtering, UMI-based error correction, and quality control. The script organizes outputs, manages dependencies, and provides enhanced error handling and logging for streamlined bioinformatics workflows.

## Features

- **FASTQ Quality Filtering**: Uses `fastp` for quality filtering based on Phred scores and low-quality base percentage thresholds.
- **UMI-based Error Correction**: Utilizes `umierrorcorrect` for precise error correction in UMI-tagged data.
- **Parallel Processing**: Supports multi-threaded processing for efficient handling of large datasets.
- **Quality Control**: Optionally generates FastQC and MultiQC reports for comprehensive quality assessment.
- **Flexible Input**: Allows various input parameters like UMI length, reference genome, and BED file for assay regions.

## Requirements

- **Dependencies**: `fastp`, `run_umierrorcorrect.py`, `bwa`, `multiqc`, `fastqc`
- **Shell**: Bash
- **Operating System**: Linux/Unix environment

The script checks for the presence of required dependencies and provides error messages if any are missing.

## Usage

The script accepts several command-line options for customization:

```bash
Usage: ./script.sh [options]

Options:
   -h, --help                  Display help message.
   -i, --input-dir             Input directory containing FASTQ files. Default is the current directory.
   -b, --bed                   Path to the assay regions BED file for umierrorcorrect.
   -r, --reference             Path to the indexed reference genome.
   -u, --umi_length            Length of the Unique Molecular Identifier or UMI. Default is 19.
   -s, --spacer_length         Spacer sequence length between reads. Default is 16.
   -t, --threads               Number of parallel jobs to run. Default is 4.
   -f, --no_filtering          Skip `fastp` filtering step.
   -q, --phred_score           Minimum Phred score for fastp quality filtering. Default is 20.
   -p, --percent_low_quality   Max percentage of low-quality bases per read. Default is 40.
   --skip-fastqc               Skip FastQC quality control.
   --skip-multiqc              Skip MultiQC report generation.
```

# Installation

Clone the repository or download the script.

Ensure dependencies (fastp, umierrorcorrect, bwa, multiqc, fastqc) are installed and accessible in your PATH.

Make the script executable:

```bash
chmod +x script.sh
```

# Example

To process FASTQ files in a specific directory, applying UMI-based error correction with default parameters:

```bash
./script.sh -i /path/to/fastq/files -r /path/to/reference/genome.fa -t 8
```

Note, that the reference genome needs bwa indexed!

# Error handling

Logs are saved in script_log.txt to record the processing of each file, error messages, and status of dependencies. The script terminates if any critical errors (like missing dependencies) occur.