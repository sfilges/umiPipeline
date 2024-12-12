# UMI processing scripts

![GitHub Release](https://img.shields.io/github/v/release/sfilges/umiPipeline?color=blue)
![CI Status](https://github.com/sfilges/umiPipeline/actions/workflows/ci.yml/badge.svg)
![Shell Check](https://img.shields.io/badge/shellcheck-passed-brightgreen)
![Language](https://img.shields.io/badge/language-bash-blue)

![License](https://img.shields.io/github/license/sfilges/umiPipeline)
![Last Commit](https://img.shields.io/github/last-commit/sfilges/umiPipeline)
![Repo Size](https://img.shields.io/github/repo-size/sfilges/umiPipeline)

## Description

This repository contains a set of scripts for processing sequencing data with UMIs.

## Installation

```bash
git clone https://github.com/sfilges/umiPipeline.git
```

## Basic Usage

- `umi_pipeline_parallel.sh`: This script processes multiple FASTQ files within nested directories in parallel.

Detailed usage instructions are provided in the docs folder.

### `umi_pipeline_parallel.sh`

This script processes multiple FASTQ files in parallel using GNU parallel, which can be supplied via the `-i`/`--input-dir` argument. The fastq files may be located in nested directories or in a single directory. If a bed file is provided via the `-b`/`--bed` argument, all samples will be processed using the provided bed file.

```bash
bash bin/umi_pipeline_parallel.sh -i data/ -r data/test_genome.fa -t 2 --skip-fastqc --skip-multiqc
```

### Reference

The reference genome should be indexed with BWA.

## Project Roadmap

- Add support for single/paired-end combinations
- Add more fine-grained adjustment of threads and number of jobs