#!/bin/bash

# A wrapper script for batch processing of paired-end FASTQ files using fastp and umierrorcorrect.
# This script takes paired-end .fastq.gz files in the specified directory, applies filtering and merging,
# and performs UMI error correction. Output files are created in the same directory as the input files.

#########################
# The command line help #
#########################
display_help() {
    echo "Usage: $(basename "$0") [options]" >&2
    echo
    echo "Options:"
    echo "   -h, --help                  Display this help message."
    echo "   -i, --input-dir             Input directory containing FASTQ files. Default is the current directory."
    echo "   -b, --bed                   Path to the assay regions BED file for umierrorcorrect. Default assumes it's in the current directory."
    echo "   -r, --reference             Path to the indexed reference genome."
    echo "   -u, --umi_length            Length of the Unique Molecular Identifier (UMI). Default is 12."
    echo "   -s, --spacer_length         Spacer sequence length between reads. Default is 16."
    echo "   -t, --threads               Number of parallel jobs to run. Default is 4."
    echo "   -f, --no_filtering          Skip fastp filtering step for FASTQ files."
    echo "   -q, --phred_score           Minimum Phred score threshold for quality filtering with fastp. Default is 20."
    echo "   -p, --percent_low_quality   Maximum percentage of low-quality bases per read. Default is 40."
    echo
    exit 1
}

##########################
# Default Parameter Values
##########################
umi_length=12
spacer_length=16
threads=4
do_filtering=true
use_bed=true
phred_score=20
percent_low_quality=40

################################
# Parse Command-Line Options
################################
while getopts ':hfi:b:r:u:s:t:q:p:' option; do
  case "$option" in
    h | --help) display_help ;;
    i | --input-dir) FILES=$OPTARG ;;
    b | --bed) BED=$OPTARG ;;
    r | --reference) REF=$OPTARG ;;
    f | --no_filtering) do_filtering=false ;;
    u | --umi_length) umi_length=$OPTARG ;;
    s | --spacer_length) spacer_length=$OPTARG ;;
    t | --threads) threads=$OPTARG ;;
    q | --phred_score) phred_score=$OPTARG ;;
    p | --percent_low_quality) percent_low_quality=$OPTARG ;;
    :) echo "Error: Missing argument for -$OPTARG" >&2; exit 1 ;;
    \?) echo "Error: Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

# Set input directory, defaulting to current directory if not provided
runDir=${FILES:-$(pwd)}

##############################
# Dependency Checks
##############################
echo "Checking for required dependencies..."
for tool in fastp run_umierrorcorrect.py bwa multiqc fastqc; do
  if ! command -v "$tool" &> /dev/null; then
    echo "Error: $tool not found. Please install $tool."
    exit 1
  fi
done
echo "All dependencies are present."

# Verify BED file and reference genome, if specified
if [[ -n $BED && ! -f $BED ]]; then
  echo "Error: Specified BED file does not exist: $BED"
  exit 1
fi
if [[ -z $REF || ! -f $REF ]]; then
  echo "Error: Reference genome file is missing or does not exist."
  exit 1
fi

##############################
# Define Processing Function #
##############################
# process_fastq_pair: Processes each pair of FASTQ files.
#                     - Runs fastp to filter and merge reads if enabled
#                     - Performs UMI error correction using umierrorcorrect

process_fastq_pair() {
    fq1="$1"
    fq2="${fq1//R1/R2}"                            # Identify the matching R2 file
    sample_name=$(basename "${fq1//.fastq.gz/}")    # Extract sample name (filename without path or extension)

    # Full paths for input files and output locations
    fq1_fullpath="$(realpath "$fq1")"
    fq2_fullpath="$(realpath "$fq2")"
    outfile="${sample_name}.merged.filtered.fastq.gz"
    outfile_fullpath="$runDir/$outfile"

    echo "Processing sample: $sample_name"

    # Run fastp for filtering if enabled
    if $do_filtering; then
        fastp --in1="$fq1" --in2="$fq2" \
              --merge --merged_out="$outfile_fullpath" \
              --qualified_quality_phred="$phred_score" \
              --unqualified_percent_limit="$percent_low_quality" \
              --thread=4

        # Wait to ensure fastp completes file creation
        sleep 2

        # Verify fastp output file exists
        if [[ ! -f "$outfile_fullpath" ]]; then
            echo "Error: fastp did not produce expected output file: $outfile_fullpath"
            echo "Skipping umierrorcorrect for $sample_name."
            return 1
        fi
    else
        # If filtering is skipped, use fq1_fullpath directly
        outfile_fullpath="$fq1"
    fi

    # Run umierrorcorrect with the generated or original output file
    if [[ -f "$outfile_fullpath" ]]; then
        if $use_bed; then
            run_umierrorcorrect.py -o "$runDir/$sample_name" -r1 "$outfile_fullpath" -r "$REF" \
                                   -mode single -ul "$umi_length" -sl "$spacer_length" \
                                   -bed "$BED" -t "$threads"
        else
            run_umierrorcorrect.py -o "$runDir/$sample_name" -r1 "$outfile_fullpath" -r "$REF" \
                                   -mode single -ul "$umi_length" -sl "$spacer_length" -t "$threads"
        fi
    else
        echo "Error: umierrorcorrect skipped as outfile is missing for sample $sample_name."
    fi
}

##################################
# Export Function and Environment
##################################
export -f process_fastq_pair
export runDir do_filtering use_bed umi_length spacer_length threads BED REF phred_score percent_low_quality

##################################
# Process all R1 FASTQ Files in Parallel
##################################
# Finds all FASTQ files with R1 in the name in the specified directory and subdirectories
# Runs the process_fastq_pair function for each file in parallel, up to the specified number of threads

find "$runDir" -type f -name "*R1*.fastq.gz" | parallel -j "$threads" process_fastq_pair {}

##############################
# Generate FastQC and MultiQC Reports
##############################
# Runs FastQC and MultiQC on processed files, generating summary reports

mkdir -p "$runDir/qc_reports"
if command -v fastqc &> /dev/null; then
  echo "Running FastQC on all FASTQ files..."
  find "$runDir" -type f -name "*.fastq.gz" | xargs fastqc -o "$runDir/qc_reports"
fi

if command -v multiqc &> /dev/null; then
  echo "Running MultiQC..."
  multiqc "$runDir" -o "$runDir/qc_reports"
fi
