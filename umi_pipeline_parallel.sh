#!/bin/bash

# Display help function
display_help() {
    echo "$(basename "$0") [-h] [-i n] [options]" >&2
    echo
    echo "   -h, --help           Display this helpful screen."
    echo "   -i, --input-dir      Input directory for fastq files. Default is current working directory."
    echo "   -b  --bed            Assay regions bed file. Default assumes file is in current working directory."
    echo "   -r  --reference      Indexed reference genome"
    echo "   -u  --umi_length     UMI length, default is 12."
    echo "   -s  --spacer_length  Spacer sequence length. Default is 16."
    echo "   -t  --threads        Number of parallel jobs to run. Default is 4."
    echo "   -f  --no_filtering   Do not use fastp to filter fastqs."
    echo "   -q  --phred_score    Min Phred score to keep with fastp filtering. Default is 20."
    echo "   -p  --percent_low_quality  Max percentage of low-quality bases per read. Default is 40."
    echo
    exit 1
}

# Default parameter values
umi_length=19
spacer_length=16
threads=4  # Adjust for parallelism, not necessarily max system cores
do_filtering=true
use_bed=true
phred_score=20
percent_low_quality=40

# Parse command line optionsq
while getopts ':hfi:b:r:u:s:t:q:p:' option; do
  case "$option" in
    h | --help) display_help; exit 0 ;;
    i | --input-dir) FILES=$OPTARG ;;
    b | --bed) BED=$OPTARG ;;
    r | --reference) REF=$OPTARG ;;
    f | --no_filtering) do_filtering=false ;;
    u | --umi_length) umi_length=$OPTARG ;;
    s | --spacer_length) spacer_length=$OPTARG ;;
    t | --threads) threads=$OPTARG ;;
    q | --phred_score) phred_score=$OPTARG ;;
    p | --percent_low_quality) percent_low_quality=$OPTARG ;;
    :) echo "Missing argument for -$OPTARG" >&2; exit 1 ;;
    \?) echo "Illegal option: -$OPTARG" >&2; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

# Set working directory
runDir=${FILES:-$(pwd)}

# Dependency check
printf "Checking dependencies...\n"
for tool in fastp run_umierrorcorrect.py bwa multiqc fastqc; do
  if ! command -v "$tool" &> /dev/null; then
    echo "$tool could not be found. Please install $tool."
    exit 1
  fi
done
echo "All dependencies are present."

# Check bed file existence if specified
if [[ $BED != "" ]] && [[ ! -f "$BED" ]]; then
  echo "Specified bed file does not exist: $BED"
  exit 1
fi

# Reference genome file check
if [[ $REF == "" ]] || [[ ! -f "$REF" ]]; then
  echo "Reference genome file is not specified or does not exist."
  exit 1
fi

# Define the processing function for parallel execution
process_fastq_pair() {
    fq1="$1"
    fq2="${fq1//R1/R2}"
    sample_name=$(basename "${fq1//.fastq.gz/}")  # Only the filename part

    # Define full paths for input and output files
    fq1_fullpath="$(realpath "$fq1")"
    fq2_fullpath="$(realpath "$fq2")"
    outfile="${sample_name}.merged.filtered.fastq.gz"
    outfile_fullpath="$runDir/$outfile"

    echo "Processing sample: $sample_name"

    # Filtering with fastp if enabled
    if $do_filtering; then
        fastp --in1="$fq1" --in2="$fq2" \
              --merge --merged_out="$outfile_fullpath" \
              --qualified_quality_phred="$phred_score" \
              --unqualified_percent_limit="$percent_low_quality" \
              --thread=4

        # Wait briefly to ensure fastp has completed
        sleep 2

        # Check if fastp completed successfully and produced the output file
        if [[ ! -f "$outfile_fullpath" ]]; then
            echo "ERROR: fastp did not produce the expected output file: $outfile_fullpath"
            echo "Skipping umierrorcorrect for $sample_name due to missing merged output."
            return 1
        fi
    else
        # If filtering is disabled, use fq1_fullpath as outfile for compatibility
        outfile_fullpath="$fq1"
    fi

    # Run umierrorcorrect if the outfile exists
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
        echo "ERROR: umierrorcorrect skipped as outfile is missing for sample $sample_name."
    fi
}

# Export function and necessary variables for parallel
export -f process_fastq_pair
export runDir do_filtering use_bed umi_length spacer_length threads BED REF phred_score percent_low_quality

# Use find and parallel to process all fastq.gz files containing R1 in subdirectories
find "$runDir" -type f -name "*R1*.fastq.gz" | parallel -j "$threads" process_fastq_pair {}

# Generate FastQC and MultiQC reports if tools are available
mkdir -p "$runDir/qc_reports"
if command -v fastqc &> /dev/null; then
  echo "Running fastqc..."
  find "$runDir" -type f -name "*.fastq.gz" | xargs fastqc -o "$runDir/qc_reports"
fi

if command -v multiqc &> /dev/null; then
  echo "Running multiqc"
