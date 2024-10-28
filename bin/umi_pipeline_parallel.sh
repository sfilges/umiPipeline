#!/bin/bash

# A batch processing script for paired-end FASTQ files using fastp and umierrorcorrect.
# Organizes outputs, manages dependencies, and provides enhanced error handling and logging.

#########################
# The command line help #
#########################
display_help() {
    echo "Usage: $(basename "$0") [options]"
    echo
    echo "Options:"
    echo "   -h, --help                  Display this help message."
    echo "   -i, --input-dir             Input directory containing FASTQ files. Default is the current directory."
    echo "   -b, --bed                   Path to the assay regions BED file for umierrorcorrect."
    echo "   -r, --reference             Path to the indexed reference genome."
    echo "   -u, --umi_length            Length of the Unique Molecular Identifier (UMI). Default is 12."
    echo "   -s, --spacer_length         Spacer sequence length between reads. Default is 16."
    echo "   -t, --threads               Number of parallel jobs to run. Default is 4."
    echo "   -f, --no_filtering          Skip fastp filtering step for FASTQ files."
    echo "   -q, --phred_score           Minimum Phred score threshold for fastp quality filtering. Default is 20."
    echo "   -p, --percent_low_quality   Maximum percentage of low-quality bases per read. Default is 40."
    echo "   --skip-fastqc               Skip FastQC quality control."
    echo "   --skip-multiqc              Skip MultiQC report generation."
    echo
    exit 1
}

##########################
# Default Parameter Values
##########################
umi_length=19
spacer_length=16
threads=4
do_filtering=true
use_bed=false
phred_score=20
percent_low_quality=40
skip_fastqc=false
skip_multiqc=false

######################
# Parse Command-Line Options
######################
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -h|--help) display_help ;;
    -i|--input-dir) FILES="$2"; shift ;;
    -b|--bed) BED="$2"; shift ;;
    -r|--reference) REF="$2"; shift ;;
    -u|--umi_length) umi_length="$2"; shift ;;
    -s|--spacer_length) spacer_length="$2"; shift ;;
    -t|--threads) threads="$2"; shift ;;
    -f|--no_filtering) do_filtering=false ;;
    -q|--phred_score) phred_score="$2"; shift ;;
    -p|--percent_low_quality) percent_low_quality="$2"; shift ;;
    --skip-fastqc) skip_fastqc=true ;;
    --skip-multiqc) skip_multiqc=true ;;
    *) echo "Unknown option: $1" ; display_help ;;
  esac
  shift
done

# Set working directory, default to current if not provided
runDir=${FILES:-$(pwd)}

# Set log file and temporary directories
LOG_FILE="$runDir/script_log.txt"
FILTERED_DIR="$runDir/filtered_fastqs"
UMI_CORRECTED_DIR="$runDir/umi_corrected_samples"
TEMP_DIR=$(mktemp -d -t umi_process_XXXXXX)
trap 'rm -rf "$TEMP_DIR"' EXIT  # Clean up temp directory on script exit

mkdir -p "$FILTERED_DIR" "$UMI_CORRECTED_DIR" "$runDir/qc_reports"
exec > >(tee -a "$LOG_FILE") 2>&1  # Log both stdout and stderr

######################
# Logging Function
######################
log_msg() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') - $1"
}

######################
# Dependency Check
######################
log_msg "Checking for required dependencies..."
missing_deps=()
for tool in fastp run_umierrorcorrect.py bwa multiqc fastqc; do
  if ! command -v "$tool" &> /dev/null; then
    missing_deps+=("$tool")
  fi
done
if [ "${#missing_deps[@]}" -ne 0 ]; then
  echo "Error: Missing dependencies: ${missing_deps[*]}"
  exit 1
fi
log_msg "All dependencies are present."

######################
# Check Paths
######################
check_paths() {
    [[ ! -d "$runDir" ]] && echo "Error: Input directory does not exist: $runDir" && exit 1
    [[ ! -f "$REF" ]] && echo "Error: Reference genome file not found: $REF" && exit 1
    if $use_bed && [[ ! -f "$BED" ]]; then
        echo "Error: BED file not found: $BED"
        exit 1
    fi
}
check_paths

##############################
# Define Processing Functions
##############################
run_fastp() {
    local fq1="$1"
    local fq2="$2"
    local outfile="$3"

    log_msg "Running fastp for files $fq1 and $fq2"
    fastp --in1="$fq1" --in2="$fq2" \
          --merge --merged_out="$outfile" \
          --qualified_quality_phred="$phred_score" \
          --unqualified_percent_limit="$percent_low_quality" \
          --thread=$((threads / 2)) || {
        log_msg "Error: fastp failed for $fq1 and $fq2"
        return 1
    }
}

run_umierrorcorrect() {
    local outfile="$1"
    local sample_name="$2"

    log_msg "Running umierrorcorrect for sample $sample_name"
    run_umierrorcorrect.py -o "$UMI_CORRECTED_DIR/$sample_name" -r1 "$outfile" -r "$REF" \
                           -mode single -ul "$umi_length" -sl "$spacer_length" \
                           ${use_bed:+-bed "$BED"} -t "$threads" || {
        log_msg "Error: umierrorcorrect failed for $sample_name"
        return 1
    }
}

process_fastq_pair() {
    fq1="$1"
    fq2="${fq1//R1/R2}"
    sample_name=$(basename "${fq1//.fastq.gz/}")

    # Define output paths
    outfile="$FILTERED_DIR/${sample_name}.merged.filtered.fastq.gz"

    log_msg "Processing sample: $sample_name"

    # Run fastp if filtering is enabled
    if $do_filtering; then
        run_fastp "$fq1" "$fq2" "$outfile" || return 1
    else
        outfile="$fq1"  # Use original file if filtering is disabled
    fi

    # Run umierrorcorrect
    run_umierrorcorrect "$outfile" "$sample_name" || return 1
}

######################
# Export and Run Parallel Processing
######################
export -f process_fastq_pair log_msg run_fastp run_umierrorcorrect
export runDir FILTERED_DIR UMI_CORRECTED_DIR TEMP_DIR do_filtering use_bed umi_length spacer_length threads BED REF phred_score percent_low_quality

find "$runDir" -type f -name "*R1*.fastq.gz" | parallel -j "$threads" process_fastq_pair {}

##############################
# Quality Control with FastQC and MultiQC
##############################
if ! $skip_fastqc && command -v fastqc &> /dev/null; then
    log_msg "Running FastQC in parallel..."
    find "$runDir" -type f -name "*.fastq.gz" | parallel -j "$threads" fastqc -o "$runDir/qc_reports" {}
fi

if ! $skip_multiqc && command -v multiqc &> /dev/null; then
    log_msg "Running MultiQC..."
    multiqc "$runDir" -o "$runDir/qc_reports"
fi

log_msg "Processing complete."
