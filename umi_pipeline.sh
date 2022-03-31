#!/bin/bash

# A wrapper script for running umierrorcorrect on multiple samples in the same
# directory. Requires '.fastq.gz' files as input and generates output folders in
# the same directory as input files. Sample names are fastq names without
# the file ending '.fastq.gz'. Requires an assay_regions.bed file in the same
# folder as the fastqs.

#########################
# The command line help #
#########################
display_help() {
    echo "$(basename "$0") [-h] [-i n] [options]" >&2
    echo
    echo "   -h, --help           Display this helpful screen."
    echo "   -i, --input-dir      Input directory for fastq files. Default is current working directory."
    echo "   -b  --bed            Assay regions bed file. Default assumes file is in current working directory."
    echo "   -r  --reference      Indexed reference genome"
    echo
    exit 1
}

##############################
# Check command line options #
#############################


while getopts 'hi:b:r:' option; do
  case "$option" in
    h | --help)
        display_help
        exit 0
        ;;
    i | --input-dir)
        FILES=$OPTARG
        ;;
    b | --bed)
        BED=$OPTARG
        ;;
    r | --reference)
        REF=$OPTARG
        ;;
    :) printf "missing argument for -%i\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

YELLOW=$(tput setaf 3)
GREEN=$(tput setaf 2)
RED=$(tput setaf 1)
NC=$(tput sgr0)

#########################
# Check file integrity and paths #
#########################

# Check dependencies
printf '%s %s %s\n' $GREEN "Checking dependencies..." $NC

# Check if fastp is installed
if ! command -v fastp &> /dev/null
  then
    printf '%s %s\n' $RED "...fastp could not be found."
    printf '%s %s\n' "Please install fastp from: " "https://github.com/OpenGene/fastp"
    exit
  else
    printf '%s %s %s\n' $GREEN "...fastp is installed." $NC
  fi

# Check if umierrorcorrect is installed
if ! command -v run_umierrorcorrect.py &> /dev/null
  then
    printf '%s %s\n' $RED "...umierrorcorrect could not be found."
    printf '%s\n' "Please install umierrorcorrect from: https://github.com/stahlberggroup/umierrorcorrect" $NC
    exit
  else
    printf '%s %s %s\n' $GREEN "...umierrorcorrect is installed." $NC
  fi

# Check if bwa is installed
if ! command -v bwa mem &> /dev/null
  then
    printf '%s %s %s\n' $RED "...bwa could not be found." $NC
    printf ' %s\n' "Please install bwa"
    exit
  else
    printf '%s %s %s\n' $GREEN "...bwa is installed." $NC
  fi

printf '%s %s %s\n' $GREEN "All dependcies are present." $NC

# Check working directory
printf '\n%s %s %s\n' $GREEN "Checking working directory..." $NC
if [[ $FILES = "" ]]
  then
    printf '%s %s %s\n' $YELLOW "...no working directory specified." $NC
    printf '%s %s %s %s\n' $GREEN "...set working directory to: " $NC `pwd`
    runDir=`pwd`
else
  printf '%s %s %s %s\n' $GREEN "...working directory is: " $NC $FILES 
  runDir=$FILES
fi

# Detect fastq files
n_files=`find $runDir -name "*.fastq.gz" | wc -l | sed 's/^ *//g'`

printf '%s %s %s %s\n' $GREEN "...detecting fastqs: "  $NC "$n_files fastqs files found."

# Check bed file
printf "\n%s %s %s\n" $GREEN "Checking bed file..." $NC
if [[ $BED = "" ]]
  then 
    printf '%s %s %s\n' $YELLOW "...no bed file specified." $NC
  else
    if test -f "$BED"
      then
        printf '%s %s %s\n' $GREEN "Using assay regions from: $BED." $NC
      else
        printf '%s %s %s %s\n' $RED "...specified bed file does not exist: " $NC $BED
        exit
    fi
fi

#Check reference genome
printf "\n%s %s %s\n" $GREEN "Checking reference genome..." $NC
if [[ $REF = "" ]]
  #Check if reference is specified
  then 
    printf "%s %s %s\n" $RED "...no reference genome file specified." $NC
    exit
else
  # Check if file exists
  if test -f "$REF"
    then
      # Check if file ending matches fasta
      if [[ "$REF" =~ .*\.(fa|fasta|fn) ]]
        then
          printf '%s %s %s %s\n' $GREEN "...using reference: " $NC "$REF."
      else
        printf '%s %s %s %s\n' $RED "...this does not seem to be a correct fasta file: " $NC $REF
        exit
      fi
  else
    printf '%s %s %s %s\n' $RED "...specified reference does not exist: " $NC $REF
    exit
  fi
fi

###########################
# Run run_umierrorcorrect #
###########################

printf '\n%s %s %s\n' $GREEN "Processing $n_files fastqs..." $NC

cd $runDir

# Processing fastq files
for fastq in *.fastq.gz ;
do
  [[ -f "$fastq" ]] || continue
  find=".fastq.gz"
  replace=""

  sample_name=${fastq//$find/$replace}

  printf "%s %s %s\n" $GREEN "$fastq => $sample_name" $NC

  printf '%s %s %s\n' $GREEN "Running fastp for fastq: $fastq " $NC

  outfile="$sample_name.filtered.fastq.gz"

  fastp \
    -i $fastq \
    -o $outfile \
    --html "$sample_name.html" \
    --trim_poly_g

  printf '%s %s %s\n' $GREEN "Running umierrorcorrect for fastq: $outfile" $NC

  run_umierrorcorrect.py \
    -o $runDir/$sample_name \
    -r1 $runDir/$outfile \
    -r $REF \
    -mode single \
    -ul 12 \
    -sl 16 \
    -bed $BED \
    -t 16 \
    --remove_large_files
done
