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
    echo "   -u --umi_length     UMI length, default is 12."
    echo "   -s --spacer_length  How long is the spacer sequence? Default is 16."
    echo "   -t  --threads        How many threads to use? Default is 16."
    echo "   -f  --no_filtering   Do not use fastp to filter fastqs."
    echo "   -q  --phred_score          Min Phread score to keep when using fastp to filter. Default is 15, typical values are 10, 15, 20, 30."
    echo "   -p  --percent_low_quality  How many percent bases in a read are allowed to be below the thrshold q value. Default is 40 (0-100)"
    echo
    exit 1
}

##############################
# Check command line options #
##############################

umi_length=12
spacer_length=16
threads=16
no_fastp=false
use_bed=true
do_filtering=true
phred_score=15
percent_low_quality=40
fastqc=false
multiqc=false

touch log.txt

while getopts ':hfi:b:r:u:s:t:q:p:' option; do
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
    f | --no_filtering)
        do_filtering=false
        ;;
    u | --umi_length)
        umi_length=$OPTARG
        ;;
    s | --spacer_length)
        spacer_length=$OPTARG
        ;;
    t | --threads)
        threads=$OPTARG
        ;;
    q | --phred_score)
        phred_score=$OPTARG
        ;;
    p | --percent_low_quality)
        percent_low_quality=$OPTARG
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

##################################
# Check file integrity and paths #
##################################

# Check dependencies
printf '%s %s %s\n' $GREEN "Checking dependencies..." $NC

# Check if fastp is installed
if ! command -v fastp &> /dev/null
  then
    printf '%s %s\n' $YELLOW "...fastp could not be found. No filtering will be performed."
    printf '%s %s\n' "Please install fastp if you want to use filtering: " "https://github.com/OpenGene/fastp"
    no_fastp=true
  else
    printf '%s %s %s\n' $YELLOW "...fastp is installed." $NC
  fi

# Check if umierrorcorrect is installed
if ! command -v run_umierrorcorrect.py &> /dev/null
  then
    printf '%s %s\n' $RED "...umierrorcorrect could not be found."
    printf '%s\n' "Please install umierrorcorrect from: https://github.com/stahlberggroup/umierrorcorrect" $NC
    exit
  else
    printf '%s %s %s\n' $YELLOW "...umierrorcorrect is installed." $NC
  fi

# Check if bwa is installed
if ! command -v bwa mem &> /dev/null
  then
    printf '%s %s %s\n' $RED "...bwa could not be found." $NC
    printf ' %s\n' "Please install bwa"
    exit
  else
    printf '%s %s %s\n' $YELLOW "...bwa is installed." $NC
  fi

# Check if multiqc is installed
if ! command -v multiqc &> /dev/null
  then
    printf '%s %s %s\n' $YELLOW "...multiqc could not be found. No merged reports will be generated." $NC
    printf ' %s\n' "Please install multqic"
  else
    printf '%s %s %s\n' $YELLOW "...multiqc is installed." $NC
    multiqc=true
  fi

# Check if fastqc is installed
if ! command -v fastqc &> /dev/null
  then
    printf '%s %s %s\n' $YELLOW "...fastqc could not be found. No merged reports will be generated." $NC
    printf ' %s\n' "Please install fastqc"
  else
    printf '%s %s %s\n' $YELLOW "...fastqc is installed." $NC
    fastqc=true
  fi

printf '%s %s %s\n' $GREEN "All dependencies are present." $NC

# Check working directory
printf '\n%s %s %s\n' $GREEN "Checking working directory..." $NC
if [[ $FILES = "" ]]
  then
    printf '%s %s %s\n' $YELLOW "...no working directory specified." $NC
    printf '%s %s %s %s\n' $GREEN "...set working directory to: " $NC `pwd`
    runDir=`pwd`
else
  printf '%s %s %s %s\n' $YELLOW "...working directory is: " $NC $FILES 
  runDir=$FILES
fi

# Detect fastq files
n_files=`find $runDir -maxdepth 1 -name "*.fastq.gz" | wc -l | sed 's/^ *//g'`

printf '%s %s %s %s\n' $YELLOW "...detecting fastqs: "  $NC "$n_files fastqs files found."

# Check bed file
printf "\n%s %s %s\n" $GREEN "Checking bed file..." $NC

if [[ $BED = "" ]]
  then 
    printf '%s %s %s\n' $YELLOW "...no bed file specified. Running without bed annotations." $NC
    use_bed=false
  else
    if test -f "$BED"
      then
        printf '%s %s %s %s\n' $YELLOW "Using assay regions from:" $NC $BED
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
          printf '%s %s %s %s\n' $YELLOW "...using reference: " $NC "$REF."
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

printf '\n%s %s %s\n' $GREEN "Processing $n_files fastq files..." $NC
printf '%s %s %s %s\n' $YELLOW "...using UMI length:" $NC $umi_length
printf '%s %s %s %s\n' $YELLOW "...using spacer length:" $NC $spacer_length
printf '%s %s %s %s\n' $YELLOW "...using threads:" $NC $threads
printf '%s %s %s %s\n' $YELLOW "...perform filtering:" $NC $do_filtering
printf '%s %s %s %s\n' $YELLOW "...using bed annotations:" $NC $use_bed

cd $runDir

# Processing fastq files
for fastq in *.fastq.gz ;
do
  [[ -f "$fastq" ]] || continue
  find=".fastq.gz"
  replace=""

  sample_name=${fastq//$find/$replace}

  printf "\n%s %s %s %s\n" $GREEN "Changing sample name:" $NC "$fastq => $sample_name" 

  if $no_fastp
    then
      outfile=$fastq
  else
    if $do_filtering
      then
        outfile="$sample_name.filtered.fastq.gz"

        printf '\n%s %s %s %s\n' $GREEN "Running fastp for fastq..." $NC $fastq
        printf '%s %s %s %s\n' $YELLOW "...using minimum Phread score:" $NC $phred_score 
        printf '%s %s %s %s\n' $YELLOW "...using max percent low quality reads:" $NC $percent_low_quality

        fastp \
          -i $fastq \
          -o $outfile \
          --html "$sample_name.html" \
          --trim_poly_g \
          -q $phred_score \
          -u $percent_low_quality
    else
      outfile=$fastq
    fi
  fi  
  
  printf '\n%s %s %s\n' $GREEN "Running umierrorcorrect for fastq: $outfile" $NC

  if $use_bed
    then
      run_umierrorcorrect.py \
      -o $runDir/$sample_name \
      -r1 $runDir/$outfile \
      -r $REF \
      -mode single \
      -ul $umi_length \
      -sl $spacer_length \
      -bed $BED \
      -t $threads
  else
    run_umierrorcorrect.py \
      -o $runDir/$sample_name \
      -r1 $runDir/$outfile \
      -r $REF \
      -mode single \
      -ul $umi_length \
      -sl $spacer_length \
      -t $threads
    fi
done

###############################
# Merging reports and cleanup #
###############################

# Generate fastqc reports
if $fastqc
  then
  printf '\n%s %s %s\n' $GREEN "Running fastqc for folder: $runDir" $NC
  fastqc $runDir/*.fastq.gz
fi  

# Move fastqc files to dedicated folder
mkdir "$runDir/qc_reports"
mv $runDir/*fastqc* qc_reports

# Generate merged reports
if $multiqc
  then
  printf '\n%s %s %s\n' $GREEN "Running multiqc for folder: $runDir" $NC
  multiqc $runDir
fi  

