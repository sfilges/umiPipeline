#!/bin/bash

# A wrapper script for running umierrorcorrect on multiple samples in the same
# directory. Requires '.fastq.gz' files as input and generates output folders in
# the same directory as input files. Sample names are fastq names without
# the file ending '.fastq.gz'. Optionally an assay_regions.bed file can be provided
# for annotation. The script handles both single-end and paired-end data, as well
# as optional paired-read-merging and quality filtering using fastp.

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
    echo "   -u  --umi_length     UMI length, default is 12."
    echo "   -s  --spacer_length  How long is the spacer sequence? Default is 16."
    echo "   -t  --threads        How many threads to use? Default is 16."
    echo "   -f  --no_filtering   Do not use fastp to filter fastqs."
    echo "   -q  --phred_score          Min Phread score to keep when using fastp to filter. Default is 15, typical values are 10, 15, 20, 30."
    echo "   -p  --percent_low_quality  How many percent bases in a read are allowed to be below the thrshold q value. Default is 40 (0-100)"
    echo "   -e  --no-paired-end        Should paired-end reads be used. Default is yes, if the flag is set, only R1 will be used."       
    echo
    exit 1
}

##############################
#   Initialize parameters    #
##############################

paired_end=true
umi_length=19
spacer_length=16
threads=16
no_fastp=false
use_bed=true
filtering=true
phred_score=20
percent_low_quality=40
fastqc=false
multiqc=false

# Initialize log file
touch log.txt

##############################
# Check command line options #
##############################

while getopts ':hfi:b:r:u:s:t:q:p:e:' option; do
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
        filtering=false
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
    e | --no-paired-end)
        paired_end=false
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

# Define colors for console output
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
#   Run umierrorcorrect   #
###########################

printf '\n%s %s %s\n' $GREEN "Processing $n_files fastq files..." $NC
printf '%s %s %s %s\n' $YELLOW "...using UMI length:" $NC $umi_length
printf '%s %s %s %s\n' $YELLOW "...using spacer length:" $NC $spacer_length
printf '%s %s %s %s\n' $YELLOW "...using threads:" $NC $threads
printf '%s %s %s %s\n' $YELLOW "...perform filtering:" $NC $do_filtering
printf '%s %s %s %s\n' $YELLOW "...using bed annotations:" $NC $use_bed

# Move to run directory
cd $runDir

# Processing fastq files
for fastq in *.fastq.gz ;
do
  # If fastq file name contains "R1"
  if [[ $fastq =~ R1 ]]
    then
    if $paired_end
      then
      # define read 1
      fq1=$fastq 
      printf '%s \n' $fq1

      # define corresponding read 2 by replacing R1 with R2
      fq2=${fastq//R1/R2} 
      printf '%s \n' $fq2
    else
      # use only read 1
      fq1=$fastq
      printf '%s \n' $fq1
    fi 
  else
    # If fastq file name does not contain "R1" continue with the next file
    continue
  fi
  
  # Define replacement
  find="_R1_001.fastq.gz"
  replace=""

  # Simplify sample name
  sample_name=${fq1//$find/$replace}

  printf "\n%s %s %s %s\n" $GREEN "Changing sample name:" $NC "$fq1 => $sample_name" 

  # Run filtering and read merging
  if $no_fastp
    then
      # If fastp is not installed, print note to console
      printf "\n%s\n" $YELLOW "Fastp is not installed."

      if $paired_end
        then
        # If data is paired-end
        printf "%s\n" $GREEN "Running UMIErrorCorrect in paired-end mode."

        if $use_bed
          then
            printf "%s\n" $GREEN "Using bed file."

            run_umierrorcorrect.py \
              -o $runDir/$sample_name \
              -r1 $runDir/$fq1 \
              -r2 $runDir/$fq2 \
              -r $REF \
              -mode paired \
              -ul $umi_length \
              -sl $spacer_length \
              -bed $BED \
              -t $threads
          else
            printf "%s\n" $YELLOW "NOT using bed file. This is not recommended."

            run_umierrorcorrect.py \
              -o $runDir/$sample_name \
              -r1 $runDir/$fq1 \
              -r2 $runDir/$fq2 \
              -r $REF \
              -mode paired \
              -ul $umi_length \
              -sl $spacer_length \
              -t $threads
          fi
        else
          # If data is not paired-end, use only R1
          printf "%s\n" $GREEN "Running UMIErrorCorrect in single-end mode."

          if $use_bed
            then
              printf "%s\n" $GREEN "Using bed file."

              run_umierrorcorrect.py \
                -o $runDir/$sample_name \
                -r1 $runDir/$fq1 \
                -r $REF \
                -mode single \
                -ul $umi_length \
                -sl $spacer_length \
                -bed $BED \
                -t $threads
          else
            printf "%s\n" $YELLOW "NOT using bed file. This is not recommended."

            run_umierrorcorrect.py \
              -o $runDir/$sample_name \
              -r1 $runDir/$fq1 \
              -r $REF \
              -mode single \
              -ul $umi_length \
              -sl $spacer_length \
              -t $threads
          fi
        fi
  else
    if $filtering
      then
        if $paired_end
        then
          printf '\n%s %s %s %s\n' $GREEN "Running fastp in paired-end mode."
          printf '%s %s %s %s\n' $GREEN "Using read 1..." $NC $fq1
          printf '%s %s %s %s\n' $GREEN "Using read 2..." $NC $fq2
          printf '%s %s %s %s\n' $YELLOW "...using minimum Phread score:" $NC $phred_score 
          printf '%s %s %s %s\n' $YELLOW "...using max percent low quality reads:" $NC $percent_low_quality

          # Define merged and filtered output file
          outfile="$sample_name.merged.filtered.fastq.gz"

          # Run fastp
          fastp \
            --in1=$fq1 \
            --in2=$fq2 \
            --merge \
            --unpaired1="${sample_name}_unpaired.fastq.gz" \
            --unpaired2="${sample_name}_unpaired.fastq.gz" \
            --merged_out=$outfile \
            --failed_out="${sample_name}_failed.fastq.gz" \
            --html "${sample_name}.html" \
            --trim_poly_g \
            --trim_poly_x \
            --qualified_quality_phred=$phred_score \
            --unqualified_percent_limit=$percent_low_quality \
            --detect_adapter_for_pe \
            --thread=$threads \
            --correction \
            --n_base_limit=3 \
            --overlap_len_require=30  \
            --length_required=100 \
            --json="${sample_name}.json" \
            --html="${sample_name}.html" \
            --report_title="${sample_name}"

          else 
            printf '\n%s %s %s %s\n' $GREEN "Running fastp in single-end mode."
            printf '%s %s %s %s\n' $GREEN "Using read 1..." $NC $fq1
            printf '%s %s %s %s\n' $YELLOW "...using minimum Phread score:" $NC $phred_score 
            printf '%s %s %s %s\n' $YELLOW "...using max percent low quality reads:" $NC $percent_low_quality

            # Define filtered output file
            outfile="$sample_name.filtered.fastq.gz"

            fastp \
            --in1=$fq1 \
            --out1=$outfile \
            --failed_out="${sample_name}_failed.fastq.gz" \
            --html "${sample_name}.html" \
            --trim_poly_g \
            --trim_poly_x \
            --qualified_quality_phred=$phred_score \
            --unqualified_percent_limit=$percent_low_quality \
            --detect_adapter_for_pe \
            --thread=$threads \
            --correction \
            --n_base_limit=3 \
            --length_required=100 \
            --json="${sample_name}.json" \
            --html="${sample_name}.html" \
            --report_title="${sample_name}"
          fi

          printf '\n%s %s %s\n' $GREEN "Running umierrorcorrect for fastq: $outfile" $NC
          printf "%s\n" $GREEN "Running UMIErrorCorrect in single-end mode."
  
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
        else
          # If fastp is not used, print note to console
          printf "\n%s\n" $YELLOW "Fastp ist not used for filtering."

          # If data is paired-end
          if $paired_end
            then
            printf "%s\n" $GREEN "Running UMIErrorCorrect in paired-end mode."

            if $use_bed
              then
                printf "%s\n" $GREEN "Using bed file."

                run_umierrorcorrect.py \
                  -o $runDir/$sample_name \
                  -r1 $runDir/$fq1 \
                  -r2 $runDir/$fq2 \
                  -r $REF \
                  -mode paired \
                  -ul $umi_length \
                  -sl $spacer_length \
                  -bed $BED \
                  -t $threads
              else
                printf "%s\n" $YELLOW "NOT using bed file. This is not recommended."

                run_umierrorcorrect.py \
                  -o $runDir/$sample_name \
                  -r1 $runDir/$fq1 \
                  -r2 $runDir/$fq2 \
                  -r $REF \
                  -mode paired \
                  -ul $umi_length \
                  -sl $spacer_length \
                  -t $threads
              fi
            else
              # If data is not paired-end, use only R1
              printf "%s\n" $GREEN "Running UMIErrorCorrect in single-end mode."

              if $use_bed
                then
                  printf "%s\n" $GREEN "Using bed file."

                  run_umierrorcorrect.py \
                    -o $runDir/$sample_name \
                    -r1 $runDir/$fq1 \
                    -r $REF \
                    -mode single \
                    -ul $umi_length \
                    -sl $spacer_length \
                    -bed $BED \
                    -t $threads
              else
                printf "%s\n" $YELLOW "NOT using bed file. This is not recommended."

                run_umierrorcorrect.py \
                  -o $runDir/$sample_name \
                  -r1 $runDir/$fq1 \
                  -r $REF \
                  -mode single \
                  -ul $umi_length \
                  -sl $spacer_length \
                  -t $threads
              fi
            fi


      fi
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

