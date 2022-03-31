# umiPipeline
 Pipeline for processing multiple fastq files and running UMIerrorcorrect.
 
 ## Dependencies
 
 - UMIErrorCorrect
 - bwa
 - (optional) fastp
 
 
 ## Usage
 
```
umi_pipeline.sh -i <path_to_fastq_dir> -r <reference_fasta> -b <optional_bed_file>
```