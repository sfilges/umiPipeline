## umiPipeline
 Pipeline for processing multiple fastq files and running UMIerrorcorrect.
 
 ### Dependencies
 
 - UMIErrorCorrect
 - bwa
 - (optional) fastp
 
The tools can be installed from pip or conda, respectively:
 
```
 conda install -c bioconda bwa
 conda install -c bioconda fastp
 pip install umierrorcorrect
```
 
 ### Requirements

- one or more fastq files
- reference genome indexed with bwa
- (recommended) bed file containing amplicon regions
 
 ### Usage
 
For basic usage only a directory containing fastq files and an indexed reference genome need o be supplied. A bed file containing amplicon annotations is optional, but recommended.

```
umi_pipeline.sh -i <path_to_fastq_dir> -r <reference_fasta> -b <optional_bed_file>
```

If fastp is installed it will automatically be used to trim adapters and perform quality filtering. An html report will be generated for each input fastq.

A number of optional flags can be set to modify the behaviour of fastp:
```
 -f  Do not perform filtering
 -q  Phred score threshold for filtering, default is 15.
 -p  Percent of reads that are allowed to have poor quality. Default is 40 (0-100).
```

For each fastq/sample a folder will be generated into which umierrorcorrect outputs will be redirected. Similarl to fastp, the arguments passed to UMIerrorocorrect can also be adjusted:
```
 -u <UMI length> Integer, default is 12.
 -s <Spacer length> Integer, default is 16.
 -t <Threads> Integer, default is 16.
```


