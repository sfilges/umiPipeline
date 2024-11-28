#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import multiprocessing
from pathlib import Path

def run_fastp(fq1, fq2, outfile, phred_score, percent_low_quality, threads):
    """Run fastp for quality filtering and merging."""
    fastp_cmd = [
        'fastp', 
        f'--in1={fq1}', 
        f'--in2={fq2}', 
        '--merge', 
        f'--merged_out={outfile}',
        f'--qualified_quality_phred={phred_score}',
        f'--unqualified_percent_limit={percent_low_quality}',
        f'--thread={threads // 2}'
    ]
    try:
        subprocess.run(fastp_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: fastp failed for {fq1} and {fq2}")
        raise e

def run_umierrorcorrect(outfile, sample_name, ref, umi_length, spacer_length, bed_file, threads):
    """Run UMI error correction."""
    umi_cmd = [
        'run_umierrorcorrect.py', 
        '-o', f'{sample_name}', 
        '-r1', outfile, 
        '-r', ref,
        '-mode', 'single', 
        f'-ul={umi_length}', 
        f'-sl={spacer_length}',
        f'-t={threads}'
    ]
    
    if bed_file:
        umi_cmd.append(f'-bed={bed_file}')
    
    try:
        subprocess.run(umi_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: umierrorcorrect failed for {sample_name}")
        raise e

def process_fastq_pair(args):
    """Process a pair of FASTQ files."""
    fq1, filtering, filtered_dir, umi_corrected_dir, ref, umi_length, spacer_length, bed_file, phred_score, percent_low_quality, threads = args
    
    fq2 = fq1.replace('R1', 'R2')
    sample_name = os.path.basename(fq1).replace('.fastq.gz', '')
    outfile = os.path.join(filtered_dir, f'{sample_name}.merged.filtered.fastq.gz')

    if filtering:
        run_fastp(fq1, fq2, outfile, phred_score, percent_low_quality, threads)
    else:
        outfile = fq1

    run_umierrorcorrect(
        outfile, 
        os.path.join(umi_corrected_dir, sample_name), 
        ref, 
        umi_length, 
        spacer_length, 
        bed_file, 
        threads
    )

def main():
    parser = argparse.ArgumentParser(description='FASTQ Processing Script')
    parser.add_argument('-i', '--input-dir', default=os.getcwd(), help='Input directory')
    parser.add_argument('-b', '--bed', help='Path to BED file')
    parser.add_argument('-r', '--reference', required=True, help='Path to reference genome')
    parser.add_argument('-u', '--umi_length', type=int, default=19, help='UMI length')
    parser.add_argument('-s', '--spacer_length', type=int, default=16, help='Spacer length')
    parser.add_argument('-t', '--threads', type=int, default=max(1, os.cpu_count() // 2), help='Parallel threads')
    parser.add_argument('-f', '--no_filtering', action='store_true', help='Skip fastp filtering')
    parser.add_argument('-q', '--phred_score', type=int, default=20, help='Minimum Phred score')
    parser.add_argument('-p', '--percent_low_quality', type=int, default=40, help='Max low-quality bases')
    parser.add_argument('--skip-fastqc', action='store_true', help='Skip FastQC')
    parser.add_argument('--skip-multiqc', action='store_true', help='Skip MultiQC')

    args = parser.parse_args()

    # Create output directories
    filtered_dir = os.path.join(args.input_dir, 'filtered_fastqs')
    umi_corrected_dir = os.path.join(args.input_dir, 'umi_corrected_samples')
    qc_dir = os.path.join(args.input_dir, 'qc_reports')
    os.makedirs(filtered_dir, exist_ok=True)
    os.makedirs(umi_corrected_dir, exist_ok=True)
    os.makedirs(qc_dir, exist_ok=True)

    # Find R1 FASTQ files
    fastq_files = [f for f in Path(args.input_dir).glob('*R1*.fastq.gz')]

    # Prepare arguments for parallel processing
    process_args = [
        (str(fq), not args.no_filtering, filtered_dir, umi_corrected_dir, 
         args.reference, args.umi_length, args.spacer_length, 
         args.bed, args.phred_score, args.percent_low_quality, args.threads) 
        for fq in fastq_files
    ]

    # Parallel processing
    with multiprocessing.Pool(processes=args.threads) as pool:
        pool.map(process_fastq_pair, process_args)

    # Quality Control
    if not args.skip_fastqc:
        subprocess.run(['fastqc', '-o', qc_dir] + [str(f) for f in fastq_files])

    if not args.skip_multiqc:
        subprocess.run(['multiqc', args.input_dir, '-o', qc_dir])

if __name__ == '__main__':
    main()