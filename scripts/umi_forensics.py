#!/usr/bin/env python3

import argparse
import os
import subprocess
from pathlib import Path
import re
import multiprocessing
from functools import partial

def find_fastq_pairs(input_dir):
    """Find all FASTQ files and pair them based on R1/R2 naming convention."""
    fastq_files = {}
    
    # Recursively search for fastq/fastq.gz files
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                filepath = os.path.join(root, file)
                
                # Extract sample name by removing R1/R2 part
                sample = re.sub(r'_R[12]_\d+\.f(ast)?q(\.gz)?$', '', file)
                
                if sample not in fastq_files:
                    fastq_files[sample] = {'R1': None, 'R2': None}
                
                if '_R1_' in file:
                    fastq_files[sample]['R1'] = filepath
                elif '_R2_' in file:
                    fastq_files[sample]['R2'] = filepath
    
    return fastq_files

def process_sample(args_dict, sample, files):
    """Process a single sample with the given parameters."""
    cmd = ['run_umierrorcorrect_forensics.py']
    
    if files['R1'] is None:
        print(f"Warning: No R1 file found for sample {sample}, skipping...")
        return
        
    cmd.extend(['-r1', files['R1']])
    
    if files['R2']:
        cmd.extend(['-r2', files['R2'], '-p'])
    
    cmd.extend([
        '-o', args_dict['output'],
        '-g', args_dict['genome'],
        '-i', args_dict['ini'],
        '-l', args_dict['library'],
        '-b', args_dict['bed']
    ])
    
    print(f"Processing sample: {sample}")
    print(f"Running command: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error processing sample {sample}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Run UMI error correction on FASTQ files')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory containing FASTQ files')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-g', '--genome', required=True, help='Reference genome path')
    parser.add_argument('-n', '--ini', required=True, help='INI file')
    parser.add_argument('-l', '--library', required=True, help='Library file')
    parser.add_argument('-b', '--bed', required=True, help='BED file with regions')
    parser.add_argument('-t', '--threads', type=int, default=multiprocessing.cpu_count(),
                      help='Number of parallel processes (default: number of CPU cores)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Find all FASTQ files and their pairs
    fastq_pairs = find_fastq_pairs(args.input_dir)
    
    # Create a dictionary of arguments to pass to process_sample
    args_dict = {
        'output': args.output,
        'genome': args.genome,
        'ini': args.ini,
        'library': args.library,
        'bed': args.bed
    }
    
    # Create a partial function with the args_dict
    process_func = partial(process_sample, args_dict)
    
    # Run the processing in parallel
    with multiprocessing.Pool(processes=args.threads) as pool:
        pool.starmap(process_func, fastq_pairs.items())

if __name__ == '__main__':
    main()