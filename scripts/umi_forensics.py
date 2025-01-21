#!/usr/bin/env python3

import argparse
import os
import subprocess
from pathlib import Path
import re

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

def main():
    parser = argparse.ArgumentParser(description='Run UMI error correction on FASTQ files')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory containing FASTQ files')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-g', '--genome', required=True, help='Reference genome path')
    parser.add_argument('-n', '--ini', required=True, help='INI file')
    parser.add_argument('-l', '--library', required=True, help='Library file')
    parser.add_argument('-b', '--bed', required=True, help='BED file with regions')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Find all FASTQ files and their pairs
    fastq_pairs = find_fastq_pairs(args.input_dir)
    
    for sample, files in fastq_pairs.items():
        cmd = ['run_umierrorcorrect_forensics.py']
        
        if files['R1'] is None:
            print(f"Warning: No R1 file found for sample {sample}, skipping...")
            continue
            
        cmd.extend(['-r1', files['R1']])
        
        # If R2 exists, add it to command with -p flag
        if files['R2']:
            cmd.extend(['-r2', files['R2'], '-p'])
        
        cmd.extend([
            '-o', args.output,
            '-g', args.genome,
            '-i', args.ini,
            '-l', args.library,
            '-b', args.bed
        ])
        
        print(f"Processing sample: {sample}")
        print(f"Running command: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing sample {sample}: {e}")

if __name__ == '__main__':
    main()