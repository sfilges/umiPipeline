#!/usr/bin/env python3

import argparse
import logging
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Tuple, Optional

def setup_logging(run_dir: Path) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        handlers=[
            logging.FileHandler(run_dir / "script_log.txt"),
            logging.StreamHandler(sys.stdout)
        ]
    )

def check_dependencies() -> bool:
    required_tools = ['fastp', 'run_umierrorcorrect.py', 'bwa', 'multiqc', 'fastqc']
    missing = [tool for tool in required_tools if not shutil.which(tool)]
    if missing:
        logging.error(f"Missing dependencies: {', '.join(missing)}")
        return False
    logging.info("All dependencies are present.")
    return True

def check_paths(run_dir: Path, ref_path: Path, bed_path: Optional[Path] = None) -> bool:
    if not run_dir.is_dir():
        logging.error(f"Input directory does not exist: {run_dir}")
        return False
    if not ref_path.is_file():
        logging.error(f"Reference genome file not found: {ref_path}")
        return False
    if bed_path and not bed_path.is_file():
        logging.error(f"BED file not found: {bed_path}")
        return False
    return True

def run_fastp(fq1: Path, fq2: Path, output_dir: Path, sample_name: str,
              merge_reads: bool, threads: int, phred_score: int,
              percent_low_quality: int) -> Tuple[Path, Optional[Path]]:
    """Run fastp on input FASTQ files."""
    try:
        if merge_reads:
            outfile = output_dir / f"{sample_name}.merged.filtered.fastq.gz"
            cmd = [
                "fastp",
                "--in1", str(fq1),
                "--in2", str(fq2),
                "--merge",
                "--merged_out", str(outfile),
                "--qualified_quality_phred", str(phred_score),
                "--unqualified_percent_limit", str(percent_low_quality),
                "--trim_poly_g",
                "--trim_poly_x",
                "--length_required", "100",
                "--thread", str(threads // 2)
            ]
            subprocess.run(cmd, check=True)
            return outfile, None
        else:
            out1 = output_dir / f"{sample_name}.filtered.R1.fastq.gz"
            out2 = output_dir / f"{sample_name}.filtered.R2.fastq.gz"
            cmd = [
                "fastp",
                "--in1", str(fq1),
                "--in2", str(fq2),
                "--out1", str(out1),
                "--out2", str(out2),
                "--qualified_quality_phred", str(phred_score),
                "--unqualified_percent_limit", str(percent_low_quality),
                "--trim_poly_g",
                "--trim_poly_x",
                "--length_required", "100",
                "--thread", str(threads // 2)
            ]
            subprocess.run(cmd, check=True)
            return out1, out2
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running fastp for {fq1} and {fq2}: {str(e)}")
        raise

def run_umierrorcorrect(r1: Path, r2: Optional[Path], output_dir: Path,
                        sample_name: str, ref_path: Path, umi_length: int,
                        spacer_length: int, bed_path: Optional[Path],
                        threads: int) -> None:
    """Run umierrorcorrect on FASTQ files."""
    try:
        cmd = [
            "run_umierrorcorrect.py",
            "-o", str(output_dir / sample_name),
            "-r1", str(r1),
            "-r", str(ref_path),
            "-ul", str(umi_length),
            "-sl", str(spacer_length),
            "-t", str(threads)
        ]
        
        if r2:
            cmd.extend(["-r2", str(r2), "-mode", "paired"])
        else:
            cmd.extend(["-mode", "single"])
            
        if bed_path:
            cmd.extend(["-bed", str(bed_path)])
            
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running umierrorcorrect for {sample_name}: {str(e)}")
        raise

def process_fastq_pair(fq1: Path, args) -> None:
    """Process a pair of FASTQ files through the pipeline."""
    fq2 = Path(str(fq1).replace("R1", "R2"))
    sample_name = fq1.stem.replace(".fastq.gz", "")
    
    logging.info(f"Processing sample: {sample_name}")
    
    filtered_dir = args.input_dir / "filtered_fastqs"
    umi_corrected_dir = args.input_dir / "umi_corrected_samples"
    
    try:
        if args.do_filtering:
            r1, r2 = run_fastp(
                fq1, fq2, filtered_dir, sample_name,
                args.merge_reads, args.threads,
                args.phred_score, args.percent_low_quality
            )
        else:
            r1, r2 = fq1, fq2
            
        run_umierrorcorrect(
            r1, r2 if not args.merge_reads else None,
            umi_corrected_dir, sample_name, args.reference,
            args.umi_length, args.spacer_length,
            args.bed, args.threads
        )
    except Exception as e:
        logging.error(f"Error processing {sample_name}: {str(e)}")
        raise

def run_qc(run_dir: Path, threads: int, skip_fastqc: bool, skip_multiqc: bool) -> None:
    """Run FastQC and MultiQC on processed files."""
    qc_dir = run_dir / "qc_reports"
    qc_dir.mkdir(exist_ok=True)
    
    if not skip_fastqc and shutil.which('fastqc'):
        logging.info("Running FastQC...")
        fastq_files = list(run_dir.rglob("*.fastq.gz"))
        for fq in fastq_files:
            try:
                subprocess.run(["fastqc", "-o", str(qc_dir), str(fq)], check=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"Error running FastQC on {fq}: {str(e)}")
    
    if not skip_multiqc and shutil.which('multiqc'):
        logging.info("Running MultiQC...")
        try:
            subprocess.run(["multiqc", str(run_dir), "-o", str(qc_dir)], check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running MultiQC: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description="Process FASTQ files with UMI error correction")
    parser.add_argument("-i", "--input-dir", type=Path, default=Path.cwd(),
                        help="Input directory containing FASTQ files")
    parser.add_argument("-b", "--bed", type=Path,
                        help="Path to the assay regions BED file")
    parser.add_argument("-r", "--reference", type=Path, required=True,
                        help="Path to the indexed reference genome")
    parser.add_argument("-u", "--umi_length", type=int, default=19,
                        help="Length of the UMI")
    parser.add_argument("-s", "--spacer_length", type=int, default=16,
                        help="Spacer sequence length")
    parser.add_argument("-t", "--threads", type=int,
                        default=mp.cpu_count() // 2,
                        help="Number of parallel jobs")
    parser.add_argument("-f", "--no_filtering", action="store_false",
                        dest="do_filtering",
                        help="Skip fastp filtering step")
    parser.add_argument("-q", "--phred_score", type=int, default=20,
                        help="Minimum Phred score threshold")
    parser.add_argument("-p", "--percent_low_quality", type=int, default=40,
                        help="Maximum percentage of low-quality bases")
    parser.add_argument("--skip-fastqc", action="store_true",
                        help="Skip FastQC quality control")
    parser.add_argument("--skip-multiqc", action="store_true",
                        help="Skip MultiQC report generation")
    parser.add_argument("--merge", action="store_true",
                        dest="merge_reads",
                        help="Merge reads before umierrorcorrect")
    
    args = parser.parse_args()
    
    # Setup directories
    args.input_dir = args.input_dir.resolve()
    for dir_name in ["filtered_fastqs", "umi_corrected_samples", "qc_reports"]:
        (args.input_dir / dir_name).mkdir(exist_ok=True)
    
    # Setup logging
    setup_logging(args.input_dir)
    
    # Check dependencies and paths
    if not check_dependencies() or not check_paths(args.input_dir, args.reference, args.bed):
        sys.exit(1)
    
    # Process FASTQ files
    try:
        fastq_files = list(args.input_dir.glob("*R1*.fastq.gz"))
        with mp.Pool(args.threads) as pool:
            pool.starmap(process_fastq_pair,
                        [(fq, args) for fq in fastq_files])
    except Exception as e:
        logging.error(f"Error during FASTQ processing: {str(e)}")
        sys.exit(1)
    
    # Run QC
    run_qc(args.input_dir, args.threads, args.skip_fastqc, args.skip_multiqc)
    logging.info("Processing complete.")

if __name__ == "__main__":
    main()