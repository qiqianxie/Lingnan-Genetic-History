#!/usr/bin/env python3
"""
Purpose: To convert phased VCF files (e.g., from Beagle) into the '.phase' format required by ChromoPainter and fineSTRUCTURE for haplotype-based analysis.

Input:
    1. A phased VCF file (gzipped) - typically containing genotypes separated by '|' or '/'.
    2. Output file path for the '.phase' format.

Output:
    - A '.phase' format file containing:
        - Line 1: Number of individuals
        - Line 2: Number of SNPs
        - Line 3: 'P' followed by physical positions of the SNPs
        - Subsequent lines: Continuous haplotype strings (e.g., '010101...') for each sample.
"""

import sys
import gzip

if len(sys.argv) != 3:
    print("Usage: python vcf2phase.py <input.vcf.gz> <output.phase>")
    sys.exit(1)

vcf_file = sys.argv[1]
phase_file = sys.argv[2]

positions = []
haps = []

print(f"Reading {vcf_file} ...")
with gzip.open(vcf_file, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            parts = line.strip().split('\t')
            N = len(parts) - 9                # Number of samples
            haps = [[] for _ in range(2 * N)] # 2 haplotypes per sample
            continue
        
        parts = line.strip().split('\t')
        positions.append(parts[1])            # Extract physical position
        
        # Extract genotypes starting from column 10 (e.g., 0|1)
        for i in range(N):
            gt = parts[9 + i].split(':')[0]
            if len(gt) >= 3:
                a1, a2 = gt[0], gt[2]
            else:
                a1, a2 = '0', '0'
            
            haps[2*i].append(a1)
            haps[2*i+1].append(a2)

print(f"Writing to {phase_file} ...")
with open(phase_file, 'w') as out:
    out.write(f"{N}\n")                  # Line 1: Number of individuals
    out.write(f"{len(positions)}\n")     # Line 2: Number of SNPs
    out.write("P " + " ".join(positions) + "\n") # Line 3: Positions
    for hap in haps:
        out.write("".join(hap) + "\n")   # Haplotype strings

print("Conversion completed successfully.")