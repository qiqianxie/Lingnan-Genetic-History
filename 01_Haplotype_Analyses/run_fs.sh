#!/bin/bash
# ==============================================================================
# Pipeline: run_fs.sh
# Author: Qiqian Xie
# Purpose: Automated, parallelized 4-stage pipeline for fineSTRUCTURE analysis
#          (Parameter estimation, ChromoPainter painting, MCMC, and Tree building).
# ==============================================================================
set -euo pipefail

WORK_DIR="./phased_data"
cd "${WORK_DIR}"

echo "========== Starting fineSTRUCTURE High-Performance Pipeline =========="

# Clean previous runs
rm -rf WGS_Project
rm -f WGS_Project.cp WGS_Project_*.out

# Assemble file lists for all 22 autosomes
PHASE_FILES=""
RECOMB_FILES=""
for i in {1..22}; do
    PHASE_FILES="${PHASE_FILES} phased_chr${i}.phase"
    RECOMB_FILES="${RECOMB_FILES} maps/chr${i}_recomb.txt"
done

# -------------------------------------------------------------
# Stage 1: Parameter Estimation (-s1indfrac 0.1)
# -------------------------------------------------------------
echo "[Stage 1] Generating parameter estimation commands..."
fs WGS_Project.cp -hpc 1 -idfile WGS.ids -phasefiles ${PHASE_FILES} -recombfiles ${RECOMB_FILES} -s1indfrac 0.1 -go

echo "[Stage 1] Executing parallel parameter estimation..."
cat WGS_Project/commandfiles/commandfile1.txt | parallel --halt now,fail=1 --progress -j 40
fs WGS_Project.cp -go

# -------------------------------------------------------------
# Stage 2: ChromoPainter Genome-Wide Chunk Painting
# -------------------------------------------------------------
echo "[Stage 2] Generating ChromoPainter painting commands..."
fs WGS_Project.cp -go 

echo "[Stage 2] Executing chunk painting across samples..."
cat WGS_Project/commandfiles/commandfile2.txt | parallel --halt now,fail=1 --progress -j 40

echo "[Stage 2] Merging painting results..."
fs WGS_Project.cp -go

# -------------------------------------------------------------
# Stage 3: MCMC Clustering (1,000,000 iterations)
# -------------------------------------------------------------
echo "[Stage 3] Initiating MCMC clustering..."
fs WGS_Project.cp -s3iters 1000000 -go

echo "[Stage 3] Sampling MCMC until convergence..."
while [ ! -s WGS_Project/commandfiles/commandfile4.txt ]; do
    if [ -f "WGS_Project/commandfiles/commandfile3.txt" ]; then
        cat WGS_Project/commandfiles/commandfile3.txt | parallel -j 4
    fi
    fs WGS_Project.cp -go
done
echo "Stage 3 MCMC sampling converged successfully."

# -------------------------------------------------------------
# Stage 4: Construct Maximum A Posteriori (MAP) Tree
# -------------------------------------------------------------
echo "[Stage 4] Building phylogenetic tree..."
cat WGS_Project/commandfiles/commandfile4.txt | parallel -j 4
fs WGS_Project.cp -go

echo "========================================================================"
echo "fineSTRUCTURE server-side pipeline completed successfully!"
echo "Core outputs generated:"
echo "  - Co-ancestry Matrix: WGS_Project_linked.chunkcounts.out"
echo "  - MCMC State:         WGS_Project_linked_mcmc.xml"
echo "  - Tree Topology:      WGS_Project_linked_tree.xml"
echo "========================================================================"