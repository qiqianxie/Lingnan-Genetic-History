#!/bin/bash
# ==============================================================================
# Original Author: joanam Joana
# Author: Qiqian Xie (To adapt this script for the 4-population pipeline)
# Purpose: Generate block-bootstrapped SFS datasets and estimate confidence 
#          intervals for the best demographic model using fastsimcoal2.
# ==============================================================================

PREFIX="Model_4pops"
VCF="WGS_1KG_HO_4pops.vcf"
POP_MAP="pop_map_4pops.txt"
PROJ="18,20,20,56"
N_BLOCKS=100

# ------------------------------------------------------------------------------
# 1. Get header and genomic data, calculate block size, and split
# ------------------------------------------------------------------------------
echo "Extracting sites and calculating block size..."
grep -v "^#" ${VCF} > ${PREFIX}.allSites
grep "^#" ${VCF} > header

# Dynamically calculate block size based on actual SNP count
TOTAL_SITES=$(wc -l < ${PREFIX}.allSites)
BLOCK_SIZE=$(( TOTAL_SITES / N_BLOCKS + 1 ))

echo "Total sites: ${TOTAL_SITES}. Splitting into ${N_BLOCKS} blocks of ~${BLOCK_SIZE} sites."

# Split sites into exactly 100 files
split -l ${BLOCK_SIZE} ${PREFIX}.allSites ${PREFIX}.sites.

# ------------------------------------------------------------------------------
# 2. Generate 100 files with randomly concatenated blocks and compute the SFS
# ------------------------------------------------------------------------------
for bs in {1..100}
do
    echo "Generating bootstrap dataset ${bs}..."
    mkdir bs_${bs}
    cd bs_${bs}

    # Add the header to our new bootstrapped vcf file
    cat ../header > ${PREFIX}.bs.${bs}.vcf
    
    # Randomly add 100 blocks (sampling with replacement)
    for r in {1..100}
    do
        cat `shuf -n1 -e ../${PREFIX}.sites.*` >> ${PREFIX}.bs.${bs}.vcf
    done
    
    # Compress the vcf file
    gzip ${PREFIX}.bs.${bs}.vcf

    # Make an SFS from the new bootstrapped file
    easySFS.py -i ${PREFIX}.bs.${bs}.vcf.gz -p ../${POP_MAP} -a -f --proj ${PROJ} -o easySFS_out > /dev/null 2>&1

    # Copy the observed SFS file into this folder and rename it to match the .tpl prefix
    cp easySFS_out/fastsimcoal2/*_MSFS.obs ./${PREFIX}_MSFS.obs
    cd ..
done

# ------------------------------------------------------------------------------
# 3. Run parameter estimation under the best model for each bootstrapped SFS
# ------------------------------------------------------------------------------
for bs in {1..100}
do
    echo "Running optimization for bootstrap ${bs}..."
    cd bs_${bs}
    
    # Run fastsimcoal2 50 times for each bootstrap dataset:
    for run in {1..50}
    do
        mkdir run${run}
        cd run${run}
        
        # Copy template, estimation, and observation files
        cp ../../${PREFIX}.tpl ./
        cp ../../${PREFIX}.est ./
        cp ../${PREFIX}_MSFS.obs ./
        
        # Execute fastsimcoal2
        fsc2709 -t ${PREFIX}.tpl -e ${PREFIX}.est -m -0 -C 10 -n 100000 -L 50 -M -q --multiSFS > /dev/null 2>&1
        
        cd ..
    done
    
    # Find the best run among the 50 iterations:
    cp ../fsc-selectbestrun.sh ./
    bash fsc-selectbestrun.sh > /dev/null 2>&1

    cd ..
done

echo "Block-bootstrapping and optimization completed!"