# Haplotype-Based Analyses (fineSTRUCTURE and IBD)

This directory contains the pipeline, custom scripts, and execution wrappers for the haplotype-based analyses described in **Methods 5.4**, including haplotype phasing, fineSTRUCTURE clustering, and high-resolution Identity-by-Descent (IBD) segment inference.

---

### Workflow Overview

#### 1. Genotype Phasing & Map Preparation
* The genome-wide VCF was split chromosome-by-chromosome using `bcftools`:
  ```bash
  for i in {1..22}; do
      bcftools view -r ${i} -O z -o chr${i}.vcf.gz input_genome.vcf.gz
  done
  ```
* Autosomal biallelic SNPs were phased using **Beagle v5.4** against the GRCh38 genetic recombination map:
  ```bash
  java -Xmx32g -jar beagle.jar gt=chr${i}.vcf.gz map=plink.chr${i}.GRCh38.map out=phased_chr${i} burnin=6 iterations=12
  ```
* Recombination maps were converted to ChromoPainter-compatible format (cM/Mb) by calculating local rates between genomic coordinates:
  ```bash
  for i in {1..22}; do
      echo -e "position\trate\tcM" > chr${i}_recomb.txt
      awk 'BEGIN{OFS="\t"} 
           {
             if(NR==1){ prev_cM=$3; prev_bp=$4; print $4, 0, $3 }
             else {
               delta_bp = $4 - prev_bp;
               delta_cM = $3 - prev_cM;
               if(delta_bp > 0) { rate = delta_cM / (delta_bp / 1000000) } else { rate = 0 };
               if(rate < 0) rate = 0;
               print $4, rate, $3;
               prev_bp = $4; prev_cM = $3;
             }
           }' plink.chr${i}.GRCh38.map >> chr${i}_recomb.txt
  done
  ```

#### 2. Format Conversion for ChromoPainter
* The phased VCFs were converted into ChromoPainter `.phase` input format using our custom matrix-transposition script:
  ```bash
  python vcf2phase.py phased_chr${i}.vcf.gz phased_chr${i}.phase
  ```

#### 3. fineSTRUCTURE Pipeline (`run_fs.sh`)
* The complete pipeline runs ChromoPainter and fineSTRUCTURE across four consecutive stages using GNU parallel for high-performance computing:
  * **Stage 1**: Parameter estimation (`-s1indfrac 0.1`).
  * **Stage 2**: Full genome ChromoPainter co-ancestry chunk painting.
  * **Stage 3**: MCMC clustering convergence (`-s3iters 1000000`).
  * **Stage 4**: Maximum a posteriori (MAP) phylogenetic tree construction.

#### 4. High-Resolution IBD Segment Extraction (`extract_ibd.py`)
* Pairwise IBD segments were detected using **hap-ibd** (`min-output=1 cM`).
* Segments were subsequently classified into length categories (**1–2 cM**, **2–5 cM**, and **>5 cM**) using `extract_ibd.py` to evaluate deep versus recent population sharing and validate WGS resolution.

---

### Files in this Directory

* `vcf2phase.py`: Custom Python script to transpose and format phased VCF genotypes into ChromoPainter phase files.
* `run_fs.sh`: Multi-threaded automation wrapper executing all four stages of the fineSTRUCTURE analysis.
* `extract_ibd.py`: Custom Python script to aggregate pairwise IBD blocks by length bins and compute mean shared chunk statistics.