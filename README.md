# Lingnan Genetic History

This repository contains the custom scripts, demographic model configurations, and parameter files used for the analyses in our study regarding the genomic diversity and population history of Lingnan populations.

Our core analytical pipelines rely heavily on mature, widely benchmarked population genetics frameworks. For standard operations (e.g., variant calling, phasing, PCA, ADMIXTURE, f-statistics, and hap-ibd), we followed standard community protocols, such as those documented in the [2019 Workshop on Computational Population Genetics](https://comppopgenworkshop2019.readthedocs.io/).

Here, we provide the specific models and custom codes to ensure full reproducibility of our advanced demographic inferences and data processing steps.

### Dependencies
Most analyses were conducted using standard tools. Ensure the following software and packages are available in your environment:
*   **Population Genetics & Phasing Tools**:
    *   [Beagle v5.4](https://faculty.washington.edu/browning/beagle/beagle.html) & [Java (>= 1.8)](https://www.java.com/)
    *   [PLINK v2.0](https://www.cog-genomics.org/plink/2.0/)
    *   [bcftools (v1.14+)](https://www.htslib.org/download/)
    *   [ADMIXTOOLS v8.0.1](https://github.com/DReichLab/AdmixTools)
    *   [ADMIXTOOLS 2 (R package)](https://uqrmaie1.github.io/admixtools/)
    *   [fastsimcoal2 (v2.7/2.8)](http://cmpg.unibe.ch/software/fastsimcoal2/)
    *   [hap-ibd](https://github.com/browning-lab/hap-ibd)
    *   [fineSTRUCTURE v4 / ChromoPainter v2](https://people.maths.bris.ac.uk/~madjl/finestructure/)
*   **Utilities & Scripting Environments**:
    *   [GNU parallel](https://www.gnu.org/software/parallel/)
    *   Python (v3.x) with `pandas`
    *   R (v4.x) with `tidyverse` and `cowplot`

### Repository Structure

*   **`01_Haplotype_Analyses/`**: Contains the pipeline and custom scripts for haplotype-based analyses (phasing, fineSTRUCTURE clustering, and IBD segment parsing).
    *   `vcf2phase.py`: Converts phased VCFs into the `.phase` format.
    *   `run_fs.sh`: Parallel 4-Stage fineSTRUCTURE Pipeline (parameter estimation, painting, MCMC, and tree construction).
    *   `extract_ibd.py`: Parses and summarizes high-resolution IBD segments into specific length bins (1–2 cM, 2–5 cM, >5 cM) to demonstrate deep ancestral connections.
*   **`02_qpGraph_Modeling/`**: Contains the R scripts used for hypothesis-driven, stepwise admixture graph modeling using the `admixtools` R package (from the basic 5-population skeleton to the final comprehensive independent-admixture model).
*   **`03_fastsimcoal2_Modeling/`**: Contains the Site Frequency Spectrum (SFS) demographic model configuration files (`.tpl` and `.est`), testing both the primary merging hypothesis and alternative topologies. It also includes the automated bash pipeline for running coalescent simulations, parameter optimization, and block-bootstrapping.

### Usage

Please navigate to each sub-directory to find specific instructions and the exact parameter files used in our study. 
*   For **Haplotype Analyses**, step-by-step documentation for phasing, ChromoPainter painting, and IBD segment extraction is available in `01_Haplotype_Analyses/README.md`.
*   For **qpGraph** analyses, you can run the R scripts directly in an environment where `admixtools 2` is installed.
*   For **fastsimcoal2**, a detailed workflow, including how to calculate AIC and 95% CIs via bootstrapping, is provided in `03_fastsimcoal2_Modeling/README.md`.

### Data Availability
The analysis-ready genomic datasets generated for this study have been deposited in Zenodo and the temporary embargo has been entirely lifted. The data is now publicly accessible at: https://doi.org/10.5281/zenodo.18289776 (DOI: 10.5281/zenodo.18289776)

### Citation
Our manuscript is currently under review. If you use the scripts, models, or data formatting pipelines from this repository, please cite our paper once published:
> **Fine-Scale Genomic Structure Unveils the Genetic History and Admixture of Han Subgroups and Indigenous Yao in Lingnan**
> *Qiqian Xie, et al.*
> *(Currently under review)*

*(Note: The citation details will be updated with the DOI and journal information upon formal publication).*

---
**Contact:** Qiqian Xie - qiqian_xie@163.com