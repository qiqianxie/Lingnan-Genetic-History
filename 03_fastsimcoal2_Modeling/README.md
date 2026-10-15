# fastsimcoal2 Demographic Modeling Pipeline

This directory contains the core configurations used to test alternative population histories, estimate divergence times, and infer effective population sizes using **fastsimcoal2** (fsc2).

### Workflow Summary

1. **SFS Generation**: The site frequency spectrum (SFS) was generated using `easySFS.py` to handle down-sampling and maximize segregating sites across target populations.
2. **Model Construction**: Various demographic models were constructed via `.tpl` (Topology) and `.est` (Estimation bounds) files to test hypotheses at different time depths:
   * `Model_A`: Evaluates the deep ancestral divergence scenario among Hakka, Guangdong Yao (GDY), and CHB populations.
   * `Model_4pops`: Represents the primary hypothesis for recent admixture and merging patterns among Lingnan populations, where Cantonese and Teochew merge first, followed by Hakka.
   * `Model_4pops_Alt`: An alternative topology explicitly testing a different merging pattern (Hakka and Cantonese merge first, joined by Teochew) to rigorously compare against the primary hypothesis.
3. **Optimization**: For each model, 50 independent runs with 100,000 coalescent simulations per cycle (`-n 100000 -L 50`) were executed to find the global maximum likelihood.
4. **Model Selection**: Akaike Information Criterion (AIC) was calculated for each optimal run to determine the best-fit scenario.
5. **Block-Bootstrapping**: To assess the robustness of the parameter estimates for the optimal model, we performed 100 replicates of block-bootstrapping. The variants were divided into 100 blocks to generate pseudo-replicate SFS datasets, and parameters were re-estimated for each replicate to calculate the 95% Confidence Intervals (CIs).

### AIC Calculation & Model Selection

We used the following standard R function to compare competing demographic topologies based on the maximum likelihoods obtained:

```R
# Extract number of parameters (k) from .est file and calculate AIC
AIC <- 2*k - 2*(MaxEstLhood / log10(exp(1)))
deltaL <- MaxObsLhood - MaxEstLhood
```

### Bootstrap 95% CI Calculation

After extracting the optimal parameters from the 100 bootstrap replicates, the 95% CIs were calculated using `R`:

```R
library(tidyverse)
# Load bootstrap summary
bs_data <- read.table("Model_4pops_bootstraps_summary_fixed.txt", header = TRUE)

# Calculate Median and 95% CI for divergence time parameters (TDIV)
ci_results <- bs_data %>%
  select(contains("TDIV")) %>% 
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Estimate") %>%
  group_by(Parameter) %>%
  summarise(
    Median = median(Estimate),
    Lower_95_CI = quantile(Estimate, probs = 0.025),
    Upper_95_CI = quantile(Estimate, probs = 0.975)
  )
print(ci_results)
```

### Acknowledgements / Third-Party Scripts

We gratefully acknowledge the community resources and scripts that facilitated this analysis:
*   **easySFS**: Used for converting VCF files into the multidimensional SFS format required by fsc2 ([isaacovercast/easySFS](https://github.com/isaacovercast/easySFS)).
*   **Utility Scripts**: The bash scripts for extracting the best run among multiple replicates (`fsc-selectbestrun.sh`) and the R scripts for visualizing demographic topologies (`plotModel.r`) were adapted from the excellent tutorials and repositories provided by **Joana Meier** (e.g., [Speciation & Population Genomics: a how-to-guide](https://speciationgenomics.github.io/fastsimcoal2/)).

### Files Included in this Repository

* `Model_A.tpl` & `.est`: Configuration for deep ancestral divergence models.
* `Model_4pops.tpl` & `.est`: Configuration for the primary sub-population topology.
* `Model_4pops_Alt.tpl` & `.est`: Configuration for the alternative topology.
* `run_block_bootstrap.sh`: The automated bash pipeline used to generate block-bootstrap SFS and perform fsc2 parameter optimization concurrently.