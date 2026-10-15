#!/usr/bin/Rscript
# Qiqian Xie 2026
# Stepwise Admixture Graph Modeling

# Load libaries
library(admixtools)
library(tidyverse)
library(cowplot)

# Define dataset prefix and directory
f2_dir <- "./f2_blocks_data"
f2_blocks <- f2_from_precomp(f2_dir)

# Define core population components
pop_outgroup <- "Mbuti"
pop_basal    <- "Tianyuan"
pop_north    <- "Wuzhuangguoliang_LN" 
pop_coastal  <- "Fujian_Qihe"
pop_south    <- "Dushan"

target_pops  <- c("GDH_Cantonese", "GDH_Hakka", "GDH_Teochew", "GDY")


# Automated Heuristic Search (0 to 8 Admixture events)
# Set seed for exact reproducibility of heuristic searches
set.seed(202605) 
max_admix <- 8
best_graphs_list <- list()

for (k in 0:max_admix) {
  cat("\nRunning automated search for numadmix =", k, "...\n")
  res_k <- find_graphs(f2_all, outpop = pop_outgroup, numadmix = k, stop_gen = 200, stop_gen2 = 50)
  
  res_k_sorted <- res_k %>% arrange(score)
  best_graphs_list[[as.character(k)]] <- res_k_sorted$edges[[1]]
}
# Output trends of Likelihood scores and Max Z-scores were plotted to determine model saturation.


# Basic 5-Population Skeleton
pops_5 <- c(pop_outgroup, pop_basal, pop_north, pop_coastal, pop_south)
f2_5 <- f2_blocks[pops_5, pops_5, ]

graph_5_base <- tibble::tribble(
  ~from, ~to,
  "root", pop_outgroup,
  "root", "node_A",
  "node_A", pop_basal,
  "node_A", "node_EastAsia",
  "node_EastAsia", "anc_North",
  "node_EastAsia", "anc_South",
  "anc_North", pop_north,
  "anc_South", "split_South",
  "anc_South", "split_Coastal",
  "split_South", pop_south,
  "split_Coastal", pop_coastal
)

fit_5 <- qpgraph(f2_5, graph_5_base, return_fstats = TRUE)
cat("STEP 1 Base Skeleton Score:", round(fit_5$score, 2), "| Max f4 |Z|:", round(max(abs(fit_5$f4$z)), 3), "\n")


# Topology Comparison
for (target in target_pops) {
  pops_6 <- c(pops_5, target)
  f2_6 <- f2_blocks[pops_6, pops_6, ]
  
  # --- Model A: Coastal Donor ---
  graph_A <- tibble::tribble(
    ~from, ~to,
    "root", pop_outgroup,
    "root", "node_A",
    "node_A", pop_basal,
    "node_A", "node_EastAsia",
    "node_EastAsia", "anc_North",
    "node_EastAsia", "anc_South",
    "anc_North", pop_north,
    "anc_North", "target_N_edge",
    "anc_South", "split_South",
    "anc_South", "split_Coastal",
    "split_South", pop_south,
    "split_Coastal", pop_coastal,
    "split_Coastal", "target_S_Coastal_edge",
    "target_N_edge", target,
    "target_S_Coastal_edge", target
  )
  fit_A <- qpgraph(f2_6, graph_A, return_fstats = TRUE)
  
  # --- Model B: Inland Donor ---
  graph_B <- tibble::tribble(
    ~from, ~to,
    "root", pop_outgroup,
    "root", "node_A",
    "node_A", pop_basal,
    "node_A", "node_EastAsia",
    "node_EastAsia", "anc_North",
    "node_EastAsia", "anc_South",
    "anc_North", pop_north,
    "anc_North", "target_N_edge",
    "anc_South", "split_South",
    "anc_South", "split_Coastal",
    "split_Coastal", pop_coastal,
    "split_South", pop_south,
    "split_South", "target_S_Inland_edge",
    "target_N_edge", target,
    "target_S_Inland_edge", target
  )
  fit_B <- qpgraph(f2_6, graph_B, return_fstats = TRUE)
}


# 3-Way Admixture Model
for (target in target_pops) {
  pops_6 <- c(pops_5, target)
  f2_6 <- f2_blocks[pops_6, pops_6, ]
  
  graph_3way <- tibble::tribble(
    ~from, ~to,
    "root", pop_outgroup,
    "root", "node_A",
    "node_A", pop_basal,
    "node_A", "node_EastAsia",
    "node_EastAsia", "anc_North",
    "node_EastAsia", "anc_South",
    "anc_North", pop_north,
    "anc_North", "target_N_edge",
    "anc_South", "split_South",
    "anc_South", "split_Coastal",
    "split_South", pop_south,
    "split_South", "target_S_Inland_edge",
    "split_Coastal", pop_coastal,
    "split_Coastal", "target_S_Coastal_edge",
    "target_S_Coastal_edge", "target_Mix_1",
    "target_S_Inland_edge", "target_Mix_1",
    "target_Mix_1", target,
    "target_N_edge", target
  )
  
  fit_3way <- qpgraph(f2_6, graph_3way, return_fstats = TRUE)
}


# Final Comprehensive Independent-Admixture Model
pops_all <- c(pop_outgroup, pop_basal, pop_north, pop_coastal, pop_south, target_pops)
f2_all <- f2_blocks[pops_all, pops_all, ]

graph_comprehensive <- tibble::tribble(
  ~from, ~to,
  "root", pop_outgroup,
  "root", "node_A",
  "node_A", pop_basal,
  "node_A", "node_EastAsia",
  "node_EastAsia", "anc_North",
  "node_EastAsia", "anc_South",
  "anc_North", pop_north,
  "anc_North", "node_North_stream",
  "anc_South", "split_South",
  "anc_South", "split_Coastal",
  "split_South", pop_south,
  "split_Coastal", pop_coastal,
  "split_Coastal", "node_Coastal_stream",
  
  # Independent Admixture for Targets Populations
  "node_North_stream", "node_canton",
  "node_Coastal_stream", "node_canton",
  "node_canton", "GDH_Cantonese",
  
  "node_North_stream", "node_hakka",
  "node_Coastal_stream", "node_hakka",
  "node_hakka", "GDH_Hakka",
  
  "node_North_stream", "node_teochew",
  "node_Coastal_stream", "node_teochew",
  "node_teochew", "GDH_Teochew",
  
  "node_North_stream", "node_yao",
  "node_Coastal_stream", "node_yao",
  "node_yao", "GDY"
)

fit_final <- qpgraph(f2_all, graph_comprehensive, return_fstats = TRUE)
max_z_final <- max(abs(fit_final$f4$z))
cat("\n【Comprehensive Model Fit】\n")
cat("Score:", round(fit_final$score, 2), "| Max |Z|:", round(max_z_final, 3), "\n")


