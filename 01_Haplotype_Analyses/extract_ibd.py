#!/usr/bin/env python3
"""
Author:  Qiqian Xie
Purpose: To parse 'hap-ibd' output files and classify Identity-by-Descent (IBD) segments into specific length bins (1-2 cM, 2-5 cM, >5 cM). 
         This demonstrates the enhanced resolution of WGS data in capturing deep/ancient ancestral connections.

Input:
    1. A directory containing hap-ibd output files (e.g., '.ibd.gz').
    2. A population mapping file (e.g., 'WGS.ids') formatted as:
        - Column 1: Individual ID
        - Column 2: Population label

Output:
    - A CSV file containing:
        - Population Pair (pop1, pop2)
        - Length Category (1-2cM, 2-5cM, >5cM)
        - Total theoretical pairs
        - Total shared length and block counts
        - Mean shared length and mean number of blocks
"""

import gzip
import os
import csv
from collections import defaultdict

def main():
    print("Starting standardized IBD extraction task...")

    # Load population map (ID to Population)
    pop_map = {}
    pop_individuals = defaultdict(set)
    
    # NOTE: Replace with relative path in actual execution
    pop_file = "./WGS.ids" 
    with open(pop_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                ind_id, pop_id = parts[0], parts[1]
                pop_map[ind_id] = pop_id
                pop_individuals[pop_id].add(ind_id)

    known_ids = sorted(pop_map.keys(), key=len, reverse=True)

    def get_clean_id(raw_id):
        if raw_id in pop_map: return raw_id
        for kid in known_ids:
            if raw_id.startswith(kid): return kid
        return None

    # Pre-calculate total theoretical pairs for each population combination
    all_pops = list(pop_individuals.keys())
    pair_counts = {}
    
    for i in range(len(all_pops)):
        for j in range(i, len(all_pops)):
            p1, p2 = sorted([all_pops[i], all_pops[j]])
            if p1 == p2:
                n = len(pop_individuals[p1])
                total_pairs = n * (n - 1) // 2
            else:
                n1 = len(pop_individuals[p1])
                n2 = len(pop_individuals[p2])
                total_pairs = n1 * n2
            pair_counts[(p1, p2)] = total_pairs

    # Scan genome-wide IBD files and bin by length
    stats = defaultdict(lambda: {'len': 0.0, 'blocks': 0})
    input_dir = "./ibd_results"

    print("Scanning IBD output files...")
    for file in os.listdir(input_dir):
        if file.endswith(".ibd.gz"):
            with gzip.open(os.path.join(input_dir, file), 'rt') as f:
                for line in f:
                    p = line.strip().split()
                    if len(p) < 8: continue
                    
                    id1 = get_clean_id(p[0])
                    id2 = get_clean_id(p[2])
                    length = float(p[7])
                    
                    if id1 and id2:
                        p1, p2 = sorted([pop_map[id1], pop_map[id2]])
                        
                        # Binning into Ancient (1-2cM), Intermediate (2-5cM), and Recent (>5cM)
                        if 1 <= length < 2: cat = "1-2cM"
                        elif 2 <= length < 5: cat = "2-5cM"
                        else: cat = ">5cM"
                        
                        stats[(p1, p2, cat)]['len'] += length
                        stats[(p1, p2, cat)]['blocks'] += 1

    # Calculate means and output
    output_file = "WGS_IBD_Network_Detailed.csv"
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pop1", "pop2", "Category", "total_pairs", "total_length", 
                         "total_blocks", "mean_shared_length", "mean_num_blocks"])
        
        for (p1, p2, cat), data in stats.items():
            t_pairs = pair_counts.get((p1, p2), 0)
            mean_len = data['len'] / t_pairs if t_pairs > 0 else 0
            mean_blocks = data['blocks'] / t_pairs if t_pairs > 0 else 0
                
            writer.writerow([
                p1, p2, cat, t_pairs, 
                round(data['len'], 4), data['blocks'], 
                round(mean_len, 6), round(mean_blocks, 6)
            ])

    print(f"Analysis complete. Results saved to {output_file}")

if __name__ == "__main__":
    main()