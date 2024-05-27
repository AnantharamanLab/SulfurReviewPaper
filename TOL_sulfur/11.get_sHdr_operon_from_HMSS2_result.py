#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    from collections import defaultdict
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from glob import glob
    from itertools import product
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: Parse the HMSS2 result to get sHdr operon 
# There are four patterns for sHdr operons:
# sHdr operon 1:	sHdrC1	sHdrB1	sHdrA	sHdrH	sHdrC2	sHdrB2
# sHdr operon 2:	sHdrC1	sHdrB1	sHdrA	sHdrH	sEtfA	sEtfB
# sHdr operon 3:	sHdrC1	sHdrB1	sHdrA	sHdrH	sHdrC2	sHdrB2	LbpA1	DsrE3B	LbpA2
# sHdr operon 4:	sHdrC1	sHdrB1	sHdrA	sHdrH	sEtfA	sEtfB	sHdrCCG1	sHdrCCG2
# The largest gene distance within a operon is 20


# Functions
def split_protein_ids(pro):
    if ',' in pro and '_' in pro: 
        return [int(x.rsplit('_', 1)[1]) for x in pro.split(',')] 
    else:
        return [int(x.rsplit('_', 1)[1]) for x in pro.split(',')] 
        
        
def is_within_distance(combination, max_distance):
    max_number = max(combination)
    min_number = min(combination)
    return (max_number - min_number) <= max_distance        


# Define sHdr operon patterns
sHdr_operon_1 = ['sHdrC1', 'sHdrB1', 'sHdrA', 'sHdrH', 'sHdrC2', 'sHdrB2']
sHdr_operon_2 = ['sHdrC1', 'sHdrB1', 'sHdrA', 'sHdrH', 'sEtfA', 'sEtfB']
sHdr_operon_3 = ['sHdrC1', 'sHdrB1', 'sHdrA', 'sHdrH', 'sHdrC2', 'sHdrB2', 'LbpA1', 'DsrE3B', 'LbpA2']
sHdr_operon_4 = ['sHdrC1', 'sHdrB1', 'sHdrA', 'sHdrH', 'sEtfA', 'sEtfB', 'sHdrCCG1', 'sHdrCCG2']


# Pre-define
patterns = [sHdr_operon_1, sHdr_operon_2, sHdr_operon_3, sHdr_operon_4]
max_distance = 20

# Step 1: Parse to get all HMSS2 annotation results
HMSS2_annotation_result_dir = '/storage1/data10/SulfurReviewPaper/TOL_sulfur/GTDBTK_DB_release202/fna_unzipped_HMSS2_result/1715149738.784964_project/per_genome'

HMSS2_annotation = defaultdict(dict)
gn2pro_list = defaultdict(list)
pro2hit = {}
annotation_result_files = glob(f"{HMSS2_annotation_result_dir}/*")

for each_HMSS2_annotation_result_file in annotation_result_files:
    gn = os.path.basename(each_HMSS2_annotation_result_file) + '.1'
    with open(each_HMSS2_annotation_result_file, 'r') as file:
        for line in file:
            cols = line.strip().split('\t')
            if len(cols) >= 3:
                pro, hit = cols[0], cols[1] # pro is like " SFZB01000003.1_20", hit is like "sHdrB1"
                pro2hit[pro] = hit
                HMSS2_annotation[gn][pro] = hit
                gn2pro_list[gn].append(pro)


# Step 2: Parse to get all gn to sHdr operon dict
def match_operon(pro_list, pattern, max_distance):
    result = False
    pattern_pro_list = ['' for _ in range(len(pattern))]
    for pro in pro_list:
        hit = pro2hit[pro]
        for i in range(len(pattern)):
            if hit == pattern[i]:
                if pattern_pro_list[i] == '':
                    pattern_pro_list[i] = pro
                else:       
                    pattern_pro_list[i] += ',' + pro
    print(pattern_pro_list)
                    
    pattern_pro_list_non_empty = [s for s in pattern_pro_list if s]    
    if len(pattern_pro_list_non_empty) == len(pattern): # It means the operon hits are full
        split_pattern_pro_list = [split_protein_ids(pro) for pro in pattern_pro_list]    
        combinations = list(product(*split_pattern_pro_list))    
        for comb in combinations:
            if is_within_distance(comb, max_distance):
                result = True  

    return result                 
        
gn2sHdr = defaultdict(list)
for gn, pro_list in gn2pro_list.items():
    for pattern in patterns:
        if match_operon(pro_list, pattern, max_distance):
            gn2sHdr[gn].append(pattern)

## Print results (or you can save them to a file)
output_file = "sHdr_operon_results.txt"  # Define the output file path
with open(output_file, "w") as f:
    for gn, operons in gn2sHdr.items():
        if operons:  # Check if operons is not empty
            operons_str = ['-'.join(pattern) for pattern in operons]  # Convert each pattern list to a string
            f.write(f"{gn}\t{','.join(operons_str)}\n")  # Join the strings with commas
