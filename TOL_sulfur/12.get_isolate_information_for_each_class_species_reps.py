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
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: Get the isolation information for each class's species representatives (131 tips in total; 26 tips of archaeal classes, 105 tips of bacterial classes)    


# Step 1 Store MAG2tax dict
MAG2tax = {} # MAG => tax
with open('GTDBTK_DB_release202/ar122_taxonomy_r202.tsv', 'r') as ar_file:
    for line in ar_file:
        tmp = line.strip().split('\t')
        MAG = re.search(r'^(.+?)\_(GC.+\_.+?)$', tmp[0]).group(2)
        tax = tmp[1]
        MAG2tax[MAG] = tax

with open('GTDBTK_DB_release202/bac120_taxonomy_r202.tsv', 'r') as bac_file:
    for line in bac_file:
        tmp = line.strip().split('\t')
        MAG = re.search(r'^(.+?)\_(GC.+\_.+?)$', tmp[0]).group(2)
        tax = tmp[1]
        MAG2tax[MAG] = tax
        
        
# Step 2 Store MAG2MAG_class dict
MAG2MAG_class = {} # MAG => MAG_class               
for MAG in sorted(MAG2tax.keys()):
    tax = MAG2tax[MAG]
    MAG_class = ""
    if re.search(r'd\_\_Archaea', tax) or re.search(r'p\_\_Proteobacteria', tax):
        MAG_class = re.search(r'^(.+?)\;o\_\_', tax).group(1)
    else:
        MAG_class = re.search(r'^(.+?)\;c\_\_', tax).group(1)
    
    MAG2MAG_class[MAG] = MAG_class


# Step 3 Store the Selected_class list
selected_class = []
with open("Sulfur_metabolic_function.for_class.mdfed.txt", "r") as f:
    for line in f:
        line = line.strip()
        if not line.startswith("Head"):
            tmp = line.split("\t")
            class_ = tmp[0]
            if re.search(r"d\_\_Archaea", class_):
                selected_class.append(class_)
            elif re.search(r"d\_\_Bacteria", class_):
                selected_class.append(class_)
f.close()
                

# Step 4 Store the species representatives
species_rep_IDs = []
arc_species_rep_IDs = [] # Store the archaeal species representative genome's ID
with open('GTDBTK_DB_release202/ar122_msa_reps_r202.faa', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('>'):
            arc_species_rep = line.split('_')[1] + '_' + line.split('_')[2]
            arc_species_rep_IDs.append(arc_species_rep)
lines.close()       

bac_species_rep_IDs = [] # Store the bacterial species representative genome's ID
with open('GTDBTK_DB_release202/bac120_msa_reps_r202.faa', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('>'):
            bac_species_rep = line.split('_')[1] + '_' + line.split('_')[2]
            bac_species_rep_IDs.append(bac_species_rep)
lines.close()    
species_rep_IDs = arc_species_rep_IDs + bac_species_rep_IDs


# Step 5 Store the species_rep to isolation dict
species_rep2isolation = {} # species_rep
with open('GTDBTK_DB_release202/ar122_metadata_r202.tsv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('accession'):
            tmp = line.split('\t')
            acc = tmp[0].split('_')[1] + '_' + tmp[0].split('_')[2]
            if acc in species_rep_IDs:
                isolation = tmp[58]
                species_rep2isolation[acc] = isolation
lines.close()

with open('GTDBTK_DB_release202/bac120_metadata_r202.tsv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('accession'):
            tmp = line.split('\t')
            acc = tmp[0].split('_')[1] + '_' + tmp[0].split('_')[2]
            if acc in species_rep_IDs:
                isolation = tmp[58]
                species_rep2isolation[acc] = isolation
lines.close()     


# Step 6 Write down the Species_rep2Isolation result
f = open('Species_rep2Isolation.txt', 'w')
header_line = 'Species representative\tClass\tncbi_isolate_source'
f.write(header_line + '\n')
for species_rep in species_rep2isolation:
    MAG_class = MAG2MAG_class[species_rep]
    if MAG_class in selected_class:
        isolation = species_rep2isolation[species_rep]
        line = species_rep + '\t' + MAG_class + '\t' + isolation
        f.write(line + '\n')
f.close()    
            
   
            