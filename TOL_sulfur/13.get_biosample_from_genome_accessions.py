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
    
# Aim: Get the BioSample IDs from genome accessions


# Step 1 Store the list of genome accessions
genome_accessions = []
with open('genome_accessions.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if len(line) > 0:
            genome_accessions.append(line)


# Step 2 Store the genome_accession2biosample dict
genome_accession2biosample = {} # genome_accession => biosample
with open('GTDBTK_DB_release202/ar122_metadata_r202.tsv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('accession'):
            tmp = line.split('\t')
            acc = tmp[0].split('_')[1] + '_' + tmp[0].split('_')[2]
            if acc in genome_accessions:
                biosample = tmp[49]
                genome_accession2biosample[acc] = biosample
lines.close()

with open('GTDBTK_DB_release202/bac120_metadata_r202.tsv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('accession'):
            tmp = line.split('\t')
            acc = tmp[0].split('_')[1] + '_' + tmp[0].split('_')[2]
            if acc in genome_accessions:
                biosample = tmp[49]
                genome_accession2biosample[acc] = biosample
lines.close()  


# Step 3 Write the result down
f = open('genome_accesssion2biosample.txt', 'w')
for genome_accession in genome_accessions:
    line = genome_accession + '\t' + genome_accession2biosample[genome_accession]
    f.write(line + '\n')
f.close()    