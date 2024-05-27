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
    
# Aim: Use the faa file to make gff files for all species rep genomes. 
#      The gff files will be used as the input for running HMSS2.
#      The format of the faa and gff files can be found at:
#      https://github.com/TSTanabe/HMSS2?tab=readme-ov-file#fasta-and-gff-file-format-examples  


def parse_faa_files(faa_dir):
    genome2all_headers_list = {}  # Dictionary to store genome name to headers list mapping
    for faa_file in os.listdir(faa_dir):
        if faa_file.endswith(".faa"):
            genome_name = os.path.splitext(faa_file)[0]  # Get the stem name of the faa file
            faa_path = os.path.join(faa_dir, faa_file)
            with open(faa_path, 'r') as faa:
                headers = []
                for line in faa:
                    if line.startswith('>'):
                        header = line.strip()[1:]
                        headers.append(header)
            genome2all_headers_list[genome_name] = headers
    return genome2all_headers_list

def generate_gff(genome2all_headers_list, gff_dir):
    for genome, headers in genome2all_headers_list.items():
        for header in headers:
            header_parts = header.split(' # ')  
            gene_id = header_parts[0]
            seq_id = gene_id.rsplit('_', 1)[0]
            start = header_parts[1]
            end = header_parts[2]  
            strand = '+'
            if header_parts[3] == '-1':
                strand = '-'             
            score = '-'
            phase = '-'
            ID = 'ID=' + gene_id + ';'     
            source = 'Prodigal_v2.6.3'   
            feature_type = 'CDS'            
            gff_file_name = f"{genome}.gff" # Start to generate the gff file
            gff_path = os.path.join(gff_dir, gff_file_name)
            with open(gff_path, 'a') as gff_file:
                gff_line = f"{seq_id}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{ID}\n"
                gff_file.write(gff_line)

def main():
    faa_dir = "GTDBTK_DB_release202/fna_unzipped"
    gff_dir = "GTDBTK_DB_release202/fna_unzipped"

    if not os.path.exists(gff_dir):
        os.makedirs(gff_dir)

    genome2all_headers_list = parse_faa_files(faa_dir)
    generate_gff(genome2all_headers_list, gff_dir)

if __name__ == "__main__":
    main()

          