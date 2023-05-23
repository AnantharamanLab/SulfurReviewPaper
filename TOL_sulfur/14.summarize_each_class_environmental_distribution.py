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
    
# Aim: Get environmental distribution result for each class (131 tips in total; 26 tips of archaeal classes, 105 tips of bacterial classes)
# Note: If any species representatives (here, abbreviated as "species") are found in a specific environmental category, we assign "1" to the corresponding class


# Step 1 Store the "Species_rep2Isolation_mdfed.txt" table 
# "Species_rep2Isolation_mdfed.txt"'s 4th column is the corresponding environment (lower level)
class2species_list = defaultdict(list)
species2env = {}
with open('Species_rep2Isolation_mdfed.txt', 'r') as lines:
    for line in lines:
        if not line.startswith('Species'): 
            line = line.rstrip('\n')
            tmp = line.split('\t')
            species, class_, env = tmp[0], tmp[1], tmp[3]
            species2env[species] = env
            class2species_list[class_].append(species)
lines.close()            
        

# Step 2 Store environmental category map
Env_cat_map = {} # env => env_cat
with open('Environmental_category_map.txt', 'r') as lines:
    for line in lines: 
        line = line.rstrip('\n')
        tmp = line.split('\t')
        Env_cat_map[tmp[0]] = tmp[1]
lines.close()     

env_cat_list = [
    "Air",
    "Saline and alkaline",
    "Soil and wetlands",
    "Terrestrial, cave, and glacial",
    "Deep subsurface and aquifer",
    "Engineered",
    "Freshwater",
    "Marine",
    "Hydrothermal vents and thermal springs",
    "Human and mammals",
    "Plants: rhizoplane",
    "other host-associated"
]


# Step 3 Get the class2env_cat2presence dict and write the result down
class2env_cat2presence = defaultdict(dict)
for class_ in class2species_list:
    for env_cat in env_cat_list:
        presence = '0'
        species_list = class2species_list[class_]
        for species in species_list:
            env = species2env[species] 
            if len(env) > 0:
                if Env_cat_map[env] == env_cat:
                    presence = '1'
        class2env_cat2presence[class_][env_cat] = presence

with open('class2env_cat2presence.txt', "w") as file:
    ## Write table header
    header = "\t".join(env_cat_list)
    file.write("Class\t" + header + "\n")

    ## Write table rows
    for class_ in class2species_list:
        row_values = [str(class2env_cat2presence.get(class_, {}).get(env_cat, '')) for env_cat in env_cat_list]
        row = class_ + "\t" + "\t".join(row_values) + "\n"
        file.write(row)        

