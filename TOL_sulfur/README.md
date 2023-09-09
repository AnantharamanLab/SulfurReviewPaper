### Tree for sulfur metabolism pathways

------

PS: The working dir in our server is "/storage1/data10/SulfurReviewPaper/TOL_sulfur"

**1** **Prepare sulfur metabolism HMMs**

Here we focus on dissimilatory sulfur metabolism HMMs. The collection contains 37 HMMs.

The custom HMM NC was copied from [METABOLIC](https://github.com/AnantharamanLab/METABOLIC).
The Kofam HMM NC was copied from ko_list. Kofam HMMs are from release Feb 2022.

NC cutoffs were all written into individual HMMs. All HMMs were provided in the folder "hmms_sulfur".

**2 Run hmmsearch for GTDB_release202_species_rep_genomes**

The GTDB release202 species rep genomes were downloaded from the [database](https://data.gtdb.ecogenomic.org/releases/release202/202.0/).

Individual genomes were annotated by Prodigal and all proteins were concatenated into "concat_gtdbtk_genomes2.faa". The concat protein file was used for hmmsearch

[script] 01.run_hmmsearch_on_GTDB_release202_species_rep_genomes.pl

**3 Parse the hmmsearch result and select DsrAB, Sdo proteins for further phylogenetic tree-based validation**

Make "HMM_result.original.txt" to store the result. It contains all HMM in the first column and all genomes in the first row. 

[script] 02.parse_hmmsearch_result_and_grasp_DsrAB_and_Sdo_proteins.pl

**4 Get DsrAB pair**

Parse to get DsrAB pair information. Manually curate the results based on the following method:

Criteria to choose DsrAB on the same scaffold:

1) Firstly pick DsrAB next to each other

2) If both proteins located next to the target, pick the one with the higher score

[script] 03.parse_to_get_DsrAB_pairs.pl

**5 Align DsrA and DsrB separately and make the phylogenetic tree based on the concatenated DsrAB alignment**

DsrA sequences (both reference and DsrA hits from this study) and DsrB sequences (both reference and DsrB hits from this study) were aligned separately. The concatenated DsrAB alignment was subjected to phylogenetic tree reconstruction.

[script] 04.make_DsrAB_tree.pl

**6 Replace the tip name of DsrAB tree**

Use both the genome ID and taxonomy to replace the genome ID in the tree generated from Step 5.

[script] 05.replace_DsrAB_tree_tips.pl

**7 Re-assign the Dsr function by Dsr operon**

Re-assign the Dsr function ("oxidative Dsr", "reductive Dsr", or "both") by Dsr operon (the assignment of DsrAB direction from the phylogenetic tree and existence of Dsr operon components). The following rule was adopted:

<img src="https://github.com/AnantharamanLab/SulfurReviewPaper/blob/main/TOL_sulfur/Dsr_direction_rule.jpg" style="zoom:80%;" />

[script] 06.re-assign_Dsr_function.pl

**8 Align Sdo sequences, and make the phylogenetic tree**

Align them with Sdo reference sequences, and make the phylogenetic tree

[script] 07.make_Sdo_tree.pl

**9 Find positive Sdo hits**

Find positive Sdo hits through two conserved sites, hETHE1 Asp196 and Asn244 (Reference: Fig. 2 in Appl Environ Microbiol. 2014 Mar;80(5):1799-806. doi: 10.1128/AEM.03281-13) 

[script] 08.find_positive_Sdo_hits.pl

**10 Conduct statistical analysis to sulfur metabolic functions**

Conduct statistical analysis to get the final sulfur metabolic function distribution across all genomes and across all classes (for Bacteria, it is phyla except for Proteobacteria).

The dependent files include "Sulfur_metabolic_function2HMM_template.txt" (Function to HMM template file), "HMM_result.original.txt" (intermediate result file from Step 3), "GTDBTK_DB_release202/ar122_taxonomy_r202.tsv" (GTDB archaeal MAG taxonomy result), "GTDBTK_DB_release202/bac120_taxonomy_r202.tsv" (GTDB bacterial MAG taxonomy result), "Dsr_function_summary.txt" (intermediate result file from Step 7), and "Sdo_positive_hits.txt" (intermediate result file from Step 9).

The output files include "Sulfur_metabolic_function.txt" and "Sulfur_metabolic_function.for_class.txt".

[script] 09.do_statistical_analysis_to_sulfur_metabolic_functions.pl

**11 Make phylogenetic tree for each class**

We picked one representative genome from each class to build the phylogenetic tree.  We used the statistical analysis result of sulfur metabolic function generated from the last step to annotate the tree.  There were 26 tips of archaeal classes, and 105 tips of bacterial classes (proteobacterial groups were classified into the class level, while other bacterial taxa were classified into the phylum level).

[script] 10.make_phylogenomic_tree_for_each_class.pl

**12 Count GTDB species number for each class**

We counted the GTDB species numbers for each class (131 taxa in total; 26 taxa of archaeal classes, 105 taxa of bacterial phyla or classes).

[script] 11.count_species_num_for_each_class.pl

**13 Get the isolation information for each class's species representatives**

We parsed the metadata to get the isolation information for each class's species representatives.

[script] 12.get_isolate_information_for_each_class_species_reps.py

**14 Get the BioSample IDs from genome accessions**   

Use the metadata to get the corresponding BioSample IDs from genome accessions

[script] 13.get_biosample_from_genome_accessions.py

**15** **Get environmental distribution result for each class** 

After we manually assigned the environmental category to each species representative from its isolation information, we summarized the environmental distribution result for each class. 

If any species representatives are found in a specific environmental category, we assign "1" to the corresponding class. Here, we only took the species representatives into consideration.

The 4th column of "Species_rep2Isolation_mdfed.txt" is the corresponding environment (lower level). Then by using "Environmental_category_map.txt", we could find the environmental category for each species representative.

[script] 14.summarize_each_class_environmental_distribution.py



