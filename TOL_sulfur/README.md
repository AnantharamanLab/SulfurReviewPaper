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

**3 Parse the hmmsearch result and grasp DsrAB, Sdo proteins for further phylogenetic tree-based validation**

Make "HMM_result.txt" to store the result. It contains all HMM in the first column and all genomes in the first row. 

[script] 02.parse_hmmsearch_result_and_grasp_DsrAB_and_Sdo_proteins.pl

**4 Get DsrAB pair**

Parse to get DsrAB pair information. Manually curate the results based on the following method:

Criteria to choose DsrAB in the same scaffold:

1) Firstly pick DsrAB next to each other

2) If both proteins located next to the target, pick the one with the higher score

[script] 03.parse_to_get_DsrAB_pairs.pl

**5 Align DsrA and DsrB separately and make the phylogenetic tree based on the concatenated DsrAB alignment**

DsrA sequences (both reference and DsrA hits from this study) and DsrB sequences (both reference and DsrB hits from this study) were aligned separately. The concatenated DsrAB alignment was subjected to phylogenetic tree reconstruction.

[script] 04.make_DsrAB_tree.pl

6 



**7 Align Sdo sequences, and make the phylogenetic tree**

Align them with Sdo reference sequences, and make the phylogenetic tree

[script] 06.make_Sdo_tree.pl









