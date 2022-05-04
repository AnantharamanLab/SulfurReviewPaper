#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get Sdo sequences and make the phylogenetic tree
# Note: This script should be run under the conda env of IQTREE_v2.1.4. To activate it, run "conda activate IQTREE_v2.1.4"

`cat Sdo_analysis/Sdo.faa Sdo_analysis/Sdo_ref.faa > Sdo_analysis/Sdo_proteins_and_ref.faa`;
`mafft Sdo_analysis/Sdo_proteins_and_ref.faa > Sdo_analysis/Sdo_proteins_and_ref.faa.mafft`;
`trimal -in Sdo_analysis/Sdo_proteins_and_ref.faa.mafft -out Sdo_analysis/Sdo_proteins_and_ref.faa.mafft.gt25perc -gt 0.25`;
`iqtree -m TESTMERGE -bb 1000 -bnni -nt AUTO -s Sdo_analysis/Sdo_proteins_and_ref.faa.mafft.gt25perc -pre Sdo_analysis/Sdo_proteins_and_ref.faa.mafft.gt25perc`;



# Subroutine
sub _store_seq{
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/>/){
			if (/\s/){
				($head) = $_ =~ /^(>.+?)\s/;
				$Seq{$head} = "";
			}else{
				($head) = $_ =~ /^(>.+?)$/;
				$Seq{$head} = "";
			}
		}else{
			$Seq{$head} .= $_;
		}
	}
	close _IN;
	return %Seq;
}