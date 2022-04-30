#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get DsrAB concatenated sequences and make the phylogenetic tree
# Note: This script should be run under the conda env of IQTREE_v2.1.4. To activate it, run "conda activate IQTREE_v2.1.4"

# Step 1 Store the Sequences of DsrA and DsrB
my %DsrA_seq = _store_seq("Dsr_analysis/DsrA.faa");
my %DsrB_seq = _store_seq("Dsr_analysis/DsrB.faa");

# Step 2 Make each DsrA and DsrB protein sequences
my %DsrA_paired_seq = (); # $mag_and_num => $dsrA_seq
my %DsrB_paired_seq = (); # $mag_and_num => $dsrB_seq
open IN, "Dsr_analysis/DsrAB_pair.mdf.txt";
while (<IN>){
	chomp;
	if (!/^DsrA/){
		my @tmp = split (/\t/);
		my $dsrA = $tmp[0];
		my $dsrB = $tmp[1];
		my $mag_and_num = $tmp[2];
		my $dsrA_seq = $DsrA_seq{">$dsrA"};
		my $dsrB_seq = $DsrB_seq{">$dsrB"};
		$DsrA_paired_seq{">$mag_and_num"} = $dsrA_seq;
		$DsrB_paired_seq{">$mag_and_num"} = $dsrB_seq;
	}
}
close IN;

## Write down the sequences
open OUT, ">Dsr_analysis/dsrA.paired.faa";
foreach my $key (sort keys %DsrA_paired_seq){
	print OUT "$key\n$DsrA_paired_seq{$key}\n";
}
close OUT;

open OUT, ">Dsr_analysis/dsrB.paired.faa";
foreach my $key (sort keys %DsrB_paired_seq){
	print OUT "$key\n$DsrB_paired_seq{$key}\n";
}
close OUT;

# Step 3 Cat each sequence and run mafft
`cat Dsr_analysis/dsrA.paired.faa Dsr_analysis/dsrA_protein_conserved_residues_alignment.degap.fasta > Dsr_analysis/dsrA.hits_and_ref.faa`;
`cat Dsr_analysis/dsrB.paired.faa Dsr_analysis/dsrB_protein_conserved_residues_alignment.degap.fasta > Dsr_analysis/dsrB.hits_and_ref.faa`;
`mafft Dsr_analysis/dsrA.hits_and_ref.faa > Dsr_analysis/dsrA.hits_and_ref.faa.mafft`;
`mafft Dsr_analysis/dsrB.hits_and_ref.faa > Dsr_analysis/dsrB.hits_and_ref.faa.mafft`;

# Step 4 Concatenate dsrAB mafft file
my %DsrA_hits_and_ref = _store_seq("Dsr_analysis/dsrA.hits_and_ref.faa.mafft");
my %DsrB_hits_and_ref = _store_seq("Dsr_analysis/dsrB.hits_and_ref.faa.mafft");

my %DsrAB_hits_and_ref = ();
foreach my $key (sort keys %DsrA_hits_and_ref){
	my $seq1 = $DsrA_hits_and_ref{$key};
	my $seq2 = $DsrB_hits_and_ref{$key};
	my $seq_concat = $seq1.$seq2;
	$DsrAB_hits_and_ref{$key} = $seq_concat;
}

open OUT, ">Dsr_analysis/dsrAB.hits_and_ref.faa.mafft";
foreach my $key (sort keys %DsrAB_hits_and_ref){
	print OUT "$key\n$DsrAB_hits_and_ref{$key}\n";
}
close OUT;

# Step 5 Run IQTREE
`trimal -in Dsr_analysis/dsrAB.hits_and_ref.faa.mafft -out Dsr_analysis/dsrAB.hits_and_ref.faa.mafft.gt25perc -gt 0.25`;
`iqtree -m TESTMERGE -bb 1000 -bnni -nt AUTO -s Dsr_analysis/dsrAB.hits_and_ref.faa.mafft.gt25perc -pre Dsr_analysis/dsrAB.hits_and_ref.faa.mafft.gt25perc`;



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