#!/usr/bin/perl

use strict;
use warnings;

# Aim: Replace DsrAB tree tips

# Step 1 Store tip map
my %Tip_map = (); # $mag => $mag_n_tax
open IN, "GTDBTK_DB_release202/ar122_taxonomy_r202.tsv";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my ($mag) = $tmp[0] =~ /^.+?\_(GC.+\_.+?)$/;
	my $tax = $tmp[1]; 
	$tax =~ s/\;/\_/g; # Replace ";" by "_"
	my $mag_n_tax = "$mag\_\_\_$tax";
	$Tip_map{$mag} = $mag_n_tax;
}
close IN;

open IN, "GTDBTK_DB_release202/bac120_taxonomy_r202.tsv";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my ($mag) = $tmp[0] =~ /^.+?\_(GC.+\_.+?)$/;
	my $tax = $tmp[1]; 
	$tax =~ s/\;/\_/g; # Replace ";" by "_"
	my $mag_n_tax = "$mag\_\_\_$tax";
	$Tip_map{$mag} = $mag_n_tax;
}
close IN;

# Step 2 Store DsrAB_pair2mag_n_tax hash
my %DsrAB_pair2mag_n_tax = (); # $dsrAB_pair => $mag_n_tax_n_num
open IN, "Dsr_analysis/dsrAB.hits_and_ref.faa.mafft.gt25perc";
while (<IN>){
	chomp;
	if (/^>GC/){ # If the header is from a dsrAB pair
		my $header = $_;
		my ($dsrAB_pair) = $header =~ /^>(.+?)$/;
		my ($mag, $num) = $dsrAB_pair =~ /^(.+?)\_\_\_(\d)$/;
		my $mag_n_tax = $Tip_map{$mag};
		my $mag_n_tax_n_num = $mag_n_tax."\_\_\_".$num;
		$DsrAB_pair2mag_n_tax{$dsrAB_pair} = $mag_n_tax_n_num;
	}
}
close IN;

# Step 3 Replace treefile
my $treefile = "";
open IN, "Dsr_analysis/dsrAB.hits_and_ref.faa.mafft.gt25perc.treefile";
while (<IN>){
	chomp;
	my $line = $_;
	$treefile = $line;
	foreach my $dsrAB_pair (sort keys %DsrAB_pair2mag_n_tax){
		my $mag_n_tax_n_num = $DsrAB_pair2mag_n_tax{$dsrAB_pair};
		$treefile =~ s/$dsrAB_pair/$mag_n_tax_n_num/g;
	}
}
close IN;

# Step 4 Write down the new treefile
open OUT, ">Dsr_analysis/dsrAB.hits_and_ref.faa.mafft.gt25perc.mdf.nwk";
print OUT "$treefile\n";
close OUT;

