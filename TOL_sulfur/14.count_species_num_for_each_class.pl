#!/usr/bin/perl

use strict;
use warnings;

# Aim: Count GTDB species numbers for each class (131 tips in total; 26 tips of archaeal classes, 105 tips of bacterial classes)


# Step 1 Store MAG2tax hash
my %MAG2tax = (); # $mag => $tax
open IN, "GTDBTK_DB_release202/ar122_taxonomy_r202.tsv";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my ($mag) = $tmp[0] =~ /^.+?\_(GC.+\_.+?)$/;
	my $tax = $tmp[1]; 
	$MAG2tax{$mag} = $tax;
}
close IN;

open IN, "GTDBTK_DB_release202/bac120_taxonomy_r202.tsv";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my ($mag) = $tmp[0] =~ /^.+?\_(GC.+\_.+?)$/;
	my $tax = $tmp[1]; 
	$MAG2tax{$mag} = $tax;
}
close IN;

# Step 2 Store %MAG_class hash
my %MAG_class = (); # $mag_class => collection of $mag, separated by "\t"
                    # For Bacteria, we classified it into phylum level, except for Proteobacteria
foreach my $mag (sort keys %MAG2tax){
	my $tax = $MAG2tax{$mag};
	my $mag_class = "";
	if ($tax =~ /d\_\_Archaea/ or $tax =~ /p\_\_Proteobacteria/){
		($mag_class) = $tax =~ /^(.+?)\;o\_\_/;
	}else{
		($mag_class) = $tax =~ /^(.+?)\;c\_\_/;
	}
	
	if (!exists $MAG_class{$mag_class}){
		$MAG_class{$mag_class} = $mag;
	}else{
		$MAG_class{$mag_class} .= "\t".$mag;
	}
}

# Step 3 Store the hash of %Species_rep_ID
my %Seq1 = _store_seq("GTDBTK_DB_release202/ar122_msa_reps_r202.faa");
my %Seq2 = _store_seq("GTDBTK_DB_release202/bac120_msa_reps_r202.faa");

my %Arc_species_rep_ID = (); # Store the archaeal species representative genome's ID
my %Bac_species_rep_ID = (); # Store the bacterial species representative genome's ID
foreach my $key (sort keys %Seq1){
	my ($arc_species_rep_id) = $key =~ /^.+?\_(GC.+\_.+?)$/;
	$Arc_species_rep_ID{$arc_species_rep_id} = 1;
}

foreach my $key (sort keys %Seq2){
	my ($bac_species_rep_id) = $key =~ /^.+?\_(GC.+\_.+?)$/;
	$Bac_species_rep_ID{$bac_species_rep_id} = 1;
}

my %Species_rep_ID = (%Arc_species_rep_ID, %Bac_species_rep_ID);

# Step 4 Count species num of each class
my %MAG_class2species_num = (); # $mag_class => $species_num
foreach my $mag_class (sort keys %MAG_class){
	my $species_num = 0;
	
	my @MAGs = split (/\t/, $MAG_class{$mag_class});
	foreach my $mag (@MAGs){
		if (exists $Species_rep_ID{$mag}){
			$species_num++;
		}
	}
	
	$MAG_class2species_num{$mag_class} = $species_num;
}

# Step 5 Write down the hash of %MAG_class2species_num
open OUT, ">MAG_class2species_num.txt";
foreach my $mag_class (sort keys %MAG_class2species_num){
	print OUT "$mag_class\t$MAG_class2species_num{$mag_class}\n";
}
close OUT;

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