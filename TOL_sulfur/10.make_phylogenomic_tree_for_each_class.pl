#!/usr/bin/perl

use strict;
use warnings;

# Aim: Make phylogenomic tree for each class (131 tips in total; 26 tips of archaeal classes, 105 tips of bacterial classes)
# Note: This script should be run under the conda env of IQTREE_v2.1.4. To activate it, run "conda activate IQTREE_v2.1.4"
=pod
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

# Step 3 Store %Arc_MAG2msa and %Bac_MAG2msa hash
my %Seq1 = _store_seq("GTDBTK_DB_release202/ar122_msa_reps_r202.faa");
my %Seq2 = _store_seq("GTDBTK_DB_release202/bac120_msa_reps_r202.faa");

my %Arc_MAG2msa = (); # Store archaeal MAG to msa alignment
my %Bac_MAG2msa = (); # Store bacterial MAG to msa alignment
foreach my $key (sort keys %Seq1){
	my ($key_new) = $key =~ /\_(GC.\_.+?)$/;
	my $key_new_arrow = ">$key_new";
	$Arc_MAG2msa{$key_new_arrow} = $Seq1{$key};
}

foreach my $key (sort keys %Seq2){
	my ($key_new) = $key =~ /\_(GC.\_.+?)$/;
	my $key_new_arrow = ">$key_new";
	$Bac_MAG2msa{$key_new_arrow} = $Seq2{$key};
}

# Step 4 Make %Arc_class2msa and %Bac_class2msa hash
my %Arc_class_list = (); # $arc_class => 1 
my %Bac_class_list = (); # $bac_class => 1 

open IN, "Sulfur_metabolic_function.for_class.mdfed.txt";
while (<IN>){
	chomp;
	if (!/^Head/){
		my @tmp = split (/\t/);
		my $class = $tmp[0];
		if ($class =~ /d\_\_Archaea/){
			$Arc_class_list{$class} = "";
		}elsif($class =~ /d\_\_Bacteria/){
			$Bac_class_list{$class} = "";
		}
	}
}
close IN;

my %Arc_class_map = (); # $arc_class => $mag (the first one in the list and has the correpsonding msa) 
my %Bac_class_map = (); # $bac_class => $mag (the first one in the list and has the correpsonding msa) 
foreach my $arc_class (sort keys %Arc_class_list){
	my @Arc_MAGs = split (/\t/, $MAG_class{$arc_class});
	foreach my $arc_mag (@Arc_MAGs){
		if (exists $Arc_MAG2msa{">$arc_mag"}){
			$Arc_class_map{$arc_class} = $arc_mag;
			last;
		}
	}
}

foreach my $bac_class (sort keys %Bac_class_list){
	my @Bac_MAGs = split (/\t/, $MAG_class{$bac_class});
	foreach my $bac_mag (@Bac_MAGs){
		if (exists $Bac_MAG2msa{">$bac_mag"}){
			$Bac_class_map{$bac_class} = $bac_mag;
			last;
		}
	}
}

my %Arc_class2msa = ();
my %Bac_class2msa = ();

foreach my $arc_class (sort keys %Arc_class_map){
	my $arc_mag = $Arc_class_map{$arc_class};
	my $arc_class_w_arrow = ">$arc_class";
	my $arc_mag_w_arrow = ">$arc_mag";
	if (exists $Arc_MAG2msa{$arc_mag_w_arrow}){ 
		$Arc_class2msa{$arc_class_w_arrow} = $Arc_MAG2msa{$arc_mag_w_arrow};
	}else{
		print "$arc_mag_w_arrow\n";
	}
}

foreach my $bac_class (sort keys %Bac_class_map){
	my $bac_mag = $Bac_class_map{$bac_class};
	my $bac_class_w_arrow = ">$bac_class";
	my $bac_mag_w_arrow = ">$bac_mag";
	$Bac_class2msa{$bac_class_w_arrow} = $Bac_MAG2msa{$bac_mag_w_arrow};
}

# Step 5 Write down %Arc_class2msa and %Bac_class2msa
open OUT, ">Arc_class2msa.faa";
foreach my $key (sort keys %Arc_class2msa){
	print OUT "$key\n$Arc_class2msa{$key}\n";
}
close OUT;

open OUT, ">Bac_class2msa.faa";
foreach my $key (sort keys %Bac_class2msa){
	print OUT "$key\n$Bac_class2msa{$key}\n";
}
close OUT;

# Step 6 Write down %Arc_class_map and %Bac_class_map
open OUT, ">Arc_class_map.txt";
foreach my $arc_mag (sort keys %Arc_class_map){
	print OUT "$arc_mag\t$Arc_class_map{$arc_mag}\n";
}

open OUT, ">Bac_class_map.txt";
foreach my $bac_mag (sort keys %Bac_class_map){
	print OUT "$bac_mag\t$Bac_class_map{$bac_mag}\n";
}

# Step 7 Run IQTREE
`trimal -in Arc_class2msa.faa -out Arc_class2msa.gt25perc -gt 0.25`;
`iqtree -m TESTMERGE -bb 1000 -bnni -nt AUTO -s Arc_class2msa.gt25perc -pre Arc_class2msa.gt25perc`;
=cut
`trimal -in Bac_class2msa.faa -out Bac_class2msa.gt25perc -gt 0.25`;
`iqtree -m TESTMERGE -bb 1000 -bnni -nt AUTO -s Bac_class2msa.gt25perc -pre Bac_class2msa.gt25perc`;

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