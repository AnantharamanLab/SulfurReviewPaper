#!/usr/bin/perl

use strict;
use warnings;

# Aim: Re-assign Dsr function based on Dsr operon

# Step 1 Store DsrAB pair information
my %DsrAB_pair = (); # $dsrAB_pair => [0] $mag_n_num [1] $scf
open IN, "Dsr_analysis/DsrAB_pair.mdf.txt";
while (<IN>){
	chomp;
	if (!/^DsrA/){
		my @tmp = split (/\t/);
		my $dsrA_gene = $tmp[0];
		my $dsrB_gene = $tmp[1];
		my $mag_n_num = $tmp[2];
		my $dsrAB_pair = "$dsrA_gene\t$dsrB_gene";
		my ($scf) = $dsrA_gene =~ /^(.+)\_/;
		
		$DsrAB_pair{$dsrAB_pair}[0] = $mag_n_num;
		$DsrAB_pair{$dsrAB_pair}[1] = $scf;
	}
}
close IN;

# Step 2 Store MAG2scfs hash
my %MAG2scfs = (); # $mag => $scfs (collection of scaffolds, separated with ",")
my %Scf2mag = (); # $scf => $mag
open IN, "GTDBTK_DB_release202/MAG2scfs.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $mag = $tmp[0];
	my $scfs = $tmp[1];
	$MAG2scfs{$mag} = $scfs;
	
	my @Scfs = split (/\,/, $scfs);
	foreach my $scf (@Scfs){
		$Scf2mag{$scf} = $mag;
	}
}
close IN;

# Step 3 Store %MAG2dsrD_presence hash
my %MAG2dsrD_presence = (); # $mag => $dsrD_presence  ("0" or "1")
open IN, "Hmmsearch_outputs/sulfur_paper_GTDBrelease202.dsrD.hmm.tblout.txt";
while (<IN>){
	chomp;
	if (!/^#/){
		my $line = $_;
		$line =~ s/(^\s+|\s+$)//g;
		$line =~ s/ +/ /g;
		my @tmp = split (/\s/, $line);
		my $protein = $tmp[0];
		my ($scf) = $protein =~ /^(.+)\_/;
		my $mag = $Scf2mag{$scf};
		
		$MAG2dsrD_presence{$mag} = 1;
	}
}
close IN;

# Step 4 Store %MAG2dsrEFH_presence hash
my %MAG2dsrEFH_presence = (); # $mag => $dsrEFH_presence  ("0" for full absence or "1" partial or full presence)
open IN, "Hmmsearch_outputs/sulfur_paper_GTDBrelease202.dsrE.hmm.tblout.txt";
while (<IN>){
	chomp;
	if (!/^#/){
		my $line = $_;
		$line =~ s/(^\s+|\s+$)//g;
		$line =~ s/ +/ /g;
		my @tmp = split (/\s/, $line);
		my $protein = $tmp[0];
		my ($scf) = $protein =~ /^(.+)\_/;
		my $mag = $Scf2mag{$scf};
		
		$MAG2dsrEFH_presence{$mag} = 1;
	}
}
close IN;

open IN, "Hmmsearch_outputs/sulfur_paper_GTDBrelease202.dsrF.hmm.tblout.txt";
while (<IN>){
	chomp;
	if (!/^#/){
		my $line = $_;
		$line =~ s/(^\s+|\s+$)//g;
		$line =~ s/ +/ /g;
		my @tmp = split (/\s/, $line);
		my $protein = $tmp[0];
		my ($scf) = $protein =~ /^(.+)\_/;
		my $mag = $Scf2mag{$scf};
		
		$MAG2dsrEFH_presence{$mag} = 1;
	}
}
close IN;

open IN, "Hmmsearch_outputs/sulfur_paper_GTDBrelease202.dsrH.hmm.tblout.txt";
while (<IN>){
	chomp;
	if (!/^#/){
		my $line = $_;
		$line =~ s/(^\s+|\s+$)//g;
		$line =~ s/ +/ /g;
		my @tmp = split (/\s/, $line);
		my $protein = $tmp[0];
		my ($scf) = $protein =~ /^(.+)\_/;
		my $mag = $Scf2mag{$scf};
		
		$MAG2dsrEFH_presence{$mag} = 1;
	}
}
close IN;


# Step 5 Store MAG2tax hash
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

# Step 6 Store MAG_n_num2dsr_direction hash
my %MAG_n_num2dsr_direction = (); # $mag_n_num => $dsr_direction (oxidative DsrAB or reductive DsrAB)
open IN, "Dsr_analysis/Phylogenetic_tree_clade_info.txt";
while (<IN>){
	chomp;
	if (!/^Tip/ and /^GC/){
		my @tmp = split (/\t/);
		my $tip_name = $tmp[0];
		$tip_name =~ s/\s/\_/g;
		my ($mag, $num) = $tip_name =~ /^(GC.\_.+?)\_.+(\d)$/;
		my $mag_n_num = "$mag\_\_\_$num";
		my $dsr_direction = "reductive DsrAB";
		if ($tmp[2] and $tmp[2] eq "oxidative DsrAB"){
			$dsr_direction = "oxidative DsrAB";
		}			
		$MAG_n_num2dsr_direction{$mag_n_num} = $dsr_direction;
	}
}
close IN;

# Step 7 Summarize the final Dsr function
my %Dsr_function_summary = (); # $dsrAB_pair => [0] $scf [1] $dsrAB_located_genome
                                              # [2] $genome_tax [3] $function [4] $dsrD_presence
											  # [5] $dsrEFH_presence [6] $final_Dsr_function
											 
foreach my $dsrAB_pair (sort keys %DsrAB_pair){
	my $scf = $DsrAB_pair{$dsrAB_pair}[1];
	my $mag_n_num = $DsrAB_pair{$dsrAB_pair}[0];
	my ($dsrAB_located_genome) = $mag_n_num =~ /^(.+?)\_\_\_/;
	my $genome_tax = $MAG2tax{$dsrAB_located_genome};
	my $function = $MAG_n_num2dsr_direction{$mag_n_num};
	
	my $dsrD_presence = 0;
	if ($MAG2dsrD_presence{$dsrAB_located_genome}){
		$dsrD_presence = 1;
	}
	
	my $dsrEFH_presence = 0;
	if ($MAG2dsrEFH_presence{$dsrAB_located_genome}){
		$dsrEFH_presence = 1;
	}	
	
	my $final_Dsr_function = "";
	# The dsrAB function (oxidative - dissimilatory sulfur oxidation or reductive - dissimilatory sulfite reduction) were determined by the following criterion:												
    # 1) If DsrAB falls into oxidative clade, DsrD is absent, and DsrEFH are present (partial presence also counted) then we assign it as oxidative DsrAB
    # 2) If DsrAB falls into reductive clade, DsrD is present, and DsrEFH are absent (full absence), then we assign it as reductive DsrAB
    # 3) If DsrD and DsrEFH both are present, and DsrAB exist, then we assign it as both oxidative and reductive DsrAB	
    # Reference: Table 2 in ISME J. 2018 Jun;12(7):1715-1728. doi: 10.1038/s41396-018-0078-0						
	if ($function eq "oxidative DsrAB" and $dsrD_presence == 0 and $dsrEFH_presence == 1){
		$final_Dsr_function = "oxidative DsrAB";
	}elsif($function eq "reductive DsrAB" and $dsrD_presence == 1 and $dsrEFH_presence == 0){
		$final_Dsr_function = "reductive DsrAB";
	}elsif($dsrD_presence == 1 and $dsrEFH_presence == 1){
		$final_Dsr_function = "both";
	}
	
	$Dsr_function_summary{$dsrAB_pair}[0] = $scf;
	$Dsr_function_summary{$dsrAB_pair}[1] = $dsrAB_located_genome;
	$Dsr_function_summary{$dsrAB_pair}[2] = $genome_tax;
	$Dsr_function_summary{$dsrAB_pair}[3] = $function;
	$Dsr_function_summary{$dsrAB_pair}[4] = $dsrD_presence;
	$Dsr_function_summary{$dsrAB_pair}[5] = $dsrEFH_presence;
	$Dsr_function_summary{$dsrAB_pair}[6] = $final_Dsr_function;
}	

# Step 8 Write down the hash of Step 7
open OUT, ">Dsr_analysis/Dsr_function_summary.txt";
print OUT "DsrA\tDsrB\tDsrAB located scaffold\tDsrAB located genome\tGenome Taxonomy\tFunction (assigned based on phylogenetic inference)\tDsrD presence\tDsrEFH presence\tOxidative or reductive Dsr\n";
foreach my $dsrAB_pair (sort keys %Dsr_function_summary){
	print OUT "$dsrAB_pair\t$Dsr_function_summary{$dsrAB_pair}[0]\t$Dsr_function_summary{$dsrAB_pair}[1]\t$Dsr_function_summary{$dsrAB_pair}[2]\t$Dsr_function_summary{$dsrAB_pair}[3]\t$Dsr_function_summary{$dsrAB_pair}[4]\t$Dsr_function_summary{$dsrAB_pair}[5]\t$Dsr_function_summary{$dsrAB_pair}[6]\n";
}
close OUT;
											  




