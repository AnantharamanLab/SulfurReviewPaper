#!/usr/bin/perl

use strict;
use warnings;

# Aim: parse hmmsearch result and grasp DsrAB and Sdo proteins
=pod
# Step 1 Store MAG to scaffold hash
my %MAG2scfs = (); # $mag => $scfs (collection of scaffolds, separated with ",")
my %Scf2mag = (); # $scf => $mag
open IN, "find /storage1/data10/SulfurReviewPaper/TOL_sulfur/GTDBTK_DB_release202/fna_unzipped -name '*.fna' | grep -v 'genes' |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($mag) = $file =~ /fna_unzipped\/(GC.+?)\_genomic/;
	my @Scfs = (); # Store all the scaffolds within this MAG
	
	# Open file and store scaffolds
	open INN, "$file";
	while (<INN>){
		chomp;
		if (/^>/){
			my $line = $_;
			my ($scf) = $line =~ /^>(.+?)\s/;
			push @Scfs, $scf;
			$Scf2mag{$scf} = $mag;
		}
	}
	close INN;
	
	$MAG2scfs{$mag} = join("\,", @Scfs);
}
close IN;

## Write down the hash
open OUT, ">GTDBTK_DB_release202/MAG2scfs.txt";
foreach my $mag (sort keys %MAG2scfs){
	print OUT "$mag\t$MAG2scfs{$mag}\n";
}
close OUT;
=cut

## Re-store the MAG2scfs hash
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

# Step 2 Parse each hmmsearch result file and generate targeted protein sequences
## Step 2.1 Store all protein sequence
my %All_pro_seq = _store_seq("/storage1/data10/SulfurReviewPaper/TOL_sulfur/GTDBTK_DB_release202/concat_gtdbtk_genomes2.faa");

## Step 2.2 Parse hmmsearch result file and store DsrAB and Sdo proteins
my %HMM_result = (); # $hmm => $mag => $hit_num
my %DsrA_faa = (); # Store all the DsrA protein sequences
my %DsrB_faa = (); # Store all the DsrB protein sequences
my %Sdo_faa = (); # Store all the Sdo protein sequences
open IN, "ls Hmmsearch_outputs/*.tblout.txt |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($hmm) = $file =~ /GTDBrelease202\.(.+?)\.hmm/;
	
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^#/){
			my $line = $_;
			$line =~ s/(^\s+|\s+$)//g;
			$line =~ s/ +/ /g;
			my @tmp = split (/\s/, $line);
			my $protein = $tmp[0];
			my ($scf) = $protein =~ /^(.+)\_/;
			my $mag = $Scf2mag{$scf};
			$HMM_result{$hmm}{$mag}++;
			
			# Store targeted protein sequences
			if ($hmm eq "dsrA"){
				my $protein_w_array = ">$protein";
				my $seq_for_this_protein = $All_pro_seq{$protein_w_array};
				my $new_protein_w_array = ">$mag\_\_\_$protein";
				$DsrA_faa{$new_protein_w_array} = $seq_for_this_protein;
			}
			
			if ($hmm eq "dsrB"){
				my $protein_w_array = ">$protein";
				my $seq_for_this_protein = $All_pro_seq{$protein_w_array};
				my $new_protein_w_array = ">$mag\_\_\_$protein";
				$DsrB_faa{$new_protein_w_array} = $seq_for_this_protein;
			}

			if ($hmm eq "K17725"){
				my $protein_w_array = ">$protein";
				my $seq_for_this_protein = $All_pro_seq{$protein_w_array};
				my $new_protein_w_array = ">$mag\_\_\_$protein";
				$Sdo_faa{$new_protein_w_array} = $seq_for_this_protein;
			}			
		}
	}
	close INN;
	
}
close IN;

## Step 2.3 Write down targeted protein sequences
open OUT, ">DsrA.faa";
foreach my $key (sort keys %DsrA_faa){
	print OUT "$key\n$DsrA_faa{$key}\n";
}
close OUT;

open OUT, ">DsrB.faa";
foreach my $key (sort keys %DsrB_faa){
	print OUT "$key\n$DsrB_faa{$key}\n";
}
close OUT;

open OUT, ">Sdo.faa";
foreach my $key (sort keys %Sdo_faa){
	print OUT "$key\n$Sdo_faa{$key}\n";
}
close OUT;

## Step 2.4 Write down the %HMM_result file
open OUT, ">HMM_result.original.txt";  # Write down the HMM result before curating DsrAB and Sdo
my $row=join("\t", sort keys %MAG2scfs);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %HMM_result){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %MAG2scfs) {       
                if (exists $HMM_result{$tmp1}{$tmp2}){
                        push @tmp, $HMM_result{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
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