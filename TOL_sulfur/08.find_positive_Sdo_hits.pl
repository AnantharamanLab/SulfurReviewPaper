#!/usr/bin/perl

use strict;
use warnings;

# Aim: Find positive Sdo hits through two conserved sites

my %Sdo_seq = ();

%Sdo_seq = _store_seq("Sdo_analysis/Sdo_proteins_and_ref.faa.mafft");

open OUT, ">Sdo_analysis/Sdo_positive_hits.txt";
foreach my $key (sort keys %Sdo_seq){
	my $seq = $Sdo_seq{$key};
	my @tmp = split ("",$seq); 
	if ($tmp[1129] eq "D" and $tmp[1310] eq "N"){
		print OUT "$key\n";
	}
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