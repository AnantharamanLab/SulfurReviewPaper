#!/usr/bin/perl

use strict;
use warnings;

# Aim: Parse to get DsrAB pairs
# Note: We required that DsrA and B should be on the same scaffold

# Step 1 Store the Sequences of DsrA and DsrB
my %DsrA_seq = _store_seq("DsrA.faa");
my %DsrB_seq = _store_seq("DsrB.faa");

# Step 2 Store the pair of DsrAB
my %DsrAB_pair = (); # $dsrA_n_dsrB => $mag "___" $n
                     # DsrA and B should be on the same scaffold

## Store the scaffold ID first
my %Scf = (); # $scf => 1
foreach my $key (sort keys %DsrA_seq){
	my ($dsrA) = $key =~ /^>(.+?)$/;
	my ($scf) = $dsrA =~ /^(.+)\_/;
	$Scf{$scf} = 1;
}

foreach my $key (sort keys %DsrB_seq){
	my ($dsrB) = $key =~ /^>(.+?)$/;
	my ($scf) = $dsrB =~ /^(.+)\_/;
	$Scf{$scf} = 1;
}

my %DsrAB_MAG_hit = (); # $mag => $dsrAB_pair_hit_num
foreach my $scf (sort keys %Scf){
	my $logic1 = 0; # Trace the positive hit for dsrA
	my $logic2 = 0; # Trace the positive hit for dsrB
	my ($mag) = $scf =~ /^(.+?)\_\_\_/;
	my @DsrA_genes = ();
	foreach my $key (sort keys %DsrA_seq){
		my ($dsrA) = $key =~ /^>(.+?)$/;
		if ($dsrA =~ /$scf\_/){
			$logic1++;
			push @DsrA_genes, $dsrA;
		}
	}

	my @DsrB_genes = ();	
	foreach my $key (sort keys %DsrB_seq){
		my ($dsrB) = $key =~ /^>(.+?)$/;
		if ($dsrB =~ /$scf\_/){
			$logic2++;
			push @DsrB_genes, $dsrB;
		}
	}
	
	my $dsrA_n_dsrB = "";
	if ($logic1 >= 1 and $logic2 >= 1){
		my $dsrA_genes = join("\,", @DsrA_genes);
		my $dsrB_genes = join("\,", @DsrB_genes);
		$dsrA_n_dsrB = "$dsrA_genes\t$dsrB_genes";
		$DsrAB_MAG_hit{$mag}++;
		my $n = $DsrAB_MAG_hit{$mag}; # The number of positive hit in this genome
		$DsrAB_pair{$dsrA_n_dsrB} = $mag."\_\_\_".$n;
	}
}

# Step 3 Write down the DsrAB pairs
`mkdir Dsr_analysis`;
open OUT, ">Dsr_analysis/DsrAB_pair.txt";
foreach my $dsrA_n_dsrB (sort keys %DsrAB_pair){
	print OUT "$dsrA_n_dsrB\t$DsrAB_pair{$dsrA_n_dsrB}\n";
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