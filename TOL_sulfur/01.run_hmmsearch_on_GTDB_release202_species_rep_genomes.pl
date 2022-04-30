#!/usr/bin/perl

use strict;
use warnings;

# Aim: Run hmmsearch on GTDB release202 species rep genomes

# Step 1 Store all the HMM 
my %HMM = (); # $hmm => 1; e.g., dsrA.hmm => 1 
open IN, "ls hmms_sulfur/*.hmm | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($hmm) = $file =~ /hmms_sulfur\/(.+?\.hmm)$/;
	$HMM{$hmm} = 1;
}
close IN;

# Step 2 Make batch run command and do parallel running
`mkdir Hmmsearch_outputs`;
my $thread = 1;
my $input_total_faa = "/storage1/data10/SulfurReviewPaper/TOL_sulfur/GTDBTK_DB_release202/concat_gtdbtk_genomes2.faa";

open OUT, ">tmp.run_hmmsearch.sh";
foreach my $hmm (sort keys %HMM){
	my $tblout_file = "Hmmsearch_outputs/sulfur_paper_GTDBrelease202.$hmm.tblout.txt";
	print OUT "hmmsearch --cut_nc --cpu $thread --tblout $tblout_file hmms_sulfur/$hmm $input_total_faa\n";
}
close OUT;

`cat tmp.run_hmmsearch.sh | parallel -j 20`;
`rm tmp.run_hmmsearch.sh`;


