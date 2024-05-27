#!/usr/bin/perl

use strict;
use warnings;

# Aim: Conduct statistical analysis to get the final sulfur metabolic function distribution across all genomes and across all classes

# Step 1 Store the original HMM result hash
my %HMM_result = (); # $hmm (e.g., K00385, w/o ".hmm") => $mag => $hit_presence
my @Header = (); # Store the header line
my %MAGs = (); # $mag => 1; Store all MAGs
open IN, "HMM_result.original.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		@Header = split (/\t/);
	}else{
		my @tmp = split (/\t/);
		my $hmm = $tmp[0];
		for(my $i=1; $i<=$#tmp; $i++){
			my $mag = $Header[$i];
			$MAGs{$mag} = 1;
			my $hit_num = $tmp[$i];
			my $hit_presence = 0;
			if ($hit_num){
				$hit_presence = 1;
			}
			$HMM_result{$hmm}{$mag} = $hit_presence;
		}
	}
}
close IN;

# Step 2 Store MAG2tax hash
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

# Step 3 Store the %MAG2Dsr_function hash
my %MAG2Dsr_function = (); # $mag => $dsr_function (e.g., "oxidative Dsr" or "reduction Dsr") => 1 
open IN, "Dsr_analysis/Dsr_function_summary.txt";
while (<IN>){
	chomp;
	if (!/^dsrA gene/){
		my @tmp = split (/\t/);
		my $mag = $tmp[3];
		my $dsr_function = "";
		if ($tmp[8]){
			$dsr_function = $tmp[8]; # Here, $dsr_function can be "oxidative Dsr", "reduction Dsr", or "both"
		}
		
		if ($dsr_function and $dsr_function ne "both"){
			$MAG2Dsr_function{$mag}{$dsr_function} = 1;
		}elsif($dsr_function and $dsr_function eq "both"){
			$MAG2Dsr_function{$mag}{"oxidative Dsr"} = 1;
			$MAG2Dsr_function{$mag}{"reduction Dsr"} = 1;
		}
	}
}
close IN;

# Step 4 Store the %MAG2Sdo hash
my %MAG2Sdo = (); # $mag => 1; indicates whether a MAG contains Sdo or not
open IN, "Sdo_analysis/Sdo_positive_hits.txt";
while (<IN>){
	chomp;
	if (/^>GC/){
		my $line = $_;
		my ($mag) = $line =~ /^>(GC.\_.+?)\_\_\_/;
		$MAG2Sdo{$mag} = 1;
	}
}
close IN;

# Step 5 Store the %MAG2sHdr hash
my %MAG2sHdr = (); # $mag => 1; indicates whether a MAG contains sHdr operon or not 
open IN, "sHdr_operon_results.txt";
while (<IN>){
	chomp;
	my $line = $_;
	my ($mag) = $line =~ /^(GC.\_.+?)\t/;
	$MAG2sHdr{$mag} = 1;
}
close IN;

# Step 6 Merge the Dsr, Sdo, sHdr hit presence info into %HMM_result
foreach my $mag (sort keys %MAGs){
	if (exists $MAG2Dsr_function{$mag}{"oxidative Dsr"}){
		$HMM_result{"oxidative Dsr"}{$mag} = 1;
	}else{
		$HMM_result{"oxidative Dsr"}{$mag} = 0;
	}
	
	if (exists $MAG2Dsr_function{$mag}{"reductive Dsr"}){
		$HMM_result{"reductive Dsr"}{$mag} = 1;
	}else{
		$HMM_result{"reductive Dsr"}{$mag} = 0;
	}

	if (exists $MAG2Sdo{$mag}){
		$HMM_result{"K17725"}{$mag} = 1;
	}else{
		$HMM_result{"K17725"}{$mag} = 0;
	}
	
	if (exists $MAG2sHdr{$mag}){
		$HMM_result{"sHdr"}{$mag} = 1;
	}else{
		$HMM_result{"sHdr"}{$mag} = 0;
	}	
}

# Step 6 Store the %Sulfur_metabolic_function2HMM_template hash
my %Sulfur_metabolic_function2HMM_template = (); # $function => $hmms
my %Sulfur_metabolic_function_HMM_ids = (); # $hmm => 1; $hmm will also include "oxidative Dsr" or "reductive Dsr"
my @Sulfur_metabolic_function_order = (); # Store the order of functions
open IN, "Sulfur_metabolic_function2HMM_template.txt";
while (<IN>){
	chomp;
	if (!/^Function/){
		my @tmp = split (/\t/);
		$Sulfur_metabolic_function2HMM_template{$tmp[0]} = $tmp[1];
		push @Sulfur_metabolic_function_order, $tmp[0];
		
		if ($tmp[1] !~ /\;/){
			my @tmp2 = split (/\,/,$tmp[1]);
			foreach my $key (@tmp2){
				$Sulfur_metabolic_function_HMM_ids{$key} = 1;
			}
		}elsif ($tmp[1] =~ /\;/){
			my @tmp2 = split (/\;/,$tmp[1]);
			foreach my $key (@tmp2){
				my @tmp3 = split (/\,/,$key);
				foreach my $key2 (@tmp3){
					if ($key2 !~ /NO/){
						$Sulfur_metabolic_function_HMM_ids{$key2} = 1;
					}
				}
			}
		}
	}
}
close IN;

# Step 7 Assign the Sulfur metabolic function for all MAGs
my %Sulfur_metabolic_function2mag2presence = (); # $function => $mag => $presence (absent or present as "0" or "1")
foreach my $function (sort keys %Sulfur_metabolic_function2HMM_template){
	foreach my $mag (sort keys %MAGs){
		my $hmms = $Sulfur_metabolic_function2HMM_template{$function};
		if ($hmms !~ /\;/){
			my @HMMs = split (/\,/, $hmms);
			my $logic = 0; # Assign the presence or not of this function
			foreach my $hmm (@HMMs){
				if ($HMM_result{$hmm}{$mag}){
					$logic = 1;
				}
			}
			if ($logic){
				$Sulfur_metabolic_function2mag2presence{$function}{$mag} = 1;
			}else{
				$Sulfur_metabolic_function2mag2presence{$function}{$mag} = 0;
			}
		}else{
			my ($hmms_1, $hmms_2) = $hmms =~ /^(.+?)\;(.+?)$/;
			if ($hmms_2 !~ /NO/){
				my $logic_1 = 0; # Assign the presence or not of function part 1 - $hmms_1
				my $logic_2 = 0; # Assign the presence or not of function part 2 - $hmms_2
				my @HMMs_1 = split (/\,/, $hmms_1);
				foreach my $hmm (@HMMs_1){
					if ($HMM_result{$hmm}{$mag}){
						$logic_1 = 1;
					}				
				}
				
				my @HMMs_2 = split (/\,/, $hmms_2);				
				foreach my $hmm (@HMMs_2){
					if ($HMM_result{$hmm}{$mag}){
						$logic_2 = 1;
					}				
				}

				if ($logic_1 and $logic_2){
					$Sulfur_metabolic_function2mag2presence{$function}{$mag} = 1;
				}else{
					$Sulfur_metabolic_function2mag2presence{$function}{$mag} = 0;
				}
			}else{
				my $logic_1 = 0; # Assign the presence or not of function part 1 - $hmms_1
				my $logic_2 = 0; # Assign the presence or not of function part 2 - $hmms_2				
				my @HMMs_1 = split (/\,/, $hmms_1);
				foreach my $hmm (@HMMs_1){
					if ($HMM_result{$hmm}{$mag}){
						$logic_1 = 1;
					}				
				}
				
				my @HMMs_2 = split (/\,/, $hmms_2);				
				foreach my $hmm (@HMMs_2){
					my ($hmm_wo_NO) = $hmm =~ /^NO\|(.+?)$/;
					if ($HMM_result{$hmm_wo_NO}{$mag}){
						$logic_2 = 1;
					}				
				}				
				
				if ($logic_1 and !$logic_2){
					$Sulfur_metabolic_function2mag2presence{$function}{$mag} = 1;
				}else{
					$Sulfur_metabolic_function2mag2presence{$function}{$mag} = 0;
				}
			}
		}
	}
}

# Step 8 Remove MAGs that have no function hits
my %Sulfur_metabolic_function2mag2presence_filtered = (); # Remove MAGs that have no function hits
my %MAGs_filtered = (); # Remove MAGs that have no function hits
foreach my $mag (sort keys %MAGs){
	my $presence_state_for_all_functions = 0;
	foreach my $function (sort keys %Sulfur_metabolic_function2mag2presence){	
		if ($Sulfur_metabolic_function2mag2presence{$function}{$mag}){
			$presence_state_for_all_functions++;
		}
	}
	
	if ($presence_state_for_all_functions){
		$MAGs_filtered{$mag} = 1;
		foreach my $function (sort keys %Sulfur_metabolic_function2mag2presence){
			$Sulfur_metabolic_function2mag2presence_filtered{$function}{$mag} = $Sulfur_metabolic_function2mag2presence{$function}{$mag};
		}
	}
}


# Step 9 Write down Sulfur metabolic function for all MAGs
open OUT, ">Sulfur_metabolic_function.txt";
my $row=join("\t", sort keys %MAGs_filtered);
print OUT "Head\t$row\n";
foreach my $tmp1 (@Sulfur_metabolic_function_order){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %MAGs_filtered) {       
                if (exists $Sulfur_metabolic_function2mag2presence_filtered{$tmp1}{$tmp2}){
                        push @tmp, $Sulfur_metabolic_function2mag2presence_filtered{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 10 Make %Sulfur_metabolic_function2mag_class2presence hash
my %MAG_class = (); # $mag_class => collection of $mag, separated by "\t"
                    # For Bacteria, we classified it into phylum level, except for Proteobacteria
foreach my $mag (sort keys %MAGs_filtered){
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

my %Sulfur_metabolic_function2mag_class2presence = (); # $function => $mag_class => $presence (absent or present as "0" or "1")
foreach my $function (sort keys %Sulfur_metabolic_function2mag2presence_filtered){
	foreach my $mag_class (sort keys %MAG_class){
		my $presence = 0;
		my @MAGs = split (/\t/, $MAG_class{$mag_class});
		foreach my $mag (@MAGs){
			if ($Sulfur_metabolic_function2mag2presence_filtered{$function}{$mag}){
				$presence = 1;
			}
		}
		$Sulfur_metabolic_function2mag_class2presence{$function}{$mag_class} = $presence;
	}
}

# Step 11 Remove $mag_class that have no function hits
my %Sulfur_metabolic_function2mag_class2presence_filtered = (); # Remove $mag_class that have no function hits
my %MAG_class_filtered = (); # Remove $mag_class that have no function hits
foreach my $mag_class (sort keys %MAG_class){
	my $presence_state_for_all_functions = 0;
	foreach my $function (sort keys %Sulfur_metabolic_function2mag_class2presence){	
		if ($Sulfur_metabolic_function2mag_class2presence{$function}{$mag_class}){
			$presence_state_for_all_functions++;
		}
	}
	
	if ($presence_state_for_all_functions){
		$MAG_class_filtered{$mag_class} = 1;
		foreach my $function (sort keys %Sulfur_metabolic_function2mag_class2presence){
			$Sulfur_metabolic_function2mag_class2presence_filtered{$function}{$mag_class} = $Sulfur_metabolic_function2mag_class2presence{$function}{$mag_class};
		}
	}
}

# Step 12 Write down Sulfur metabolic function for all Class
open OUT, ">Sulfur_metabolic_function.for_class.txt";
my $row2=join("\t", sort keys %MAG_class_filtered);
print OUT "Head\t$row2\n";
foreach my $tmp1 (@Sulfur_metabolic_function_order){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %MAG_class_filtered) {       
                if (exists $Sulfur_metabolic_function2mag_class2presence_filtered{$tmp1}{$tmp2}){
                        push @tmp, $Sulfur_metabolic_function2mag_class2presence_filtered{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

