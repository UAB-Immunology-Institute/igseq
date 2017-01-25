#!/usr/bin/env perl

###################################################################################################
#
# AbIsotyper.pl 
#
# Reads in FASTA formatted files and determines the Ab isotype that the  
#   sequence came from based on constant primer regions.  Will search for exact matches to labels
#	specified in AbIsotyper_labels.inf as well as truncated versions of these labels.
#
#
# by Christopher Fucile
#
######################################
#
# Division of Allergy, Immunology and Rheumatology
# University of Rochester Medical Center
# Copyright (c) 2013, All Rights Reserved
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
#
# USAGE:
#
#    AbIsotyper.pl <input directory> -diff:(1-X)
#    AbIsotyper.pl <input directory> -diff:(1-X)  < -size:(1-X) OR -percent:(1-100) >
#    AbIsotyper.pl <input directory> -diff:(1-X)  < -size:(1-X) OR -percent:(1-100) > -seed:(1-X)
#
###################################################################################################

use strict;
use warnings;

use Text::LevenshteinXS qw(distance);
use Carp;
use File::Basename;
use File::Path;
use POSIX;
use Data::Dumper;
use Parallel::ForkManager;
use Sys::CPU qw(cpu_count);
use Cwd;
use FindBin;

our $VERSION = '2013.03.29.1';

use 5.012;

####################################
##           Variables            ##
####################################
my $working_dir = $FindBin::RealBin;
$working_dir =~ s/\/lib//;
my %data;
my %labels;
my $difference;
my $biggestLabelSize=0;
my $inputDir;
my $outDir;
my $cutoff;
my $seed;
my $percent_diff;
my $sampleReduce = 0;
####################################

# get input directory

    $inputDir = $ARGV[0];
    
    foreach my $arg(@ARGV){    
    	$difference = $1   if  $arg =~ /\-diff:(\d+)/;
		$cutoff = $1       if ($arg =~ /\-size:(\d+)/ && !$percent_diff);
		$seed = $1		   if  $arg =~ /\-seed:(\d+)/;
		$percent_diff = $1 if ($arg =~ /\-percent:(\d+)/ && !$cutoff);
	}
	
	if($cutoff || $percent_diff){$sampleReduce = 1;}
	
    my @files;

	chomp $inputDir;

	croak "Error: Please remove whitespaces from directory name. $inputDir\n" if $inputDir =~ /\s+/;

    opendir DIR, "$inputDir" or croak "ERROR: can not open directory $inputDir: $!";

    my @tmpfiles = readdir DIR;
    my $numOfJunk = grep{ /^\./} @tmpfiles;
    splice @tmpfiles, 0, $numOfJunk;
    unshift @tmpfiles, '';
    push @files, split ' ', join " $inputDir", @tmpfiles;
	@files = grep{ /\.fasta$/} @files;

    closedir DIR;
	my $base = basename($inputDir);
	
	
# time stamp for running script and log file

    my @MONTH = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
    my $LOG_FILE = "$base\_LOG.txt";
    my $RUN_TIME = POSIX::strftime("%Y-%m-%d, %H:%M:%S", localtime);

# grab primers from file
    open my $settingsINF, "<". "$working_dir/config/primers.inf" or die "Error opening settings file $working_dir/config/primers.inf $!\n";
    my %settings;
    my %fwdPrimers;
    my %revPrimers;
    my $fwdSection=0;
    my $revSection=0;
    while (<$settingsINF>){
    	my $line = $_;
		chomp $line;
    	if ($line =~ /^# Forward/){
    		$fwdSection = $fwdSection ? 0 : 1; # if	 fwdSection flag is on turn it off, if it is off turn it on
    		next;
    	}elsif($line =~ /^# Reverse/){
    		$revSection = $revSection ? 0 : 1; # if	 revSection flag is on turn it off, if it is off turn it on
    		next;
		}
    	
    	if ($fwdSection){
    		if ($line =~ /R|Y|M|K|S|W|H|B|V|D|N/){
				$line = replaceOdd($line);
    		}
    		$fwdPrimers{$line}=1;
    		my $reverse = reverse $line;
			$reverse =~ tr/ACGTacgt()/TGCAtgca)(/;
			$fwdPrimers{$reverse}=1;
    		next;
    	}
    	
    	if ($revSection){
    		if ($line =~ /R|Y|M|K|S|W|H|B|V|D|N/){
				$line = replaceOdd($line);
    		}
    		$revPrimers{$line}=1;
    		my $reverse = reverse $line;
			$reverse =~ tr/ACGTacgt()/TGCAtgca)(/;
			$revPrimers{$reverse}=1;
    		next;
    	}
    	
		my ($data) = split('#', $line);
		my ($parameter, $value) = split ':', $data;
		$settings{$parameter} = $value;
    }

# output directory

    $outDir = "$working_dir/AbIsotyper\_Out/$base\_".POSIX::strftime("%Y_%m_%d_%H_%M_%S",localtime);    
    mkpath($outDir) or die "Error Making output directory $!\n";
    
# append to log file whatever shows up on terminal (errors only show up on STDERR on terminal)

    open STDOUT, "| tee -a $outDir'/'$LOG_FILE" or croak "Tee stdout error";
    print "*********************************************\n";
    print "> script start: $RUN_TIME\n";


# read labels into %labels variable

    open INFSETTINGS, "< $working_dir/config/AbIsotyper_labels.inf" or croak "ERROR: can't open AbIsotyper_labels.inf\nMake sure file is in same directory as AbIsotyper.pl $!\n";
    while (<INFSETTINGS>){
        chop;
        my @row = split ':';
        if(exists $labels{$row[1]}){
            print "Warning possible identical sequences in label file.\n";
            print "Ensure all sequences are unique for labeling.\n";
        }
        $labels{uc($row[1])} = uc($row[0]);
        $biggestLabelSize = length($row[1]) if length($row[1]) > $biggestLabelSize;
    }
    croak "ERROR: Problem with labels file (AbIsotyper_labels.inf). Check usage.\n" if scalar keys %labels == 0;

    close INFSETTINGS;

# read and process data 
	print "> Processing files...\n";
	my $pm = new Parallel::ForkManager(cpu_count()-1);
	my $fwdPrimer_count=0;
	my $revPrimer_count=0;
    foreach my $file(@files){                                                   # for each file
    	my $pid = $pm->start and next;                                          # fork a new process
    	
    	srand($seed) if $seed;                                                  # set rng seed
        open INF, "<$file" or croak "ERROR: can not open input file $file\n";
		
		$file =~ s/.*\/(.*)/$1/;
		open OPF, ">> $outDir/$file" or croak "Error opening: $outDir/$file\n$!\n";

        my $i=0;
        my %counts;
        my $flag;
        
        my %barcodes;
        
        while (<INF>){
        	my $line1 = $_;
        	my $line2 = <INF>;
        	
			chop($line1);
        	chop($line2);
        	
        	$i++;
        	
            if($line1 =~ /^\>/){
            	my ($fwdBarcode, $revBarcode);
            	($line2, $fwdBarcode) = removePrimer($line2, \%fwdPrimers, \$fwdPrimer_count, \%settings);
            	($line2, $revBarcode) = constantRegionPrimer($line2, \%revPrimers, \$revPrimer_count, \%settings);
            	next if ! $line2;

            	$barcodes{$fwdBarcode.'|'.$revBarcode}++;
            	$line1 = type($line1, $line2, \%labels, \%counts);

            	if($line1){
					my ($c) = $line1 =~ />(\w+?)\_/;
					$counts{$c}++;
            		print OPF $line1."\n";
            		print OPF lc $line2."\n";
            	}else{
            		die "ERROR: NOTHING ASSINGED TO A SEQUENCE";
                }
            }else{
				die "ERROR: Problem with lines in input file. ".
                    "Check formatting\n";
            }
        }
        print "> $file\tNumber of sequences with barcodes identified: ".scalar(keys(%barcodes))."\n";
        close INF;
        close OPF;
        
        print "> Isotype proportions for $file\n";
        foreach my $count(keys %counts){
        	print ">\t$count: $counts{$count}\n";
    	}
    	
    	open my $TMP, '>', "barcode_output.txt" or die "Error opening file $!\n";
    	map {print $TMP "$_\t$barcodes{$_}\n"} keys %barcodes;
    	close $TMP;
    	
        $pm->finish(0);
    }
    $pm->wait_all_children;

############## Subroutines #####################



sub replaceOdd{
	my $line = shift;
	my %primerSub = (
	'R' => '(G|A)',
	'Y' => '(T|C)',
	'M' => '(M|C)',
	'K' => '(G|T)',
	'S' => '(G|C)',
	'W' => '(A|T)',
	'H' => '(A|C|T)',
	'B' => '(G|T|C)',
	'V' => '(G|C|A)',
	'D' => '(G|T|A)',
	'N' => '(G|T|A|C)'
	);
	$line =~ s/R/$primerSub{R}/gi;
	$line =~ s/Y/$primerSub{Y}/gi;
	$line =~ s/M/$primerSub{M}/gi;
	$line =~ s/K/$primerSub{K}/gi;
	$line =~ s/S/$primerSub{S}/gi;
	$line =~ s/W/$primerSub{W}/gi;
	$line =~ s/H/$primerSub{H}/gi;
	$line =~ s/B/$primerSub{B}/gi;
	$line =~ s/V/$primerSub{V}/gi;
	$line =~ s/D/$primerSub{D}/gi;
	$line =~ s/N/$primerSub{N}/gi;
	return $line;
}

sub removePrimer{
	my $line2 = shift;
	my $primers = shift;
	my $count = shift;
	my $settings = shift;
	
	my $regex = join ')|(', keys %{$primers};
	$regex = '('.$regex.')';

	if ($line2 =~ s/(?'barcode'(?:.{8}|.{4}))($regex)//i){
		$$count++;
		return ($line2,$+{barcode});
	}
	return 0;
}

sub constantRegionPrimer{
	my $line2 = shift;
	my $constantRegionPrimers = shift;
	my $count = shift;
	my $settings = shift;
	
	my $regex = join ')|(', keys %{$constantRegionPrimers};
	$regex = '('.$regex.')';

	if ($line2 =~ /(?:$regex)(?'barcode'(?:.{8}|.{4}))/i){
		$$count++;
		return ($line2,$+{barcode});
	}
}

sub type{ #determines isotype for a sequence
	my $line1  = shift;
	my $line2  = shift;
	my $labels = shift;
	my $counts = shift;

	$line1 =~ s/(.*?:.*?:).*?-/$1/;
	my $verbose = 0;
	print "Seq:$line2\n" if $verbose;
	my $newlabel;
	my $dist;
	my @matches;
	foreach my $label(keys %{$labels}){
		print "Label:$label\n" if $verbose;
		if($line2 =~ /$label/){                                                # if label match is found directly
			my $matchtmp = $line1;
			my @tmpSplit = split ('', $line1);
			splice(@tmpSplit, 0, 1, ('>', uc"$labels{$label}", '_'));
			$line1 = join( '',  @tmpSplit);

			return ($line1);                                                    # add label, go to next seq
		}else{                                                                  # otherwise see if theres a match within the cut off range
			if(length($line2)<= length($label)){next;}
			my @tmpSplit  = split( '', $line2);
			my $sequence = join( '', splice( @tmpSplit, -(length($label))));
			
			my $tmpdist = distance($sequence, $label);
			if($tmpdist <= $difference){                   #levenshtein distance matched within cut off add label
				@tmpSplit = split ('', $line1);
				splice(@tmpSplit, 0, 1, ('>', uc"$labels->{$label}", '_'));
				if (!defined($dist) || $dist > $tmpdist){
					$newlabel = uc $labels->{$label};
					$dist = $tmpdist;
				}elsif($dist == $tmpdist){
					if (uc $newlabel ne uc $labels->{$label}){
 						print "Warning ambigous isotype found!\n $line1 \n $newlabel \t $labels->{$label}\n";
 					}
				}
			}
		}
		my $hold = <STDIN> if $verbose;
	}


	if ($newlabel){
		my @tmpSplit = split ('', $line1);
		splice(@tmpSplit, 0, 1, ('>', $newlabel, '_'));
		$line1 = join ('', @tmpSplit);
		return ($line1);                                              # label added go to next seq
	}
	
	my $trunc;# = truncated( [$line1,$line2], $labels, $biggestLabelSize, $difference, $counts); # check truncated labels against it.

	if(!$trunc){                      # otherwise no match and not a truncated seq, add U instead       
		my @tmpSplit = split ('', $line1);
		splice(@tmpSplit, 0, 1, ('>', 'U', '_'));
		$line1 = join( '',  @tmpSplit);   
		return ($line1);
	}else{
		return ($trunc);
	}
	
	return 0;
}

# sub type{ #determines isotype for a sequence
# 	my $line1  = shift;
# 	my $line2  = shift;
# 	my $labels = shift;
# 	my $counts = shift;
# 	$line1 =~ s/(.*?:.*?:).*?-/$1/;
# 	my $verbose = 0;
# 	print "Seq:$line2\n" if $verbose;
# 	my $newlabel;
# 	my $dist;
# 	my @matches;
# 	
# 	foreach my $label(keys %{$labels}){
# 		print "Label:$label\n" if $verbose;
# 		if($line2 =~ /$label/){                                                # if label match is found directly
# 			print "MATCH: $labels{$label}\n" if $verbose;
# 			my $matchtmp = $line1;
# 			my @tmpSplit = split ('', $line1);
# 			splice(@tmpSplit, 0, 1, ('>', uc"$labels{$label}", '_'));
# 			$line1 = join( '',  @tmpSplit);
# 
# 			$counts->{uc"$labels->{$label}"}++;
# 			return ($line1);                                                    # add label, go to next seq
# 		}else{                                                                  # otherwise see if theres a match within the cut off range
# 			if(length($line2)<= length($label)){next;}
# 			my @tmpSplit  = split( '', $line2);									
# 			my $sequence = join( '', splice( @tmpSplit, -(length($label))));	
# 
# 			my $tmpdist = distance($sequence, $label);
# 			print "STRCMP: $tmpdist\n" if $verbose;
# 			if($tmpdist <= $difference){                   #levenshtein distance matched within cut off add label
# 				@tmpSplit = split ('', $line1);
# 				splice(@tmpSplit, 0, 1, ('>', uc"$labels->{$label}", '_'));
# 				$line1 = join ('', @tmpSplit);
# 				return ($line1);
# 				
# 				
# # 				if (!defined($dist) || $dist > $tmpdist){
# # 					$newlabel = uc $labels->{$label};
# # 					$counts->{uc"$labels->{$label}"}++ if !defined($dist);
# # 					
# # 					$dist = $tmpdist;
# # 				}elsif($dist == $tmpdist){
# # 					if (uc $newlabel ne uc $labels->{$label}){
# # 						print "Warning ambigous isotype found!\n $line1 \n $newlabel \t $labels->{$label}\n";
# # 					}
# # 				}
# 			}
# 		}
# #		my $hold = <STDIN> if $verbose;
# 	}
# 
# 
# # 	if ($newlabel){
# # 		my @tmpSplit = split ('', $line1);
# # 		splice(@tmpSplit, 0, 1, ('>', $newlabel, '_'));
# # 		$line1 = join ('', @tmpSplit);
# # 		return ($line1);                                              # label added go to next seq
# # 	}
# 	
#   # otherwise no match and not a truncated seq, add U instead       
# 
# 	my @tmpSplit = split ('', $line1);
# 	splice(@tmpSplit, 0, 1, ('>', 'U', '_'));
# 	$line1 = join( '',  @tmpSplit);   
# 
# 	$counts->{U}++;
# 	return ($line1);
# 	
# }
# 


sub truncated{ #similar to above algorithm except truncating the sequence labels incase of early sequencing termination
    my $seq         = shift;
    my $labels      = shift;
    my $largest     = shift;
    my $difference  = shift;
#    my $counts      = shift;

    my %truncatedLabels;
    
    foreach my $label(keys %{$labels}){
        my @tmpSplit  = split ('', $label);
        $truncatedLabels{lc join('', splice(@tmpSplit, 0, 30) )} = $labels{$label};
    }
    if(length($seq->[1]) <= $largest){return 0;}
    my @tmpSplit  = split ('', lc$seq->[1]);
    my $tmpSeq = join('', splice(@tmpSplit, -$largest));
    
    foreach my $label (keys %truncatedLabels){
        if($tmpSeq =~ /$label/){
            @tmpSplit = split ('', $seq->[0]);
            splice(@tmpSplit, 0, 1, ('>', uc"$truncatedLabels{$label}", '_'));
            $seq->[0] = join( '',  @tmpSplit);

            #$counts->{uc $truncatedLabels{$label} }++;
            return $seq->[0];
        }else{
			my @tmpSplit  = split( '', $tmpSeq);
			my $sequence = join( '', splice( @tmpSplit, -(length($label))));

			if(distance($sequence, $label)<= $difference){ # levenshtein distance metric
				@tmpSplit = split ('', $seq->[0]);
				splice(@tmpSplit, 0, 1, ('>', uc"$truncatedLabels{$label}", '_'));
				$seq->[0] = join( '',  @tmpSplit);

				#$counts->{uc $truncatedLabels{$label} }++;
				return $seq->[0];
			}
        }
    }
    
    return 0;
}





print "> script finished: ".POSIX::strftime("%Y-%m-%d, %H:%M:%S", localtime)."\n";
close STDOUT;
1;
__END__

=head1 NAME
 
 AbIsotyper - Reads in FNA files which contain sequences in FASTA format and
 determines the Ab isotype that the sequence came from based on constant primer
 regions. Will search for exact matches to labels specified in AbIsotyper_labels.inf
 as well as truncated versions of these labels.  Also includes options for 
 sample reduction of sequences either percentage based or by a set number.

=head1 VERSION

 2013.03.29.1

=head1 DEPENDENCIES
 
 Text::LevenshteinXS
 Carp
 File::Path
 Parallel::ForkManager

=head1 SYNOPSIS
 
 perl AbIsotyper.pl <input directory> -diff:(1-X)
 perl AbIsotyper.pl <input directory> -diff:(1-X)  < -size:(1-X) OR -percent:(1-100) >
 perl AbIsotyper.pl <input directory> -diff:(1-X) -seed:(1-X)  < -size:(1-X) OR -percent:(1-100) >

=head1 DESCRIPTION
 
 This script reads a directory consisting of the FASTA formatted files and compares 
  the sequences to a list of labels in the AbIsotyper_labes.inf file. Matched labels 
  will add an identifier to the beginning of the FASTA comment line for that sequence.
  If no match is found a U will be used.  It will declare near matches up to the 
  specified cut off value.  That is if a sequence and a label are dissimilar by 1 bp 
  and the cut off is 2 it will label that sequence.
 
=head1 SETTINGS
 
 -diff:(1-X)
    Takes a number which will be used as a cutoff value for variance in similarity.
    Units is bp. (2 is suggested)

 -seed:(1-X)  OPTIONAL
 	Seed for a random number generator. Set to a specific value for 
    reproducibility.
    
 -size:(1-X)  OPTIONAL
    Takes a number which will be used as a cutoff value for sample reduction.
    
 -percent:(1-100)  OPTIONAL
    Takes a percentage which will be used as a cutoff value for sample reduction.
 
=head1 OUTPUT
 
 Generates identical files with the same respective names and the input modified
  as outlined above.  As well as additional files based on options representing a
  portion of the sequences as specified in the settings.
 
=head1 CAVETS
 
 Don't try to specify a sequence label for U.
 
=head1 REVISIONS
 
 2012.10.22
  - Integrated ability for sample size reduction.
 
 2012.10.18
  - Fixed a bug causing truncated sequences not directly matching to be	labeled as U 
   (~1% of sequences)
  - Parallelized files, should see significant speed up while doing multiple files
  - Memory utilization improved

=head1 AUTHOR
 
 Christopher Fucile
 
 Division of Allergy, Immunology and Rheumatology
 University of Rochester Medical Center
 Copyright (c) 2012, All Rights Reserved
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
 subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all copies or substantial
 portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
=cut
