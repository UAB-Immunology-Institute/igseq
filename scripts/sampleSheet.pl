#!/usr/bin/env perl

###################################################################################################
#
# sampleSheet.pl 
#
#
#
# by Christopher Fucile
#
######################################
#
# Informatics Institute
# University of Alabama at Birmingham
# Copyright (c) 2017, All Rights Reserved
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
#
###################################################################################################

use strict;
use warnings;

use File::Find;
use File::Spec;
use Data::Dumper;
use FindBin;

our $VERSION = '2017.01.20.1';

use 5.012;

my $path = $FindBin::RealBin;
my ($vol, $dir, $files) = File::Spec->splitpath($path);

my $inputDIR = $ARGV[0];

my @files;
find(\&wanted, $inputDIR);
$dir = File::Spec->catdir($dir, 'config');
open my $OPF, ">", File::Spec->catpath($vol, $dir, 'SampleSheet.csv') or die "Error opening output file $!\n";

my $first = 1;
foreach my $file (@files){
	open my $INF, "<", $file or die "Error opening file $!\n";
	my $line1 = <$INF>;
	my @rest = <$INF>;
	print $OPF $line1 if $first;
	print $OPF @rest;
	$first = 0;
}


sub wanted{
	push @files, $File::Find::name if /SampleSheet\.csv/;
}

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
