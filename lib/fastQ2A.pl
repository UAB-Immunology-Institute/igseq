#!/usr/bin/perl

###################################################################################################
#
# fastQ2A.pl 
#
# Converts FASTQ to FASTA format for join files from fastq-join part of ea-utils
#
# by Christopher Fucile
#
######################################################
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
# USAGE:
#
#    perl fastQ2A.pl <input file>
#
###################################################################################################

use strict;
use warnings;

use Carp;
use File::Basename;
use File::Path qw(make_path);
use FindBin;
use Data::Dumper;

our $VERSION = '2012.10.23.1';

if(!(@ARGV)){ croak "ERROR: no input directory\n";}

my $input = $ARGV[0];
my $base = basename($input);
my $dir = $FindBin::RealBin;
$dir =~ s/\/lib/\/fastaOutput/;

make_path($dir);

print "> Converting file $input\n";

open INF, "<". $input or die "Can not open input file $input: $!\n";

$base =~ s/\.fastq$/\.fasta/;
	

open OPF, "> $dir/$base" or die "Can not open output file $dir/$base: $!\n";
	
my $lineNumber = 1;
	
while(<INF>){
	my $line1 = $_;
	my $line2 = <INF>;
	my $line3 = <INF>;
	my $line4 = <INF>;

	if ($line1 !~ /^\@/){die "Error in fastq file: New sequence entry not starting in \@\n";}
	$line1 =~ s/^\@/\>/;
	print OPF $line1;
	print OPF $line2;
}

close INF;
close OPF;


1;

=head1 NAME
 
 fastQ2A.pl
 
=head1 SYNOPSIS
 
 perl fastQ2A.pl <input file>
 
=head1 DESCRIPTION
 
 Converts FASTQ to FASTA format for join files from fastq-join part of ea-utils
 
=head1 SETTINGS
 
 None
 
=head1 OUTPUT
 
 FASTA files with the same name as the input.
 
=head1 CAVEATS
 
 
 
=head1 AUTHOR
 
 Christopher Fucile
 University of Rochester Medical Center
 Division of Allergy, Immunology & Rheumatology
 Copyright (c) 2013, All Rights Reserved
 
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