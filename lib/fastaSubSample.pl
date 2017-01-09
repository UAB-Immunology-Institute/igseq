#!/usr/bin/perl

###############################################################################
#
# fastaSubSample.pl 
#
# Randomly selects an approximate number of sequences from a FASTA file. Can
#	also break the file into chunks.
#
# by Christopher Fucile
#
######################################
#
# Division of Allergy, Immunology and Rheumatology
# University of Rochester Medical Center
# Copyright (c) 2012, All Rights Reserved
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
#    perl fastaSubSample.pl <input dir>
#    perl fastaSubSample.pl <input dir> < -size:(1-X) > -subsample:(y|n)
#    perl fastaSubSample.pl <input dir> < -size:(1-X) > -subsample:(y|n) -seed:(1-X)
#
# DEPENDS ON:
#	Parallel::ForkManager
#
###############################################################################
our $VERSION = '2013.03.29.1';
use 5.012;

use strict;
use warnings;

use File::Path qw(make_path);
use POSIX;
use Carp;
use Parallel::ForkManager;
use Data::Dumper;
use FindBin;

my $dir = $FindBin::RealBin;
$dir =~ s/\/lib//;

my $inputDir = $ARGV[0];
my $seed;
my $cutoff;
my $subsample;

foreach my $arg(@ARGV){    
	$cutoff = $1       if  $arg =~ /\-size:(\d+)/;
	$seed = $1		   if  $arg =~ /\-seed:(\d+)/;
	$subsample = $1	   if  $arg =~ /\-subsample:(y|n)/;
}

$inputDir =~ s|(\/$)||;
croak "Error: Please remove whitespaces from directory name. $inputDir\n" if $inputDir =~ /\s+/;

opendir DIR, "$inputDir" or croak "ERROR: can not open directory $inputDir: $!";

my @files;
my @tmpfiles = readdir DIR;
my $numOfJunk = grep{ /^\./} @tmpfiles;
splice @tmpfiles, 0, $numOfJunk;
unshift @tmpfiles, '';
push @files, split ' ', join " $inputDir/", @tmpfiles;
@files = grep{ /\.fasta$/} @files;
closedir DIR;

#my $outDir = ($inputDir =~ /\/(.*)/) ? $1 : $inputDir;
my $outname;
my @MONTH = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my $LOG_FILE = "fastaSubSample_LOG.txt";
my $RUN_TIME = POSIX::strftime("%Y-%m-%d, %H:%M:%S", localtime);

my $outDir = "$dir/fastaSubSample\_Out/";
make_path($outDir);

open STDOUT, "| tee -a $outDir'/'$LOG_FILE" or croak "Tee stdout error";
#print "*********************************************\n";
print "> script start: $RUN_TIME\n";

my $pm = new Parallel::ForkManager(scalar(@files));
print "> Processing files...\n";
foreach my $file(@files){
    srand($seed) if $seed;
	#my $pid = $pm->start and next;
	open INF, "<$file" or die "Error opening input file: $file\n $!\n";
	if(<INF> !~ /^\>/){die "File does not appear to be FASTA.\n";}
	seek INF, 0, 0;
	my $size = grep {/^\>/} <INF>;

	if($cutoff >= $size){warn "Warning size of file too small. Specify smaller amount of sequences in order to sub-set. $file\n";}
	
	my $numFiles += int($size/$cutoff);
	$numFiles++ if $size % $cutoff;
    
    my %file_assignments;
    for my $file_number(1..$numFiles-1){
        for (1..$cutoff){
            my $num = int(rand($size*2));
            while (($num % 2) == 0 || exists($file_assignments{$num})){
                $num = int(rand($size*2));
            }
            #print "$_\t$num\n";
            $file_assignments{$num} = $file_number;
        }
    }
    
    for (my $line_num=1; $line_num < $size*2; $line_num += 2){
        $file_assignments{$line_num} = $numFiles if !exists($file_assignments{$line_num});
    }
    
	for (1..$numFiles){
        my $i = 1;
		my $fileNum = $_;
		my $output = ($file =~ /.*\/(.*)\./) ? $1 : "output.txt";
		open OPF, ">", "$outDir$output\_FILE_$fileNum.fasta" or die "Error opening output file:$outDir/$output\_FILE_$fileNum.fasta\n  $!\n";
		
		seek INF, 0, 0;
		my $flag = 0;
        
		while(<INF>){
			my $line1 = $_;
			my $line2 = <INF>;
			if($fileNum == $file_assignments{$i}){
				print OPF $line1;
				print OPF $line2;
			}
            $i += 2;
		}

		close OPF;
		print "> Output written to $outDir$output\_FILE_$fileNum.fasta\n";
	}
	close INF;
	$pm->finish(0);
}
$pm->wait_all_children;



print "> script finished: ".POSIX::strftime("%Y-%m-%d, %H:%M:%S", localtime)."\n";
close STDOUT;
1;

__END__

=head1 NAME
 
 fastaSubSample.pl 
   - Randomly selects a number of sequences from a FASTA file.

=head1 VERSION

 2013.03.29.1
 
=head1 DEPENDENCIES
 
 Parallel::ForkManager

=head1 SYNOPSIS
 
 perl fastaSubSample.pl <input directory>  < -size:(1-X) >
 perl fastaSubSample.pl <input directory>  < -size:(1-X) > -seed:(1-X)
 
=head1 DESCRIPTION
 
 This script reads a directory consisting of the FASTA formatted files and will
 randomly select a subset 
 
=head1 SETTINGS
 
 -seed:(1-X)  OPTIONAL
 	Seed for a random number generator. Set to a specific value for 
    reproducibility.
    
 -size:(1-X)  OPTIONAL
    Takes a number which will be used as a cutoff value for sample reduction.
    
 -percent:(1-100)  OPTIONAL
    Takes a percentage which will be used as a cutoff value for sample reduction.
 	 
=head1 OUTPUT
 
 Generates identical files with the the number of sequences added to the name.
 
=head1 CAVETS
 
 
 
=head1 REVISIONS
 

=head1 AUTHOR
 
 Christopher Fucile
 
 Division of Allergy, Immunology and Rheumatology
 University of Rochester Medical Center
 Copyright (c) 2012, All Rights Reserved
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute, license,
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

