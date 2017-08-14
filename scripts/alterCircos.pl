#!/usr/bin/env perl

###################################################################################################
#
#
#
###################################################################################################

use strict;
use warnings;
use Data::Dumper;
use DBI;
use DBD::SQLite;

my @files = @ARGV;
my $verbose = 1;

print Dumper @files if $verbose;


my $pop = 'AN1';

usage("Error, can not read file $files[0] $!\n") if ! -f $files[0];
usage("Error, can not read directory $files[1] $!\n") if ! -d $files[1];
usage("Error, no isotype information supplied $!\n") if ! $files[2];

$pop = $files[2];

my $dbh = DBI->connect("DBI:SQLite:dbname=$files[0]", "", "", {RaiseError=>1, AutoCommit =>0}) or die $DBI::errstr;

my $query = "select lineageID, population, count(*) from sequence where population == ? group by lineageID, population";
print Dumper $query if $verbose;

print "Querying Database...\n";
my $sth = $dbh->prepare($query);
$sth->execute($pop);
my $dbSeq = $sth->fetchall_hashref(['lineageID', 'population'] );
print Dumper scalar keys %{$dbSeq} if $verbose;

open my $INF, "<". "$files[1]/circos_lineage_links.txt" or die "Error $!\n";
my @links = <$INF>;
close $INF;
chomp @links;

mkdir "circos_output";
open my $OPF, ">". "circos_output/circos_lineage_links_OUTPUT.txt" or die "Error $!\n";
print "Filtering Circos data...\n";
my @matchingSeqGroups;
LINEAGE:
	foreach my $lineageID (keys %{$dbSeq}){
		#next if scalar(keys $dbSeq->{$lineageID}) < 2;
		foreach my $population (keys %{$dbSeq->{$lineageID}}){
			#next if ! grep {/E/} keys %{$dbSeq->{$lineageID}{$population}};
			push @matchingSeqGroups, [$lineageID, $population];
		}
	}

foreach my $group (@matchingSeqGroups){
	my @matches = grep {/$group->[1].*?#$group->[0]$/} @links;
	if (@matches){
		print $OPF join "\n", @matches;
		print $OPF "\n";
	}
}

close $OPF;
$dbh->disconnect;

print "Done.\n";

sub usage{
	my $message = shift;
	print "Usage:\nperl alterCircos <SQLite database file> <Directory of circos files> <Isotype>\n";
	die "$message";
}


=head1 NAME
 

 
=head1 SYNOPSIS
 

 
=head1 DESCRIPTION
 
 
=head1 SETTINGS
 
 
=head1 OUTPUT
 
 
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