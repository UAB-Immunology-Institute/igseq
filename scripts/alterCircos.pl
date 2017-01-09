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

usage("Error with file $files[0] $!\n") if ! -f $files[0];
usage("Error with directory $files[1] $!\n") if ! -d $files[1];

my $dbh = DBI->connect("DBI:SQLite:dbname=$files[0]", "", "", {RaiseError=>1, AutoCommit =>0}) or die $DBI::errstr;

my $query = "select lineageID, isotype, count(*) from sequence group by lineageID, isotype";

my $dbSeq = $dbh->selectall_hashref($query, ['lineageID', 'isotype'] );
print Dumper $dbSeq;


sub usage{
	my $message = shift;
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