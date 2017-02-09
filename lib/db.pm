package db;
###################################################################################################
#
# summary.pm
#
# summarizes IgSeq analysis data
#
# by Christopher Fucile
########################
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
#    use summary;
#
###################################################################################################
use strict;
use warnings;

use File::Path qw(make_path);
use Data::Dumper;




my %statements;
open my $SQL, "<". File::Spec->catfile($FindBin::Bin,'lib','sql.txt') or die "Error opening SQL statement file $!\n";
my $statement;
while (<$SQL>){
	chomp;
	if ($_ =~ '^\/\*'){
		($statement) = $_ =~ '\/\*(.*?)\*\/';
	}else{
		$statements{$statement} .= $_;
	}
}
close $SQL;


sub createSequenceTable {
	my $dbh = shift;
	
	my $sql = "DROP TABLE IF EXISTS sequence";
	$dbh->do($sql) or die $DBI::errstr;

	$sql = $statements{'createSequenceTable'};

	$dbh->do($sql) or die $DBI::errstr;
	
	return 0;
}

sub insertSequence {
	my $dbh = shift;
	my $data = shift;
	
	my $sql = $statements{'insertSequence'};

	$dbh->do($sql, @{$data}) or return $DBI::errstr;
	
	return 0;
}

sub createIndicies {
	my $dbh = shift;
	
	my $sql1 = "CREATE INDEX idx1 ON sequence(isotype)";
	my $sql2 = "CREATE INDEX idx2 ON sequence(VGene)";
	my $sql3 = "CREATE INDEX idx3 ON sequence(JGene)";
	my $sql4 = "CREATE INDEX idx4 ON sequence(lineageID)";
	my $sql5 = "CREATE INDEX idx5 ON sequence(population)";
	my $sql6 = "CREATE INDEX idx6 ON sequence(population, lineageID)";
	my $sql7 = "CREATE INDEX idx7 ON sequence(population, isotype)";
	my $sql8 = "CREATE INDEX idx8 ON sequence(population, isotype, lineageID)";
		
	$dbh->do($sql1) or return $DBI::errstr;
	$dbh->do($sql2) or return $DBI::errstr;
	$dbh->do($sql3) or return $DBI::errstr;
	$dbh->do($sql4) or return $DBI::errstr;
	$dbh->do($sql5) or return $DBI::errstr;
	$dbh->do($sql6) or return $DBI::errstr;
	$dbh->do($sql7) or return $DBI::errstr;
	$dbh->do($sql8) or return $DBI::errstr;
}

sub updateSequence {
	my $dbh = shift;
	my $sql = "";
	
	$dbh->do($sql) or return $DBI::errstr;
	
	return 0;
}


sub selectSequences {
	my $dbh = shift;
	my $where = shift;
	my $sql = $statements{'selectSequences'} . $where;

	my $sth = $dbh->prepare($sql);
	my $results = $sth->execute() or die $DBI::errstr;

	return $sth->fetchall_arrayref();
}


sub deleteSequence {
	my $dbh = shift;
	my $sql = "";
	
	$dbh->do($sql) or return $DBI::errstr;
	
	return 0;
}

1;