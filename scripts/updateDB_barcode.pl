#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use DBI;
use DBD::SQLite;
use File::Spec;
use Scalar::Util qw(looks_like_number);

my ($db,$barcodeFolder) = @ARGV;

print "Please backup any database before proceeding!!!\n";
print "Proceed (y/n): ";
my $input = <STDIN>;
chomp $input;
exit if $input ne 'y';

if (-f $db){
	my $dbh = DBI->connect("DBI:SQLite:dbname=$db", "", "", {RaiseError =>1, AutoCommit => 0}) or die $DBI::errstr;
	print "DB connected, updating...\n";

	opendir my $IND, $barcodeFolder or die "Error reading directory\n$!\n" if -d $barcodeFolder;
	my @files = grep {! /^\./ } readdir $IND;

	my $sql = "PRAGMA user_version";
	my $sth = $dbh->prepare($sql) or die $DBI::errstr;
	$sth->execute() or die $DBI::errstr;
	my $version = $sth->fetchrow_array;
	$sth->finish;

	if( !looks_like_number($version) || $version < 1){
		$sql = "ALTER TABLE sequence ADD COLUMN Forward_Barcode TEXT";
		$dbh->do($sql) or die $DBI::errstr;

		$sql = "ALTER TABLE sequence ADD COLUMN Reverse_Barcode TEXT";
		$dbh->do($sql) or die $DBI::errstr;

		$sql = "PRAGMA user_version = 1";
		$dbh->do($sql) or die $DBI::errstr;
	}

	$sql = "UPDATE sequence SET Forward_Barcode = ?, Reverse_Barcode = ? WHERE sequenceID = ?";
	$sth = $dbh->prepare($sql) or die $DBI::errstr;

	foreach my $file (@files){
		open my $INF, '<', File::Spec->catfile($barcodeFolder,$file) or die "Error opening file\n$!\n" if -f File::Spec->catfile($barcodeFolder,$file);
		my @lines = <$INF>;
		chomp @lines;

		foreach my $line (@lines){
			my ($seqID,$barcodes) = split('\t',$line);
			my ($fwd,$rev) = split('\|',$barcodes);
			$sth->execute($fwd,$rev,$seqID);
		}
	}

	$dbh->do('commit');
	$dbh->disconnect();
}



sub usage {

}