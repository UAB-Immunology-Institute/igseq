#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use File::Spec;
use Cwd;
use DBI;
use DBD::SQLite;

my $dir = cwd;
my $query = pop @ARGV;
my @dbs = @ARGV;

usage() if $query !~ /select.*/i;

print "Working...\n";
foreach my $db (@dbs){
	my $dbh = DBI->connect("DBI:SQLite:dbname=$db", "", "", {RaiseError =>1, AutoCommit => 0}) or die $DBI::errstr;
	my $db_name = fileparse($db);
	$db_name =~ s/\.db/\.fasta/;
	my $sth = $dbh->prepare($query) or return $DBI::errstr;
	$sth->execute();
	my $results = $sth->fetchall_arrayref();

	open my $OPF, '>', File::Spec->catfile($dir,"$db_name");
	foreach my $line (@{$results}){
		my $seq = pop @{$line};
		my $id = join('|', @{$line});
		print $OPF ">". $id  ."\n";
		print $OPF uc $seq . "\n";
	}
	close $OPF;

	$dbh->disconnect();
}

sub usage{
	print "perl extractSequences.pl <path to database>... '<sql query>'\n";
	exit;
}