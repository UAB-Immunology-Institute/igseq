#!/usr/bin/env perl

#use strict;
use warnings;

use Data::Dumper;
use DBI;
use DBD::SQLite;
use File::Basename;
use DBI qw(:sql_types);

my $settingsfile = pop @ARGV;
my $queryFile = pop @ARGV;
my @dbs = @ARGV;

my $query;
my $type;

if ($queryFile =~ /2_.*?/){
	$type = 'nt';
	$query  = "select sequenceID, lineageID, population, Vgene, Dgene, Jgene, NT_CDR3, NT_VDJ from sequence where Vgene like ? and Dgene like ? and Jgene like ? and length(NT_CDR3) = ?";
	print "Recognized nucleotide sequences...\n";
}elsif($queryFile =~ /4_.*?/){
	$type = 'aa';
	$query = "select sequenceID, lineageID, population, Vgene, Dgene, Jgene, AA_CDR3, AA_VDJ from sequence where Vgene like ? and Dgene like ? and Jgene like ? and length(AA_CDR3) = ?";
	print "Recognized amino acid sequences...\n";
}else{
	usage();
}

my %settings;
open my $fh, '<', $settingsfile or die "ERROR: can't open lineage.inf\n$!\n";
while (<$fh>){
    chop;
    my ($non_comment) = split '#';
    my @row = split ':', $non_comment;
    $settings{$row[0]} = $row[1];
}
die "ERROR: Problem with settings file (lineage.inf). Check usage.\n" if scalar keys %settings !=13;
close $fh;

open $fh, '<', $queryFile or die "Error $!\n";
my $header = <$fh>;
my @querySeqs = <$fh>;
close $fh;


open my $OPF, '>', 'sequenceQuery_output.txt' or die "Error $!\n";
print $OPF "QueryID\tDB\tHit_SeqID\tLineageID\tPopulation\tVgene\tDgene\tJgene\tCDR3\tVDJ\tSimilarity\n";
my @results;
#"select sequenceID, lineageID, population, Vgene, Dgene, Jgene, NT_CDR3 from sequence where Vgene = ? and Jgene = ? and length(NT_CDR3) = ?"

foreach my $querySeq (@querySeqs){
	my @line = split(/\t/, $querySeq);
	my ($Vgene) = $line[3] =~ /(IG\wV.*?)\*.*?\s/;
	my ($Jgene) = $line[4] =~ /(IG\wJ.*?)\*.*?\s/;
	my ($Dgene) = $line[5] =~ /(IG\wD.*?)\*.*?\s/;
	my $CDR3 	= $line[14];
	my $count 	= 0;

	foreach my $db (@dbs){
		my $db_name = basename($db);
		#print "Searching $db_name...\n";
		my $dbh = DBI->connect("DBI:SQLite:dbname=$db", "", "", {RaiseError =>1, AutoCommit => 0}) or die $DBI::errstr;

		my $sth = $dbh->prepare($query);

		if( 	$settings{'useVDJ'}       eq 'n' ||
			   ($settings{'v_Variations'} eq 'y' &&
			    $settings{'d_Variations'} eq 'y' &&
			    $settings{'j_Variations'} eq 'y' )){

				$sth->bind_param( 1, '%', SQL_VARCHAR );
				$sth->bind_param( 2, '%', SQL_VARCHAR );
				$sth->bind_param( 3, '%', SQL_VARCHAR );
				$sth->bind_param( 4, length($CDR3), SQL_INTEGER );

		}elsif( $settings{'useVDJ'} 	  eq 'y' &&
			   ($settings{'v_Variations'} eq 'y' &&
			    $settings{'d_Variations'} eq 'y' &&
			    $settings{'j_Variations'} eq 'y' )){

				die "Error in settings.  useVDJ can not be set to 'y' while v, d and j variations are allowed.\n";
			
		}elsif( $settings{'useVDJ'} 	  eq 'y' &&
			   ($settings{'v_Variations'} eq 'n' &&
			    $settings{'d_Variations'} eq 'n' &&
			    $settings{'j_Variations'} eq 'n' )){

				#$query = "select sequenceID, lineageID, population, Vgene, Dgene, Jgene, NT_CDR3 from sequence where Vgene = ? and Dgene = ? and Jgene = ? and length(NT_CDR3) = ?";
				$sth->bind_param( 1, $Vgene, SQL_VARCHAR );
				$sth->bind_param( 2, $Dgene, SQL_VARCHAR );
				$sth->bind_param( 3, $Jgene, SQL_VARCHAR );
				$sth->bind_param( 4, length($CDR3), SQL_INTEGER );

		}elsif(	$settings{'useVDJ'} 	  eq 'y' &&
		   	    $settings{'v_Variations'} eq 'y' ){
			if( $settings{'d_Variations'} eq 'y' ){

				#$query = "select sequenceID, lineageID, population, Vgene, Dgene, Jgene, NT_CDR3 from sequence where Jgene = ? and length(NT_CDR3) = ?";
				$sth->bind_param( 1, '%', SQL_VARCHAR );
				$sth->bind_param( 2, '%', SQL_VARCHAR );
				$sth->bind_param( 3, $Jgene, SQL_VARCHAR );
				$sth->bind_param( 4, length($CDR3), SQL_INTEGER );

			}elsif($settings{'j_Variations'} eq 'y'){

				#$query = "select sequenceID, lineageID, population, Vgene, Dgene, Jgene, NT_CDR3 from sequence where Dgene = ? and length(NT_CDR3) = ?";
				$sth->bind_param( 1, '%', SQL_VARCHAR );
				$sth->bind_param( 2, $Dgene, SQL_VARCHAR );
				$sth->bind_param( 3, '%', SQL_VARCHAR );
				$sth->bind_param( 4, length($CDR3), SQL_INTEGER );

			}else{

				#$query = "select sequenceID, lineageID, population, Vgene, Dgene, Jgene, NT_CDR3 from sequence where Dgene = ? and Jgene = ? and length(NT_CDR3) = ?";
				$sth->bind_param( 1, '%', SQL_VARCHAR );
				$sth->bind_param( 2, $Dgene, SQL_VARCHAR );
				$sth->bind_param( 3, $Jgene, SQL_VARCHAR );
				$sth->bind_param( 4, length($CDR3), SQL_INTEGER );
			}

		}elsif(	$settings{'useVDJ'} 	  eq 'y' &&
		   	    $settings{'d_Variations'} eq 'y' ){
			if( $settings{'j_Variations'} eq 'y'){

				#$query = "select sequenceID, lineageID, population, Vgene, Dgene, Jgene, NT_CDR3 from sequence where Vgene = ? and length(NT_CDR3) = ?";
				$sth->bind_param( 1, $Vgene, SQL_VARCHAR );
				$sth->bind_param( 2, '%', SQL_VARCHAR );
				$sth->bind_param( 3, '%', SQL_VARCHAR );
				$sth->bind_param( 4, length($CDR3), SQL_INTEGER );

			}else{

				#$query = "select sequenceID, lineageID, population, Vgene, Dgene, Jgene, NT_CDR3 from sequence where Vgene = ? and Jgene = ? and length(NT_CDR3) = ?";
				$sth->bind_param( 1, $Vgene, SQL_VARCHAR );
				$sth->bind_param( 2, '%', SQL_VARCHAR );
				$sth->bind_param( 3, $Jgene, SQL_VARCHAR );
				$sth->bind_param( 4, length($CDR3), SQL_INTEGER );
			}

		}elsif(	$settings{'useVDJ'} 	  eq 'y' &&
		   	    $settings{'j_Variations'} eq 'y' ){
				#$query = "select sequenceID, lineageID, population, Vgene, Dgene, Jgene, NT_CDR3 from sequence where Vgene = ? and Dgene = ? and length(NT_CDR3) = ?";
				$sth->bind_param( 1, $Vgene, SQL_VARCHAR );
				$sth->bind_param( 2, $Dgene, SQL_VARCHAR );
				$sth->bind_param( 3, '%', SQL_VARCHAR );
				$sth->bind_param( 4, length($CDR3), SQL_INTEGER );

		}else{

		 die "Error. Problem with settings file $!\n";
		}

		$sth->execute();
		my $r = $sth->fetchall_arrayref();

		for (0..$#{$r}){
			my @tmp = @{ $r->[$_] };
			my $score = hamming($tmp[6], $CDR3);
			my $sim = 100 * (length($tmp[6]) - $score) / length($CDR3);
			if ($sim >= $settings{'cdr3_homology'}){
					$count++;
	                print $OPF join("\t", $line[1], $db_name, @tmp, $sim)."\n";
	                #map { push @{$results->{$i}}, $_ } @{$preclusters{$sorted_seq[$j]}{'indices'}};
	                #splice @sorted_seq, $j, 1;
	        }
	    }
	    $dbh->disconnect;
	}
	print "$count hits\n";
	#print $OPF join("\t", @results)."\n";

}


sub hamming {# subroutine returns number of characters that are different in two strings of equal length
	my $a = shift;
	my $b = shift;
	die "ERROR: strings must be same length" if length($a) != length($b);
	my $n = length($a);
	my $h = 0;
	for my $k (0 .. $n - 1) { $h++ if substr($a, $k, 1) ne substr($b, $k, 1); }
	return $h;
}

sub usage{
	print "perl sequenceQuery.pl <'db'>... <'Vquest 2 or 4 file'> <'lineage.inf'>\n";
	exit;
}