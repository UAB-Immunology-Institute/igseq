package lineage;
require Exporter;

@ISA = qw(Exporter);
@EXPORT_OK = qw(lineage);


###################################################################################################
#
# lineage.pm
#
# takes a merged sequence file (output of merge454.pl) and computes lineages
#
#
# by Alex Rosenberg
#
# June 23, 2011
#
######################################
# Modified December 5, 2012:
#
# Added ability to ignore one wrong V, D or J sequence_id mismatch
# Added settings file lineage.inf containing user settings (see below)
# New output: histogram output of sequence length vs number of sequences 	X_SEQvsCOUNT_NO_SINGLETONS.txt  
#                                                                           X_SEQvsCOUNT_WITH_SINGLETONS.txt
# New output: list of lineage sizes for each lineage individually           X_HISTOGRAM_Second.txt
# New output: list of representative sequences from each lineage            X_LINEAGE_REPS.txt
# New output: list of settings from the run                                 X_RUN_SETTINGS.txt
# New output: list of unique V id's and how many sequences are present      X_UNIQUE_V_COUNTS.txt
# New output: list of mutation rates per lineage                            X_LINEAGE_MUTATION.txt
#
# Parallelized using Parallel::ForkManager for determining similarieties between sequences
#
# More memory efficient.  Eliminated Graph.pm dependency which resulted in high memory usage.
#
# by Christopher Fucile
#
# Modified March 27, 2013:
#
# Change of core algorithm from graph based approach to agglomerative clustering
#  approach.  Resulting in significant increases of speed and memory efficiencies.
# Results are slightly affected by algorithm changes mostly due to tenuous
#  relationships constructed in the graph based approach being broken, resulting
#  in a small number of the larger clusters to be split down into smaller pieces.
#
# by Christopher Fucile
######################################
#
# Laboratory of Ignacio Sanz
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
#    lineage(<input file>, <output file>)
#
###################################################################################################
our $VERSION = '2012.11.7.1';

# std libs
    use strict;
    use POSIX;
    use Carp;
    use File::Path qw(make_path);
    use File::Basename;
    use Tie::File;
    use Data::Dumper;
    use FindBin;
    use File::Find;
    use File::Spec;
    use DBI;
	use lib "$FindBin::Bin/";

# dependencies
    use Statistics::Basic qw(stddev);
    use Parallel::ForkManager;
    use Sys::CPU qw(cpu_count);
    use Devel::Size qw(total_size size);
    use Devel::Symdump;
    use List::MoreUtils qw(indexes);
    use DBD::SQLite;
    use Algorithm::Cluster qw(treecluster);
    use circos_config;
    use expanded_clones;
    

# global constants
    
    my $verbose = 0;
    my $sampling = 50000;
    my $singletons = 1;
    
    
sub lineage { 
	
	my $path = $_[2];
    my $outDir = File::Spec->catdir($path,'lineageOutput');
    make_path($outDir);
    
# time stamp for running script and log file

	my @MONTH = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	my $LOG_FILE = 'lineage454_log.txt';
	my $RUN_TIME = POSIX::strftime("%Y-%m-%d, %H:%M:%S", localtime);

# get file names
	my $db_name = $_[4];
	my $outfile = $_[3];
	my $input_file = $_[0];
	my @infileTie;
	my $header;

	my $inDir = File::Spec->catdir($path,'mergeOutput');
	opendir my $DIR, $inDir or die "Error opening folder $!\n";
	my @pops = grep {!/^\./} readdir $DIR;

	my @paths;
	foreach my $pop (@pops){
		find( sub{push @paths,$File::Find::name if /$pop\_$input_file/}, $path);
	}

	foreach my $file (@paths){
		if (-f $file){
			open my $INF, "<". $file or die "Error opening $file $!\nSample may not contain this isotype\n";
			my ($pop) = $file =~ /$inDir\/(.*?)\//;
			$header = <$INF>;
			while (<$INF>){
				push @infileTie, ($pop."\t".$_);
			}
			close $INF; 
		}else{
			print "Error opening $file\nSample may not contain this isotype\n";
			next;
		}
	}

	my $infileSize = scalar @infileTie;
	unshift @infileTie, $header;
	chomp @infileTie;

# settings file for user options
    my %settings;
    my $settingsfile = $_[1];
    open my $fh, "<". $settingsfile or croak "ERROR: can't open lineage.inf\nMake sure file is in same directory as lineage454.pl\n";
    while (<$fh>){
        chop;
        my ($non_comment) = split '#';
        my @row = split ':', $non_comment;
        $settings{$row[0]} = $row[1];
    }
    croak "ERROR: Problem with settings file (lineage.inf). Check usage.\n" if scalar keys %settings !=13;
    close $fh;


	print "*********************************************\n";
	print POSIX::strftime("%H:%M:%S", localtime)."> script start: $RUN_TIME\n";
	print POSIX::strftime("%H:%M:%S", localtime)."> Initializing data structures\n";
    $| = 0;


# initialize some data structures and determine size of file
# this will clump together groups based on variation settings
# i.e. allowing d variations creates the following hash:
# $data->{v_gene_ids}{'d'}{j_gene_ids}
# where gene_ids are each unique id present 
my $data;
my @allnodes;

for (1..$infileSize){
	my @tmp = split "\t", $infileTie[$_];

	my $junction = $settings{'useJunction'} eq 'n' ? 'junction' : junctionAnalysis($infileTie[$_]);

	if( 	$settings{'useVDJ'}       eq 'n' ||
		   ($settings{'v_Variations'} eq 'y' &&
		    $settings{'d_Variations'} eq 'y' &&
		    $settings{'j_Variations'} eq 'y' )){
		push @{$data->{'v'}{'d'}{'j'}{length($tmp[17])}{$junction} }, $_;
		push @allnodes, $_;

	}elsif( $settings{'useVDJ'} 	  eq 'y' &&
		   ($settings{'v_Variations'} eq 'y' &&
		    $settings{'d_Variations'} eq 'y' &&
		    $settings{'j_Variations'} eq 'y' )){
		die "Error in settings.  useVDJ can not be set to 'y' while v, d and j variations are allowed.\n";
		
	}elsif( $settings{'useVDJ'} 	  eq 'y' &&
		   ($settings{'v_Variations'} eq 'n' &&
		    $settings{'d_Variations'} eq 'n' &&
		    $settings{'j_Variations'} eq 'n' )){
		push @{$data->{$tmp[3]}{$tmp[4]}{$tmp[5]}{length($tmp[17])}{$junction} }, $_;
		push @allnodes, $_;

	}elsif(	$settings{'useVDJ'} 	  eq 'y' &&
	   	    $settings{'v_Variations'} eq 'y' ){
		
		if( $settings{'d_Variations'} eq 'y' ){
			push @{$data->{'v'}{'d'}{$tmp[5]}{length($tmp[17])}{$junction} }, $_;
			push @allnodes, $_;
		}elsif($settings{'j_Variations'} eq 'y'){
			push @{$data->{'v'}{$tmp[5]}{'j'}{length($tmp[17])}{$junction} }, $_;
			push @allnodes, $_;
		}else{
			push @{$data->{'v'}{$tmp[5]}{$tmp[4]}{length($tmp[17])}{$junction} }, $_;
			push @allnodes, $_;
		}


	}elsif(	$settings{'useVDJ'} 	  eq 'y' &&
	   	    $settings{'d_Variations'} eq 'y' ){
		
		if( $settings{'v_Variations'} eq 'y' ){								#
			push @{$data->{'v'}{'d'}{$tmp[4]}{length($tmp[17])}{$junction} }, $_;		# can probably be removed
			push @allnodes, $_;												#
		}elsif($settings{'j_Variations'} eq 'y'){
			push @{$data->{$tmp[3]}{'d'}{'j'}{length($tmp[17])}{$junction} }, $_;
			push @allnodes, $_;
		}else{
			push @{$data->{$tmp[3]}{'d'}{$tmp[4]}{length($tmp[17])}{$junction} }, $_;
			push @allnodes, $_;
		}


	}elsif(	$settings{'useVDJ'} 	  eq 'y' &&
	   	    $settings{'j_Variations'} eq 'y' ){
		
		if( $settings{'d_Variations'} eq 'y' ){								#
			push @{$data->{$tmp[3]}{'d'}{'j'}{length($tmp[17])}{$junction} }, $_;		# can probably be removed
			push @allnodes, $_;												#
		}elsif($settings{'v_Variations'} eq 'y'){							#
			push @{$data->{'v'}{$tmp[5]}{'j'}{length($tmp[17])}{$junction} }, $_;		#
			push @allnodes, $_;												#
		}else{
			push @{$data->{$tmp[3]}{$tmp[5]}{'j'}{length($tmp[17])}{$junction} }, $_;
			push @allnodes, $_;
		}


	}else{
	 die "Error. Problem with settings file $!\n";
	}
}

print POSIX::strftime("%H:%M:%S", localtime)."> Number of sequences read: ".scalar(@allnodes)."\n";

# more data structure initialization
my $numCC;
my $numSC;
my $numEdges;
my %singletons = map {$_ =>1} @allnodes;
my %hist;
my %seqLength;
my $num = cpu_count();
my $popmut;
my $total_results;
$total_results->{'settings'} = \%settings;
@{$total_results->{'pops'}} = grep {!/^\./} @pops;

# A.) loop through each array which is grouped by v gene id, j gene id and cdr3
#			length.
#	1. first level of hash is v gene
#		2. second level is d gene
#			3. third level is j gene
#               4. fourth level is cdr3 length
#                   5. fifth level is the array filled with the line indicies of
#                       sequences matching that criteria
#

print POSIX::strftime("%H:%M:%S", localtime)."> Determining connections...\n";
foreach my $vGene (sort keys %{$data}){ # A1.
	my $Vcount=0;
	
	print POSIX::strftime("%H:%M:%S", localtime)."> Analyzying $vGene\n" if $verbose;
	
	foreach my $dGene (keys %{$data->{$vGene}}){# A2.
	
		print POSIX::strftime("%H:%M:%S", localtime)."> Analyzing $dGene\n" if $verbose;
	
		foreach my $jGene (keys %{$data->{$vGene}{$dGene}}){ # A3.
			
			print POSIX::strftime("%H:%M:%S", localtime)."> Analyzing $jGene\n" if $verbose;
			
			foreach my $length (keys %{$data->{$vGene}{$dGene}{$jGene}}){ # A4.
				print POSIX::strftime("%H:%M:%S", localtime)."> CDR3 length: $length\n" if $verbose;


				
				#print Dumper @seqs;
				#my $junctionGroups = junctionAnalysis(\@seqs, \%settings);
				
				foreach my $jgroup (keys %{$data->{$vGene}{$dGene}{$jGene}{$length}}){
				
					my @tmp;
					my $adjList;
					foreach (@{$data->{$vGene}{$dGene}{$jGene}{$length}{$jgroup}}){ # A5.
	
						#my @row = split "\t", $infileTie[$_];
						#push @tmp, $row[30];
						push @tmp, $infileTie[$_];
						$Vcount++;
					}
				
					my $results;
					my $mut_results;

					if(scalar(@tmp)>0){

						###########################
						#$mut_results = mut_cluster(\@tmp, \%settings);
						$results = cluster(\@tmp, \%settings);
						#$results = cluster_2(\@tmp, \%settings);

						###########################

					}else{
						print Dumper @tmp;
						croak "Error: nothing in array $!\n";
					}
				
					#if ($mut_results && scalar(@{$mut_results})>0){
					#	map {print $TMP_OUT $_->[0]."\t".$_->[1]."\n"; } @{$mut_results};
					#}

					if ($results){ 
						foreach my $cluster (keys %{$results}){
							if(scalar(@{$results->{$cluster}}) > 1){ #non-singletons
								$total_results->{'numCC'}++;
								my $mut;
								my @tmpMut;
								my $v_gene;
								foreach my $tmpnode (@{$results->{$cluster}}){
									my $node = $data->{$vGene}{$dGene}{$jGene}{$length}{$jgroup}[$tmpnode];
									my @tmpline = split '\t', $infileTie[$node];

									chomp($infileTie[$node]);
									$infileTie[$node] .= "\t$jgroup";
									push @{$total_results->{'clusters'}{$total_results->{'numCC'}}}, $node;
								
									#print $fh_info "$numCC\t$infileTie[$node]\n";
								
									my @seqID = $tmpline[1] =~ /(.*?):.*:(\d+:\d+)$/;
									$seqID[1] =~ s/:/\-/ if $seqID[1];

									#print $fh_fasta ">$numCC|" . join('|',@seqID) . "\n";
									(my $seqtmp = $tmpline[11]) =~ s/\.//g;
									$total_results->{'seqLength'}{'lineages'}{length($seqtmp)}++;
									#print $fh_fasta $seqtmp . "\n";
								
									my $vMut = ( $tmpline[38] =~ /(\d+)(\()?/ ) ? $1 : undef;
									my $vLength = ($tmpline[36] =~ /(\d+)(\()?/ ) ? $1 : undef;
									if ($vLength){
										my $seqMut = $vMut / $vLength * 100;
										$mut += $seqMut; #v-gene mutation rate
										$popmut += $seqMut; #overall v-gene mutation rate
										push (@tmpMut, $seqMut);
										$v_gene = $tmpline[3];
									}
								}
								#print $fh_mutation "$numCC\t".scalar(@{$results->{$cluster}})."\t".$v_gene."\t".( $mut / scalar(@{$results->{$cluster}}) ). "\t".stddev(\@tmpMut)."\n";
								$total_results->{'lineage_size_count'}{ scalar @{ $results->{$cluster} } }++;
							}else{#singletons
								$total_results->{'numSC'}++;
								my $tmpnode = $results->{$cluster}[0];
								my $node = $data->{$vGene}{$dGene}{$jGene}{$length}{$jgroup}[$tmpnode];
							
								push @{$total_results->{'singletons'}}, $node;
							
								my @tmpline = split '\t', $infileTie[$node];
								(my $seqtmp = $tmpline[11]) =~ s/\.//g;
								my $seqMut = ( $tmpline[38] =~ /(\d+)(\()?/ ) ? $1 : undef;
								$popmut += $seqMut;
								#print $fh_sinfo "$numSC\t$infileTie[$node]\n";
							
								my @seqID = $tmpline[1] =~ /(.*?):.*:(\d+:\d+)$/;
								$seqID[1] =~ s/:/\-/ if $seqID[1];
							
								#print $fh_sfasta ">0|".join('|',@seqID)."\n";
								#print $fh_sfasta $seqtmp."\n";
							
								$total_results->{'seqLength'}{'singletons'}{length($seqtmp)}++;
							}
						}
					}
				}
			}
		}
	}
	#print $fh_v "$vGene\t$Vcount\n";
}
$total_results->{'header'} = $header;

#sub_sample($total_results, $sampling);

output($outfile, $total_results, \@infileTie, $path, \%settings, $db_name);



# close script nicely

	print POSIX::strftime("%H:%M:%S", localtime)."> script finished\n";
	my $END_TIME = POSIX::strftime("%H:%M:%S", localtime);
	print POSIX::strftime("%H:%M:%S", localtime)."> script end time: $END_TIME\n";



}

####################################
####################################
#           Subroutines            #
####################################
####################################

sub junctionAnalysis {
	my $seq = shift;
	my $settings = shift;
	my $verbose =0;
	
	my %results;
	
	my $vd5 = 3;  #v-d junction length to consider at 5' end
	my $vd3 = 3;  #v-d junction length to consider at 3' end
	my $dj5 = 3;  #d-j junction length to consider at 5' end
	my $dj3 = 3;  #d-j junction length to consider at 3' end
	

		my @row = split "\t", $seq;
		my $vd5_chunk = $row[11];
		my $dj3_chunk = $row[19];
		my $cdr3 = $row[33];
		
		my $junction1_index;
		my $junction2_index;
		
		my $junction1_seq = substr $row[11], ($vd5* -1);
		my $junction2_seq = substr $row[19], 0, $dj3;
		
		print "vd5:\t$vd5_chunk\ncdr3:\t$cdr3\ndj3:\t$dj3_chunk\n" if $verbose;
		print "J1-left: $junction1_seq\n" if $verbose;
		
		if ($cdr3 =~ /$vd5_chunk/gi){
			$junction1_index = pos($cdr3);
		}else{
			for (0..length($vd5_chunk)){
				$vd5_chunk = substr $vd5_chunk, 1;
				if ($cdr3 =~ /^$vd5_chunk/gi){
					$junction1_index = pos($cdr3);
					last;
				}
			}
		}

		$junction1_seq .= substr $cdr3, $junction1_index, $vd3;
		print "J1 full: $junction1_seq\n" if $verbose;
		print "J2-right: $junction2_seq\n" if $verbose;
		
		if ($cdr3 =~ /.*$dj3_chunk/gi){
			$junction2_index = pos($cdr3) - length($dj3_chunk);
		}else{
			for (0..length($dj3_chunk)){
				print "$dj3_chunk\n" if $verbose;
				$dj3_chunk = substr $dj3_chunk, 0, -1;
				print "$dj3_chunk\n" if $verbose;
				if ($cdr3 =~ /$dj3_chunk$/gi){
					$junction2_index = pos($cdr3) - length($dj3_chunk);
					last;
				}
			}
		}
		
		$junction2_seq = substr($cdr3, ($junction2_index - $vd3), $vd3) . $junction2_seq;
		
		print "J2 full: $junction2_seq\n" if $verbose;
		print "Junction 1 match end at $junction1_index\n" if $verbose;
		print "Junction 2 match end at $junction2_index\n" if $verbose;
		
		#my $hold = <STDIN>;
		#print "J2: $junction2_seq\n\n";
		
	#	push @{$results{"$junction1_seq\_$junction2_seq"}}, join "\t", @row;
		
		#print "$row[1]\n";
		#print "V:\t$vd5_chunk\nCDR3:\t$cdr3\nJ:\t$dj3_chunk\n\n";

	
	
	return "$junction1_seq\_$junction2_seq";
}

sub hamming {# subroutine returns number of characters that are different in two strings of equal length
	my $a = shift;
	my $b = shift;
	croak "ERROR: strings must be same length" if length($a) != length($b);
	my $n = length($a);
	my $h = 0;
	for my $k (0 .. $n - 1) { $h++ if substr($a, $k, 1) ne substr($b, $k, 1); }
	return $h;
}

sub germlineSeq {


}

sub germlineAnalysis{
	my $table = shift;
	my $settings = shift;
	
	
	
		
	my $results;
	my %preclusters;
    my $index=0;
    foreach (@{$table}){
    	my @row = split "\t", $_;
        $preclusters{$row[17]}{'size'}++;
        push @{$preclusters{$row[17]}{'indices'}}, $index;
        $index++;
    }
	
    my @sorted_seq = sort { $preclusters{$b}{'size'} <=> $preclusters{$a}{'size'} } keys %preclusters;
    
    for my $i (0..$#sorted_seq ){
        map { push @{$results->{$i}}, $_ } @{$preclusters{$sorted_seq[$i]}{'indices'}};
        for (my $j = $#sorted_seq; $j >= $i+1; $j--){
            my $score = hamming($sorted_seq[$i], $sorted_seq[$j]);
            my $hom = 100 * (length($sorted_seq[$i]) - $score) / length($sorted_seq[$j]);
            if ($hom >= $settings->{'cdr3_homology'}){
                map { push @{$results->{$i}}, $_ } @{$preclusters{$sorted_seq[$j]}{'indices'}};
                splice @sorted_seq, $j, 1;
            }
        }
        last if $i >= $#sorted_seq;
    }

    my $sum;
    foreach my $cluster (keys %{$results}){
        $sum += scalar(@{$results->{$cluster}});
    }
    
}

sub cluster{
    my $table = shift;
    my $settings = shift;

    my $results;
    my %preclusters;
    my $index=0;
    foreach (@{$table}){
    	my @row = split "\t", $_;
        $preclusters{$row[17]}{'size'}++;
        push @{$preclusters{$row[17]}{'indices'}}, $index;
        $index++;
    }

    my @sorted_seq = sort { $preclusters{$b}{'size'} <=> $preclusters{$a}{'size'} } keys %preclusters;
    
    for my $i (0..$#sorted_seq ){
        map { push @{$results->{$i}}, $_ } @{$preclusters{$sorted_seq[$i]}{'indices'}};
        for (my $j = $#sorted_seq; $j >= $i+1; $j--){
            my $score = hamming($sorted_seq[$i], $sorted_seq[$j]);
            my $hom = 100 * (length($sorted_seq[$i]) - $score) / length($sorted_seq[$j]);
            if ($hom >= $settings->{'cdr3_homology'}){
                map { push @{$results->{$i}}, $_ } @{$preclusters{$sorted_seq[$j]}{'indices'}};
                splice @sorted_seq, $j, 1;
            }
        }
        last if $i >= $#sorted_seq;
    }

    my $sum;
    foreach my $cluster (keys %{$results}){
        $sum += scalar(@{$results->{$cluster}});
    }
    croak "error, not all sequences assigned to a cluster $!\n" if $sum != scalar(@{$table});
    return $results;
}

sub output{
	my $outfile = $_[0];
	my $results = $_[1];
	my $data = $_[2];
	my $path = $_[3];
	my $settings = $_[4];
	my $db_name = $_[5];
	my @mapping = (2,undef,undef,undef,1,0,4,6,5,7,9,8,10,11,12,13,14,15,16,17,18,19,20,21,162,23,24,25,26,27,28,29,30,31,32,33,34,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161);

	print Dumper $settings;
	my $fileoutput = 1 if $settings->{'fileoutput'} eq 'y';

	make_path(File::Spec->catdir($path,'lineageOutput',$outfile,'GROUPED_LINEAGES'));
	make_path(File::Spec->catdir($path,'db'));
	my $dbh;

	if ($db_name){
		$dbh = DBI->connect("DBI:SQLite:dbname=$path/db/$db_name.db", "", "", {RaiseError =>1, AutoCommit => 0}) or die $DBI::errstr;
		print "DB connected, generating output...\n";
		$dbh->do(qq{PRAGMA synchronous=OFF;});
		$dbh->do(qq{PRAGMA temp_store=MEMORY;});
		$dbh->do(qq{PRAGMA cache_size=800000;});
		$dbh->do(qq{PRAGMA journal_mode = TRUNCATE;});
		db::createSequenceTable($dbh);
	}
	
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
	
	my $sql = $statements{'insertSequence'} if $db_name;
	my $sth = $dbh->prepare($sql) or return $DBI::errstr if $db_name;

	my %pop_conn;
	my %uniq_iso;

	open my $COMP_COUNTS, ">". "$path/lineageOutput/$outfile/GROUPED_LINEAGE_COUNTS.txt" or die "Error opening file GROUPED_LINEAGE_COUNTS.txt $!\n";
	print $COMP_COUNTS "Grouped Lineage\tIGHV\tTotal Number of Populations\tTotal Number of Sequences";
	foreach my $pop (sort @{$results->{'pops'}}){
		$pop =~ s|.*\/(\S+)$|$1|;
		print $COMP_COUNTS "\tNumber of Sequences $pop";
	}
	print $COMP_COUNTS "\n";
	
	my %counts;
	my $inserts = 0;
	foreach my $lineage (sort {$a <=> $b} keys %{$results->{'clusters'}} ){
		my $v_gene;		
		open my $COMP_INFO,  ">". "$path/lineageOutput/$outfile/GROUPED_LINEAGES/GROUPED_LINEAGES_$lineage\_INFO.txt"  or die "Error opening info file $!\n" if $fileoutput;
		open my $COMP_FASTA, ">". "$path/lineageOutput/$outfile/GROUPED_LINEAGES/GROUPED_LINEAGES_$lineage\_FASTA.txt" or die "Error opening info file $!\n" if $fileoutput;
		print $COMP_INFO "LineageID\tPopulation\t$results->{'header'}" if $fileoutput;
		
		foreach my $index ( @{$results->{'clusters'}{$lineage}} ){
			print $COMP_INFO "$lineage\t$data->[$index]\n" if $fileoutput;

			my @tmp_row = split '\t', $data->[$index];

			my $seq = $tmp_row[9] =~ s/\.//gr;
			print $COMP_FASTA ">$tmp_row[0]|$lineage|$tmp_row[1]\n" if $fileoutput;
			print $COMP_FASTA $seq."\n" if $fileoutput;
			
			my $iso = 'U';
			if ($tmp_row[1] =~ /^(.+?)\_/){
				$iso = $1;
			}
			$counts{$lineage}{'counts'}{$tmp_row[0]}++;
			$counts{$lineage}{'total'}++;
			$counts{$lineage}{'isotypes'}{$iso}++;
			$uniq_iso{$iso} = 1;
			$v_gene = $tmp_row[3];
			
			push @{$pop_conn{'lin_by_pop'}{$tmp_row[0]}{$lineage}}, $index;
			
			my @data = ($lineage, split("\t", $data->[$index]));
			chomp @data;
			my @insert;
			for (0..$#mapping){
				if ($_ == 3){
					$insert[$_] = "$iso";
					next;
				}elsif (! defined($mapping[$_])){
					$insert[$_] = 'null';
					next;
				}elsif ($data[$mapping[$_]] =~ /^(\d+)\s*$/) {
					$insert[$_] = $1;
				}else{
					$insert[$_] = "$data[$mapping[$_]]";
				}
			}
			$sth->execute(@insert) if $db_name;
			$inserts++ if $db_name;
			if ($db_name && ($inserts % 10000 == 0)){ $dbh->do('commit'); $dbh->do('begin');}
		}
		close $COMP_INFO if $fileoutput;
		
		map {$pop_conn{'pop_count'}{$_} += $counts{$lineage}{'counts'}{$_} } keys %{$counts{$lineage}{'counts'}};
		$pop_conn{'lineages'}{join '.', sort keys %{$counts{$lineage}{'counts'}}}++;
		$pop_conn{'seq'}{join '.', sort keys %{$counts{$lineage}{'counts'}}} += $counts{$lineage}{'total'};



		print $COMP_COUNTS "$lineage\t$v_gene\t".scalar(keys %{$counts{$lineage}{'counts'}})."\t$counts{$lineage}{'total'}";
		foreach my $pop(sort @{$results->{'pops'}}){
			$pop =~ s|.*\/(\S+)$|$1|;
			if ($counts{$lineage}{'counts'}{$pop}){
				print $COMP_COUNTS "\t$counts{$lineage}{'counts'}{$pop}";
			}else{
				print $COMP_COUNTS "\t0";
			}
		}
		print $COMP_COUNTS "\n";
				
	}

	open my $COMP_ISOTYPE, ">". "$path/lineageOutput/$outfile/GROUPED_LINEAGE_ISOTYPE_DISTRIBUTIONS.txt" or die "Error opening file GROUPED_LINEAGE_ISOTYPE_DISTRIBUTIONS.txt $!\n";
	print $COMP_ISOTYPE "Lineage ID\t";
	
	foreach my $isotype (sort keys %uniq_iso){
		$isotype = "Unidentified" if uc $isotype eq 'U';
		print $COMP_ISOTYPE "$isotype Count\tPercent $isotype\t";
	}
	print $COMP_ISOTYPE "Total Number\n";

	foreach my $lineage (sort {$a <=> $b} keys %{$results->{'clusters'}} ){
		print $COMP_ISOTYPE "$lineage\t";
		foreach my $isotype (sort keys %uniq_iso){
			$counts{$lineage}{'isotypes'}{$isotype} = 0 if ! $counts{$lineage}{'isotypes'}{$isotype};
			print $COMP_ISOTYPE "$counts{$lineage}{'isotypes'}{$isotype}\t".
			100 * ($counts{$lineage}{'isotypes'}{$isotype} / $counts{$lineage}{'total'})."\t";
		}
		print $COMP_ISOTYPE "$counts{$lineage}{'total'}\n";
	}
	
	close $COMP_COUNTS;
	close $COMP_ISOTYPE;

if($singletons){
		open my $COMP_SINGLETON_INFO,  ">". "$path/lineageOutput/$outfile/SINGLETONS_INFO.txt"  or die "Error opening singleton file $!\n" if $fileoutput;
		open my $COMP_SINGLETON_FASTA, ">". "$path/lineageOutput/$outfile/SINGLETONS_FASTA.txt" or die "Error opening singleton file $!\n" if $fileoutput;
	
		print $COMP_SINGLETON_INFO "SingletonID\tPopulation\t$results->{'header'}" if $fileoutput;

		my $s_id = 0;

		foreach my $singleton (@{$results->{'singletons'}}){
			$s_id++;
			my @tmp_row = split '\t', $data->[$singleton];
			my $seq = $tmp_row[9] =~ s/\.//gr;
			print $COMP_SINGLETON_FASTA ">$tmp_row[0]|$s_id.1|$tmp_row[1]\n" if $fileoutput;
			print $COMP_SINGLETON_FASTA $seq."\n" if $fileoutput;
			print $COMP_SINGLETON_INFO "$s_id.1\t$data->[$singleton]\n" if $fileoutput;
		
			push @{$pop_conn{'lin_by_pop'}{$tmp_row[0]}{"$s_id.1"}}, $singleton;
			
			my $iso = 'U';
			if ($tmp_row[1] =~ /^(.+?)\_/){
				$iso = $1;
			}
			
			my @data = ("$s_id.1", split("\t", $data->[$singleton]));
			chomp @data;
			my @insert;
			for (0..$#mapping){
				if ($_ == 3){
					$insert[$_] = "$iso";
					next;
				}elsif (! defined($mapping[$_])){
					$insert[$_] = 'null';
					next;
				}elsif ($data[$mapping[$_]] =~ /^(\d+)\s*$/) {
					$insert[$_] = $1;
				}else{
					$insert[$_] = "$data[$mapping[$_]]";
				}
			}
			$sth->execute(@insert) if $db_name;
			$inserts++ if $db_name;
			if ($db_name && ($inserts % 10000 == 0)){ $dbh->do('commit'); $dbh->do('begin');}

			$pop_conn{'pop_count'}{$tmp_row[0]}++;
			$pop_conn{'lineages'}{$tmp_row[0]."singletons"}++;
			$pop_conn{'seq'}{$tmp_row[0].".singletons"}++;
		}
		close $COMP_SINGLETON_INFO if $fileoutput;
		close $COMP_SINGLETON_FASTA if $fileoutput;
}
	$dbh->commit if $db_name;
	
	open my $COMP_POP_CONNEC, ">". "$path/lineageOutput/$outfile/POPULATION_CONNECTIONS.txt" or die "Error opening POPULATION_CONNECTIONS.txt file $!\n";
	print $COMP_POP_CONNEC "Population Combination\tNumber of Lineages\tNumber of Sequences\n";
	foreach my $combo (sort {length($a) <=> length($b)}keys %{$pop_conn{'lineages'}}){
		print $COMP_POP_CONNEC "$combo\t$pop_conn{'lineages'}{$combo}\t$pop_conn{'seq'}{$combo}\n";
	}
	close $COMP_POP_CONNEC;
	
	open my $COMP_seq_matches, ">". "$path/lineageOutput/$outfile/SEQUENCE_MATCHES.txt" or die "Error opening sequence_matches file $!\n";
	

	my ($cutoff_type, $cutoff_value) = split '\|', $settings->{"circos_cutoffs"};


	my $cutoffs = expanded_clones::expanded($cutoff_type, \%{$pop_conn{'lin_by_pop'}}, \%{$pop_conn{'pop_count'}}, $cutoff_value);
	make_path("$path/lineageOutput/$outfile/circos/");
	
	
	my $lin_cutoffs;
	#$lin_cutoffs = expanded_clones::expanded('d', \%{$pop_conn{'lin_by_pop'}}, \%{$pop_conn{'pop_count'}}, 60);
	
	circos_config::print_config("$path","$path/lineageOutput/$outfile/circos/", $pop_conn{'pop_count'}, $cutoffs, $lin_cutoffs, \%uniq_iso);
	
	open my $CIRCOS, ">". "$path/lineageOutput/$outfile/circos/circos_lineage_populations.txt" or die "Error opening circos_lineage_populations file $!\n";
	my $chr = 1;
	foreach my $pop(sort keys %{$pop_conn{'pop_count'}}){
		print $CIRCOS "chr - $pop $pop 0 ".$pop_conn{'pop_count'}{$pop}." chr$chr\n";
		$chr++;
	}
	close $CIRCOS;
	db::createIndicies($dbh) if $db_name;
	$dbh->commit if $db_name;

	print_lineage_bands(\%pop_conn, $outfile, $path);
	print_lineage_links(\%pop_conn, $results, $data, $outfile, $path);
	print_connections_data(\%pop_conn, $outfile, $path);
	print_variability_data(\%pop_conn, $results, $data, $outfile, 'd', $path);
	print_variability_data(\%pop_conn, $results, $data, $outfile, 'j', $path);
	print_isotype_data(\%pop_conn, $results, $data, $outfile, $cutoffs, $path, \%uniq_iso);
	print_mutation_data(\%pop_conn, $results, $data, $outfile, $path);
	
}

sub print_mutation_data{
	my $pop_conn = shift;
	my $results = shift;
	my $data = shift;
	my $outfile = shift;
	my $path = shift;

	my %avg;


	open my $CIRCOS,  	 ">". "$path/lineageOutput/$outfile/circos/circos_V_mutations_data.txt" or die "Error opening file circos_V_mutations_data.txt $!\n";
	open my $CIRCOS_avg, ">". "$path/lineageOutput/$outfile/circos/circos_V_avg_mutation_data.txt" or die "Error opening file circos_V_avg_mutation_data.txt $!\n";
	foreach my $lineage(keys %{$pop_conn->{'circos_boundries'}}){
		foreach my $pop (keys %{$pop_conn->{'circos_boundries'}{$lineage}}){
			if ($lineage !~ /^s/){
				my $start = $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'start'};
				foreach my $seq (@{$results->{'clusters'}{$lineage}}){
					my @row = split '\t', $data->[$seq];
					if ($row[0] eq $pop){
						my $rate = $row[36] == 0 ? 0 : $row[38] / $row[36] * 100;
						print $CIRCOS "$pop $start $start $rate #$lineage\n";
						$start++;
						$avg{$pop} += $rate;
					}
				}

			}else{
				my $start = $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'start'};
				my ($index) = $lineage =~ /^s(.*)/;
				my @row = split '\t', $data->[$results->{'singletons'}[$index-1]];
				my $rate = $row[36] == 0 ? 0 : $row[38] / $row[36] * 100;
				print $CIRCOS "$pop $start $start $rate #$lineage\n";
				$avg{$pop} += $rate;
			}
		}
	}

	foreach my $pop(keys %avg){
		for (0..$pop_conn->{'pop_count'}{$pop}){
			print $CIRCOS_avg "$pop $_ $_ ". ($avg{$pop} / $pop_conn->{'pop_count'}{$pop}) ."\n";
		}
	}	
	close $CIRCOS;
	close $CIRCOS_avg;
	return 1;
}

sub print_isotype_data{
	my $pop_conn = shift;
	my $results = shift;
	my $data = shift;
	my $outfile = shift;
	my $cutoffs = shift;
	my $path = shift;
	my $uniqIsotypes = shift;
	my $verbose = 0;
	
	my %isotypes;

	open my $CIRCOS, ">". "$path/lineageOutput/$outfile/circos/circos_lineage_isotype_data.txt" or die "Error opening file $path/lineageOutput/$outfile/circos/circos_lineage_isotype_data.txt $!\n";
	foreach my $lineage(keys %{$pop_conn->{'circos_boundries'}}){

		foreach my $pop (keys %{$pop_conn->{'circos_boundries'}{$lineage}}){
			if ($pop_conn->{'circos_boundries'}{$lineage}{$pop}{'stop'} > $cutoffs->{$pop}){
				my %vars;
				map{$vars{$_}=0}keys %{$uniqIsotypes};
				my $sum = 0;

				foreach my $seq (@{$results->{'clusters'}{$lineage}}){
					my @row = split '\t', $data->[$seq];
					if ($row[0] eq $pop){
						my ($isotype) = ($row[1] =~ /^(.+?)\_/) ? $1 : 'U';
						$vars{$isotype}++;
						$sum++;
					}
				}
				
				print $CIRCOS "$pop $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'start'} $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'stop'} ";
				my @percents;
				map {push @percents, ( $sum ? ($vars{$_}/$sum * 100) : 0);} sort keys %{$uniqIsotypes};

				print $CIRCOS join ',', @percents;
				print $CIRCOS " #$lineage\n";

			}else{
				my $size = $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'stop'} - $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'start'};
				
				if ($size > 1){
					foreach my $seq (@{$results->{'clusters'}{$lineage}}){
						my @row = split '\t', $data->[$seq];
						
						if ($row[0] eq $pop){
							my ($isotype) = ($row[1] =~ /^(.+?)\_/) ? $1 : 'U';
							$isotypes{$pop}{$size}{$isotype}++;
						}
					}
				}else{
					my ($index) = $lineage =~ /^s(.*)/;
					my @row = split '\t', $data->[$results->{'singletons'}[$index-1]];
					my ($isotype) = ($row[1] =~ /^(.+?)\_/) ? $1 : 'U';
					$isotypes{$pop}{$size}{$isotype}++;
				}
 			}
		}
		print "\n" if $verbose;
		my $hold = <STDIN> if $verbose;
	}

	foreach my $pop (keys %{$pop_conn->{'group_boundries'}}){
		foreach my $size (keys %{$isotypes{$pop}}){
			print $CIRCOS "$pop $pop_conn->{'group_boundries'}{$pop}{$size}{'start'} $pop_conn->{'group_boundries'}{$pop}{$size}{'stop'} ";
			my @percents;
			foreach (sort keys %{$uniqIsotypes}){
				my $total;
				map {$total += $_} values $isotypes{$pop}{$size};
				my $value = $isotypes{$pop}{$size}{$_} ? ($isotypes{$pop}{$size}{$_} / $total * 100) : 0;
				push @percents, $value;
			}
			print $CIRCOS join ',', @percents;
			print $CIRCOS " #size:$size\n";
		}
	}

	close $CIRCOS;

	return 1;
}

sub print_variability_data{
	my $pop_conn = shift;
	my $results = shift;
	my $data = shift;
	my $outfile = shift;
	my $v_d_j = shift;
	my $path = shift;
	
	open my $CIRCOS, ">". "$path/lineageOutput/$outfile/circos/circos_$v_d_j\_variability_data.txt" or die "Error opening file circos_$v_d_j\_variablility_data.txt $!\n";
	foreach my $lineage(keys %{$pop_conn->{'circos_boundries'}}){
		foreach my $pop (keys %{$pop_conn->{'circos_boundries'}{$lineage}}){
			my %vars;
			if ($lineage !~ /^s/){
				foreach my $seq (@{$results->{'clusters'}{$lineage}}){
					my @row = split '\t', $data->[$seq];
					$vars{$row[4]}=1 if $v_d_j eq 'v';
					$vars{$row[6]}=1 if $v_d_j eq 'd';
					$vars{$row[5]}=1 if $v_d_j eq 'j';
				}
				print $CIRCOS "$pop $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'start'} $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'stop'} ". scalar(keys %vars) ." #$lineage\n";
			}else{
				print $CIRCOS "$pop $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'start'} $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'stop'} 1 #$lineage\n";
			}
		}
	}
	close $CIRCOS;
	return 1;
}

sub print_connections_data{
	my $pop_conn = shift;
	my $outfile = shift;
	my $path = shift;
	
	open my $CIRCOS, ">". "$path/lineageOutput/$outfile/circos/circos_connections_data.txt" or die "Error opening file circos_connections_data.txt $!\n";
	foreach my $lineage(keys %{$pop_conn->{'circos_boundries'}}){
		my $size = scalar(keys %{$pop_conn->{'circos_boundries'}{$lineage}})-1;
		foreach my $pop (keys %{$pop_conn->{'circos_boundries'}{$lineage}}){
			print $CIRCOS "$pop $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'start'} $pop_conn->{'circos_boundries'}{$lineage}{$pop}{'stop'} $size #$lineage\n";
		}
	}
	close $CIRCOS;
	return 1;
}

sub print_lineage_links{
	my $pop_conn = shift;
	my $results = shift;
	my $data = shift;
	my $outfile = shift;
	my $path = shift;
	
	open my $CIRCOS, ">". "$path/lineageOutput/$outfile/circos/circos_lineage_links.txt" or die "Error opening file circos_lineage_links.txt $!\n";
	
	foreach my $lineage(keys %{$results->{'clusters'}}){
		my @pops = keys %{$pop_conn->{'circos_boundries'}{$lineage}};
		my $color = int(rand(256)).",".int(rand(256)).",".int(rand(256)).",0.5";
		for my $j (0..$#pops){
			for my $i ($j+1..$#pops){
				print $CIRCOS 	"$pops[$j] ".
								"$pop_conn->{'circos_boundries'}{$lineage}{$pops[$j]}{'start'} ".
								"$pop_conn->{'circos_boundries'}{$lineage}{$pops[$j]}{'stop'} ".
								"$pops[$i] ".
								"$pop_conn->{'circos_boundries'}{$lineage}{$pops[$i]}{'start'} ".
								"$pop_conn->{'circos_boundries'}{$lineage}{$pops[$i]}{'stop'} ".
								"color=($color) #$lineage\n";
			}
		}
	}
	
	close $CIRCOS;
	return 1;
}

sub print_lineage_bands{
	my $pop_conn = shift;
	my $outfile = shift;
	my $path = shift;

	open my $CIRCOS, ">". "$path/lineageOutput/$outfile/circos/circos_lineage_bands.txt" or die "Error opening file circos_lineage_bands.txt $!\n";
	
	foreach my $pop (keys %{$pop_conn->{'lin_by_pop'}}){
		my $count = 0;
		foreach my $lineage (sort {scalar @{$pop_conn->{'lin_by_pop'}{$pop}{$a}} <=> scalar @{$pop_conn->{'lin_by_pop'}{$pop}{$b}}} keys %{$pop_conn->{'lin_by_pop'}{$pop}}){
			$pop_conn->{'circos_boundries'}{$lineage}{$pop}{'start'} = $count;
			$pop_conn->{'group_boundries'}{$pop}{scalar(@{$pop_conn->{'lin_by_pop'}{$pop}{$lineage}})}{'start'} = $count if (!$pop_conn->{'group_boundries'}{$pop}{scalar(@{$pop_conn->{'lin_by_pop'}{$pop}{$lineage}})}{'start'} || $count < $pop_conn->{'group_boundries'}{$pop}{scalar(@{$pop_conn->{'lin_by_pop'}{$pop}{$lineage}})}{'start'});
			
			print $CIRCOS "$pop $count ". ($count+scalar(@{$pop_conn->{'lin_by_pop'}{$pop}{$lineage}})) ." #$lineage\n";
			$count += scalar(@{$pop_conn->{'lin_by_pop'}{$pop}{$lineage}});
			$pop_conn->{'circos_boundries'}{$lineage}{$pop}{'stop'} = $count;
			$pop_conn->{'group_boundries'}{$pop}{scalar(@{$pop_conn->{'lin_by_pop'}{$pop}{$lineage}})}{'stop'} = $count;

		}
	}

	close $CIRCOS;
	return 1;
}
1;
__END__

=head1 NAME

 lineage - Computes clonal lineages between sequences and groups them
			 together

=head1 SYNOPSIS

 perl lineage.pl <input file> <output file>

=head1 DESCRIPTION

 This script takes input from the merge.pl script and computes clonal
 lineages from the sequences presented. Clonal lineages are determined based on
 V, D and J similarity, allowing for one variation in the three, as well as
 CDR3 length and homology. CDR3 homology cut off levels can be adjusted in the
 lineage.inf settings file which should be contained in the same directory as
 lineage.pl.

 It contains the following settings:

    lineage_rep_size_cutoff:(0-X)
    cdr3_homology:(0-100)
    v_Variations:(y|n)
    d_Variations:(y|n)
    j_Variations:(y|n)
    useVDJ:(y|n)
    runIsotypes:(y|n)
    runWithIdenticals:(y|n)
    runWithoutIdenticals:(y|n)
    runV_genes:(y|n)

 Cdr3_homology can be a number ranging from 0 - 100 representing the percent
    homology of the cdr3 segments between sequences.

 Lineage_rep_size_cutoff can be a whole number ranging from 0 to the size of
    the largest lineage created by the program.  It controls what lineages are
    output to the '_LINEAGE_REPS.txt' file.
 
 vdjVariations can be either y or n indicating whether to allow variations in
    matching of v,d and j gene id's.
 
 useVDJ can be used to specify whether to consider v d and j gene assignments
    while clustering sequences.
 
 runIsotypes can be used to specify whether to run lineage on the individual
    isotypes as created by merge
 
 runWithIdenticals can be used to specify whether to run lineage on the
    identical output as created by merge
 
 runWithoutIdenticals can be used to specify whether to run lineage on the
    non-identical output as created by merge
 
 runV_genes can be used to specify whether to run lineage on the individual
    v_gene output as created by merge

=head1 AUTHOR

 Christopher Fucile & Alex Rosenberg
 
 Laboratory of Ignacio Sanz
 Division of Allergy, Immunology and Rheumatology
 University of Rochester Medical Center
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
