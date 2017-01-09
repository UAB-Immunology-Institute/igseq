package merge;
require Exporter;

@ISA = qw(Exporter);
@EXPORT_OK = qw(merge);

###################################################################################################
#
# merge.pm
#
# merges several primary data files useful for sequence analysis
#
# by Alex Rosenberg
# 
# June 23, 2011
#
########################
# Modified February 17, 2015:
#
# Added ability to include sequence alignments that were missing D alignments from IHHMune.
# Parses out allele identifiers on gene names ( *01_(L1) )
# Added output capability of a file consisting of sequences that are productive but not in both IHHMune and VQuest. 
# Added ability to remove identical and nearly identical sequences.
# Re-wrote outputTable function to be more efficient ( ~99.4% speed increase, not noticed at small file sizes)
# Included individual outputs for isotypes as well as v-gene filter output.
#
# by Christopher Fucile
########################
#
# Division of Allergy, Immunology and Rheumatology
# University of Rochester Medical Center
#
# copyright (c) 2015, All Rights Reserved
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
#    merge(<Vquest 2 file>, <Vquest 4 file>, <Vquest 8 file>)
#
###################################################################################################

use strict;
use POSIX;
use Carp;
use Class::Struct;
use Scalar::Util qw(looks_like_number);
use Tie::File;
use Data::Dumper;
use File::Path qw(make_path);
use lib "$FindBin::RealBin/";

my $path = $FindBin::RealBin;

sub merge{
	my $f2   = shift;
	my $f4   = shift;
	my $f8   = shift;
	my $fout = shift;
	
# expected number of columns in each file type
	my $IHMMUNE_NCOL  = 17;
	my $VQUEST2_NCOL  = 18;
	my $VQUEST4_NCOL  = 18;
	my $VQUEST8_NCOL  = 130;
	
# indices of columns to keep in each file type (note - column numbers, not perl indices)	
	my $VQUEST2_ICOL  = [ 1 .. 17 ];
	my $VQUEST4_ICOL  = [ 1, 6 .. 17 ];
	my $VQUEST8_ICOL  = [ 1, 4 .. 13 ];
	
# time stamp for running script and log file
	my $LOG_FILE = 'merge_log.txt';
	my $RUN_TIME = POSIX::strftime("%Y-%m-%d, %H:%M:%S", localtime);

# get settings from merge.inf
    my %settings;
    my $mergeConfig = File::Spec->catfile(($path, 'config'), 'merge.inf');
    open my $INFSETTINGS, "<" . $mergeConfig or croak "ERROR: can't open merge configuration files $!\n";
	my $flag = 0;
	my @idfilter;
    while (<$INFSETTINGS>){
        chomp;
		if($_ =~/\#/){
			$flag = 1;
		}elsif($flag==1){
			push(@idfilter, $_);
		}else{
			my @row = split ':';
			$settings{$row[0]} = $row[1];
		}
    }
    close $INFSETTINGS;
    croak "ERROR: Problem with settings file (merge.inf). Too many inputs. Check usage. \n" if scalar keys %settings !=7;
	
# make output directories
	my $filePrefix = $fout;
	my $outDir = File::Spec->catfile(($path,'mergeOutput'),$filePrefix);
    eval { make_path($outDir) };
    if ($@) {
        print "Couldn't create directory: $@ $!\n";
    }

# user option for ignoring no d gene alignments
    if($settings{include_NoDgene} eq 'y' ){
            print "> Including alignments missing D gene matches.\n";}else{
            print "> Excluding alignments missing D gene matches.\n";}


# user option for removing identical sequences up to X bp
    if(!(looks_like_number($settings{identical_Cutoff}))){
        print "> Invalid setting, only removing exactly identical sequences.\n";
        $settings{identical_Cutoff}=0;
    }

# append to log file whatever shows up on terminal (errors only show up on STDERR on terminal)
	$LOG_FILE = File::Spec->catfile(($path,'mergeOutput',"$filePrefix"),$LOG_FILE);
	open STDOUT, "| tee -a $LOG_FILE" or croak "Tee stdout error";
	print "*********************************************\n";
	print "> script start: $RUN_TIME\n";

# write file names to log file
	print "> Vquest file 2 = " . $f2 . "\n";
	print "> Vquest file 4 = " . $f4 . "\n";
	print "> Vquest file 8 = " . $f8 . "\n";

# open files
	open my $INF_v2, $f2 or die "Error opening file $f2 $1 $!\n";
	open my $INF_v4, $f4 or die "Error opening file $f4 $1 $!\n";
	open my $INF_v8, $f8 or die "Error opening file $f8 $1 $!\n";
	
# check file lengths
	while(<$INF_v2>){};
	my $v2_rowCount = $.;
	seek $INF_v2, 0, 0;
	
	while(<$INF_v4>){};
	my $v4_rowCount = $.;
	seek $INF_v4, 0, 0;
	
	while(<$INF_v8>){};
	my $v8_rowCount = $.;
	seek $INF_v8, 0, 0;
	
	die "Error, not all V-quest files have same length.\n" if ($v2_rowCount != $v4_rowCount && $v2_rowCount != $v8_rowCount);
	
# read header
	my $v2_header = <$INF_v2>;
	my $v4_header = <$INF_v4>;
	my $v8_header = <$INF_v8>;
	
	chomp($v2_header, $v4_header, $v8_header);

	my @v2_header = split "\t", $v2_header;
	splice @v2_header, 3, 0, @v2_header[3,4,5];
	map {s/ and allele//} @v2_header[3,4,5];
	my @v4_header = split "\t", $v4_header;
	my @v8_header = split "\t", $v8_header;
	
	map {$_ = 'Vquest_2_'.$_;} @v2_header;
	map {$_ = 'Vquest_4_'.$_;} @v4_header;
	map {$_ = 'Vquest_8_'.$_;} @v8_header;

# create output files
	
	my $outFile = File::Spec->catfile($outDir, "$filePrefix\_withIdenticals.txt");
	my @output;
	#tie my @output, 'Tie::File', $outFile or die "Error making output file $outFile $!\n";
	$output[0] =  join("\t", @v2_header[ 1 ..20 ], @v4_header[ 1, 6 .. 17 ], @v8_header[ 1, 4 .. 129 ])."\n";

# process v2
	print "> Processing V-quest 2 file...\n";
	my %idkeepv2;
	my $index = 0;
	my %isotypes;
	my %sequenceCopies;
	while(<$INF_v2>){
		chomp;
		my @row = split "\t", $_;
		my ($seq_nospace) = $row[6] =~ s/\.//gr;
		if ((($row[2] !~ /(No results|unproductive|unknown)/)                      ||
            ($row[2] =~ /No results/   && $settings{include_NoResults}    eq 'y')  ||
            ($row[2] =~ /unproductive/ && $settings{include_Unproductive} eq 'y')  ||
            ($row[2] =~ /unknown/      && $settings{include_Unknown}      eq 'y') )&&
            ($row[4] =~ /[IGH|TRB]J/)												     &&
            (length($seq_nospace) >= $settings{'length_threshold'})					 ){
            	$index++; # start at 1 since header is at 0
                                
                my $v = $row[3] =~ s/.*?\_?([IG|TRB].*?)(\*.*)/$1/r;
                my $j = $row[4] =~ s/.*?\_?([IG|TRB].*?)(\*.*)/$1/r;
                my $d = $row[5] =~ s/.*?\_?([IG|TRB].*?)(\*.*)/$1/r;

                push @{$sequenceCopies{$row[6]}{'indicies'}}, $index;
                $sequenceCopies{$row[6]}{'nogap'} = $row[6] =~ s/\.//gr;
                $sequenceCopies{$row[6]}{'count'}++;
				
				$row[3] =~ s/.*?\_?([IG|TRB].*?)(\s.*)/$1/;
                $row[4] =~ s/.*?\_?([IG|TRB].*?)(\s.*)/$1/;
                $row[5] =~ s/.*?\_?([IG|TRB].*?)(\s.*)/$1/;

				splice @row, 3, 0, ($v, $j, $d);
				
				$idkeepv2{$row[1]} = $index;
                push @output, join "\t", @row[ 1 .. 20];
                my ($isotype) = $row[1] =~ /^(.w?)\_/;
                $isotype = 'U' if ! $isotype;
                $isotypes{$isotype}{$index} = $index;
        }	
	}

my $totalCount = 0;
my %identified;
#find single copy sequences
print "Looking for copies of sequences...\n";
	my @seqs = keys %sequenceCopies;
# OUTER:
# 	for (my $i=0; $i <= $#seqs; $i++){
# 		if ($sequenceCopies{$seqs[$i]}{'count'} > 1){
# 			$totalCount += $sequenceCopies{$seqs[$i]}{'count'};
# 		}
# 		next OUTER if $identified{$i};
# 		my $lenI = length($seqs[$i]);
# 		for (my $j=$i+1; $j <= $#seqs; $j++){
# 			next OUTER if $identified{$j};
# 			my $lenJ = length($seqs[$j]);
# 			if ($lenI == $lenJ){
# 				my $sequencedPartI = $sequenceCopies{$seqs[$i]}{'nogap'};
# 				my $sequencedPartJ = $sequenceCopies{$seqs[$j]}{'nogap'};
# 				my $distance = hamming($seqs[$i], $seqs[$j]);
# 				if ($distance == abs(length($sequencedPartI) - length($sequencedPartJ)) ){
# 					#print "Found $distance\t". abs(length($sequencedPartI) - length($sequencedPartJ))."\n";
# 					$totalCount += $sequenceCopies{$seqs[$i]}{'count'} if $sequenceCopies{$seqs[$i]}{'count'} == 1;
# 					$totalCount += $sequenceCopies{$seqs[$j]}{'count'} if $sequenceCopies{$seqs[$j]}{'count'} == 1;
# 					$identified{$i}=1;
# 					$identified{$j}=1;
# 				}
# 			}
# 		}
# 	}

# print number of seq
	print "> Found $#output sequences meeting initial criteria...\n";
# process v4
	while (<$INF_v4>){
		chomp;
		my @row = split "\t", $_;
		if ($idkeepv2{$row[1]}){
			$output[$idkeepv2{$row[1]}] .= "\t". join "\t", @row[1,6..17];
		}
	}
	#exit;
# process v8
	while (<$INF_v8>){
		chomp;
		my ($line) = $_ =~ s/\(.*?\)//gr;
		my @row = split "\t", $line;
		if ($idkeepv2{$row[1]}){
			$output[$idkeepv2{$row[1]}] .= "\t". join("\t", @row[1,4..129]) . "\n";
		}
	}

# print results with identicals
	print "> Generating result files...\n";
	open my $OPF, ">".$outFile or die "Error opening output file $outFile $!\n";
	print $OPF @output;
	close $OPF;
	
# genearte isotype outputs
	foreach my $isotype (keys %isotypes){
		#next if $isotype eq 'U';
		$outFile = File::Spec->catfile($outDir, "$filePrefix\_$isotype\_withIdenticals.txt");
		open $OPF, ">".$outFile or die "Error opening output file $outFile $!\n";
		print $OPF $output[0];
		foreach my $row (keys %{$isotypes{$isotype}}){
			print $OPF $output[$isotypes{$isotype}{$row}];
		}
		close $OPF;
	}
	
	
# find identical sequences
	if ($settings{'separate_Identicals'} eq 'y'){
		my %identID;
		for(my $i = 1; $i <= $#output; $i++){
			my @row = split "\t", $output[$i];
			my ($seq) = lc $row[7] =~ s/\.//gr;
			$identID{$seq} = $i;
		}
		$outFile = File::Spec->catfile($outDir, "$filePrefix\_withOutIdenticals.txt");
		open $OPF, ">". $outFile or die "Error opening output file $outFile $!\n";
		print $OPF $output[0];
		foreach my $seq (keys %identID){
			print $OPF $output[$identID{$seq}];
		}
		close $OPF;
		foreach my $isotype ( keys %isotypes ){
			next if $isotype eq 'U';
			$outFile = File::Spec->catfile($outDir, "$filePrefix\_$isotype\_withOutIdenticals.txt");
			open $OPF, ">".$outFile or die "Error opening output file $outFile $!\n";
			print $OPF $output[0];
			foreach my $seq (%identID){
				my @row = split "\t", $output[$identID{$seq}];
				my ($seqIsotype) = $row[0] =~ /^(.*?)\_/;
				if ($isotype eq $seqIsotype){
					print $OPF $output[$identID{$seq}];
				}else{
					next;
				}
			}
			close $OPF;
		}
	}
	print "> Done.\n";
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

1;