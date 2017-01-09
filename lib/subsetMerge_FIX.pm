package subsetMerge;
###################################################################################################
#
# subsetMerge.pm
#
# merges v-quest files that were subsetted for 150k sequence batch processing.
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
#    use subsetMerge;
#
###################################################################################################
use strict;
use warnings;

use File::Path qw(make_path remove_tree);
use Data::Dumper;
#use String::LCSS_XS qw(lcss);
use Parallel::ForkManager;

sub mergeSubset{
    my $vquest_files = shift;
    my $path = shift;
    
    open my $SETTINGS, "<". "$path/config/subsetMerge.inf" or die "Error opening file $!\n";

    my %sample_sets;
print "Merging subsets...\n";
	my $group="";
    while (<$SETTINGS>){
        chomp;

        next if /##subsets##|run:y/;
        if (/run:(\w)/){
        	return 0 if $1 ne 'y';
        }
		
		if (/^#(.*)/){
			$group = $1;
			next;
		}
		push @{$sample_sets{$group}}, $_;

    }
	close $SETTINGS;

    remove_tree('$path/subsetMerged') if -e '$path/subsetMerged';
    my $pm = new Parallel::ForkManager(scalar(keys %sample_sets));
my %total_results;
my @filePrefix = (	"2_IMGT-gapped-nt-sequences_",
					"4_IMGT-gapped-AA-sequences_",
					"8_V-REGION-nt-mutation-statistics_");

########## Callback ##############

$pm->run_on_finish(
		sub{
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
			
			if (defined($data_structure_reference)){

				@total_results{keys %{$data_structure_reference}} = values %{$data_structure_reference};
			}
		}
	);
########## Callback ##############


    foreach my $sample (keys %sample_sets){
       my $pid = $pm -> start and next;
        my %results;
        print "Starting $sample\n";
        make_path("$path/subsetMerged/$sample") if (! -e "$path/subsetMerged/$sample");
        my @v_quest_sets = @{$sample_sets{$sample}};

        for my $j (0..2){
            my $first = 1;
            my $count =0;
            print "\tWorking on file $filePrefix[$j]$sample.txt\n";
            push @{$results{$sample}}, "$filePrefix[$j]$sample.txt";
            foreach my $subset (@v_quest_sets){
                print "\t\tSubset $subset\n";
                open my $INF, "<". "$path/data/$subset/$vquest_files->{$subset}[$j]" if -e "$path/data/$subset/$vquest_files->{$subset}[$j]" or die "Error opening input file $!\n";
                open my $OPF, ">>". "$path/subsetMerged/$sample/$filePrefix[$j]$sample.txt" or die "Error creating file $!\n";
                
                my $header = <$INF>;
                print $OPF $header if $first;
                $first = 0;
                
                while(<$INF>){
                    #need to adjust the counts for adding additional files
                    my $line = $_;
                    my @split_line = split '\t', $line;
                    $count++;
                    print $OPF "$count\t".join("\t",@split_line[1..$#split_line]);
                }
                close $INF;
                close $OPF;
            }
        }
        $pm->finish(0, \%results);
    }
    $pm->wait_all_children;

    return \%total_results;
}

1;