package subsetMerge;
###################################################################################################
#
# subsetMerge.pm
#
# merges v-quest files that were subsetted for 150k sequence batch processing.
#
#########
# 
# Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License
# 
# By exercising the Licensed Rights (defined below), You accept and agree to be 
# bound by the terms and conditions of this Creative Commons 
# Attribution-NonCommercial-ShareAlike 4.0 International Public License ("Public 
# License"). To the extent this Public License may be interpreted as a contract, 
# You are granted the Licensed Rights in consideration of Your acceptance of these 
# terms and conditions, and the Licensor grants You such rights in consideration of 
# benefits the Licensor receives from making the Licensed Material available under these 
# terms and conditions.
# 
# See full licence and definitions in file IGSEQ_LICENSE.PDF
# 
# This software was used to process raw sequence data for the manuscript:
# "Diversity, cellular origin and autoreactivity of antibody-secreting
#     cell expansions in acute systemic lupus erythematosus"
# Tipton CM, Fucile CF, Darce J, Chida A, Ichikawa T, Gregoretti I,
#     Schieferl S, Hom J, Jenks S, Feldman RJ, Mehr R, Wei C, Lee FE,
#     Cheung WC, Rosenberg AF, Sanz I.
# Nature Immunology 2015
# 
# 
# Required attribution:
# 
# By Christopher Fucile (1)
# also
# Christopher Tipton (2)
# Inaki Sanz (2)
# Alex Rosenberg (1)
# 
# (1) Division of Allergy, Immunology and Rheumatology, University of
#     Rochester Medical Center, Rochester NY
# (2) Division of Rheumatology, Department of Medicine, Emory Univeristy,
#     Atlanta GA
# 
# Copyright (c) 2015
# University of Rochester
# Box URV
# 601 Elmwood Avenue
# Rochester, NY 14642
# 
#########
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
use Parallel::ForkManager;

sub mergeSubset{
    my $vquest_files = shift;
    my $path = shift;
    
    open my $SETTINGS, "<". "$path/config/subsetMerge.inf" or die "Error opening file $!\n";

    my %sample_sets;

    while (<$SETTINGS>){
        chomp;
        next if /##subsets##|run:y/;
        return 0 if /run:n/;
        my ($group) = split '\s';
        push @{$sample_sets{$group}}, "2_IMGT-gapped-nt-sequences.txt";
        push @{$sample_sets{$group}}, "4_IMGT-gapped-AA-sequences.txt";
        push @{$sample_sets{$group}}, "8_V-REGION-nt-mutation-statistics.txt";
    }

    remove_tree('$path/subsetMerged') if -e '$path/subsetMerged';
    my $pm = new Parallel::ForkManager(scalar(keys %sample_sets));
    foreach my $sample (keys %sample_sets){
        my $pid = $pm -> start and next;
        
        print "Starting $sample\n";
        make_path("$path/subsetMerged/$sample") if (! -e "$path/subsetMerged/$sample");
        my @v_quest_sets = grep {/$sample(\_f\d+$|$)/i} keys %{$vquest_files};
        for my $j (0..2){
            my $first = 1;
            my $count =0;
            print "\tWorking on file $sample_sets{$sample}->[$j]\n";
            foreach my $subset (@v_quest_sets){
                print "\t\tSubset $subset\n";
                open my $INF, "<". "$path/data/$subset/$vquest_files->{$subset}[$j]" if -e "$path/data/$subset/$vquest_files->{$subset}[$j]" or die "Error opening input file $!\n";
                
                open my $OPF, ">>". "$path/subsetMerged/$sample/$sample_sets{$sample}->[$j]" or die "Error creating file $!\n";
                
                my $header = <$INF>;
                print $OPF $header if $first;
                $first = 0;
                
                while(<$INF>){
                    my $line = $_;
                    my @split_line = split '\t', $line;
                    my ($line_id) = split ':_', $split_line[1];
                    $count++;
                    my $info = join("\t",@split_line[2..$#split_line]);

                    print $OPF "$count\t$line_id\t".$info;
                }
                close $INF;
                close $OPF;
            }
        }
        $pm->finish(0);
    }
    $pm->wait_all_children;

    return \%sample_sets;
}

1;