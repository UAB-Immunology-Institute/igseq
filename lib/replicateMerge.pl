#!/usr/bin/perl

###################################################################################################
#
# replicateMerge.pl
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
#    perl replicateMerge.pl <input directory>
#
###################################################################################################

our $VERSION = '2013.03.29.1';

use strict;
use warnings;

use Math::Round qw(round);
use Data::Dumper;
use Carp;
use Tie::File;
use File::Path;
use Cwd;

use Sys::CPU qw(cpu_count);
use Parallel::ForkManager;
use Text::LevenshteinXS qw(distance);

srand(2013);
use FindBin;
my $working_dir = $FindBin::RealBin;
$working_dir =~ s/\/lib//;
my $output_dir = "$working_dir/replicateMergeOutput";

File::Path::make_path("$working_dir/replicateMergeOutput");

my $dir = $ARGV[0];

opendir my $DIR, $dir or die "Can not open directory $dir $!\n";
my @files = grep{/\_join\.fastq/} readdir $DIR;
closedir $DIR;


open STDOUT, "| tee $output_dir/output.txt" or croak "Tee stdout error";

my %settings;
my %replicates;
my $replicate_flag = 0;
my $threshold_flag = 0;
my $threshold_setting_count = 0;

#read settings

print "> Reading settings\n";
open my $INF, "<"."$working_dir/config/replicateConfig.inf" or die "Error opening $working_dir/config/replicateConfig.inf $!\n";
while (<$INF>){
    chomp;
    if( $_ =~/^\##/ ){
        $replicate_flag = 1;
    }elsif( $_ =~/^\#/ || $_ =~/^$/){
        next;
    }elsif($replicate_flag == 1){
        my ($replicate_id, $replicates) = split ':', $_;
        push(@{$replicates{$replicate_id}}, split '\s', $replicates);
    }elsif( $_ =~/\%/ ){
        if ($threshold_flag == 0){
            $threshold_flag = 1;
            $threshold_setting_count++;
        }elsif($threshold_flag == 1){
            $threshold_flag = 0;
        }
    }elsif($threshold_flag){
        my @row = split ':';
        $settings{"$row[0].$threshold_setting_count"} = $row[1];
    }else{
        my @row = split ':';
        $settings{$row[0]} = $row[1];
    }
}

print "> Reading input\n";
#for each set of replicates specified in settings
if (scalar(keys %replicates) > 1){
    my $pm = new Parallel::ForkManager(cpu_count()-1);
    my $childCount = 0;
    foreach my $set (keys %replicates){
        $childCount++;
        print "Starting checks on replicate set $childCount...\n";
        my $pid = $pm->start and next;
        my @file_set;
        foreach my $replicate ( @{$replicates{$set}} ){
            my ($file) = grep {/$replicate\_/i} @files;
            $file = "$dir".$file;
            push @file_set, $file;
        }

        my ($results, $print, $stats) = readData(\@file_set, \%settings);
        my ($merge_print, $merge_stats) = mergeReplicates(\%settings, $results, scalar(@file_set));
        @{$stats}{keys %{$merge_stats}} = values %{$merge_stats};
        print $print.$merge_print."\n";
        printStats($set, $stats);
        printOutput($results, "$set.fastq");
        
        $pm->finish(0);
    }
    $pm->wait_all_children;
}else{
    map { $_ = $dir.$_ } @files;
    my $pm = new Parallel::ForkManager(cpu_count()-1);
#    my $pm = new Parallel::ForkManager(1);
    foreach my $file (@files){
        print "Starting checks on $file\n";
        my $pid= $pm->start and next;
        my @file_set;
        
        push @file_set, $file;
        my ($results, $print, $stats) = readData(\@file_set, \%settings);
        my ($merge_print, $merge_stats) = mergeReplicates(\%settings, $results, scalar(@file_set));
        @{$stats}{keys %{$merge_stats}} = values %{$merge_stats};
        print $print.$merge_print."\n";
        
        $file =~ s|.*\/||g;
        $file =~ s|\.fastq||g;
        printStats($file, $stats);
        printOutput($results, "$file.fastq");
        $pm->finish(0);
    }
    $pm->wait_all_children;
}



############
############
sub printStats{
    my ($set, $stats) = @_;
    open my $OPF, ">" . "$output_dir/$set"."_stats.txt" or die "error opening file $!\n";
    print $OPF "\t$set\n";
    print $OPF "Initial\t$stats->{'total_initial'}\n".
          "After Length Threshold\t$stats->{'after_length'}\n".
          "After Merging\t$stats->{'after_merge'}\n".
          "Unique left After Quality Threshold\t$stats->{'uniq_left'}\n".
          "Total After Quality Threshold\t$stats->{'total_left'}\n";
    close $OPF;
}



############
############
sub highQualSelection{
    my ($length, $seq, $amount, $results ) = @_;

    my %set;
    my @ids;
    foreach my $file( keys %{$results->{$length}{$seq}}){
        for (0..$#{$results->{$length}{$seq}{$file}{'id'}}){
            $set{$results->{$length}{$seq}{$file}{'id'}[$_]}{'qual'} = $results->{$length}{$seq}{$file}{'qual'}[$_];
            $set{$results->{$length}{$seq}{$file}{'id'}[$_]}{'qual_value'} = round( sequenceQuality($results->{$length}{$seq}{$file}{'qual'}[$_]));
            $set{$results->{$length}{$seq}{$file}{'id'}[$_]}{'file'} = $file;
        }
    }
    
    @ids = sort {$set{$a}{'qual_value'} <=> $set{$b}{'qual_value'} } keys %set;

    die "Error, something odd is going on $!\n" if $amount > scalar(@ids);
    
    @ids = @ids[0..($#ids-$amount)];
SEQ_TO_REMOVE:
    foreach my $id (@ids){
        for (0..$#{$results->{$length}{$seq}{ $set{$id}{'file'} }{'id'}}){
            if($results->{$length}{$seq}{ $set{$id}{'file'} }{'id'}[$_] eq $id){
                splice( @{$results->{$length}{$seq}{ $set{$id}{'file'} }{'id'}}, $_, 1);
                splice( @{$results->{$length}{$seq}{ $set{$id}{'file'} }{'qual'}}, $_, 1);
                if (scalar(@{$results->{$length}{$seq}{ $set{$id}{'file'} }{'qual'}}) == 0){
                    delete $results->{$length}{$seq}{ $set{$id}{'file'} };
                }
                next SEQ_TO_REMOVE;
            }
        }
    }
}

############
############
sub lowQualSelection{
    my ($length, $seq, $amount, $results ) = @_;
    
    my %set;
    my @ids;
    foreach my $file( keys %{$results->{$length}{$seq}}){
        for (0..$#{$results->{$length}{$seq}{$file}{'id'}}){
            $set{$results->{$length}{$seq}{$file}{'id'}[$_]}{'qual'} = $results->{$length}{$seq}{$file}{'qual'}[$_];
            $set{$results->{$length}{$seq}{$file}{'id'}[$_]}{'qual_value'} = round( sequenceQuality($results->{$length}{$seq}{$file}{'qual'}[$_]));
            $set{$results->{$length}{$seq}{$file}{'id'}[$_]}{'file'} = $file;
        }
    }
    
    @ids = sort {$set{$b}{'qual_value'} <=> $set{$a}{'qual_value'} } keys %set;
    
    die "Error, something odd is going on $!\n" if $amount > scalar(@ids);

    @ids = @ids[0..($#ids-$amount)];
SEQ_TO_REMOVE:
    foreach my $id (@ids){
        for (0..$#{$results->{$length}{$seq}{ $set{$id}{'file'} }{'id'}}){
            if($results->{$length}{$seq}{ $set{$id}{'file'} }{'id'}[$_] eq $id){
                splice( @{$results->{$length}{$seq}{ $set{$id}{'file'} }{'id'}}, $_, 1);
                splice( @{$results->{$length}{$seq}{ $set{$id}{'file'} }{'qual'}}, $_, 1);
                if (scalar(@{$results->{$length}{$seq}{ $set{$id}{'file'} }{'qual'}}) == 0){
                    delete $results->{$length}{$seq}{ $set{$id}{'file'} };
                }
                next SEQ_TO_REMOVE;
            }
        }
    }
}


############
############
sub randomSelection{
    my ($length, $seq, $amount, $set, $results ) = @_;
RANDSEQ:
    while((scalar(@$set)-$amount) > 0){
        my $index = int rand scalar(@$set);
        foreach my $file(keys %{$results->{$length}{$seq}}){
            for (0..$#{$results->{$length}{$seq}{$file}{'id'}}){
                if($results->{$length}{$seq}{$file}{'id'}->[$_] eq $set->[$index]){
                    splice( @{$results->{$length}{$seq}{$file}{'id'}}, $_, 1);
                    splice( @{$results->{$length}{$seq}{$file}{'qual'}}, $_, 1);
                    splice @$set, $index, 1;
                    if (scalar(@{$results->{$length}{$seq}{ $file }{'qual'}}) == 0){
                        delete $results->{$length}{$seq}{ $file };
                    }
                    next RANDSEQ;
                }
            }
        }
    }
}



############
############
sub readData{
    my ($file_set, $settings) = @_;

    my %results;
    my $total;
    my $total_count_length_f;
    my $file_counter = 0;
    my $print;
    my %stats;
    
    foreach my $file(@{$file_set}){
        $file_counter++;
        open my $INF, "<". $file or die "Can not open input file $file:$!\n";
        my @tmp_fastq = <$INF>;
#		tie my @tmp_fastq, 'Tie::File', $file or die "Error tie-ing file $file $!\n";  added 2014.10.30, removed? 2016.02.17
        close $INF;
        chomp(@tmp_fastq);
        carp "File does not appear to be fastq format $!\n" if $tmp_fastq[0] !~ /^\@/;
        
        my $count_initial  = 0;
        my $count_length_f = 0;
        my $count_length_p = 0;
        
        
        for (my $j = 0; $j <= $#tmp_fastq; $j += 4){
            
            $count_initial++;
            
            if (length($tmp_fastq[$j+1]) <= $settings->{'lengthThreshold'}){
                $count_length_f++;
                next;
            }else{
                if(defined($results{length($tmp_fastq[$j+1])}{$tmp_fastq[$j+1]})){
                    $count_length_p++;
                    push @{$results{length($tmp_fastq[$j+1])}{$tmp_fastq[$j+1]}{$file}{'id'}}, $tmp_fastq[$j];
                    push @{$results{length($tmp_fastq[$j+1])}{$tmp_fastq[$j+1]}{$file}{'qual'}}, $tmp_fastq[$j+3];
                }else{
                    $count_length_p++;
                    push @{$results{length($tmp_fastq[$j+1])}{$tmp_fastq[$j+1]}{$file}{'id'}}, $tmp_fastq[$j];
                    push @{$results{length($tmp_fastq[$j+1])}{$tmp_fastq[$j+1]}{$file}{'qual'}}, $tmp_fastq[$j+3];
                }
            }
        }
        $print .= "> File: $file\n\t".
        "Initial Number of sequences: $count_initial\n\t".
        "Number of sequences not passing length threshold $settings->{'lengthThreshold'}: $count_length_f\n";
        $total += $count_initial;
        $total_count_length_f += $count_length_f;
    }
    $print .= "Total number of initial squences: $total\n".
            "Total number of sequences after length threshold: ".($total - $total_count_length_f)."\n";
    $stats{'total_initial'} = $total;
    $stats{'after_length'} =($total - $total_count_length_f);
    return \%results, $print, \%stats;
}


############
############
sub mergeReplicates{
    my ($settings, $results, $numOfReplicates) = @_;
    my $uniq_count_in_all   = 0;
    my $total_count_in_all  = 0;
    my $sequences_in_one_file = 0;
    my $total_count_after_selection = 0;
    my $threshold_count_before = 0;
    my $threshold_unique_sequences_removed = 0;
    my $threshold_count_removed = 0;
    my $print;
    my %stats;

#TODO: move this loop into a sub function and have calls from within the if statements
LENGTH:
    foreach my $length(keys %{$results}){
        my @seqs = keys %{$results->{$length}};
    SEQS:
        for my $j(0..$#seqs){

            if (scalar(keys %{$results->{$length}{$seqs[$j]}}) < $settings{'numReplicatesSeqPresentIn'}){
                foreach my $file(keys %{$results->{$length}{$seqs[$j]}}){
                    $sequences_in_one_file += scalar(@{$results->{$length}{$seqs[$j]}{$file}{'id'}});
                }
                delete $results->{$length}{$seqs[$j]};
            }else{
            FILE:
                foreach my $file(keys %{$results->{$length}{$seqs[$j]}}){
                    my @quality = @{$results->{$length}{$seqs[$j]}{$file}{'qual'}};
                    my $index = 0;
                    for (0..$#quality){
                        $threshold_count_before++;
                        if (!threshold($quality[$_], $settings)){
                            $threshold_count_removed++;
                            splice( @{$results->{$length}{$seqs[$j]}{$file}{'id'}},   $index, 1);
                            splice( @{$results->{$length}{$seqs[$j]}{$file}{'qual'}}, $index, 1);
                            if (scalar(@{$results->{$length}{$seqs[$j]}{$file}{'qual'}}) == 0){
                                delete $results->{$length}{$seqs[$j]}{$file};
                                if (scalar(keys %{$results->{$length}{$seqs[$j]}}) == 0){
                                    delete $results->{$length}{$seqs[$j]};
                                    $threshold_unique_sequences_removed++;
                                    next SEQS;
                                    if(scalar(keys %{$results->{$length}}) == 0){
                                        delete $results->{$length};
                                        next LENGTH;
                                    }
                                }
                                next FILE;
                            }
                        }else{
                            $index++;
                        }
                    }
                }
                if ($settings->{'mergeReplicatesCopyNumber'} eq 'smallest'){
                    my $smallest = 0;
                    my @set;
                    foreach my $file(keys %{$results->{$length}{$seqs[$j]}}){
                        my $test_size = scalar @{$results->{$length}{$seqs[$j]}{$file}{'id'}};
                        $total_count_in_all += $test_size;
                        $smallest = ($smallest) ? (($test_size < $smallest) ? $test_size : $smallest) : $test_size;
                        push @set, @{$results->{$length}{$seqs[$j]}{$file}{'id'}};
                    }
                    
                    randomSelection  ($length, $seqs[$j], $smallest, \@set, $results) if $settings->{'mergeReplicatesCopyQuality'} eq 'random';
                    highQualSelection($length, $seqs[$j], $smallest, $results)        if $settings->{'mergeReplicatesCopyQuality'} eq 'high_quality';
                    lowQualSelection ($length, $seqs[$j], $smallest, $results)        if $settings->{'mergeReplicatesCopyQuality'} eq 'low_quality';
                    
                    map {$total_count_after_selection += scalar(@{$results->{$length}{$seqs[$j]}{$_}{'id'}})} keys %{$results->{$length}{$seqs[$j]}};
                    
                }elsif($settings->{'mergeReplicatesCopyNumber'} eq 'largest'){
                    my $largest = 0;
                    my @set;
                    foreach my $file(keys %{$results->{$length}{$seqs[$j]}}){
                        my $test_size = scalar @{$results->{$length}{$seqs[$j]}{$file}{'id'}};
                        $total_count_in_all += $test_size;
                        $largest = ($largest) ? (($test_size > $largest) ? $test_size : $largest) : $test_size;
                        push @set, @{$results->{$length}{$seqs[$j]}{$file}{'id'}};
                    }


                    randomSelection  ($length, $seqs[$j], $largest, \@set, $results) if $settings->{'mergeReplicatesCopyQuality'} eq 'random';
                    highQualSelection($length, $seqs[$j], $largest, $results)        if $settings->{'mergeReplicatesCopyQuality'} eq 'high_quality';
                    lowQualSelection ($length, $seqs[$j], $largest, $results)        if $settings->{'mergeReplicatesCopyQuality'} eq 'low_quality';

                    map {$total_count_after_selection += scalar(@{$results->{$length}{$seqs[$j]}{$_}{'id'}})} keys %{$results->{$length}{$seqs[$j]}};

                }elsif($settings->{'mergeReplicatesCopyNumber'} eq 'average'){
                    my $total;
                    my @set;

                    foreach my $file(keys %{$results->{$length}{$seqs[$j]}}) {
                        $total += scalar @{$results->{$length}{$seqs[$j]}{$file}{'id'}};
                        push @set, @{$results->{$length}{$seqs[$j]}{$file}{'id'}};
                    }

                    $total_count_in_all += $total;
                    my $average = round($total / $numOfReplicates);

                    randomSelection  ($length, $seqs[$j], $average, \@set, $results) if $settings->{'mergeReplicatesCopyQuality'} eq 'random';
                    highQualSelection($length, $seqs[$j], $average, $results)        if $settings->{'mergeReplicatesCopyQuality'} eq 'high_quality';
                    lowQualSelection ($length, $seqs[$j], $average, $results)        if $settings->{'mergeReplicatesCopyQuality'} eq 'low_quality';
                    

                    map {$total_count_after_selection += scalar(@{$results->{$length}{$seqs[$j]}{$_}{'id'}})} keys %{$results->{$length}{$seqs[$j]}};
                }elsif($settings->{'mergeReplicatesCopyNumber'} eq 'sum'){
                    my $total;
                    my @set;
                    
                    foreach my $file(keys %{$results->{$length}{$seqs[$j]}}) {
                        $total += scalar @{$results->{$length}{$seqs[$j]}{$file}{'id'}};
                        push @set, @{$results->{$length}{$seqs[$j]}{$file}{'id'}};
                    }
                    
                    $total_count_in_all += $total;

                    randomSelection  ($length, $seqs[$j], $total, \@set, $results) if $settings->{'mergeReplicatesCopyQuality'} eq 'random';
                    highQualSelection($length, $seqs[$j], $total, $results)        if $settings->{'mergeReplicatesCopyQuality'} eq 'high_quality';
                    lowQualSelection ($length, $seqs[$j], $total, $results)        if $settings->{'mergeReplicatesCopyQuality'} eq 'low_quality';
                    

                    map {$total_count_after_selection += scalar(@{$results->{$length}{$seqs[$j]}{$_}{'id'}})} keys %{$results->{$length}{$seqs[$j]}};
                }
            }
        }
        $uniq_count_in_all += scalar(keys %{$results->{$length}});
        if(!keys %{$results->{$length}}){
            delete $results->{$length};
        }
    }
    $print = "Number of sequences not in $settings{'numReplicatesSeqPresentIn'} of the replicates: $sequences_in_one_file\n".
        "Number of total sequences after merging replicates: $threshold_count_before\n".
        "Number of total sequences removed through thresholding: $threshold_count_removed\n".
        "Number of total sequences left after thresholding: ".($threshold_count_before - $threshold_count_removed)."\n".
        "Number of unique sequences removed by thresholding: $threshold_unique_sequences_removed\n".
        "Number of unique sequences represented across all replicates after thresholding: $uniq_count_in_all\n";

    $stats{'after_merge'}   = $threshold_count_before;
    $stats{'total_left'}    = $total_count_after_selection;
    $stats{'uniq_left'}     = $uniq_count_in_all;
    return $print, \%stats;

}



sub printOutput{
    my $results = shift;
    my $name = shift;
    open my $OUTPUT, ">". "$output_dir/$name" or die "Error $!\n";
    foreach my $length (keys %{$results}){
        foreach my $seq (keys %{$results->{$length}}){
            foreach my $file (keys %{$results->{$length}{$seq}}){
                for (0..$#{$results->{$length}{$seq}{$file}{'id'}}){
                    print $OUTPUT "$results->{$length}{$seq}{$file}{'id'}->[$_]\n";
                    print $OUTPUT "$seq\n";
                    print $OUTPUT "+\n";
                    print $OUTPUT "$results->{$length}{$seq}{$file}{'qual'}->[$_]\n";
                }
            }
        }
    }
    close $OUTPUT;
}

################
################
sub mergeInDelMut{
    my $results    = shift;
    my $settings   = shift;
    my $file   = shift;
    my $numOfReplicates = shift;

    my $num_of_gapped_only = 0;
    my $numToCheck = 0;
    my $numToCheck_last6 = 0;
    my $numToCheck_first6 = 0;
    my %count;
    $|=1;
    foreach my $length(sort {$b <=> $a} keys %{$results}){ #sorted decending order, only need to compare to sequences the same length and the next length
        my @basepairs_j = sort keys %{$results->{$length}};
        my @basepairs_k = sort keys %{$results->{$length}};
            
        for my $j(0..$#basepairs_j){
            if(scalar(@basepairs_j)>1){
                for my $k($j+1..$#basepairs_k-1){           #sequence comparison of same length
                    
                    my $dist = Text::LevenshteinXS::distance($basepairs_j[$j], $basepairs_k[$k]);
                    
                    if ($dist <= $settings->{'InDelMutNumberBasepairDifference'}){
                        my $bp_j = join( '', @{[splice( @{[split('', $basepairs_j[$j])]}, -6)]});
                        my $bp_k = join( '', @{[splice( @{[split('', $basepairs_k[$k])]}, -6)]});
                        my $dist2 = Text::LevenshteinXS::distance($bp_j, $bp_k);
                        $bp_j = join( '', @{[splice( @{[split('', $basepairs_j[$j])]}, 0, 6)]});
                        $bp_k = join( '', @{[splice( @{[split('', $basepairs_k[$k])]}, 0, 6)]});
                        my $dist3 = Text::LevenshteinXS::distance($bp_j, $bp_k);
                        if ($dist2 == $dist){
                            $numToCheck_last6++;
                        }elsif($dist3 == $dist){
                            $numToCheck_first6++;
                        }
                        
                        $numToCheck++;
                    }
                }
            }
            for my $l(1..$settings->{'InDelMutNumberBasepairDifference'}){
                

                next if !exists $results->{$length - $l};
                my @basepairs_minusDiff = keys %{$results->{$length - $l}};
                for my $k(0..$#basepairs_minusDiff){   #seq comparison with length - $l
                    my $dist = Text::LevenshteinXS::distance($basepairs_j[$j], $basepairs_minusDiff[$k]);
                    
                    if ($dist <= $settings->{'InDelMutNumberBasepairDifference'}){
                        my $bp_j = join( '', @{[splice( @{[split('', $basepairs_j[$j])]}, -6)]});
                        my $bp_k = join( '', @{[splice( @{[split('', $basepairs_minusDiff[$k])]}, -6)]});
                        my $dist2 = Text::LevenshteinXS::distance($bp_j, $bp_k);
                        $bp_j = join( '', @{[splice( @{[split('', $basepairs_j[$j])]}, 0, 6)]});
                        $bp_k = join( '', @{[splice( @{[split('', $basepairs_minusDiff[$k])]}, 0, 6)]});
                        my $dist3 = Text::LevenshteinXS::distance($bp_j, $bp_k);
                        if ($dist2 == $dist){
                            $numToCheck_last6++;
                        }elsif($dist3 == $dist){
                            $numToCheck_first6++;
                        }

                        
                        $numToCheck++;
                        if ( $l == $dist){ #difference between sequences is only gaps
                            $num_of_gapped_only++;
                        }
                    }
                }
            }
        }
    }
    print "Number to check $numToCheck\n";
    print "Number in last 6: $numToCheck_last6\n";
    print "Number in first 6: $numToCheck_first6\n";
    print "Number of sequence compairisons with only gaps $num_of_gapped_only\n";

}


################
################
sub sequenceQuality{
    my $qual_seq = shift;
    my $quality;
    my @qual_seq = split '', $qual_seq;
    for (0..$#qual_seq){
        $quality += ord($qual_seq[$_])-33;
    }
    return ($quality/scalar(@qual_seq));
}


################
################
sub checkQualityDifference($$){
    my ($q1, $q2, $seq1, $seq2, $threshold) = @_;

    die "error $!\n" if length($q1) != length($q2);
    
    my @q1 = split '', $q1;
    my @q2 = split '', $q2;
    
    my @seq1 = split '', $seq1;
    my @seq2 = split '', $seq2;
    
    for (0..$#q1){
        if( abs((ord($q1[$_])-33) - (ord($q2[$_])-33)) >= $threshold && $seq1[$_] ne $seq2[$_]){
            my $diff = (ord($q1[$_])-33) - (ord($q2[$_])-33);
            return ($seq1[$_],$seq2[$_],$q1[$_],$q2[$_],$diff);
        }
        return 0;
    }
}


################
################
sub threshold($$){
    my ($qual, $settings) = @_;

    
    $qual = [split '', $qual];
    my %count;

    my @qualityThresholds = grep {/qualityThreshold/} sort keys %{$settings};
    my $numberOfSettings = scalar @qualityThresholds;

    foreach my $base (@{$qual}){
        for (0..$numberOfSettings-1){
            my $settingNumber = $_ +1;
            $count{$settingNumber} = 0 if !exists $count{$settingNumber};
            $count{$settingNumber}++   if (ord($base)-33) < $settings->{"qualityThreshold.$settingNumber"};
            if($settings->{"basepairThreshold.$settingNumber"} != 0){
                return 0 if $count{$settingNumber} >= $settings->{"basepairThreshold.$settingNumber"} ;
            }else{
                return 0 if $count{$settingNumber} > $settings->{"basepairThreshold.$settingNumber"};
            }
        }
    }
    return 1;
}


################
################
sub mergeIdenticals($$){
    my ($file_set, $settings) = @_;
    
    my %sequences;
    
    foreach my $file (@{$file_set}){
        open my $INF, "<". $file or die "Can not open input file $file: $!\n";
        while(<$INF>){
            my $line1 = $_;
            my $line2 = <$INF>;
            my $line3 = <$INF>;
            my $line4 = <$INF>;
            
            chomp($line1, $line2, $line3, $line4);
        }
        close $INF;
    }
}


1;


=head1 NAME
 
 replicateMerge.pl
 
=head1 SYNOPSIS
 
 perl replicateMerge.pl  
 
=head1 DESCRIPTION
 
 Performs quality filtering as well as the ability to merge technical replicates
    of samples.
 
=head1 SETTINGS
 
 See the replicateConfig.inf file
 
=head1 OUTPUT
 
 FASTQ files filtered based on settings.
 
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