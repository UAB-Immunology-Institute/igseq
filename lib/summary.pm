package summary;
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
use Parallel::ForkManager;

use constant V_GENES  => ['IGHV1', 'IGHV2', 'IGHV3', 'IGHV4', 'IGHV5', 'IGHV6', 'IGHV7'];

sub IgSeq_summary{
    my $datetime  = shift;
    my $filegroup = shift;

    my ($pop, $other) = $filegroup =~ /lineageOutput\/(.*?)\/(.*$)/;
    
    $other =~ s|\/|\.|g;
    
    opendir my $DIR, $filegroup if -e $filegroup or die "Can not open directory $filegroup $!\n";
    my @files = readdir $DIR;
    closedir $DIR;
    
    make_path('/Applications/IgSeq/Summary');
    my $OPF;

    open $OPF, ">". "/Applications/IgSeq/Summary/$pop\_IgSeq.summary.txt" or die
        "Error opening output file $pop\_IgSeq.summary.txt $!\n";
    
    print $OPF "Population\t$pop\n";
    mutations($filegroup, $OPF, \@files);
    D_calculations($filegroup, $OPF, \@files);
    isotype_info($filegroup, $OPF, \@files);
    my $v_434_amount = v_gene_info($filegroup, $OPF, \@files);
    v_434($filegroup, $OPF, \@files, $v_434_amount);
    

    close $OPF;
    return 0;
}

sub mutations{
    my $filegroup = shift;
    my $OPF = shift;
    my $tmp = shift;
    my $v_434_amount = shift;
    
    my @files = grep {/_INFO\./} @{$tmp};
    my $count = 0;
    my $mutations = 0;
    
    foreach my $file (@files){
        open my $INF, "<". "$filegroup/$file" if -f "$filegroup/$file" or die "Error opening file $file $!\n";
        my $header = <$INF>;
        my $index = ($file =~ /SINGLETONS/) ? 50 : 51;
        
        while (<$INF>){
            chomp;
            $count++;
            my @line = split '\t';
            $mutations += $line[$index];
        }
        close $INF;
    }
    my $tmp_output = ( $count != 0) ? ($mutations / $count ) : 0;
    print $OPF "Average_Population_Mutation_Rate\t$tmp_output\n";
    
    @files = grep {/_LINEAGE_MUTATION/} @{$tmp};
    $count = 0;
    $mutations = 0;
    
    foreach my $file (@files){
        open my $INF, "<". "$filegroup/$file" if -f "$filegroup/$file" or die "Error opening file $file $!\n";
        my $header = <$INF>;
        my $index = 3;

        while (<$INF>){
            chomp;
            $count++;
            my @line = split '\t';
            $mutations += $line[$index];
        }
        close $INF;
    }
    $tmp_output = ( $count != 0) ? ($mutations / $count ) : 0;
    print $OPF "Average_Mutation_Rate_of_Lineage_Averages\t$tmp_output\n";
    
    return 0;
}

sub v_434{
    my $filegroup = shift;
    my $OPF = shift;
    my $tmp = shift;
    my $v_434_amount = shift;
    
    my @files = grep {/_INFO\./} @{$tmp};
    my $amount = 0;
    
    foreach my $file (@files){
        open my $INF, "<". "$filegroup/$file" if -f "$filegroup/$file" or die "Error opening file $file $!\n";
        my $header = <$INF>;
        my $index = ($file =~ /SINGLETONS/) ? 37 : 38;

        while (<$INF>){
            next if !/IGHV4-34/;
            chomp;
            
            my @line = split '\t';
            
            if($line[$index] =~ /AVY$/i){
                $amount++;
            }
        }
        close $INF;
    }
    my $avy = ($v_434_amount != 0) ? ($amount / $v_434_amount * 100) : 0;
    print $OPF "AVY_IGHV4-34_percent.of.V4-34\t$avy\n";
    
    
    return 0;
}

sub v_gene_info{
    my $filegroup = shift;
    my $OPF = shift;
    my $tmp = shift;
    
    my ($file) = grep {/UNIQUE_V_COUNTS/} @{$tmp};
    
    my %v_gene_families;
    my $size;
    my $v_434_amount = 0;
    
    open my $INF, "<". "$filegroup/$file" if -f "$filegroup/$file" or die "Error opening file $file $!\n";
    my $header = <$INF>;
    
    while(<$INF>){
        chomp;
        my ($id, $num) = split '\s';
        $v_434_amount = $num if $id =~ /IGHV4-34/i;
        ($id) = $id =~ /(IGHV\d+)\-/;
        
        $v_gene_families{$id} += $num;
        $size += $num;
    }
    

    
    for (0..6){
        my $tmp_output_1 = exists($v_gene_families{V_GENES->[$_]}) ? $v_gene_families{V_GENES->[$_]} : 0;
        print $OPF V_GENES->[$_]."\_amount\t$tmp_output_1\n";
        
        my $tmp_output = ($size > 0 && exists($v_gene_families{V_GENES->[$_]}) ) ? ($v_gene_families{V_GENES->[$_]}/$size *100) : 0;
        print $OPF V_GENES->[$_]."\_percent\t$tmp_output\n";
    }
    
    print $OPF "IGHV4-34_amount\t$v_434_amount\n";
    
    my $tmp_output = ($size != 0) ? ($v_434_amount / $size * 100) : 0;
    print $OPF "IGHV4-34_percent.of.total\t$tmp_output\n";
    
    $tmp_output = ($v_gene_families{'IGHV4'} != 0) ? ($v_434_amount / $v_gene_families{'IGHV4'} * 100) : 0;
    print $OPF "IGHV4-34_percent.of.IGHV4\t$tmp_output\n";
    
    close $INF;
    return (0, $v_434_amount);
}


sub isotype_info{
    my $filegroup = shift;
    my $OPF = shift;
    my $tmp = shift;
    
    my @files = grep {/\_FASTA\./} @{$tmp};
    
    my $size;
    my %isotypes;
    
    foreach my $file (@files){
        open my $INF, "<". "$filegroup/$file" if -f "$filegroup/$file" or die "Error opening file $file $!\n";
        while (<$INF>){
            if (/^\>/){
                chomp;
                $size++;
                my $line = $_;
                
                my ($isotype) = $line =~ /^\>\d+\|(\w)\_/;
                $isotypes{$isotype}++;
                
            }
        }
        close $INF;
    }
    
    foreach my $isotype (sort keys %isotypes){
        print $OPF "$isotype\_Isotype_amount\t$isotypes{$isotype}\n";
        my $tmp_output = ($size != 0) ? ($isotypes{$isotype}/$size * 100) : 0;
        print $OPF "$isotype\_Isotype_percent\t$tmp_output\n";
    }
    
    
    return 0;
}



sub D_calculations{
    my $filegroup = shift;
    my $OPF = shift;
    my $tmp = shift;

    my @files = @{$tmp};
    my $size = 0;
    my %lineageHistogram;
    my ($seq_d80, $seq_d50, $seq_d20);
    my ($lin_d80, $lin_d50, $lin_d20);
    
    my ($file) = grep {/LINEAGE_HISTOGRAM/} @files;
    open my $INF, "<". "$filegroup/$file" if -f "$filegroup/$file" or die "Error opening file $file $!\n";
    my $header = <$INF>;
    while(<$INF>){
        chomp;
        my ($linSize,$linNum) = split '\s';
        $lineageHistogram{$linSize} = $linNum;
        $size += $linSize * $linNum;
    }
    
    $seq_d80 = 0.8 * $size;
    $seq_d50 = 0.5 * $size;
    $seq_d20 = 0.2 * $size;
    
    $size = 0;
    my $count = 0;
    foreach my $linSize (sort {$b<=>$a} keys %lineageHistogram){
        for (1..$lineageHistogram{$linSize}){
            $count++;
            $size += $linSize;
            my $fraction = 0;
            next if $linSize == 0;
            $fraction = ( ($linSize - ($size - $seq_d80)) / $linSize) if ($size > $seq_d80 && !defined($lin_d80));
            $lin_d80 = ($count-1 + $fraction) if ($size > $seq_d80 && !defined($lin_d80));
            
            $fraction = ( ($linSize - ($size - $seq_d50)) / $linSize) if $size > $seq_d50 && !defined($lin_d50);
            $lin_d50 = ($count-1 + $fraction) if $size > $seq_d50 && !defined($lin_d50);
            
            $fraction = ( ($linSize - ($size - $seq_d20)) / $linSize) if $size > $seq_d20 && !defined($lin_d20);
            $lin_d20 = ($count-1 + $fraction) if $size > $seq_d20 && !defined($lin_d20);
        }
    }
    print $OPF "D80\t$lin_d80\n";
    print $OPF "D50\t$lin_d50\n";
    print $OPF "D20\t$lin_d20\n";

    close $INF;
    return 0;
}


sub summary_join{
    
    opendir my $ODH, "/Applications/IgSeq/Summary" or die "Error opening directory /Applications/IgSeq/Summary $!\n";
    my @summary_pops = grep {!/^\./} readdir $ODH;
    closedir $ODH;
    
    my @output = ("Population\t",
        "Average_Population_MutationRate\t",
        "Average_Mutation_Rate_of_Lineage_Averages\t",
        "D80\t",
        "D50\t",
        "D20\t",
        "A_Isotype_amount\t",
        "A_Isotype_percent\t",
        "G_Isotype_amount\t",
        "G_Isotype_percent\t",
        "M_Isotype_amount\t",
        "M_Isotype_percent\t",
        "U_Isotype_amount\t",
        "U_Isotype_percent\t");
    
    for (0..6){
        push @output, (V_GENES->[$_]."\_amount\t",
            V_GENES->[$_]."\_percent\t");
    }
    
    push @output, "IGHV4-34_amount\t",
        "IGHV4-34_percent.of.total\t",
        "IGHV4-34_percent.of.IGHV4\t",
        "AVY_IGHV4-34_percent.of.V4-34\t";
    
    foreach my $file (@summary_pops){
        my $index = 0;
        open my $IFH, "<". "/Applications/IgSeq/Summary/$file" or die "Error opening file $file $!\n";
        while (<$IFH>){
            chomp;
            my ($junk, $value) = split "\t";
            
            $output[$index] .= $value."\t";
            $index++;
        }
        close $IFH;
    }
    
    map {s/\t$/\n/ig} @output;
    
    open my $OFH, ">". "/Applications/IgSeq/Summary/IgSeq.summary.txt" or die "Error opening file IgSeq.summary.txt $!\n";
    print $OFH @output;
    close $OFH;
    return 0;
}

1;