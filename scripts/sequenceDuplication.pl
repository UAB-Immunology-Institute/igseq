#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use DBI;
use DBD::SQLite;
use FindBin;

my $cwd = $FindBin::Bin;
my $db = $ARGV[0];

my $dbh = DBI->connect("DBI:SQLite:dbname=$cwd/$db", "", "", {RaiseError =>1, AutoCommit => 0}) or die $DBI::errstr;

my $sql = "select population,sequenceID,lineageID,NT_VDJ from sequence";

my $results = $dbh->selectall_arrayref($sql) or return $DBI::errstr;

my %seqs; 
my @overlapLineageID;
my %pops;
foreach my $seq (@{$results}){
	my ($noBarcodeSeq) = $seq->[3] =~ s/(\.+)\w{6}/$1\.\.\.\.\.\./ir;
	my ($sub) = $seq->[0] =~ s/.*?((1095)|(1325)).*$/$1/r;
	$seqs{$seq->[2]}{$noBarcodeSeq}{'subject'}{$sub} = 1;

	push @{ $seqs{$seq->[2]}{$noBarcodeSeq}{'pop'} }, $seq->[0];
	push @{ $seqs{$seq->[2]}{$noBarcodeSeq}{'seqID'} }, $seq->[1];
	push @{ $seqs{$seq->[2]}{$noBarcodeSeq}{'barcodeSeq'} }, $seq->[3];

	$pops{$seq->[0]} = 1;
}
# my @pops = ("S1_35D7_SWM_Ca09_pos",
# 			"S2_35D7_PB",
# 			"S3_35D7_SWM_Ca09_neg",
# 			"S4_35D7_DNM_Ca09_pos",
# 			"S5_35D7_Naive",
# 			"S6_35D7_DNM_Ca09_neg",
# 			"S7_35D14_SWM_Ca09_pos",
# 			"S8_35D14_Plasmablasts",
# 			"S9_35D14_SWM_Ca09_neg",
# 			"S10_35D14_DNM_Ca09_pos",
# 			"S11_35D14_Naive",
# 			"S12_35D14_DNM_Ca09_neg",
# 			"S13_40D7_SWM_Ca09_pos_Tbet_low",
# 			"S14_40D7_SWM_Ca09_neg_Tbet_low",
# 			"S15_40D7_SWM_Ca09_pos_Tbet_high",
# 			"S16_40D7_Plasmablasts",
# 			"S17_40D7_Naive",
# 			"S18_40D14_SWM_Ca09_pos_Tbet_low",
# 			"S19_40D14_SWM_Ca09_neg_Tbet_low",
# 			"S20_40D14_SWM_Ca09_pos_Tbet_high",
# 			"S21_40D14_Plasmablasts",
# 			"S22_40D14_Naive",
# 			"S23_40D7_SWM_Ca09_pos_Tbet_high",
# 			"S24_40D14_SWM_Ca09_neg_Tbet_low");
my @pops = ("1095_PBL_DN",
			"1095_PBL_DP",
			"1095_PBL_NAIVE",
			"1095_PBL_POP2",
			"1095_PBL_POP3",
			"1095_PBL_SM",
			"1095_np_naive",
			"1095_np_pop1",
			"1095_np_pop2",
			"1095_np_pop3",
			"1095_np_pop4",
			"1095_np_pop5",
			"1325_NP_NAIVE",
			"1325_NP_POP2",
			"1325_NP_POP3",
			"1325_NP_POP4",
			"1325_PBL_DN",
			"1325_PBL_DP",
			"1325_PBL_NAIVE",
			"1325_PBL_POP2",
			"1325_PBL_POP3",
			"1325_PBL_SM");
undef $results;
my $count = 0;
my %uniqByPop;
my $verbose = 0;
open my $OPF2, '>', 'sublineages.txt' or die "Error opening file $!\n";
print $OPF2 "LineageID\tSub-lineage\tPop1\tPop2\t# Pop1\t# Pop2\n";
foreach my $lineageID (keys %seqs){
	foreach my $seq (keys %{$seqs{$lineageID}}){
		if (1 < scalar keys %{ $seqs{$lineageID}{$seq}{'subject'} }){
			my %uniquePops;
			my $total = 0;
			$count++;
			map { $uniquePops{$_} = 0 if ! defined $uniquePops{$_};
				  $uniquePops{$_}++;
				  $total++; } @{$seqs{$lineageID}{$seq}{'pop'}} ;
			#print Dumper @{$seqs{$lineageID}{$seq}{'pop'}};
			#print Dumper %uniquePops;
			foreach my $i1 (0..$#pops){
				next if ! defined $uniquePops{$pops[$i1]};
				next if $pops[$i1] !~ /1095/;
				#print "$pop1\t";
				foreach my $i2 ($i1..$#pops){
					next if ! defined $uniquePops{$pops[$i2]};
					next if $pops[$i2] !~ /1325/;
					print $OPF2 "$lineageID\t$count\t$pops[$i1]\t$pops[$i2]\t$uniquePops{$pops[$i1]}\t$uniquePops{$pops[$i2]}\n";
					#print "$pop2\n";
					#print Dumper %uniqByPop;
					
				  	$uniqByPop{$pops[$i1]}{$pops[$i2]} = defined $uniqByPop{$pops[$i1]}{$pops[$i2]} ? $uniqByPop{$pops[$i1]}{$pops[$i2]} + $uniquePops{$pops[$i1]} : $uniquePops{$pops[$i1]}; 
				}
			}
			#print Dumper %uniqByPop;
			#my $hold = <STDIN>;
			
		}
	}
}
close $OPF2;

open my $OPF, '>', 'uniqueByPopulation.txt' or die "Error opening file $!\n";
my $header = join("\t", @pops);
print $OPF "\t$header\n";
my @pops2 = @pops;
foreach my $pop1 (@pops){
	print $OPF "$pop1\t";
	foreach my $pop2 (@pops2){
		my $value = $uniqByPop{$pop1}{$pop2} ? $uniqByPop{$pop1}{$pop2} : 0;
		print $OPF "$value\t";
	}
	print $OPF "\n";
}
close $OPF;
