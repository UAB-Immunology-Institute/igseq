#!/usr/bin/env perl



use strict;
use warnings;

use Cwd;
use File::Path qw(make_path);
use POSIX;
use Carp;
use Parallel::ForkManager;
use Data::Dumper;

my $size = $ARGV[3];
my @files = @ARGV[0..2]];
my $seed = 2013;
#my @sizes = (300,600,1000);


open my $INF_2, "<". $files[0] or die usage("Error opening file $files[0] $!\n") if -f $files[0];
open my $INF_4, "<". $files[1] or die usage("Error opening file $files[1] $!\n") if -f $files[1];
open my $INF_8, "<". $files[2] or die usage("Error opening file $files[2] $!\n") if -f $files[2];

my $header_2 = <$INF_2>;
my $header_4 = <$INF_4>;
my $header_8 = <$INF_8>;

my @data_2 = <$INF_2>;
my @data_4 = <$INF_4>;
my @data_8 = <$INF_8>;


my $file_count;
#foreach my $size(@sizes){
	my %used;
	srand($size);
	my ($output_2) = $files[0] =~ /(.*?)\.txt$/;
	my ($output_4) = $files[1] =~ /(.*?)\.txt$/;
	my ($output_8) = $files[2] =~ /(.*?)\.txt$/;
	
	$file_count++;
	
	$output_2 .= "\_$size\_subsample\.txt";
	$output_4 .= "\_$size\_subsample\.txt";
	$output_8 .= "\_$size\_subsample\.txt";

	open my $OPF_2, ">". $output_2 or die "Error opening file \n";
	open my $OPF_4, ">". $output_4 or die "Error opening file \n";
	open my $OPF_8, ">". $output_8 or die "Error opening file \n";

	print $OPF_2 $header_2;
	print $OPF_4 $header_4;
	print $OPF_8 $header_8;
	my $cur_size=0;
	while ($cur_size != $size){
		my $index = int(rand($#data_2));
		if (! defined($used{$index})){
			print $OPF_2 $data_2[$index];
			print $OPF_4 $data_4[$index];
			print $OPF_8 $data_8[$index];
			
			$used{$index}=1;
			$cur_size++;

		}
	}
	close $OPF_2;
	close $OPF_4;
	close $OPF_8;
#}
close $INF_2;
close $INF_4;
close $INF_8;
print "Done\n";


sub usage{
	my $message = shift;
	print "$message";
	print "perl vquest_subsample.pl <V-quest 2 file> <V-quest 4 file> <V-quest 8 file> <Sample Size>\n";
}