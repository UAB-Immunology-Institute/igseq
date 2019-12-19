#!/usr/bin/env perl
use strict;
use warnings;

use Cwd;
use File::Path qw(make_path);
use POSIX;
use Carp;
use Parallel::ForkManager;
use List::MoreUtils qw(indexes);
use Data::Dumper;

#my $size = $ARGV[3];
my @folders = @ARGV;
my $seed = 2013;
my @sizes = (100000);
my $productive = 1;

foreach my $folder (@folders){
	opendir my $DIR, $folder or die "Error opening input folder $!\n";
	my @files = sort grep {/^(2|4|8)\_/} readdir $DIR;
	print Dumper @files;

	open my $INF_2, '<', "$folder/".$files[0] or die usage("Error opening file $files[0] $!\n") if -f "$folder/".$files[0];
	open my $INF_4, '<', "$folder/".$files[1] or die usage("Error opening file $files[1] $!\n") if -f "$folder/".$files[1];
	open my $INF_8, '<', "$folder/".$files[2] or die usage("Error opening file $files[2] $!\n") if -f "$folder/".$files[2];

	my $header_2 = <$INF_2>;
	my $header_4 = <$INF_4>;
	my $header_8 = <$INF_8>;

	my @data_2 = <$INF_2>;
	my @data_4 = <$INF_4>;
	my @data_8 = <$INF_8>;
my $prodCount = grep {/\tproductive\t/} @data_2;

	my $file_count;
	foreach my $size(@sizes){
		my $outfolder = $folder."_$size";
		my %used;
		srand($seed);
		make_path($folder."_$size");
		

		my ($output_2) = $files[0] =~ /(.*?)\.txt$/;
		my ($output_4) = $files[1] =~ /(.*?)\.txt$/;
		my ($output_8) = $files[2] =~ /(.*?)\.txt$/;
		
		$file_count++;
		
		$output_2 .= "\_$size\_subsample\.txt";
		$output_4 .= "\_$size\_subsample\.txt";
		$output_8 .= "\_$size\_subsample\.txt";

		open my $OPF_2, '>', "$outfolder/".$output_2 or die "Error opening file \n";
		open my $OPF_4, '>', "$outfolder/".$output_4 or die "Error opening file \n";
		open my $OPF_8, '>', "$outfolder/".$output_8 or die "Error opening file \n";

		print $OPF_2 $header_2;
		print $OPF_4 $header_4;
		print $OPF_8 $header_8;
		
		my @indices;
		if ($prodCount <= $size){
			@indices = indexes {/\tproductive\t/} @data_2;
			print $OPF_2 @data_2[@indices];
			print $OPF_4 @data_4[@indices];
			print $OPF_8 @data_8[@indices];
			last;
						
		}
		
		my $cur_size=0;
		while ($cur_size != $size){
			my $index = int(rand($#data_2));
			if (! defined($used{$index}) && ($data_2[$index] =~ /\tproductive\t/)){
				print $OPF_2 $data_2[$index];
				print $OPF_4 $data_4[$index];
				print $OPF_8 $data_8[$index];
				
				$used{$index}=1;
				$cur_size++;

			}
			last if ( ($cur_size == $#data_2) || ($cur_size == $prodCount) ); #subset size is larger than data set size
		}
		close $OPF_2;
		close $OPF_4;
		close $OPF_8;
	}
	close $INF_2;
	close $INF_4;
	close $INF_8;
	print "Done\n";
}

sub usage{
	my $message = shift;
	print "$message";
	print "perl vquest_subsample.pl <V-quest 2 file> <V-quest 4 file> <V-quest 8 file> <Sample Size>\n";
}