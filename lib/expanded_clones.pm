package expanded_clones;

use warnings;
use strict;
use Data::Dumper;
use File::Find;


sub expanded{
	my $type = shift;
	
	return lin_measurement(shift, shift, shift) if lc($type) eq "lin";
	return seq_measurement(shift, shift, shift) if lc($type) eq "seq";
	return delta_measurement(shift, shift, shift) if lc($type) eq "delta";
	return d_measurement(shift, shift, shift) if lc($type) eq "d";
	
}

sub lin_meausrement{
	my $data = shift;
	my $total = shift;
	my $cutoff = shift;
	
	my %cut_offs;
	
	foreach my $pop (keys %{$data}){
		my $previous = 0;
		$cut_offs{$pop} = 0;
		foreach my $group ( sort { scalar(@{$data->{$pop}{$a}}) <=> scalar(@{$data->{$pop}{$b}}) } keys %{$data->{$pop}} ){
			$previous = (scalar(@{$data->{$pop}{$group}}) / $total->{$pop} * 100) if !$previous;
			last if( (scalar(@{$data->{$pop}{$group}}) / $total->{$pop} * 100) - $previous >= $cutoff);
			$cut_offs{$pop} += scalar(@{$data->{$pop}{$group}});
			$previous = (scalar(@{$data->{$pop}{$group}}) / $total->{$pop} * 100);
		} 
	}
	return \%cut_offs;
}

sub seq_measurement{
	my $data = shift;
	my $total = shift;
	my $cutoff = shift;


	my %cut_offs;
	
	foreach my $pop (keys %{$total}){
		$cut_offs{$pop} = $total->{$pop} - $cutoff;
	}
	return \%cut_offs;
}

sub delta_measurement{
	my $data = shift;
	my $total = shift;
	my $cutoff = shift;

	my %cut_offs;

	foreach my $pop (keys %{$data}){
		my $previous = 0;
		$cut_offs{$pop} = 0;
		foreach my $group ( sort { scalar(@{$data->{$pop}{$a}}) <=> scalar(@{$data->{$pop}{$b}}) } keys %{$data->{$pop}} ){
			$previous = (scalar(@{$data->{$pop}{$group}}) / $total->{$pop} * 100) if !$previous;
			last if( (scalar(@{$data->{$pop}{$group}}) / $total->{$pop} * 100) - $previous >= $cutoff);
			$cut_offs{$pop} += scalar(@{$data->{$pop}{$group}});
			$previous = (scalar(@{$data->{$pop}{$group}}) / $total->{$pop} * 100);
		} 
	}
	return \%cut_offs;
}

sub d_measurement{
	my $data = shift;
	my $total = shift;
	my $cutoff = (shift) / 100;

	my %cut_offs;
	
	foreach my $pop (keys %{$total}){
		my $size;
		my $previous = 0;
		my $seq_cut = int($cutoff * $total->{$pop});
		foreach my $group (sort { scalar(@{$data->{$pop}{$a}}) <=> scalar(@{$data->{$pop}{$b}}) } keys %{$data->{$pop}} ){
			$size += scalar(@{$data->{$pop}{$group}});
			$previous = $size if ! $previous;
			if ($size > $seq_cut){ last }
			if ($size == $seq_cut){ $previous = $size; last;}
			$previous = $size;
		}
		$cut_offs{$pop} = $previous;
	}
	return \%cut_offs;
}

1;
__END__
