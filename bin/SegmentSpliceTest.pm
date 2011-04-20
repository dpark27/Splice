#!/usr/bin/perl -w

package SegmentSpliceTest;

use strict;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Statistics::Basic qw(:all nofill);
use Statistics::Distributions;

sub new {
	my $class = shift;
	my %params = @_;

	unless($params{ 'readlength' }) {
		die "Required Parameter readlength not provided for SegmentSpliceTest";
	}

	my $self = { read_length => $params{ 'readlength' } };

	bless $self, $class;
	return $self;
}

sub GetPairedPoints {
	my $self = shift;
	my ($window_start, $forward) = @_;

	my @paired_points = ();	
	my $half_read_length = int(($self -> { read_length }) / 2);
	for(my $i = 0; $i < $half_read_length; $i++) {
		my $first;
		my $second;
		if($forward) {
			$first = $window_start + $i;
			$second = $first + 46; #$half_read_length;
		}
		else {
			$first = $window_start - $half_read_length;
			$second = $first - $half_read_length;
		}
		
		my $pair = "$first,$second";	
		push(@paired_points, $pair);
	}	

	return \@paired_points;
}

sub TestDifferences {
	my $self = shift;
	my ($diffs_ref, $gene_position) = @_;

	my @diffs = @$diffs_ref;
	my $dbar = 0;
	if($#diffs > 0) {
		$dbar = sum(@diffs) / ($#diffs + 1);
	}

	my $sd = stddev(@diffs); 
	if($sd eq 'n/a' || $sd == 0) { 
		return 0;
	} 
	
	my $t = 0;
	if($#diffs && $sd) {
		$t = $dbar / ($sd / sqrt($#diffs + 1)); 
	}

	return $t;
}

sub GetCriticalTValue {
	my $self = shift;
	my ($deg_freedom, $pvalue) = @_;

	#print "P-Value: $pvalue\nDeg: $deg_freedom\n";
	my $t = Statistics::Distributions::tdistr($deg_freedom, $pvalue);
	#print "Critical Value: $t\n"; 

	return $t;
}

sub ScanForSplice {
	my $self = shift;
	my ($indexes_ref, $pvalue, $gene_start, $chr) = @_;

	my $critical_t = $self -> GetCriticalTValue(($self -> { read_length } / 2) - 1, $pvalue);	

	my @indexes = @$indexes_ref;
	my $indexes_length = @indexes;	

	my $count = 0;	
	my $prev_slice = 0;
	my @return_vals = ();
	for(my $i = 0; ($i + $self -> { read_length }) < $indexes_length; $i++) {
		my $random_paired_points_ref = $self -> GetPairedPoints($i, 1);
		my @random_paired_points = @$random_paired_points_ref;
		
		my @paired_value_diffs = ();
		my @firsts = ();
		my @seconds = ();
		foreach my $pair(@random_paired_points) {
			my @data = split(/,/, $pair);
			my $first = $data[0];
			my $second = $data[1];
			if(!defined($indexes[$first]) || !defined($indexes[$second])) {
				last;
			}
			
			my $paired_value_difference = $indexes[$first] - $indexes[$second]; 
			push(@paired_value_diffs, $paired_value_difference);
		}

		my $test_result = $self -> TestDifferences(\@paired_value_diffs, $gene_start - 200 + $i);
		if(abs($test_result) > $critical_t) {
			my $point = $gene_start - 200 + $i;
			$point = "$chr\:" . $point;		
	
			if($test_result > 0) {
				my $marker = 3;
				push(@return_vals, "$point,$marker");
			}
			else {
				my $marker = 5;
				push(@return_vals, "$point,$marker");
			}

			#if($prev_slice != $gene_start - 200 + $i - 1) {
			#	print "\n\n";
			#	print $gene_start - 200 + $i . ",$i ,$test_result\n";$count++;
			#}
			#else {
			#	print $gene_start - 200 + $i . ",$test_result\n";
			#}

			#$prev_slice = $gene_start + $i - 200;
		}
	}

	return \@return_vals;
}

1;
