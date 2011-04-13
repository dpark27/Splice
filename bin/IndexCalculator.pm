#!/usr/bin/perl -w

package IndexCalculator;

use strict;

use FileHandle;
use Getopt::Long;
use Bio::DB::Sam;

require BEDParser;
require RPKMCalculator;

sub new {
	my $class = shift;
	my %params = @_;
	
	for my $required(qw{ bam genome }) {
		if($params{ $required } eq '') {
			die "Required parameter $required not provided";
		}
	}

	my $self = { bam => $params{ 'bam' },
     		     genome => $params{ 'genome' }};

	bless $self, $class;	
	return $self;	
}

sub calculate_indexes {
	my $self = shift;
	my ($bed_definition) = @_;

	my $bed_parser = BEDParser -> new;
	my %bed_values = $bed_parser -> parse_gene_definition($bed_definition);
	
	my $sam = Bio::DB::Sam -> new(-bam => $self -> { bam },
				      -fasta => $self -> { genome },
				      -autoindex => 0);

	print "Calculating read depths across exons\n";
	my $chr = $bed_values{ -chr };
	my $exons_ref = $bed_values{ -exons };
	my @exons = @$exons_ref;
	my @exon_depths = ();
	foreach my $exon(@exons) {
		my @coords = split(/,/, $exon); 
		my $depths = $self -> read_depths($sam, $chr, $coords[0], $coords[1]);	
		@exon_depths = (@exon_depths, @$depths);
	}	
	
	print "Calculating median\n";
	my $median_exon_read_depth = $self -> median(\@exon_depths);

	print "Calculating all depths\n";	
	my $gene_start = $bed_values{ -genebegin };
	my $gene_stop = $bed_values{ -geneend };
	my $all_depths_ref = $self -> read_depths($sam, $chr, $gene_start, $gene_stop);
	my @all_depths = @$all_depths_ref;

	my @indexes = ();
	foreach my $depth(@all_depths) {
		my $index = $depth;
		if($median_exon_read_depth) {
			$index = ($depth / $median_exon_read_depth);
		}

		push(@indexes, $index);
	}

	return \@indexes;
}

sub median {
	my $self = shift;
	my ($array_ref) = @_;

	my @array = @$array_ref;

	my $median;
	my $mid = int @array/2;
	if(@array % 2) {
		$median = $array[ $mid ];
	}
	else {
		$median = ($array[ $mid - 1 ] + $array[ $mid ])/2;
	}

	return $median;
}

sub unique {
	my $self = shift;
	my ($array_ref) = @_;

	my @array = @$array_ref;
	my %observed = ();
	my @uniques = ();
	foreach my $value(@array) {
		my $boolean = $observed{ $value };
		unless ($boolean) {
			push(@uniques, $value);
			$observed{ $value } = 1;
		}
	}
	
	return \@uniques;
}

sub read_depths {
	my $self = shift;
	my ($sam, $chr, $start, $end) = @_;
	
	my ($coverage) = $sam -> features(-type => 'coverage', -seq_id => $chr, -start => $start,  -end => $end); 
	my @data_points = $coverage -> coverage;
	
	my $size = @data_points;
	if($size != $end - $start + 1) {
		my $index = 0;
		while($index < $end - $start + 1) {
			push(@data_points, 0);
			
			$index++;
		}
	}

	return \@data_points;
}

1;
