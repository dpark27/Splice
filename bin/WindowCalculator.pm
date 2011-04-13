#!/usr/bin/perl -w

package WindowCalculator;

use strict;

use Getopt::Long;
use Bio::DB::Sam;

sub new {
	my $class = shift;
	my %params = @_;
	
	for my $required(qw{ bam genome size }) {
		if($params{ $required } eq '') {
			die "Required parameter $required not provided";
		}
	}

	my $self = { bam => $params{ 'bam' },
     		     genome => $params{ 'genome' },
		     size => $params{ 'size' }};

	bless $self, $class;	
	return $self;	
}

sub calculate_gene_window {
	my $self = shift;
	my ($bed_definition) = @_;

	my %bed_values = $self -> parse_bed_definition($bed_definition);
	
	my $sam = Bio::DB::Sam -> new(-bam => $self -> { bam },
				      -fasta => $self -> { genome });

	my $depths = $self -> read_depth($sam, "chr12", "113344738", "113345024");	
	my $unique_depths = $self -> unique($depths);
	print "@$unique_depths\n";
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

sub read_depth {
	my $self = shift;
	my ($sam, $chr, $start, $end) = @_;
	
	my @depths = ();	
	my $callback = sub {
		my ($seqid, $pos, $pileup) = @_;
	
		if($pos >= $start && $pos <= $end) {
			my $depth = @$pileup;
			push(@depths, $depth);
		}
	};

	my $segment = "$chr:$start\-$end";	
	$sam -> pileup($segment, $callback);
	
	return \@depths;	
}

sub parse_bed_definition {
	my $self = shift;
	my ($bed_definition) = @_;

	my @data = split(/\t/, $bed_definition);

	my @exon_starts = split(/,/, $data[8]);
	my @exon_ends = split(/,/, $data[9]);
	my @exon_coords = ();
	
	my $i = 0;
	while($i < @exon_starts) {
		my $coord = $exon_starts[$i] . "," . $exon_ends[$i];
		$i++;
		push(@exon_coords, $coord);
	}	
	
	my %return_values = (-chr => $data[1],
			     -start => $data[3],
			     -end => $data[4],
			     -exons => \@exon_coords);
	
	return %return_values;
}

1;
