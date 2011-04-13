#!/usr/bin/perl -w

package BEDParser;

use strict;

sub new {
	my $class = shift;
	my $self = {};
	
	bless $self, $class;
	return $self;
}

sub parse_gene_definition {
        my $self = shift;
        my ($bed_definition) = @_;

        my @data = split(/\t/, $bed_definition);

        my @exon_starts = split(/,/, $data[9]);
        my @exon_ends = split(/,/, $data[10]);
        my @exon_coords = ();

        my $i = 0;
        while($i < @exon_starts) {
                my $coord = $exon_starts[$i] . "," . $exon_ends[$i];
                $i++;
                push(@exon_coords, $coord);
        }

        my %return_values = (-id => $data[1],
                             -chr => $data[2],
                             -genebegin => $data[4],
                             -geneend => $data[5],
                             -exons => \@exon_coords);

        return %return_values;
}

1;
