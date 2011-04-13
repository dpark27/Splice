#!/usr/bin/perl -w

package RPKMCalculator;

use strict;

use Bio::DB::Sam;
use Getopt::Long;

require BEDParser;

sub new {
	my $class = shift;
	my $self = { total_mapped => 0 };

	bless $self, $class;
	return $self;
}

sub GetTotalMappedCount {
	my $self = shift;
	my($sam) = @_;

	my @reads = $sam -> features(-type => 'match');
	
	my $mapped_count = 0;
	foreach my $read(@reads) {
		unless($read -> unmapped) {
			$mapped_count++;	
		}
	}	

	$self -> { total_mapped } = $mapped_count;
}

sub GetTranscriptLength {
	my $self = shift;
	my ($exons_ref) = @_;

	my $transcript_length = 0;
	my @exons = @$exons_ref;
	foreach my $exon(@exons) {
		my @coords = split(/,/, $exon);
		
		$transcript_length = $transcript_length + ($coords[1] - $coords[0]);
	}

	return $transcript_length;
}

sub GetMappedToExonCount {
	my $self = shift;
	my ($chr, $sam, $exons_ref) = @_;

	my $mapped_to_exon_count = 0;

	my @exons = @$exons_ref;
	foreach my $exon(@exons) {
		my @coords = split(/,/, $exon);
	
		my $segment = $sam -> segment($chr, $coords[0], $coords[1]);
		my @reads = $segment -> features;
		my $reads_count = @reads;
		
		$mapped_to_exon_count = $mapped_to_exon_count + $reads_count;
	}

	return $mapped_to_exon_count;
}

sub CalculateRPKM {
	my $self = shift;
	my($bed_definition, $bam, $genome) = @_;

	my $bed_parser = BEDParser -> new;
	my %bed_values = $bed_parser -> parse_gene_definition($bed_definition);

	my $sam = Bio::DB::Sam -> new(-bam => $bam,
				      -fasta => $genome,
				      -autoindex => 0); 

	my $N = $self -> { total_mapped }; 
	if($N == 0) {
		$self -> GetTotalMappedCount($sam);
		$N = $self -> { total_mapped };
	}

	my $C = $self -> GetMappedToExonCount($bed_values{ -chr }, $sam, $bed_values{ -exons });
	my $L = $self -> GetTranscriptLength($bed_values{ -exons });

	my $rpkm_val = (10^9 * $C)/($N * $L);
	
	return $rpkm_val;
}

1;
