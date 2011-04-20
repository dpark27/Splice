#!/usr/bin/perl -w

package SpliceConfirmer;

use FileHandle;
use Bio::DB::Sam;

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

sub CreateBlastInput {
	my $self = shift;
	my ($slice_points_ref) = @_;

	my $seqs_ref = $self -> ConvertToSequences($splice_points_ref);

	# TODO: base this off bam file using basename
	my $output_file = "test.fa";	
	$self -> PermutePossibleJunctions($output_file, $seqs_ref);

	
}

sub RunBlast {

}

sub ParseBast {

}

sub PermutePossibleJunctions {
	my $self = shift;
	my ($output_file, $sequences_ref) = @_;

	my $out = FileHandle -> new;
	if($out -> open("> $output_file")) {
		
		my @seqs = @$sequences_ref;
		my $seqs_count = @seqs;
		my $i = 0;
		while($i < $seqs_count) {
			my $current_seq = pop(@seqs);
			my @current_data = split(/,/, $current_seq);
			my $current_marker = $current_data[0];
			my $current_chr = $current_data[1];
			my $current_start = $current_data[2];
			my $current_end = $current_data[3];
			my $current_dna = $current_data[4];

			foreach my $seq(@seqs) {
				my @data = split(/,/, $seq);
				
				my $marker = $data[0];
				if($marker != $current_marker) {
					my $chr = $data[1];
					my $start = $data[2];
					my $end = $data[3];
					my $dna = $data[4];

					my $fasta;
					if($current_marker == 3) {
						$fasta = ">$chr:$start\-$end , $chr:$current_start\-$current_end\n";
						$fasta = $fasta . $dna . $current_dna . "\n";
					}
					else {
						$fasta = ">$chr:$current_start\-$current_end , $chr:$start\-$end\n";
						$fasta = $fasta . $current_dna . $dna . "\n";
					}

					print $fasta;
				}	
			}

			push(@seqs, $current_seq);	
			$i++;
		}
	}
}

sub ConvertToSequence {
	my $self = shift;
	my ($splices_ref) = @_;
	
	my $sam = Bio::DB::Sam -> new(-bam => $self -> { bam },
				      -fasta => $self -> { genome },
				      -autoindex => 0);

	my @results = ();
	my @splices = @$splices_ref;
	foreach my $splice(@splices) {
		my @data = split(/,/, $splice);
		
		my $chr = $data[0];
		my $start = $data[1];
		my $test_result = $data[2];

		my $dna;
		my $end;
		my $marker;
		if($test_result > 0) {
			$end = $start + 20;
			$marker = 3;
			$dna = $sam -> seq($chr, $start, $end);
		}
		else {
			$end = $start;
			$start = $start - 20;
			$marker = 5;
			$dna = $sam -> seq($chr, $start, $end);
		}

		my $info = "$marker,$chr,$start,$end,$dna";
		push(@results, $info);
	}	

	return \@results;
}

1;
