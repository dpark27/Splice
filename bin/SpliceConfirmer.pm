#!/usr/bin/perl -w

package SpliceConfirmer;

use File::Basename;
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

sub CreateBlastableFile {
	my $self = shift;
	my ($splices_ref, $genotype) = @_;
	
	my ($fives_ref, $threes_ref) = $self -> CreateSequences($splices_ref);
	my $combos_ref = $self -> PermuteCombinations($fives_ref, $threes_ref);
	my @combos = @$combos_ref; 
	
	my $out = FileHandle -> new;
	if($out -> open("> /Volumes/MacintoshHD/Splice/ReadyToBlast/" . $genotype . ".fasta")) { print "here\n";
		foreach my $combo(@combos) {
			$out -> print($combo);	
		}

		$out -> close;
	}


}

sub PermuteCombinations {
	my $self = shift;
	my ($fives_ref, $threes_ref) = @_;

	my @fives = @$fives_ref; #print "@fives\n"; 
	my @threes = @$threes_ref; #print "@threes\n";

	my @combos = ();
	foreach my $three(@threes) {
		my @three_data = split(/,/, $three);
		
		my $three_marker = $three_data[0];
		my $three_point = $three_data[1];
		my $three_dna = $three_data[2];

		foreach my $five(@fives) {
			my @five_data = split(/,/, $five);
			
			my $five_marker = $five_data[0];
			my $five_point = $five_data[1];
			my $five_dna = $five_data[2];

			my $combo = ">$three_point | $five_point\n" . "$three_dna" . "$five_dna\n"; #print "COMBO: $combo\n";
			push(@combos, $combo);
		}
	}

	return \@combos;
}

sub CreateSequences {
	my $self = shift;
	my ($splices_ref) = @_;

	my $sam = Bio::DB::Sam -> new(-bam => $self -> { bam },
				      -fasta => $self -> { genome },
				      -autoindex => 0);
	my @fives = ();
	my @threes = ();	

	my %splices = %$splices_ref;
	while(my ($point, $marker) = each %splices) {
		my @data = split(/\:/, $point);
		my $chr = $data[0];
		my $site = $data[1]; print "$marker\n"; 

		my $dna;
		if($marker == 5) {
			$dna = $sam -> seq($chr, $site, $site + 20); #print "push 5\n";
			
			my $entry = "$marker,$point,$dna"; #print $entry . "\n";
			push(@fives, $entry);
		}
		else {
			$dna = $sam ->seq($chr, $site - 20, $site); #print "push 3\n";
		
			my $entry = "$marker,$point,$dna"; #print $entry . "\n";
			push(@threes, $entry);
		}
	}

	
	return (\@fives, \@threes);
}

1;
