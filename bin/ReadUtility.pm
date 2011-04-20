#!/usr/bin/perl -w

package ReadUtility;

use strict;

use FileHandle;
use File::Basename;
use File::Path;
use Bio::DB::Sam;
use Getopt::Long;

require BEDParser;

sub new {
	my $class = shift;
	my $self = { total_mapped => 0 };

	bless $self, $class;
	return $self;
}

sub CreateBLASTDB {
	my $self = shift;
	my ($file) = @_;
	
	# TODO: make storage location variable
	my $blast_storage = "/Volumes/Raid-1/BlastStorage";

	if(-e $file) {
		system("mkblastdb -in $file -out $blast_storage -dbtype nucl");
	}
	else {
		die "No such fasta file for unmapped reads\n";
	}
}

sub MoveFiles {
	my $self = shift;
	my ($dir, $mvdir, $files_ref) = @_;

	my @files = @$files_ref;
	foreach my $file(@files) {
		system("mv $dir\/$file $mvdir");
	}
}

sub GetUnmappedReads {
	my $self = shift;
	my ($bam, $genome, $append_to) = @_;

	my $sam = Bio::DB::Sam -> new(-bam => $bam,
			      	      -fasta => $genome,
				      -autoindex => 0);

	my @reads = $sam -> features(-type => 'match');

	my $dir = dirname($bam);
	$dir = $dir . "\/Unmapped\/";
	unless(-d $dir) {
		mkpath($dir);	
	}
	
	my $file_name = basename($bam);
	$file_name = $file_name . "\.unmapped\.fasta";

	my $out = FileHandle -> new;
	if($append_to ne '') {
		if($out -> open(">> $append_to")) {
			foreach my $read(@reads) {
				if($read -> unmapped) {
					$out -> print(">" . $read -> query -> name . " | $file_name\n");
					$out -> print($read -> query -> dna . "\n");
				}
			}
		}

		return $append_to;
	}
	else {
		if($out -> open("> " . $dir . $file_name)) { print $dir . $file_name . "\n";
			foreach my $read(@reads) {
				if($read -> unmapped) {
					$out -> print(">" . $read -> query -> name . " | $file_name\n");
					$out -> print($read -> query -> dna . "\n");
				}
			}
		}
	}
	
	$out -> close;

	return $dir. $file_name;
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

	return 1;
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

	print "C:$C\nL:$L\nN:$N\n";
	
	my $rpkm_val = (10 ** 9 * $C)/($N * $L);
	
	return $rpkm_val;
}

1;
