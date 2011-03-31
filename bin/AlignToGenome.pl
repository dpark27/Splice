#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use File::Basename;

GetOptions('reference=s' => \(my $reference = ''),
	   'data=s' => \(my $data = ''),
	   'ibase=s' => \(my $index_base = ''),
	   'threads=s' => \(my $threads = 1));

&PerformAlignment($reference, $data, $index_base, $threads);

sub PerformAlignment {
	my ($refernce, $data, $index_base, $threads) = @_;

	my $bowtie_index_dir = $ENV{ 'BOWTIE_INDEXES' };
	my $index_arg = "$bowtie_index_dir$index_base";

	my $threads_arg = "-p $threads";
	my $sam_arg = "-S";

	my $data_arg = $data;
	$data_arg =~ s/\s/\\ /g;

	my $output_dir = $ENV{ 'SPLICE' } . "Alignments/";
	my $output_base = basename($reference);
	my @output_base_parts = split(/\./, $output_base);
	$output_base = $output_base_parts[0];
	my $output_arg = $output_dir . $output_base . "\." . time . "\.sam";

	system("bowtie --best -l 28 -e 200 -n 2 $index_arg $data_arg $threads_arg $sam_arg $output_arg");

	my $sam_input = $output_arg;
	my $bam_output = $output_arg;
	$bam_output =~ s/\.sam$/\.bam/;
	system("samtools view -bSo $bam_output $sam_input");

	my $sorted_bam_output = $bam_output;
	$sorted_bam_output =~ s/\.bam$/\.sorted/gi;
	system("samtools sort $bam_output $sorted_bam_output");
	
	$sorted_bam_output = $sorted_bam_output . ".bam";
	system("samtools index $sorted_bam_output");
}

