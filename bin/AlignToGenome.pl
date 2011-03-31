#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use File::Basename;
use Bio::Tools::Run::Bowtie;

GetOptions('reference=s' => \(my $reference = ''),
	   'data=s' => \(my $data = ''),
	   'ibase=s' => \(my $index_base = ''),
	   'threads=s' => \(my $threads = 1));

&PerformAlignment($reference, $data, $index_base, $threads);

sub PerformAlignment {
	my ($refernce, $data, $index_base, $threads) = @_;

	# TODO: Enable in future when BioPerl has more support for Bowtie and out of beta
	#my $bowtie_factory = Bio::Tools::Run::Bowtie -> new;
	#$bowtie_factory -> run(-seq => $data, -ind => $index_base);

	my $bowtie_index_dir = $ENV{ 'BOWTIE_INDEXES' };
	my $index_arg = "$bowtie_index_dir$index_base";
	$index_arg =~ s/ /\\ /g;
	
	my $threads_arg = "-p $threads";
	my $sam_arg = "-S --sam-nohead";

	my $data_arg = $data;
	$data_arg =~ s/ /\\ /g;

	my $output_dir = $ENV{ 'SPLICE' } . "Alignments/";
	my $output_base = basename($reference);
	my @output_base_parts = split(/\./, $output_base);
	$output_base = $output_base_parts[0];
	my $output_arg = $output_dir . $output_base . "\." . time . "\.sam";
	$output_arg =~ s/ /\\ /g;

	system("bowtie $index_arg $data_arg $threads_arg $sam_arg $output_arg");
}
