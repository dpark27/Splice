#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use File::Path qw(make_path);
use File::Basename;

GetOptions('reference=s' => \(my $reference = ''),
	   'data=s' => \(my $data = ''),
	   'outputdir=s' => \(my $output_dir = ''),
	   'datadir=s' => \(my $data_dir = ''),
	   'ibase=s' => \(my $index_base = ''),
	   'threads=s' => \(my $threads = 1));

if($data_dir ne '') { 
	opendir(my $dh, $data_dir) || die "Cannot open $data_dir";
	my @files = grep{ /\.fastq$/ } readdir($dh);
	closedir($dh);

	foreach my $data(@files) {
		&PerformAlignment($reference, $output_dir, "$data_dir$data", $index_base, $threads);
	}
}
else {
	if($data ne '') {
		&PerformAlignment($reference, $output_dir, $data, $index_base, $threads);
	}
	else {
		die "No Data";
	}
}

sub PerformAlignment {
	my ($refernce, $output_dir, $data, $index_base, $threads) = @_;

	my $bowtie_index_dir = $ENV{ 'BOWTIE_INDEXES' };
	my $index_arg = "$bowtie_index_dir$index_base";

	my $threads_arg = "-p $threads";
	my $sam_arg = "-S";

	my $data_arg = $data;
	$data_arg =~ s/\s/\\ /g;

	unless(-d "$output_dir") {
		make_path("$output_dir");
	} 

	my $output_base = basename($data);
	my @output_base_parts = split(/\./, $output_base);
	$output_base = $output_base_parts[0];
	
	my $output_arg = $output_dir . $output_base . "\." . time . "\.sam";
	$output_arg =~ s/\s/\\ /g;
	
	print "bowtie --best -l 24 -e 170 $index_arg $data_arg $threads_arg $sam_arg $output_arg\n";
	system("bowtie --best -l 24 -e 170 $index_arg $data_arg $threads_arg $sam_arg $output_arg");

	my $sam_input = $output_arg;
	my $bam_output = $output_arg;
	$bam_output =~ s/\.sam$/\.bam/;
	
	print "samtools view -bSo $bam_output $sam_input\n";
	system("samtools view -bSo $bam_output $sam_input");

	my $sorted_bam_output = $bam_output;
	$sorted_bam_output =~ s/\.bam$/\.sorted/gi;
	
	print "samtools sort $bam_output $sorted_bam_output\n";
	system("samtools sort $bam_output $sorted_bam_output");
	
	$sorted_bam_output = $sorted_bam_output . ".bam";
	
	print "samtools index $sorted_bam_output\n";
	system("samtools index $sorted_bam_output");
}

