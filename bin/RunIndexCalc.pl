#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;
use File::Basename;

require IndexCalculator;
require RPKMCalculator;
require BEDParser;

GetOptions('bam=s' => \(my $bam = ''),
	   'genome=s' => \(my $genome = ''),
	   'bedfile=s' => \(my $bed_file = ''),
	   'outputdir=s' => \(my $output_dir = ''));

&CreateIndexFiles($bam, $genome, $bed_file, $output_dir);

sub CreateIndexFiles {
	my ($bam, $genome, $bed_file, $output_dir) = @_;

	my $index_calc = IndexCalculator -> new('bam' => $bam,
					  	'genome' => $genome);

	my $rpkm_calc = RPKMCalculator -> new;
	
	my $in = FileHandle -> new;
	if($in -> open("< $bed_file")) { 
		my $output_file = basename($bam);	

		my $out = FileHandle -> new;
		if($out -> open("> $output_dir\\$output_file\.index")) { 
			my $bed_parser = BEDParser -> new;
			while(my $bed_definition = $in -> getline) {
				my %bed_values = $bed_parser -> parse_gene_definition($bed_definition);
				print "$bed_values{ -id }\n";
				
				my $rpkm = $rpkm_calc -> CalculateRPKM($bed_definition, $bam, $genome);
				print "RPKM: $rpkm\n";	
	
				my $gene_def = "$bed_values{ -id };$bed_values{ -chr };$bed_values{ -genebegin };$bed_values{ -geneend };";
				
				my $exon_defs = "";
				my $exon_coords_ref = $bed_values{ -exons };
				my @exon_coords = @$exon_coords_ref;
				foreach my $coord(@exon_coords) {
					$coord =~ s/,/;/g;
					$exon_defs = $exon_defs . ";" .  $coord;
				}
				
				$out -> print($gene_def);
				$out -> print($exon_defs);
				$out -> print(";");

				my $indexes_ref = $index_calc -> calculate_indexes($bed_definition);
				my @indexes = @$indexes_ref;
				foreach my $index(@$indexes_ref) {
					$out -> print($index . ',');	
				}
			
				$out -> print("\n");
			}
		
			$out -> close;
		}

		$in -> close;
	}
}
