#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;
use File::Basename;

require IndexCalculator;
require ReadUtility;
require BEDParser;
require SegmentSpliceTest;
require SpliceConfirmer;

GetOptions('bam=s' => \(my $bam = ''),
	   'bams=s' => \(my $bams_dir = ''),
	   'genome=s' => \(my $genome = ''),
	   'bedfile=s' => \(my $bed_file = ''),
	   'outputdir=s' => \(my $output_dir = ''));

if($bams_dir ne '') {
	&RunOnMultipleBams($bams_dir, $bed_file, $genome);
}
else {
	&CreateIndexFiles($bam, $genome, $bed_file, $output_dir);
}

sub RunOnMultipleBams {
	my ($bams_dir, $bed_file, $genome) = @_;

	opendir(my $dh, $bams_dir);
	my @bams = grep { /\.sorted\.bam$/ } readdir($dh);
	closedir($dh);

	my $bed_in = FileHandle -> new;
	if($bed_in -> open("< $bed_file")) {
		my $bed_parser = BEDParser -> new;
		my $read_utility = ReadUtility -> new;
	
		while(my $bed_definition = $bed_in -> getline) {
			my %splice_set = ();
			my $fasta_file = '';			

			my %bed_values = $bed_parser -> parse_gene_definition($bed_definition);
			print "$bed_values{ -id }\n";

			my $bam_file;
			foreach my $bam(@bams) { 
				# TODO: only uncomment if you need to create blastdb fasta
				#if($fasta_file eq '') {
				#	$fasta_file = $read_utility -> GetUnmappedReads($bams_dir . $bam, $genome, "/Volumes/Raid-1/BlastStorage/all.fa");
				#}
				#else {
			#		$read_utility ->  GetUnmappedReads($bams_dir . $bam, $genome, $fasta_file);
			#	}
				
				my $rpkm = $read_utility -> CalculateRPKM($bed_definition, $bams_dir . $bam, $genome);
				print "RPKM: $rpkm\n";
				if($rpkm > 100) { 
					$bam_file = $bams_dir . $bam;
					my $index_calc = IndexCalculator -> new('bam' => $bam_file,
										'genome' => $genome);

					my $indexes_ref = $index_calc -> calculate_indexes($bed_definition);
					my $scanner = SegmentSpliceTest -> new('readlength' => 46);

					my $possible_splices_ref = $scanner -> ScanForSplice($indexes_ref, .0001, $bed_values{ -genebegin }, $bed_values{ -chr });
					if(!keys(%splice_set)) {
						my @poss_splices = @$possible_splices_ref;
						foreach my $poss_splice(@poss_splices) {
							my @data = split(/,/, $poss_splice);
							
							my $point = $data[0];
							my $type = $data[1];

							$splice_set{ $point } = $type;
						}	
					}
					else {
						my %intersect = ();
						my @poss_splices = @$possible_splices_ref;
						foreach my $poss_splice(@poss_splices) {
							my @data = split(/,/, $poss_splice);
							
							my $point = $data[0];
							my $type = $data[1];
							if($splice_set{ $point }) {
								$intersect{ $point } = $splice_set{ $point };
							}
													
						}

						%splice_set = %intersect;
					}
				}
			}
	
			print "\n\n";
			my @sortorder = sort keys %splice_set;
			my $count = 0;
			foreach my $key(@sortorder) {
				print $count . ": " . $key . " " . $splice_set{ $key } . "\n";
				$count++;
			}
		
#			$splice_set{ 'chr12:111838909' } = 3;
#			$splice_set{ 'chr12:111839000' } = 5;
			#my $confirmer = SpliceConfirmer -> new('bam' => $bam_file,
							   #    'genome' => $genome);
	
			#$confirmer -> CreateBlastableFile(\%splice_set, $bed_values{ -id });	
		}
	}
}

sub CreateIndexFiles {
	my ($bam, $genome, $bed_file, $output_dir) = @_;

	my $index_calc = IndexCalculator -> new('bam' => $bam,
					  	'genome' => $genome);

	my $read_utility = ReadUtility -> new;
	#my $fasta_file = $read_utility -> GetUnmappedReads($bam, $genome);
	#$read_utility -> CreateBLASTDB($fasta_file);
	
	my $in = FileHandle -> new;
	if($in -> open("< $bed_file")) { 
		my $output_file = basename($bam);	

		my $out = FileHandle -> new;
		if($out -> open("> $output_dir\$output_file\.index")) { 
			my $bed_parser = BEDParser -> new;
			while(my $bed_definition = $in -> getline) {
				my %bed_values = $bed_parser -> parse_gene_definition($bed_definition);
				print "$bed_values{ -id }\n";
				
				# TODO: determine minimum RPKM for saying gene is expressed
				#my $rpkm = $read_utility -> CalculateRPKM($bed_definition, $bam, $genome);
				#print "RPKM: $rpkm\n";

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

				my $scanner = SegmentSpliceTest -> new( 'readlength' => 46 );
				$scanner -> ScanForSplice($indexes_ref, .0001, $bed_values{ -genebegin });

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
