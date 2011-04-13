#!/usr/bin/perl -w

use strict;

use Getopt::Long;
require WindowCalculator;

GetOptions('bam=s' => \(my $bam = ''),
	   'genome=s' => \(my $genome = ''),
	   'size=s' => \(my $size = 10));

my $calc = WindowCalculator -> new('bam' => $bam,
				   'genome' => $genome,
				   'size' => $size);

$calc -> calculate_gene_window("uc009zwf.2	chr12	+	113344738	113369988	113344847	113369727	6	113344738,113346340,113348855,113354313,113355351,113369682,	113345024,113346629,113349040,113354543,113355505,113369988,	P00973-4	uc009zwf.2");
