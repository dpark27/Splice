#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use File::Basename;
use File::Copy;

my $reference_file = '';
my $index_base = '';
GetOptions('reference=s' => \$reference_file, 'ibase=s' => \$index_base);

my $genome_dirname = dirname($reference_file);
my $bowtie_indexes_dir = $ENV{ 'BOWTIE_INDEXES' };

system("bowtie-build \"$reference_file\" \"$bowtie_indexes_dir$index_base\"");
