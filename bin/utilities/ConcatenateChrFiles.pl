#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;

my $chr_files_path = '';
my $output_file = '';
GetOptions('path=s' => \$chr_files_path, 'output=s' => \$output_file);

#print $chr_files_path . "\n" . $output_file . "\n";

opendir(my $dh, $chr_files_path) || die "No such directory $chr_files_path";
my @chr_files = grep{ /\.fa$/ } readdir($dh);
close $dh;

my $out = FileHandle -> new;
if($out -> open("> $chr_files_path/$output_file")) {
	foreach my $chr_file(@chr_files) {
		print "$chr_file\n";

		my $in = FileHandle -> new;
		if($in -> open("< $chr_files_path/$chr_file")) {
			while(my $line = $in -> getline) {
				$out -> print("$line");
			}
		}
	
		$in -> close;
	}
}

