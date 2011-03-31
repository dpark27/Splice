#!/usr/bin/perl -w

use strict;

use FileHandle;
use File::Basename;
use Getopt::Long;

GetOptions('dir=s' => \(my $data_dir = ''));

opendir(my $dh, $data_dir) || die "No such directory $data_dir";
my @data_files = grep{ /.fastq$/ } readdir($dh);
close $dh;

my %file_names = ();
foreach my $file(@data_files) { 
	my $file_basename = basename($file);
	my @base_parts = split(/_/, $file_basename);
	if(!$file_names{ $base_parts[0] }) { 
		$file_names{ $base_parts[0] } = $file_basename;
	}
	else {
		$file_names{ $base_parts[0] } = $file_names{ $base_parts[0] } . ";$file_basename";
	}
}

my @sort = sort keys %file_names;
foreach my $key(@sort) {
	my $out = FileHandle -> new;
	
	print "Now processing $key\n";	
	if($out -> open("> $data_dir$key" . "\.fastq")) {
		my $in_list = $file_names{ $key };
		my @in_files = split(/;/, $in_list);
		foreach my $in_file(@in_files) {
			print "File: $in_file\n";		
	
			my $in = FileHandle -> new;
			if($in -> open("< $data_dir$in_file")) {
				while(my $line = $in -> getline) {
					$out -> print($line);
				}
			}

			$in -> close;
		}	
	}

	$out -> close;
}
