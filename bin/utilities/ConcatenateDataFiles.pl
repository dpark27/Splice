#!/usr/bin/perl -w

use strict;

use FileHandle;
use File::Basename;
use Getopt::Long;

GetOptions('dir=s' => \(my $data_dir = ''),
	   'all=s' => \(my $concat_all = 0));

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

if($concat_all) {
	my $out = FileHandle -> new;
	if($out -> open("> $data_dir" . "all\.fastq")) {
		print "Now processing all.fastq\n";

		foreach my $key(@sort) {
			my $in = FileHandle -> new;
			if($in -> open("< $data_dir$key\.fastq")) {
				print "File: $key.fastq\n";

				while(my $line = $in -> getline) {
					$out -> print($line);
				}
			}

			$in -> close;
		}
	}

	$out -> close;
}
