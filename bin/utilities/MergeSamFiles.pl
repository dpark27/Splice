#!/usr/bin/perl -w

use strict;

use Getopt::Long;

GetOptions('dir=s' => \(my $dir = ''),
	   'out=s' => \(my $out = ''));

if($dir ne '' && $out ne '') {
	opendir(my $dh, $dir);
	my @files = grep{ /\.bam$/ } readdir($dh);
	closedir($dh);
	
	my $list_in = '';
	$dir =~ s/\s/\\ /g;
	foreach my $file(@files) {
		$list_in = $list_in . " $dir". $file;
	}
	
	$out =~ s/\s/\\ /g;

	print "samtools merge $out $list_in\n";
	system("samtools merge $out $list_in");
}
else {
	die "No directory or output name provided\n";
}
