opendir(my $dh, "/Volumes/MacintoshHD/Splice/Alignments/hg18/Argonne/Unmapped/") || die; 
my @db=grep{ /\.[a-z]{3}$/ } readdir($dh); 
print @db;
