my $usage = "perl check-cmsearch-v-cmscan.pl <cmsearch non-deoverlapped tblout output> <cmscan --fmt 2 --oskip tblout output>";
if(scalar(@ARGV) != 2) { 
  die $usage;
}

my($cmsearch_file, $cmscan_file) = (@ARGV);

# deoverlap $cmsearch_file
system("perl ./cmsearch-deoverlap.pl $cmsearch_file");
my $deoverlap_file = $cmsearch_file . ".deoverlapped";

# strip cmsearch
my $strip_cmsearch_file = $deoverlap_file . ".strip";
system("perl ./strip-cmsearch-tblout.pl $deoverlap_file | sort > $strip_cmsearch_file");

# strip/convert cmscan
my $strip_cmscan_file = $cmscan_file . ".strip";
system("perl ./convert-cmscan-to-cmsearch-tblout.pl $cmscan_file | sort > $strip_cmscan_file");

$diff = `diff $strip_cmsearch_file $strip_cmscan_file`;

if($diff ne "") { 
  die "ERROR diff failed:\n$diff\n";
}
else { 
  printf("Files are identical\n");
}

