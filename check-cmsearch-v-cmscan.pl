use Getopt::Long;

my $usage;
$usage  = "cmsearch-deoverlap v0.05 [Dec 2019]\n\n";
$usage .= "Usage:\n\n";
$usage .= "perl check-cmsearch-v-cmscan.pl [OPTIONS] <cmsearch non-deoverlapped tblout output> <cmscan --fmt 2 --oskip tblout output>\n\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t--maxkeep    : only remove hits that have a higher scoring overlap that is not removed\n";
$usage .= "\t\t--clanin <s> : only remove overlaps within clans, read clan info from file <s> (cmscan run with --oclan --clanin)\n\n";

my $in_clanin  = undef; # defined if --clanin option used
my $do_maxkeep = 0;
&GetOptions( "maxkeep"  => \$do_maxkeep,
             "clanin=s" => \$in_clanin );

if(scalar(@ARGV) != 2) { 
  die $usage;
}

my($cmsearch_file, $cmscan_file) = (@ARGV);

# deoverlap $cmsearch_file
my $opts = ""; 
if(defined $in_clan) { $opts .= " --clanin $in_clanin "; }
if($do_maxkeep)      { $opts .= " --maxkeep "; }

system("perl ./cmsearch-deoverlap.pl $opts $cmsearch_file");
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

