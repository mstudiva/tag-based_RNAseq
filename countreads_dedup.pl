
#!/usr/bin/perl

print "
countreads_dedup.pl : counts the number of Illumina reads in a bunch of deduplicated files
argument - glob to fastq files, default \.dedup
";

my $glob="\.dedup";
if ($ARGV[0]) { $glob=$ARGV[0];}

opendir THIS, ".";
my @fqs=grep /$glob/,readdir THIS;
my $f;
my $nrd;
foreach $f(@fqs){
	$nrd=`cat $f | wc -l`;
	$nrd=$nrd/4;
	print "$f\t$nrd\n";
}
