#!/usr/bin/perl

print "

countreads_host.pl : counts the number of host-mapped reads in a bunch of fastq files
argument - glob to fastq files, default \.host.clean

";

my $glob="\.host.clean";
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
