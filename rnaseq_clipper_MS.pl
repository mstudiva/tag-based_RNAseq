#!/usr/bin/perl

$usage= "

rnaseq_clipper.pl  : 

Clips 5'-leader off Illumina fastq reads in RNA-seq
Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (mstudiva@fau.edu) 

Removes duplicated reads sharing the same degenerate header and 
the first 6 bases of the sequence (reads containing N bases in this
region are discarded, too)

prints to STDOUT

arguments:
1 : fastq file name
2 : string to define the leading sequence, default '[ATGC]{4}G{3}'
'keep' : optional flag to say whether the sequences without leader should be kept. 
		 By default, they are discarded.

Example:
rnaseq_clipper.pl D6.fq

					 
";

my $fq=shift or die $usage; # shift pulls the filename into variable fq, die kills script if no filename is given, then loops back to beginning of script# shift pulls the filename into variable fq, die kills script if no filename is given, then loops back to beginning of script
my $lead=""; # adaptor variable is identified but left undefined
my $keep=1; # keep variable is identified but left undefined, default is to remove sequences without headers
if ($ARGV[0]) { $lead=$ARGV[0];} # if a first argument is given, use it for the adaptor sequence
else { $lead="[ATGC]{4}G{3}";} # otherwise use the default
if ($ARGV[1]) {$keep=1;} # if second argument is given, sequences missing adaptor are kept
# specifying variable names to prevent autovivification
my $trim=0; # trim variable is identified but left undefined
my $name=""; # name variable is identified but left undefined
my $name2=""; # name2 variable is identified but left undefined
my $seq=""; # sequence variable is identified but left undefined
my $qua=""; # quality score variable is identified but left undefined
my %seen={}; # seen variable is identified but left undefined
open INP, $fq or die "cannot open file $fq\n"; # open input of given filename, or kills script with message
my $ll=3; # ll is initially set to three (option 1)
my $nohead=0; # no header variable is identified but left undefined
my $dups=0; # duplicates variable is identified but left undefined
my $ntag; # number of tag variable is identified but left undefined
my $tot=0; # total variable is identified but left undefined
my $goods=0; # good variable is identified but left undefined
while (<INP>) { # while reading input files
	if ($ll==3 && $_=~/^(\@.+)$/ ) { # if option 1 and variable _ is something (.+) between the beginning (^) and end ($) of a string starting with @ (header of each read)
		$tot++; # returns the next consecutive value of total, aka counts total reads
		$name2=$1; # name2 is set to the read header starting with @
		if ($seq=~/^($lead)(.+)/) {	# if the sequence starts with the adaptor, followed by something (read sequence)
			my $start=substr($2,0,6); # start variable is the first 6 characters of the read sequence
			my $idtag=$1.$start; # defines idtag as the adaptor.first 6 bases of sequence
			if (!$seen{$idtag} and $idtag!~/N/) { # if the particular idtag is seen and does not include an N
				$seen{$idtag}=1; # sets the count for that particular idtag equal to one
				$trim=length($1); # and trim is equal to the length of the adaptor
				print "$name\n$2\n+\n",substr($qua,$trim),"\n"; # prints read header, sequence, +, quality scores minus the corresponding adaptor regions, separated by new lines
				$goods++; # adds +1 for every new adaptor seen
			}
			elsif ($seen{$idtag}) { $dups++; } # or if non-unique adaptor seen, add +1 to duplicate count
			else { $ntag++; } # and adds +1 to number of duplicate tags found
		}
		elsif ($keep and $name) { print "$name\n$seq\n+\n",$qua,"\n";} # if keep and name are defined, prints name, sequence, +, quality scores with new lines in between
		else {$nohead++;} # if no adaptor is found, add +1 to nohead count
		$seq=""; # sequence variable is identified but left undefined
		$ll=0; # ll variable is identified but left undefined
		$qua=""; # quality score variable is identified but left undefined
		@sites=(); # sites array variable is identified but left undefined
		$name=$name2; # sets the something string to be the original sequence
	}
	elsif ($ll==0){ # after the above, if ll is equal to zero (noheader, option 2)
		chomp; # remove any newline characters from the end of each read
		$seq=$_; # the sequence is kept the same without trimming (noheader kept)
		$ll=1; # ll set to one
	}
	elsif ($ll==2) { # or if ll is equal to two (option 3)
		chomp; # remove any newline characters from the end of each read
		$qua=$_; # quality scores kept the same without trimming (full sequence)
		$ll=3; # ll set to three
	}
	else { $ll=2;} # otherwise ll is set to two
}

                if ($seq=~/^($lead)(.+)/) {
                        my $start=substr($2,0,6);
                        my $idtag=$1.$start;
                        if (!$seen{$idtag} and $idtag!~/N/) {
                                $seen{$idtag}=1;
                                $trim=length($1);
                                print "$name\n$2\n+\n",substr($qua,$trim),"\n";
                                $goods++;
                        }
                        elsif ($seen{$idtag}) { $dups++; }
                        else { $ntag++; }
                }
elsif ($keep and $name) { print "$name\n$seq\n+\n",$qua,"\n";}
# outputs the following strings with tabs in between:
# filename, total reads, good reads, duplicates, missing header, number of N nucleotides in header
warn "$fq\ttotal:$tot\tgoods:$goods\tdups:$dups\tnoheader:$nohead\tN.in.header:$ntag\n";
