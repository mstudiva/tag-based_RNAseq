# Tag-based RNA-seq reads processing pipeline, version November 4, 2021
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com)
# for use on the FAU KoKo HPC


#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions – please make sure to read them before copy-pasting.

# log onto cluster
ssh mstudiva@koko-login.hpc.fau.edu


#------------------------------
# installing RNA-seq scripts and setting up the workspace

# switch to home directory
cd

# unless you have done it in the past, make directory called bin,
# all your scripts should go in there:
mkdir bin

# switch to bin:
cd bin

# clone github repository
git clone https://github.com/mstudiva/tag-based_RNAseq.git

# move files from subdir tag-based_RNAseq-master to bin/:
mv tag-based_RNAseq/* .

# remove the tag-based_RNAseq directory
rm -rf tag-based_RNAseq

# remove the TACC version of launcher_creator.py from the bin directory (preinstalled on KoKo)
rm launcher_creator.py

# If you have not previously, download BaseSpaceCLI
wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O $HOME/bin/bs

chmod +x ~/bin/bs

# go to the website and confirm authorization by logging in to your BaseSpace acct.
bs auth

#------------------------------
## Download and concatenate reads with a launcher_creator script

echo '#!/bin/bash' > downloadReads.sh
echo 'bs download project --concurrency=high -q -n ####### -o .' >> downloadReads.sh
# -n is the BaseSpace project name and -o is the output directory

echo "find . -name '*.gz' -exec mv {} . \;" >> downloadReads.sh
echo 'rmdir SA*' >>downloadReads.sh
echo 'mkdir ../concatReads' >> downloadReads.sh
echo 'cp *.gz ../concatReads' >> downloadReads.sh
echo 'cd ../concatReads' >> downloadReads.sh
echo 'mergeReads.sh -o mergeTemp' >> downloadReads.sh
# -o is the directory to put output files in

echo 'rm *L00*' >> downloadReads.sh
echo "find . -name '*.gz' -exec mv {} . \;" >> downloadReads.sh
echo 'gunzip *.gz' >> downloadReads.sh
echo 'rmdir mergeTemp' >> downloadReads.sh

chmod +x downloadReads.sh

launcher_creator.py -b 'srun downloadReads.sh' -n downloadReads -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch --mem=200GB downloadReads.slurm

#-------------------------------
# Concatenating sequence files

# double check you have the number of files you should
ll *.fastq | wc -l

# If your samples are split across multiple files from different lanes,
# concatenating the corresponding fastq files by sample:
# ngs_concat.pl commonTextInFastqFilenames  "FilenameTextImmediatelyBeforeSampleID(.+)FilenameTextImmediatelyAfterSampleID"

echo "ngs_concat.pl 'Text-' '(.+)-\d'" > concat
launcher_creator.py -j concat -n concat -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch concat.slurm

# this one-liner replaces any sample numbers with the three digit version
for f in *.fq; do x="${f##*_}"; mv "$f" "${f%_*}$(printf '_%03d.fq' "${x%.fq}")"; done

# double check you have the correct number of files as samples
ll *.fq | wc -l

# look at the reads:
# head -50 SampleName.fq
# note that every read has four lines, the ID line starts with @HWI

# this little one-liner will show sequence-only in file:
# head -100 SampleName.fq | grep -E '^[NACGT]+$'

# to count the number of reads in all samples
echo "countreads_raw.pl > countreads_raw.txt" > count
launcher_creator.py -j count -n count -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count.slurm

#------------------------------
# Create conda environment
# Uncomment and run below if you don't have a conda env. set up
# module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

conda create -n sctld cutadapt

# Removing adaptors and low quality reads
echo '#!/bin/bash' > trim.sh
echo 'conda activate sctld' >> trim.sh
for F in *.fq; do
echo "tagseq_clipper.pl $F | cutadapt - -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o ${F/.fq/}.trim" >>trim.sh;
done

# Does not work with launcher_creator, consider breaking up script and running multiple jobs
sbatch -o trim.o%j -e trim.e%j --mem=200GB trim.sh

# how the job is doing?
squeue -u mstudiva

# double check you have the same number of files as samples
ll *.trim | wc -l

# but did the trimming really work?
# Use the same one-liner as before on the trimmed file to see if it is different
# from the raw one that you looked at before:

# head -100 SampleName.fq | grep -E '^[NACGT]+$'

# head -100 SampleName.trim | grep -E '^[NACGT]+$'
# the long runs of base A should be gone

# double-check that the rnaseq_clipper did not filter out too many reads by looking at the trim.e####### file
# make sure you're not losing too many reads to duplicates
# rename as a txt file
mv trim.e####### trim.txt

# to save time in case of issues, move the concatenated fq files to backup directory
mv *.fq ~/rawReads/

# to count the number of reads in trimmed samples
echo "countreads_trim.pl > countreads_trim.txt" > count_trim
launcher_creator.py -j count_trim -n count_trim -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count_trim.slurm

#------------------------------
# download and format reference transcriptome:

mkdir ~/db/
cd ~/db/
# copy your transcriptome .fasta file(s) to db/

# creating bowtie2 index for transcriptome(s)
# replace 'Host' and 'Symbiont' with your respective filenames
echo 'bowtie2-build Host.fasta Host' > btb
echo 'bowtie2-build Symbiont.fasta Symbiont' >> btb
launcher_creator.py -j btb -n btb -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch btb.slurm


#------------------------------
## Mapping reads to host/symbiont transcriptomes with bowtie2

# if working with split host/symbiont transcriptomes, need to run FIVE rounds of mapping: 1) all reads to symbiont, 2) remaining reads to coral, 3) mapped coral reads to coral (again, for alignment rate counts), 4) symbiont reads to coral, and 5) remaining symbiont reads to symbiont

mkdir symbionts
mkdir junk

# replace 'Symbionts' with your symbiont filename, and 'Host' with your host filename
# maps reads to symbiont reference first, and splits symbiont reads (.sym) to alternate subdirectory
# outputs .host files of remaining reads for coral mapping
tagseq_bowtie2_launcher.py -g ~/db/Symbiont -f .trim -n maps --split -u host -a sym --aldir symbionts --launcher -e studivanms@gmail.com
sbatch maps.slurm

# delete the .sam files since the symbiont reads need further mapping to clean up conserved genes
rm *.sam

# conduct mapping on coral reads (.host files) to coral reference
# outputs .host.sam files for coral gene counts, and .host.clean files for coral alignment rates
tagseq_bowtie2_launcher.py -g ~/db/Host -f .host -n maps2 --split -u un -a clean --undir junk --launcher -e studivanms@gmail.com
sbatch maps2.slurm

# delete the .sam files since you will generate them later from clean coral reads
rm *.host.sam

# conduct second round of mapping on coral reads (.host.clean files) to coral reference, just for alignment rate calculations
tagseq_bowtie2_launcher.py -g ~/db/Host -f .host.clean -n maps3 --launcher -e studivanms@gmail.com
sbatch maps3.slurm

cd symbionts/
mkdir junk

# removing genes from symbiont reads that align to both coral/symbiont references by conducting another round of mapping on .sym files
# outputs .clean files of true symbiont reads for one final rounds of symbiont mapping
tagseq_bowtie2_launcher.py -g ~/db/Host -f .sym -n maps4 --split -u clean -a host --aldir junk --launcher -e studivanms@gmail.com
sbatch maps4.slurm

# delete the .sam files from the zoox mapping to coral reference
rm *.sam

# conduct final round of mapping on true symbiont reads (.sym.clean files) to symbiont reference
# outputs .sym.clean.sam files for symbiont counts
tagseq_bowtie2_launcher.py -g ~/db/Symbiont -f .clean -n maps5 --launcher -e studivanms@gmail.com
sbatch maps5.slurm

mv *.sym.clean.sam ..
mv *.sym.clean ..
cd ..

# double check you have the same number of host and symbiont .sam files as samples
ll *.host.clean.sam | wc -l
ll *.sym.clean.sam | wc -l

# to count the number of mapped host and symbiont reads for mapping efficiency
# calculate mapping efficiency from these values compared to trimmed reads in Excel
# remember to delete the results from the .sam files
echo "countreads_host.pl > countreads_host.txt" > count_align
echo "countreads_sym.pl > countreads_sym.txt" >> count_align
launcher_creator.py -j count_align -n count_align -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count_align.slurm


#------------------------------
## OPTIONAL: For mapping paired-end reads, including to separate host/symbiont transcriptomes with bowtie2, follow the script bowtie2_pairedend_README.txt


#------------------------------
# generating read-counts-per gene: (again, creating a job file to do it simultaneously for all)

# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file and genes
cp ~/annotate/Host_seq2iso.tab ~/db/
cp ~/annotate/Symbiont_seq2iso.tab ~/db/

samcount_launch_bt2.pl '\.host.clean.sam$' /home/mstudiva/db/Host_seq2iso.tab > sc
samcount_launch_bt2.pl '\.sym.clean.sam$' /home/mstudiva/db/Symbiont_seq2iso.tab >> sc
launcher_creator.py -j sc -n sc -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch sc.slurm

# OPTIONAL: For single-reference alignments
samcount_launch_bt2.pl '\.sam$' /home/mstudiva/db/Host_seq2iso.tab > sc
launcher_creator.py -j sc -n sc -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch sc.slurm

# double check you have the same number of files as samples (double if you have split coral/symbiont reads)
ll *.counts | wc -l

# assembling them all into a single table:
srun expression_compiler.pl *.host.clean.sam.counts > allc_host.txt
srun expression_compiler.pl *.sym.clean.sam.counts > allc_sym.txt

# OPTIONAL: For single-reference alignments
srun expression_compiler.pl *.sam.counts > allc_host.txt

# how do the files look?
head allc_host.txt
head allc_sym.txt

# let's remove those annoying chains of extensions from sample names:
cat allc_host.txt | perl -pe 's/\.trim\.host\.clean\.sam\.counts//g'>allcounts_host.txt
cat allc_sym.txt | perl -pe 's/\.trim\.sym\.clean\.sam\.counts//g' >allcounts_sym.txt

# OPTIONAL: For single-reference alignments
cat allc_host.txt | perl -pe 's/\.sam\.counts//g'>allcounts_host.txt

# OPTIONAL: For paired-end alignments
cat allc_host.txt | perl -pe 's/\.fastq\.host\.clean\.sam\.counts//g'>allcounts_host.txt
cat allc_sym.txt | perl -pe 's/\.fastq\.sym\.clean\.sam\.counts//g' >allcounts_sym.txt

head allcounts_host.txt
head allcounts_sym.txt

# display full path to where you were doing all this:
pwd
# copy the path

#------------------------------
# open new terminal window on Mac, or WinSCP on Windows
# navigate to the directory you want the file to be
# (on Mac, you can just drag the destination folder from Finder into command line to get the full path).

# this copies all .txt files, including allcounts, allc, alignrate, read counts
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/*.txt .

# copy the file from KoKo using scp (in WinSCP, just paste the path you just copied into an appropriate slot (should be self-evident) and drag the allcounts.txt file to your local directory)

# DONE! Next, we will be using R to make sense of the counts...
# BUT FIRST, read tagSeq_analysis_README.txt for the next steps
