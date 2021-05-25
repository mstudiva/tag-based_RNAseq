username# Tag-based RNA-seq reads processing pipeline, version May 10, 2021
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com)
# for use on the FAU KoKo HPC

#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- email@gmail.com by your actual email;
#	- username with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions – please make sure to read them before copy-pasting.

# log onto cluster
ssh username@koko-login.hpc.fau.edu
# enter FAU password and accept the Duo Mobile prompt on your phone

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
git clone https://github.com/z0on/tag-based_RNAseq.git

# move files from subdir tag-based_RNAseq-master to bin/:
mv tag-based_RNAseq/* .

# remove the tag-based_RNAseq directory
rm -rf tag-based_RNAseq

# remove the TACC version of launcher_creator.py from the bin directory (preinstalled on KoKo)
rm launcher_creator.py

# clone github repository with modified scripts for use with M. cavernosa/Cladocopium spp. and Eli Meyer's library prep
git clone https://github.com/mstudiva/Transcriptional-plasticity-mesophotic-Mcav.git

# move files from subdirectory and delete empty subdirectory
mv tagseq-modified-scripts/* .
rm -rf tagseq-modified-scripts

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

launcher_creator.py -b 'srun downloadReads.sh' -n downloadReads -q shortq7 -t 06:00:00 -e email@gmail.com
sbatch --mem=200GB downloadReads.slurm

#-------------------------------
# Concatenating sequence files

# double check you have the number of files you should
ll *.fastq | wc -l

# If your samples are split across multiple files from different lanes,
# concatenating the corresponding fastq files by sample:
# ngs_concat.pl commonTextInFastqFilenames  "FilenameTextImmediatelyBeforeSampleID(.+)FilenameTextImmediatelyAfterSampleID"

echo "ngs_concat.pl 'Text-' '(.+)-\d'" > concat
launcher_creator.py -j concat -n concat -q shortq7 -t 6:00:00 -e email@gmail.com
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
launcher_creator.py -j count -n count -q shortq7 -t 6:00:00 -e email@gmail.com
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
squeue -u username

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
launcher_creator.py -j count_trim -n count_trim -q shortq7 -t 6:00:00 -e email@gmail.com
sbatch count_trim.slurm

#------------------------------
# download and format reference transcriptome:

# there are already-annotated M. cavernosa and Cladocopium transcriptomes on github and Dropbox
mkdir annotate
cd annotate
git clone https://github.com/mstudiva/Mcav-Cladocopium-Annotated-Transcriptome.git
mv Mcav-Annotated-Transcriptome/* .
rm -rf Mcav-Annotated-Transcriptome

wget -O Mcavernosa_Cladocopium.fasta https://www.dropbox.com/s/4s093k2wzimavp8/Mcavernosa_Cladocopium.fasta

cp ~/annotate/Mcavernosa_Cladocopium.fasta ~/db/
cd db

# creating bowtie2 index for your transcriptome:
echo 'bowtie2-build Mcavernosa_Cladocopium.fasta Mcavernosa_Cladocopium' > btb
launcher_creator.py -j btb -n btb -q shortq7 -t 6:00:00 -e email@gmail.com
sbatch btb.slurm

#------------------------------
# mapping reads to the transcriptome with bowtie2

# creating a list of mapping commands, one per reads file:
srun tagseq_bowtie2map.pl "trim$" ~/db/Mcavernosa_Cladocopium > maps
# If you have a file named count_trim, you need to edit maps to remove a line
nano maps
# find the line with count_trim, then Ctrl+K to cut, then save
launcher_creator.py -j maps -n maps -q shortq7 -t 6:00:00 -e email@gmail.com
sbatch maps.slurm

# how is the job?
squeue -u username

# complete! I got a bunch of large .sam files.
ll

# double check you have the same number of files as samples
ll *.trim.sam | wc -l

# what is the mapping efficiency? This will find relevant lines in the "job output" file
# that was created while the mapping was running
grep "overall alignment rate" maps.e####### > alignrate.txt
nano alignrate.txt

#------------------------------
# generating read-counts-per gene: (again, creating a job file to do it simultaneously for all)

# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file
# and genes. Typically, each gene is represented by several contigs in the transcriptome.
# Mcavernosa_Cladocopium_seq2iso.tab is in your annotate directory
cp ~/annotate/Mcavernosa_Cladocopium_seq2iso.tab ~/db/

samcount_launch_bt2.pl '\.sam$' /home/username/db/Mcavernosa_Cladocopium_seq2iso.tab > sc
launcher_creator.py -j sc -n sc -q shortq7 -t 6:00:00 -e email@gmail.com
sbatch sc.slurm

# check on the job
squeue -u username

# done! a bunch of .counts files were produced.
ll

# double check you have the same number of files as samples
ll *.trim.sam.counts | wc -l

# assembling them all into a single table:
srun expression_compiler.pl *.sam.counts > allc.txt

# how does the allc.txt look?
head allc.txt

# let's remove those annoying chains of extensions from sample names:
cat allc.txt | perl -pe 's/\.trim\.sam\.counts//g' >allcounts.txt
head allcounts.txt

# display full path to where you were doing all this:
pwd
# copy the path

#------------------------------
# open new terminal window on Mac, or WinSCP on Windows
# navigate to the directory you want the file to be
# (on Mac, you can just drag the destination folder from Finder into command line to get the full path).

# this copies all .txt files, including allcounts, allc, alignrate, read counts
cd /path/to/local/directory
scp username@koko-login.fau.edu:~/path/to/HPC/directory/*.txt .

# copy the file from KoKo using scp (in WinSCP, just paste the path you just copied into an appropriate slot (should be self-evident) and drag the allcounts.txt file to your local directory)

# DONE! Next, we will be using R to make sense of the counts...
# BUT FIRST, read tagSeq_analysis_README.txt for the next steps
