## Tag-based RNA-seq (Tag-Seq) reads processing pipeline, version August 16, 2023
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com) for use on the FAU KoKo HPC

## BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# log onto cluster
ssh mstudiva@koko-login.hpc.fau.edu


#------------------------------
## Installing RNA-seq scripts and setting up the workspace

# all scripts will go in a bin directory
mkdir bin
cd bin

# clone github repository
git clone https://github.com/mstudiva/tag-based_RNAseq.git

# move files from subdirectory tag-based_RNAseq-master to bin
mv tag-based_RNAseq/* .
rm -rf tag-based_RNAseq
# remove the TACC version of launcher_creator.py from the bin directory (preinstalled on KoKo)
rm launcher_creator.py

# if you have not previously, download BaseSpaceCLI
wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O $HOME/bin/bs
chmod +x ~/bin/bs

# authorization of your BaseSpace account to KoKo
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
echo 'for file in *.gz; do mv "${file}" "${file/-2/}"; done' >> downloadReads.sh
# this removes the '-2' in filenames from duplicate samples for downstream merging

echo "for file in *.fastq.gz; do echo $file | awk -F_ '{ printf("%03d_%s\n", $1, substr($0, index($0, $2))); }' | xargs -I{} mv $file {}; done" >> downloadReads.sh
# this replaces any sample numbers with the three digit version

echo 'mergeReads.sh -o mergeTemp' >> downloadReads.sh
# -o is the directory to put output files in

echo 'rm *L00*' >> downloadReads.sh
echo "find . -name '*.gz' -exec mv {} . \;" >> downloadReads.sh
echo 'gunzip *.gz' >> downloadReads.sh
echo 'rmdir mergeTemp' >> downloadReads.sh

chmod +x downloadReads.sh

launcher_creator.py -b 'srun downloadReads.sh' -n downloadReads -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch --mem=200GB downloadReads.slurm

# double check you have the correct number of files as samples
ll *.fastq | wc -l

# look at the reads
head -50 SampleName.fastq
# note that every read has four lines, the ID line starts with @HWI

# shows only nucleotide sequences in file
head -100 SampleName.fastq | grep -E '^[NACGT]+$'

# to count the number of reads in all samples
echo "countreads.pl > countreads_raw.txt" > count
launcher_creator.py -j count -n count -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count.slurm


#------------------------------
## Creating conda environments for specialized modules

# uncomment and run below if you don't have conda set up
# module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

conda create -n cutadaptenv cutadapt
# this is specifically for cutadapt, which doesn't play well with other KoKo modules


#------------------------------
## Removing adaptors and low quality reads

echo '#!/bin/bash' > trim.sh
echo 'conda activate cutadaptenv' >> trim.sh
for F in *.fastq; do
echo "tagseq_clipper.pl $F | cutadapt - -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o ${F/.fastq/}.trim" >>trim.sh;
done

# does not work with launcher_creator, consider breaking up script and running multiple jobs
sbatch -o trim.o%j -e trim.e%j --mem=200GB trim.sh

# checking the status of the job
squeue -u mstudiva

# double check you have the same number of files as samples
ll *.trim | wc -l

# did the trimming work?
head -100 SampleName.fastq | grep -E '^[NACGT]+$'
head -100 SampleName.trim | grep -E '^[NACGT]+$'
# the long runs of base A should be gone

# double-check that the rnaseq_clipper did not filter out too many reads by looking at the trim.e####### file
# make sure you're not losing too many reads to duplicates
# rename as a txt file
mv trim.e####### trim.txt

# to save time in case of issues, move the concatenated fastq files to backup directory
mv *.fastq ~/rawReads/

# to count the number of reads in trimmed samples
echo "countreads_trim.pl > countreads_trim.txt" > count_trim
launcher_creator.py -j count_trim -n count_trim -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count_trim.slurm


#------------------------------
## Download and format reference transcriptome

mkdir ~/db/
cd ~/db/
# copy your transcriptome fasta file(s) to db/

# creating bowtie2 index for multiple transcriptomes
module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5
# replace 'Host' and 'Symbiont' with your respective filenames
echo 'bowtie2-build Host.fasta,Symbiont.fasta,Symbiont2.fasta Host_concat' > btb
launcher_creator.py -j btb -n btb -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch btb.slurm


#------------------------------
## Mapping reads to host/symbiont transcriptomes with bowtie2

module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5

# if there are multiple available reference transcriptomes for your species, create some test directories, copy over a handful of trim files, and test out alignment to each reference
cp {001..009}.trim # to copy a numbered sequence of filenames
tagseq_bowtie2_launcher.py -g ~/db/Host -f .trim -n mapstest --launcher -e studivanms@gmail.com
sbatch mapstest.slurm
# look at each maps.e####### file for overall mapping efficiency and % of >1 matches, then compare among transcriptomes to find the best alignment rates

# now running the full set of files
mkdir unaligned

# bowtie will separate out aligned and unaligned reads for counts below
# it will also pick the best alignment from the concatenated transcriptomes
tagseq_bowtie2_launcher_concat.py -g ~/db/Host_concat/Host_concat -f .trim -n maps --split -u un -a al --undir unaligned --launcher -e studivanms@gmail.com
sbatch maps.slurm

## OPTIONAL: For mapping paired-end reads, including to separate host/symbiont transcriptomes with bowtie2, follow the script bowtie2_pairedend_README.txt


#------------------------------
## Mapping efficiency

# double check you have the same number of host and symbiont .sam files as samples
ll *.sam | wc -l

# to count the number of mapped reads for mapping efficiency
# calculate mapping efficiency from these values compared to trimmed reads in Excel
echo "countreads_align.pl > countreads_align.txt" > count_align
launcher_creator.py -j count_align -n count_align -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count_align.slurm


#------------------------------
## Generating read counts per gene

# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file and genes
cp ~/annotate/Host_concat_seq2iso.tab ~/db/
# if working with multiple reference transcriptomes, concatenate the seq2iso tables
cat Host_seq2iso.tab Sym_seq2iso.tab > Host_concat_seq2iso.tab

module load samtools-1.10-gcc-8.3.0-khgksad
samcount_launch_bt2.pl '\.sam$' /home/mstudiva/db/Host_concat_seq2iso.tab > sc
launcher_creator.py -j sc -n sc -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch sc.slurm

# double check you have the same number of files as samples
ll *.counts | wc -l

# assembling them all into a single table:
srun expression_compiler.pl *.counts > allc.txt

# how do the counts look?
head allc.txt

# let's remove those annoying chains of extensions from sample names
cat allc.txt | perl -pe 's/\.trim\.sam\.counts//g'>allcounts.txt

# OPTIONAL: for paired-end alignments
cat allc.txt | perl -pe 's/\.fastq\.sam\.counts//g'>allcounts.txt

head allcounts.txt


#------------------------------
## Downloading gene count files

# open new terminal window on Mac, or WinSCP on Windows, and navigate to desired directory
cd /path/to/local/directory

# this copies all .txt files, including allcounts, allc, alignrate, read counts
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*.txt .

# DONE! Next, we will be using R to make sense of the counts...
# BUT FIRST, read tagSeq_analysis_README.txt for the next steps
