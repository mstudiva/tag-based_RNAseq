# Tag-based RNA-seq (Tag-Seq) reads processing pipeline, version January 18, 2023
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com) for use on the FAU KoKo HPC


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

# goes to the website to confirm authorization by logging in to your BaseSpace acct.
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

# look at the reads:
# head -50 SampleName.fastq
# note that every read has four lines, the ID line starts with @HWI

# this little one-liner will show sequence-only in file:
# head -100 SampleName.fastq | grep -E '^[NACGT]+$'

# to count the number of reads in all samples
echo "countreads.pl > countreads_raw.txt" > count
launcher_creator.py -j count -n count -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count.slurm

#------------------------------
# Create conda environment
# Uncomment and run below if you don't have a conda env. set up
# module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

conda create -n cutadapt cutadapt

# Removing adaptors and low quality reads
echo '#!/bin/bash' > trim.sh
echo 'conda activate condaenv' >> trim.sh
for F in *.fastq; do
echo "tagseq_clipper.pl $F | cutadapt - -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o ${F/.fastq/}.trim" >>trim.sh;
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

# head -100 SampleName.fastq | grep -E '^[NACGT]+$'

# head -100 SampleName.trim | grep -E '^[NACGT]+$'
# the long runs of base A should be gone

# double-check that the rnaseq_clipper did not filter out too many reads by looking at the trim.e####### file
# make sure you're not losing too many reads to duplicates
# rename as a txt file
mv trim.e####### trim.txt

# to save time in case of issues, move the concatenated fq files to backup directory
mv *.fastq ~/rawReads/

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
module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5
# replace 'Host' and 'Symbiont' with your respective filenames
echo 'bowtie2-build Host.fasta Host' > btb
echo 'bowtie2-build Symbiont.fasta Symbiont' >> btb
launcher_creator.py -j btb -n btb -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch btb.slurm


#------------------------------
## Mapping reads to host/symbiont transcriptomes with bowtie2

# if there are multiple available reference transcriptomes for your species, create some test directories, copy over a handful of trim files, and test out alignment to each reference
cp {001..009}.trim # to copy a numbered sequence of filenames
tagseq_bowtie2_launcher.py -g ~/db/Host -f .trim -n mapstest --launcher -e studivanms@gmail.com
sbatch mapstest.slurm
# look at each maps.e####### file for overall mapping efficiency and % of >1 matches, then compare among transcriptomes to find the best alignment rates

# if working with split host/symbiont transcriptomes, need to run FOUR rounds of mapping: 1) all reads to symbiont, 2) remaining reads to host, 3) symbiont reads to host to remove co-occurring genes, and remaining symbiont reads to symbiont for final .sam files
mkdir symbionts
mkdir junk

# replace 'Symbionts' with your symbiont filename, and 'Host' with your host filename
# maps reads to symbiont reference first, and splits symbiont reads (.sym) to alternate subdirectory
# outputs .host files of remaining reads for host mapping
tagseq_bowtie2_launcher.py -g ~/db/Symbiont -f .trim -n maps --split -u host -a sym --aldir symbionts --launcher -e studivanms@gmail.com
sbatch maps.slurm

# delete the .sam files since the symbiont reads need further mapping to clean up co-occurring genes
rm *.sam

# OPTIONAL: if using two symbiont transcriptomes (multiple genera), run these two rounds of mapping instead before proceeding, one for each symbiont transcriptome
# start with the more abundant symbiont genus
tagseq_bowtie2_launcher.py -g ~/db/Symbiont -f .trim -n maps --split -u sym2 -a sym --aldir symbionts --launcher -e studivanms@gmail.com
sbatch maps.slurm

rm *.sam

tagseq_bowtie2_launcher.py -g ~/db/Symbiont2 -f .sym2 -n maps1b --split -u host -a sym2 --aldir symbionts --launcher -e studivanms@gmail.com
sbatch maps1b.slurm
# END OPTIONAL

# conduct mapping on host reads (.host files) to host reference
# outputs .host.sam files for host gene counts, and .host.clean files for host alignment rates
tagseq_bowtie2_launcher.py -g ~/db/Host -f .host -n maps2 --split -u un -a clean --undir junk --launcher -e studivanms@gmail.com
sbatch maps2.slurm

cd symbionts/
mkdir junk

# removing genes from symbiont reads that align to both host/symbiont references by conducting another round of mapping on .sym files
# outputs sym.clean files for symbiont alignment rates
tagseq_bowtie2_launcher.py -g ~/db/Host -f .sym -n maps3 --split -u clean -a host --aldir junk --launcher -e studivanms@gmail.com
sbatch maps3.slurm

# OPTIONAL
# running on the second symbiont reference files
tagseq_bowtie2_launcher.py -g ~/db/Host -f .sym2 -n maps3b --split -u clean -a host --aldir junk --launcher -e studivanms@gmail.com
sbatch maps3b.slurm
# END OPTIONAL

# delete the .sam files from the symbiont mapping to host
rm *.sam

# conduct final round of mapping on true symbiont reads (.sym.clean files) to symbiont reference
# outputs .sym.clean.sam files for symbiont gene counts
tagseq_bowtie2_launcher.py -g ~/db/Symbiont -f .sym.clean -n maps4 --launcher -e studivanms@gmail.com
sbatch maps4.slurm

mv *.sym.clean.sam ..
mv *.sym.clean ..
cd ..

# OPTIONAL
# running on the second symbiont reference files
tagseq_bowtie2_launcher.py -g ~/db/Symbiont2 -f .sym2.clean -n maps4b --launcher -e studivanms@gmail.com
sbatch maps4b.slurm

mv *.sym2.clean.sam ..
mv *.sym2.clean ..
cd ..
# END OPTIONAL

# double check you have the same number of host and symbiont .sam files as samples
ll *.host.sam | wc -l
ll *.sym.clean.sam | wc -l
ll *.sym2.clean.sam | wc -l # OPTIONAL

# to count the number of mapped host and symbiont reads for mapping efficiency
# calculate mapping efficiency from these values compared to trimmed reads in Excel
# remember to delete the results from the .sam files
echo "countreads_host.pl > countreads_host.txt" > count_align
echo "countreads_sym.pl > countreads_sym.txt" >> count_align
echo "countreads_sym.pl .sym2.clean > countreads_sym2.txt" >> count_align # OPTIONAL
launcher_creator.py -j count_align -n count_align -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count_align.slurm


#------------------------------
## OPTIONAL: For mapping paired-end reads, including to separate host/symbiont transcriptomes with bowtie2, follow the script bowtie2_pairedend_README.txt
# END OPTIONAL


#------------------------------
# generating read-counts-per gene: (again, creating a job file to do it simultaneously for all)

# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file and genes
cp ~/annotate/Host_seq2iso.tab ~/db/
cp ~/annotate/Symbiont_seq2iso.tab ~/db/

module load samtools-1.10-gcc-8.3.0-khgksad
samcount_launch_bt2.pl '\.host.sam$' /home/mstudiva/db/Host_seq2iso.tab > sc
samcount_launch_bt2.pl '\.sym.clean.sam$' /home/mstudiva/db/Symbiont_seq2iso.tab >> sc
launcher_creator.py -j sc -n sc -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch sc.slurm

# OPTIONAL: For single-reference alignments
samcount_launch_bt2.pl '\.sam$' /home/mstudiva/db/Host_seq2iso.tab > sc
launcher_creator.py -j sc -n sc -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch sc.slurm
# END OPTIONAL

# double check you have the same number of files as samples (double if you have split host/symbiont reads)
ll *.counts | wc -l

# assembling them all into a single table:
srun expression_compiler.pl *.host.sam.counts > allc_host.txt
srun expression_compiler.pl *.sym.clean.sam.counts > allc_sym.txt

# OPTIONAL: For single-reference alignments
srun expression_compiler.pl *.sam.counts > allc_host.txt
# END OPTIONAL

# how do the files look?
head allc_host.txt
head allc_sym.txt

# let's remove those annoying chains of extensions from sample names:
cat allc_host.txt | perl -pe 's/\.trim\.host\.sam\.counts//g'>allcounts_host.txt
cat allc_sym.txt | perl -pe 's/\.trim\.sym\.clean\.sam\.counts//g' >allcounts_sym.txt

# OPTIONAL: For single-reference alignments
cat allc_host.txt | perl -pe 's/\.sam\.counts//g'>allcounts_host.txt

# OPTIONAL: For paired-end alignments
cat allc_host.txt | perl -pe 's/\.fastq\.host\.clean\.sam\.counts//g'>allcounts_host.txt
cat allc_sym.txt | perl -pe 's/\.fastq\.sym\.clean\.sam\.counts//g' >allcounts_sym.txt
# END OPTIONAL

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
