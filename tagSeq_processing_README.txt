# Tag-based RNA-seq reads processing pipeline, version September 11, 2020
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com)
# for use on the FAU KoKo HPC

# log onto cluster
ssh mstudiva@koko-login.fau.edu
# enter FAU password

#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions – please make sure to read them before copy-pasting.

#------------------------------
# installing RNA-seq scripts and setting up the workspace

# switch to home directory
cd

# unless you have done it in the past, make directory called bin,
# all your scripts should go in there:
mkdir bin
mkdir backup
mkdir tagseq
mkdir db

# switch to bin:
cd bin

# clone github repository
git clone https://github.com/z0on/tag-based_RNAseq.git

# move files from subdir tag-based_RNAseq-master to bin/:
mv tag-based_RNAseq/* .

# remove the tag-based_RNAseq directory
rm -rf tag-based_RNAseq

# remove the TACC version of launcher_creator.py from the bin directory
rm launcher_creator.py

# clone github repository with modified scripts for use with M. cavernosa/Cladocopium spp. and Eli Meyer's library prep
git clone https://github.com/mstudiva/Transcriptional-plasticity-mesophotic-Mcav.git

# move files from subdirectory and delete empty subdirectory
mv tagseq-modified-scripts/* .
rm -rf tagseq-modified-scripts

# to use Google Drive for file storage, you have to install gdrive
# download the file from the following link to your personal computer
# https://docs.google.com/uc?id=0B3X9GlR6EmbnQ0FtZmJJUXEyRTA&export=download
# scp the file to your bin directory
scp gdrive-linux-x64 mstudiva@koko-login.fau.edu:~/bin/

mv gdrive-linux-x64 gdrive
# make the file executable
chmod +x gdrive
srun gdrive about

# copy the resulting link and log in to your Google Drive account
# allow gdrive access to your account
# copy the resulting token code and paste it into your ssh client window

#------------------------------
# setting up working environment:

# switch to home directory
cd

# start nano editor on the file that will hold your environmental settings
nano .bashrc

# 11 Sept 2020 Update: Many of these modules are now deprecated. Use 'module avail' on KoKo to determine the most up-to-date version to include in your .bashrc
# paste these lines in SECTION 1:
module load gcc
module load slurm
module load openmpi/gcc
module load launcher
module load fastx_toolkit
module load bowtie2/2.3.0
module load blast

# paste this line in SECTION 2:
export PATH="$HOME/bin/:$PATH"

# press these keys to save it and exit (note what happens on the screen):
	ctrl-O - enter
	ctrl-X - enter

# make the environment changes take effect (will happen automatically next time you log in):
source .bashrc

#------------------------------
# let's check if the environment is set correctly: the programs should run from any location now

# switch to home directory
cd

# try a few commands: you should see messages about their usage and/or error messages
# about missing arguments, but not "command not found"
samcount.pl
launcher_creator.py

#------------------------------
# downloading sequence data:

# once each run .zip file is downloaded, upload to a web hosting service like Dropbox or Google Drive
# FAU has free Google Drive storage up to 10TB

cd backup

# Google Drive instructions
# first, use the following line to display all .zip files in your Google Drive account
gdrive list --query "name contains '.zip'"

# copy the Ids for each of your .zip run files and replace below
echo 'gdrive download 0B40lh4qYLpPVYmd6RUF5emVaRWc' > wget
echo 'gdrive download 0B40lh4qYLpPVRUc0cmhsTjdJSk0' >> wget
echo 'gdrive download 0B40lh4qYLpPVcm00RVp3SlVwT0E' >> wget
echo 'gdrive download 0B40lh4qYLpPVVm5OdHdSMTdBX0U' >> wget

launcher_creator.py -j wget -n wget -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch wget.slurm

echo 'unzip ACB9MJANXX.zip' > unzip
echo 'unzip ACB32LANXX.zip' >> unzip
echo 'unzip BCB5VAANXX.zip' >> unzip
echo 'unzip BCB9MKANXX.zip' >> unzip

launcher_creator.py -j unzip -n unzip -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch unzip.slurm

# move any fastq files from subdirectories, then delete the subdirectories
mv BCB9MKANXX/* .
rm -rf BCB9MKANXX

#-------------------------------
# unzipping and concatenating sequence files

# double check you have the number of files you should
ll *.fastq.gz | wc -l

mv *.fastq.gz ~/tagseq/
cd ~/tagseq/

# creating and launching a cluster job to unzip all files:
ls *.gz | perl -pe 's/(\S+)/gunzip $1/' >gunz
launcher_creator.py -j gunz -n gunz -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch gunz.slurm

# check status of your job (PD : pending, not enough nodes; R : running; nothing printed on the screen - complete)
squeue -u mstudiva

# double check you have the number of files you should
ll *.fastq | wc -l

# If your samples are split across multiple files from different lanes,
# concatenating the corresponding fastq files by sample:
# ngs_concat.pl commonTextInFastqFilenames  "FilenameTextImmediatelyBeforeSampleID(.+)FilenameTextImmediatelyAfterSampleID"

echo "ngs_concat.pl 'MS_' '(.+)_[ATCG]{6}'" > concat
launcher_creator.py -j concat -n concat -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch concat.slurm

# this one-liner replaces any sample numbers with the three digit version
for f in MS_*.fq; do x="${f##*_}"; mv "$f" "${f%_*}$(printf '_%03d.fq' "${x%.fq}")"; done

# double check you have the correct number of files as samples
ll *.fq | wc -l

# look at the reads:
# head -50 SampleName.fq
head -50 MS_009.fq
# note that every read has four lines, the ID line starts with @HWI

# this little one-liner will show sequence-only in file MS_009.fq:
# head -100 SampleName.fq | grep -E '^[NACGT]+$'
head -100 MS_009.fq | grep -E '^[NACGT]+$'

# move all the .fastq files to the backup directory
mv *.fastq ~/backup/

# to count the number of reads in all samples
echo "countreads_raw.pl > countreads_raw.txt" > count
launcher_creator.py -j count -n count -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count.slurm

#------------------------------
# adaptor and quality trimming:

# creating and launching the cleaning process for all files in the same time:
ls *fq | perl -pe 's/(\S+)/rnaseq_clipper_MS\.pl $1 \| fastx_clipper -a AAAAAAAA -l 20 -Q33 \| fastx_clipper -a AGATCGGAAG -l 20 -Q33 \| fastq_quality_filter -q 20 -p 90 -Q33 >$1.trim/' >clean
launcher_creator.py -j clean -n clean -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch clean.slurm

# how the job is doing?
squeue -u mstudiva

# It is complete! I got a bunch of .trim files that are non-empty!
ll
# double check you have the same number of files as samples
ll *.fq.trim | wc -l

# but did the trimming really work?
# Use the same one-liner as before on the trimmed file to see if it is different
# from the raw one that you looked at before:

# head -100 SampleName.fq | grep -E '^[NACGT]+$'
head -100 MS_009.fq | grep -E '^[NACGT]+$'

# head -100 SampleName.fq.trim | grep -E '^[NACGT]+$'
head -100 MS_009.fq.trim | grep -E '^[NACGT]+$'
# the long runs of base A should be gone

# double-check that the rnaseq_clipper did not filter out too many reads by looking at the clean.e####### file
nano clean.e2569132
# make sure you're not losing too many reads to duplicates
# rename as a txt file
mv clean.e2569132 clean.txt

# to save time in case of issues, move the concatenated fq files to backup directory
mv *.fq ~/backup/

# to count the number of reads in trimmed samples
echo "countreads_trim.pl > countreads_trim.txt" > count_trim
launcher_creator.py -j count_trim -n count_trim -q shortq7 -t 6:00:00 -e studivanms@gmail.com
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
launcher_creator.py -j btb -n btb -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch btb.slurm

#------------------------------
# mapping reads to the transcriptome with bowtie2

# creating a list of mapping commands, one per reads file:
cd ~/tagseq/

srun iRNAseq_bowtie2map.pl "trim$" ~/db/Mcavernosa_Cladocopium > maps
# If you have a file named count_trim, you need to edit maps to remove a line
nano maps
# find the line with count_trim, then Ctrl+K to cut, then save
launcher_creator.py -j maps -n maps -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch maps.slurm

# how is the job?
squeue -u mstudiva

# complete! I got a bunch of large .sam files.
ll

# double check you have the same number of files as samples
ll *.fq.trim.sam | wc -l

# what is the mapping efficiency? This will find relevant lines in the "job output" file
# that was created while the mapping was running
grep "overall alignment rate" maps.e2571269 > alignrate.txt
nano alignrate.txt

#------------------------------
# almost done! Just two small things left:
# generating read-counts-per gene: (again, creating a job file to do it simultaneously for all)

# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file
# and genes. Typically, each gene is represented by several contigs in the transcriptome.
# Mcavernosa_Cladocopium_seq2iso.tab is in your annotate directory
cp ~/annotate/Mcavernosa_Cladocopium_seq2iso.tab ~/db/

samcount_launch_bt2.pl '\.sam$' /home/mstudiva/db/Mcavernosa_Cladocopium_seq2iso.tab > sc
launcher_creator.py -j sc -n sc -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch sc.slurm

# check on the job
squeue -u mstudiva

# done! a bunch of .counts files were produced.
ll

# double check you have the same number of files as samples
ll *.fq.trim.sam.counts | wc -l

# assembling them all into a single table:
srun expression_compiler.pl *.sam.counts > allc.txt

# how does the allc.txt look?
head allc.txt

# let's remove those annoying chains of extensions from sample names:
cat allc.txt | perl -pe 's/\.fq\.trim\.sam\.counts//g' >allcounts.txt
head allcounts.txt

# display full path to where you were doing all this:
pwd
# copy the path!

#------------------------------
# whew. Now just need to copy the result to your laptop!

# open new terminal window on Mac, or WinSCP on Windows
# navigate to the directory you want the file to be
# (on Mac, you can just drag the destination folder from Finder into command line to get the full path).

# this copies all .txt files, including allcounts, allc, alignrate, read counts
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/*.txt .

# copy the file from KoKo using scp (in WinSCP, just paste the path you just copied into an appropriate slot (should be self-evident) and drag the allcounts.txt file to your local directory)

# DONE! Next, we will be using R to make sense of the counts...
# BUT FIRST, read tagSeq_analysis_README.txt for the next steps
