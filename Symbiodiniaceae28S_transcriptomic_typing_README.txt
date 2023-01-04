# Quick n' dirty Symbiodiniaceae typing from transcriptomic datasets, version September 11, 2020
# Created by Michael Studivan (studivanms@gmail.com)
# for use on the FAU KoKo HPC

# log onto cluster
ssh username@koko-login.hpc.fau.edu
# enter FAU password

#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- username with your KoKo user name

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions – please make sure to read them before copy-pasting.

#------------------------------
# initial setup and file transfers
# create a directory in the desired location and temporarily move all your .fastq files there
mkdir zoox
cd zoox
mv ~/path/to/directory/*.fastq .

# from your local machine, scp the combined Symbiodiniaceae 28S genes fasta to the HPC
scp 28S.fasta username@koko-login.hpc.fau.edu:~/path/to/directory/.

#------------------------------
# then login and navigate to the directory
# build a bowtie2 index of Symbiodiniaceae 28S genes
echo 'bowtie2-build 28S.fasta Symbiodiniaceae_28S' > btb
launcher_creator.py -j btb -n btb -q shortq7 -t 6:00:00
sbatch btb.slurm

#------------------------------
# lists all the .fastq files in the directory, then applies the bowtie2 and samtools commands to align and count 28S sequences to sample transcripts
ls *fastq | perl -pe 's/(\S+)/bowtie2 -U $1 -x Symbiodiniaceae_28S --threads 20 -q --score-min L,0,0 --very-sensitive --end-to-end \| samtools view -Sb -o $1.bam/' >align
launcher_creator.py -j align -n align -t 6:00:00 -q shortq7
sbatch align.slurm

# check on the status of the job
squeue -u username

# double check you have the same number of files as samples
ll *.bam | wc -l

# once finished, rename the job error output containing alignment rates to a .txt file
mv align.e####### alignrate.txt
nano alignrate.txt

# move the .fq files back to the original directory
mv *.fq ~/path/to/directory/

# lists all the .bam files, then applies the sort command
ls *bam | perl -pe 's/(\S+)/samtools sort -o $1.sort $1/' >sort
launcher_creator.py -j sort -n sort -t 6:00:00 -q shortq7
sbatch sort.slurm

# check on the status of the job
squeue -u username

# double check you have the same number of files as samples
ll *.sort | wc -l

# once done, lists all .sort files for index and count commands
ls *.sort | perl -pe 's/(\S+)/samtools index $1 \| samtools idxstats $1 > $1.counts/' >count
launcher_creator.py -j count -n count -t 6:00:00 -q shortq7
sbatch count.slurm

# check on the status of the job
squeue -u username

# double check you have the same number of files as samples
ll *.counts | wc -l

# combines all results into a single horizontal .txt
paste *.counts > counts.txt

#------------------------------
# scp .txt files (align rate and counts) to local machine, then build spreadsheet with all samples as columns
scp username@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*.txt .
