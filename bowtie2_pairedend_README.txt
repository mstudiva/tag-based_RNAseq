# Mapping paired-end reads to transcriptomes with bowtie2, version January 22, 2022
# Also includes sequential mapping to host/symbiont transcriptomes
# Created by Michael Studivan (studivanms@gmail.com)
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
## For single-reference transcriptomes (aposymbiotic species)

# First need to create a template script containing the correct sequence of left, right, and output names (along with standard flags)
# ex. "bowtie2 --local -x /mnt/beegfs/home/mstudiva/db/Host -1 R1_Host_1.fastq -2 R2_Host_1.fastq -S Host_1.fastq.sam --no-hd --no-sq --no-unal"
launcher_creator.py -j bowtie.sh -n bowtie -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch bowtie.slurm

# Easy peasy!


#------------------------------
## Mapping paired-end reads to host/symbiont transcriptomes with bowtie2

# Buckle up.

mkdir symbionts
mkdir junk

# replace 'Symbiont' with your symbiont filename, and 'Host' with your host filename
# You need to create template scripts for each of the bowtie steps containing the correct sequence of left, right, and output names (along with standard flags)

# maps reads to symbiont reference first, and splits symbiont reads (.sym) to alternate subdirectory
# outputs .host files of remaining reads for host mapping
# ex. "bowtie2 --local -x /mnt/beegfs/home/mstudiva/db/Symbiont -1 R1_Host_1.fastq -2 R2_Host_1.fastq -S Host_1.fastq.sam --no-hd --no-sq --no-unal --al-conc symbionts/Host_1.fastq.sym --un-conc ./Host_1.fastq.host"
launcher_creator.py -j bowtie.sh -n maps -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch maps.slurm

# delete the .sam files since the symbiont reads need further mapping to clean up conserved genes
rm *.sam

# conduct mapping on host reads (.host files) to host reference
# outputs .host.sam files for host gene counts, and .host.clean files for host alignment rates
# ex. "bowtie2 --local -x /mnt/beegfs/home/mstudiva/db/Host -1 Host_1.fastq.1.host -2 Host_1.fastq.2.host -S Host_1.fastq.host.sam --no-hd --no-sq --no-unal --al-conc ./Host_1.fastq.host.clean --un-conc junk/Host_1.fastq.host.un"
launcher_creator.py -j bowtie2.sh -n maps2 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch maps2.slurm

# delete the .sam files since you will generate them later from clean host reads
rm *.host.sam

# conduct second round of mapping on host reads (.host.clean files) to host reference, just for alignment rate calculations
# ex. "bowtie2 --local -x /mnt/beegfs/home/mstudiva/db/Host -1 Host_1.fastq.1.host.clean -2 Host_1.fastq.2.host.clean -S Host_1.fastq.host.clean.sam --no-hd --no-sq --no-unal"
launcher_creator.py -j bowtie3.sh -n maps3 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch maps3.slurm

cd symbionts/
mkdir junk

# removing genes from symbiont reads that align to both host/symbiont references by conducting another round of mapping on .sym files
# outputs .clean files of true symbiont reads for one final rounds of symbiont mapping
# ex. "bowtie2 --local -x /mnt/beegfs/home/mstudiva/db/Host -1 Host_1.fastq.1.sym -2 Host_1.fastq.2.sym -S Host_1.fastq.sym.sam --no-hd --no-sq --no-unal --al-conc ./Host_1.fastq.sym.clean --un-conc junk/Host_1.fastq.sym.host"
launcher_creator.py -j bowtie4.sh -n maps4 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch maps4.slurm

# delete the .sam files from the zoox mapping to host reference
rm *.sam

# conduct final round of mapping on true symbiont reads (.sym.clean files) to symbiont reference
# outputs .sym.clean.sam files for symbiont counts
# ex. "bowtie2 --local -x /mnt/beegfs/home/mstudiva/db/Symbiont -1 Host_1.fastq.1.sym.clean -2 Host_1.fastq.2.sym.clean -S Host_1.fastq.sym.clean.sam --no-hd --no-sq --no-unal"
launcher_creator.py -j bowtie5.sh -n maps5 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
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
echo "countreads_host.pl '.host.1.clean' > countreads_host.txt" > count_align
echo "countreads_host.pl '.host.2.clean' >> countreads_host.txt" >> count_align
echo "countreads_sym.pl '.sym.1.clean' > countreads_sym.txt" >> count_align
echo "countreads_sym.pl '.sym.2.clean' >> countreads_sym.txt" >> count_align
launcher_creator.py -j count_align -n count_align -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count_align.slurm
