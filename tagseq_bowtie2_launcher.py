#!/usr/bin/env python3
################################################################################
#Script Name	: tagseq_bowtie_launcher.py
#Description	: Aligns trimmed reads to reference transcriptome, and separates symbiont reads from coral
#Author       	: Ryan J Eckert, modified by Michael S Studivan
#Email         	: ryan.j.eckert@gmail.com, studivanms@gmail.com

#Modified from '2bRAD_bowtie_launcher.py' by Ryan Eckert
#(https://github.com/RyanEckert/Stephanocoenia_FKNMS_PopGen)
################################################################################
import os, argparse, sys, glob

parser = argparse.ArgumentParser(description='Prints a list of bowtie2 commands to be used with launcher_creator.py')
parser.add_argument("-f", type=str, help="File type to add, don't include . [e.g. fastq]", action = "store", default='trim')
parser.add_argument("-g", type=str, help="Path to reference genome", action = "store", default='')
parser.add_argument("--keep_unal", help="When used unaligned reads will be included in the sam alignment files", action = "store_true")
parser.add_argument("--split", help="If used, split reads into unaligned/aligned files", action = "store_true")
parser.add_argument("-a", type=str, help="Optional appended name for aligned reads", action = "store", default='aligned')
parser.add_argument("-u", type=str, help="Optional appended name for unaligned reads", action = "store", default='unaligned')
parser.add_argument("--undir", type=str, help="Optional directory for unaligned reads", action = "store", default='.')
parser.add_argument("--aldir", type=str, help="Optional directory for unaligned reads", action = "store", default='.')
parser.add_argument("-n", type=str, help="Name of output file to write launcher commands to", action = "store", default='maps')
parser.add_argument("--launcher", help="If used, create .slurm file with launcher_creator.py", action = "store_true")
parser.add_argument("-e", help="If used, email to send SLURM info to", action = "store")

args = parser.parse_args()
split = args.split
launcher = args.launcher
e = args.e
f = args.f
g = args.g
n = args.n
a = args.a
u = args.u
ad = args.aldir
ud = args.undir
keep = args.keep_unal

outfile = open(n, 'w')

for file in glob.iglob(r'*'+f):
    list = 'bowtie2 --local -x '+g+' -U '+file+ ' -S '+file+'.sam'
    if split:
        if keep: outfile.write(list+' --al '+ad+'/'+file+'.'+a+' --un '+ud+'/'+file+'.'+u+'\n')
        else: outfile.write(list+' --no-hd --no-sq --no-unal --al '+ad+'/'+file+'.'+a+' --un '+ud+'/'+file+'.'+u+'\n')
    else:
        if keep: outfile.write(list+'\n')
        else: outfile.write(list+' --no-hd --no-sq --no-unal'+'\n')

outfile.close()

if launcher:
    os.system("module load launcher/3.5");
    if e:
        os.system('launcher_creator.py -j '+n+' -n '+n+' -q shortq7 -t 6:00:00 -e '+e)
    else:
        os.system('launcher_creator.py -j '+n+' -n '+n+' -q shortq7 -t 6:00:00')
