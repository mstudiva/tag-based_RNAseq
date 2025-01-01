#### PACKAGES ####

library(tidyverse)
library(VennDiagram)
library(pheatmap)
library(reshape2)
library(RColorBrewer)


#### ORTHOFINDER ####

# Install orthofinder on your local machine using the tutorials (https://davidemms.github.io/menu/tutorials.html)

# Copy your translated protein fasta files (_out_PRO.fas) that you want to compare into a directory called 'orthofinder'
# If you have not already filtered by the longest contig per isogroup (by using fasta2SBH.pl during transcriptome annotation), follow step 7 of tutorial 2 above

# Run the following command in Terminal: 'orthofinder -f orthofinder/'
# Check the number of genes assigned to orthogroups (e.g., 'OrthoFinder assigned 33019 genes (90.8% of total) to 10278 orthogroups')
# Ideally, it should be >80%


#### ORTHOLOGS ####

# load("orthofinder_DEGs.RData") # if previously run
orthologs <- read.table(file = "orthofinder/OrthoFinder/Results_Dec31/Orthologues/Orthologues_Durusdinium_out_PRO/Durusdinium_out_PRO__v__Breviolum_out_PRO.tsv", sep = "\t", header = TRUE, quote="", fill=FALSE)

orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  mutate(Protein_symD = trimws(Protein_symD)) %>%
  mutate(Protein_symB = trimws(Protein_symB)) %>%
  unique() -> orthologs_unique
write.table(orthologs_unique, file = "orthologs_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)