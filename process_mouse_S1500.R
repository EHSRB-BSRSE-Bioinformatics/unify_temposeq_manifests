# script to fill in blank IDs in TempO-seq manifests (as much as possible) and unify columns across probe designs

# Script written by Lauren but mostly stolen from Kate

# final design will be: 
# Probe ID, Gene Symbol, Gene Ensembl ID, Entrez ID, Ligated sequence, Aligned transcript IDs

####################
# Mouse S1500+ 1.2 #
####################

# Using the updated manifest available at https://ntp.niehs.nih.gov/ntp/research/areas/tox21/mouse_s1500_manifest_10182018.xlsx
# Which have Entrez IDs for each probe

# The main issue with this set is missing Ensembl IDs

###############
#Load packages#
###############

library(tidyverse)
library(biomaRt)
library(annotate)
library(seqinr)
library(rentrez)
library(stringr)

#################
#Getting started#
#################

# read in manifest file
original_manifests_directory <- "original_manifests"

mouse_s1500_1.2_file <- file.path(original_manifests_directory, "mouse_s1500_manifest_10182018.csv")

manifest <- read.csv(mouse_s1500_1.2_file, sep=',')

# First, do some data conversion

manifest[manifest == 'NULL'] <- NA

# Rename and select columns according to standardized human S1500 CSVs
manifest <- manifest %>%
  dplyr::rename("Probe_Name" = "prelim_mouse_name",
                "Gene_Symbol" = "prelim_mouse_gene",
                "Probe_ID" = "prelim_mouse_id",
                "Probe_Sequence" = "prelim_mouse_ligseq",
                "Entrez_ID" = "Entrez_Gene_ID..Sciome.generated.",
                "Transcripts_Targeted" = "prelim_mouse_targeted_enst") %>%
  dplyr::select("Probe_ID","Probe_Name","Gene_Symbol","Entrez_ID","Probe_Sequence","Transcripts_Targeted")

# Set up output file
output_directory <- "output_manifests"
output_filename <- "Mouse_S1500_1-2_standardized.csv"

#####################################################################################################################

# Get Ensembl IDs, since all rows have Entrez IDs

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

ensembl_for_missing <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'),
                             filters = 'entrezgene_id',
                             values = manifest$Entrez_ID,
                             mart = ensembl)

# Find rows where manifest gene symbol and biomart gene name match
# Mutate to title case, sometimes mismatch is just case sensitivity
entrez_matched <- manifest %>%
  inner_join(ensembl_for_missing, by=c("Entrez_ID" = "entrezgene_id")) %>%
  mutate(Gene_Symbol = str_to_title(Gene_Symbol),
         external_gene_name =  str_to_title(external_gene_name)) %>%
  filter(Gene_Symbol == external_gene_name)

# Find rows where manifest gene symbol and biomart gene name don't match
symbol_mismatch <- manifest %>%
  inner_join(ensembl_for_missing, by=c("Entrez_ID" = "entrezgene_id")) %>%
  mutate(Gene_Symbol = str_to_title(Gene_Symbol),
         external_gene_name =  str_to_title(external_gene_name)) %>%
  filter(Gene_Symbol != external_gene_name)

nrow(symbol_mismatch)
# 211 rows with mismatch

check_n_matches <- entrez_matched %>% group_by(Probe_ID) %>% mutate(n_copies = n())
check_n_matches %>% filter(n_copies > 1)
