# script to fill in blank IDs in TempO-seq manifests (as much as possible) and unify columns across probe designs

# Script written by Lauren but mostly stolen from Kate

# final design will be: 
# Probe ID, Gene Symbol, Gene Ensembl ID, Entrez ID, Ligated sequence, Aligned transcript IDs

####################
# Mouse S1500+ 1.2 #
####################

# Using the updated manifest available at https://ntp.niehs.nih.gov/ntp/research/areas/tox21/mouse_s1500_manifest_10182018.xlsx
# Which have Entrez IDs for most(!) probes (not quite all)

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
library(XML)
library(BiocManager)
# Install if needed BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

#################
#Getting started#
#################

# read in manifest file
original_manifests_directory <- "original_manifests"

mouse_s1500_1.2_file <- file.path(original_manifests_directory, "mouse_s1500_manifest_10182018.csv")

manifest <- read.csv(mouse_s1500_1.2_file, sep=',')

# Set up output file
output_directory <- "output_manifests"
output_filename <- "Mouse_S1500_1-2_standardized.csv"


# A bit of data conversion

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

# Set up functions (stolen from Kate)
check_data_unique <- function(x){
  n <- x %>% nrow()
  n_unique <- x %>% unique() %>% nrow()
  stopifnot(n==n_unique)
  n_unique_probeid <- x %>% dplyr::select(Probe_ID) %>% unique() %>% nrow()
  stopifnot(n==n_unique_probeid)
}

get_gene_info_from_ncbi_using_entrez <- function(entrez_id){
  Sys.sleep(4) # added to stop the HTTP 300 timeout errors from NCBI
  
  fetch_result <- entrez_fetch(db="gene",id=entrez_id, rettype="xml", parsed=TRUE)
  fetch_list <- xmlToList(fetch_result)
  dbtags <- fetch_list$Entrezgene$Entrezgene_gene[["Gene-ref"]][["Gene-ref_db"]]
  external_array <- as.vector(as.matrix(as.data.frame(dbtags)))
  if(length(external_array)>2){
    external_dbs <- external_array[seq(1, length(external_array), 2)]
    external_ids <- external_array[seq(2, length(external_array), 2)]
    external_data <- data.frame(db=external_dbs, id=external_ids)
    ensembl_id <- external_data %>% filter(db=="Ensembl") %>% dplyr::select(id)
  } else {
    ensembl_id = NA
  }
  entrez_chromosome <- fetch_list$Entrezgene$Entrezgene_source$BioSource$BioSource_subtype$SubSource$SubSource_name
  gene_symbol <- fetch_list$Entrezgene$Entrezgene_gene[["Gene-ref"]][["Gene-ref_locus"]]
  
  gene.df <- data.frame(entrez_id="",ensembl_id="",gene_symbol="",entrez_chromosome="") 
  gene.df$entrez_id <- entrez_id
  gene.df$ensembl_id <- ifelse(length(ensembl_id)>0, ensembl_id, NA)
  gene.df$gene_symbol <- gene_symbol
  gene.df$entrez_chromosome <- entrez_chromosome
  return(gene.df)
}
#####################################################################################################################

# Search for ensembl gene IDs using ensembl transcript IDs

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

ensembl_from_transcript <- getBM(attributes = c('entrezgene_id',
                                            'ensembl_gene_id',
                                            'external_gene_name',
                                            'mgi_symbol',
                                            'chromosome_name',
                                            'start_position',
                                            'end_position'),
                             filters = 'ensembl_transcript_id',
                             values = manifest$Transcripts_Targeted,
                             mart = ensembl)


embltr_matched <- manifest %>%
  inner_join(ensembl_from_transcript, by=c("Gene_Symbol" = "mgi_symbol")) %>%
  filter(Entrez_ID == entrezgene_id)

# Without filtering for Entrez_ID == entrezgene_id, some probes get duplicate rows
# All duplicate rows are identical, except that biomart found 2 Entrez IDs.
# From looking online, both Entrez do indeed match to the same ensembl ID.
# Which is a pain!
# So I'm taking the entrezID that matches the one from NIH's fixed up manifest

# Find duplicated rows
check_n_matches <- embltr_matched %>% group_by(Probe_ID) %>% mutate(n_copies = n())
duplicate_rows1 <- check_n_matches %>% filter(n_copies > 1)

stopifnot(check_n_matches %>% filter(n_copies > 1) %>% nrow() == 0)

output_manifest1 <- embltr_matched %>%
  dplyr::select("Probe_ID",
                "Probe_Name",
                "Gene_Symbol",
                "ensembl_gene_id",
                "Entrez_ID",
                "Probe_Sequence",
                "Transcripts_Targeted")

###############################################################################################################
# Deal with remaining probes

remaining1 <- manifest %>%
  dplyr::filter(! Probe_ID %in% output_manifest1$Probe_ID)

remaining_noentrez <- remaining %>%
  dplyr::filter(is.na(Entrez_ID))

#############################################################################################################
# Get Ensembl IDs based on Entrez IDs

# Discussed with Matt 08 Feb 2024: don't use external_gene_name, it's not applicable for mouse
# Try using mgi_symbol

ensembl_from_entrez <- getBM(attributes = c('entrezgene_id',
                                            'ensembl_gene_id',
                                            'external_gene_name',
                                            'mgi_symbol',
                                            'chromosome_name',
                                            'start_position',
                                            'end_position'),
                             filters = 'entrezgene_id',
                             values = remaining1$Entrez_ID,
                             mart = ensembl)


# Filter for rows where manifest gene symbol and biomart gene name match
# Mutate to title case, sometimes mismatch is just case sensitivity
entrez_matched <- remaining1 %>%
  inner_join(ensembl_from_entrez, by=c("Entrez_ID" = "entrezgene_id")) %>%
  mutate(Gene_Symbol = str_to_title(Gene_Symbol),
         mgi_symbol =  str_to_title(mgi_symbol)) %>%
  dplyr::filter(Gene_Symbol == mgi_symbol)

# Find rows where manifest gene symbol and biomart gene name don't match
# Will have to deal with these somehow
symbol_mismatch <- remaining1 %>%
  inner_join(ensembl_from_entrez, by=c("Entrez_ID" = "entrezgene_id")) %>%
  dplyr::mutate(Gene_Symbol = str_to_title(Gene_Symbol),
                mgi_symbol =  str_to_title(mgi_symbol)) %>%
  dplyr::filter(Gene_Symbol != mgi_symbol)

nrow(symbol_mismatch)
# 212 rows with mismatch

# Find duplicate rows
check_n_matches2 <- entrez_matched %>% group_by(Probe_ID) %>% mutate(n_copies = n())
duplicate_rows2 <- check_n_matches2 %>% filter(n_copies > 1)

# Deal with probes matching to multiple ensembl IDs
# Just looked up the transcripts from manifest online, found ensembl ID from that
ensembles_to_remove <- c("ENSMUSG00000107877")

entrez_matched <- entrez_matched %>%
  filter(!ensembl_gene_id %in% ensembles_to_remove)

# Checks
check_n_matches2 <- entrez_matched %>% group_by(Probe_ID) %>% mutate(n_copies = n())
duplicate_rows2 <- check_n_matches2 %>% filter(n_copies > 1)
stopifnot(check_n_matches2 %>% filter(n_copies > 1) %>% nrow() == 0)


output_manifest2 <- entrez_matched %>%
  dplyr::select("Probe_ID",
                "Probe_Name",
                "Gene_Symbol",
                "ensembl_gene_id",
                "Entrez_ID",
                "Probe_Sequence",
                "Transcripts_Targeted") %>%
  bind_rows(output_manifest1)

check_data_unique(output_manifest2)

###########################################################################################################

remaining2 <- manifest %>%
  dplyr::filter(! Probe_ID %in% output_manifest2$Probe_ID)

# Search biomart by gene symbols

by_gene_name <- getBM(attributes = c('entrezgene_id',
                                     'ensembl_gene_id',
                                     'external_gene_name',
                                     'mgi_symbol',
                                     'chromosome_name',
                                     'start_position',
                                     'end_position'),
                      filters = 'mgi_symbol',
                      values = remaining2$Gene_Symbol,
                      mart = ensembl)

by_gene_name_singles <- by_gene_name %>%
  group_by(external_gene_name) %>%
  mutate(num_results = n())
by_gene_name_singles <- by_gene_name_singles %>% filter(num_results == 1)

matched_by_gene_name <- remaining2 %>%
  inner_join(by_gene_name_singles, by=c('Gene_Symbol' = 'mgi_symbol'))

# Entrez IDs don't match, but I'm suspicious of some of the entrez IDs that came from NIH anyway
# And in a few spot look-ups of the transcripts targeted, they match the ensembl ID found by biomart

output_manifest3 <- matched_by_gene_name %>%
  dplyr::select("Probe_ID",
                "Probe_Name",
                "Gene_Symbol",
                "ensembl_gene_id",
                "Entrez_ID",
                "Probe_Sequence",
                "Transcripts_Targeted") %>%
  bind_rows(output_manifest2)

check_data_unique(output_manifest3)

###############################################################################################################

remaining3 <- manifest %>%
  dplyr::filter(! Probe_ID %in% output_manifest3$Probe_ID)

# I don't trust the entrez ID for some of them
# So what if I look for the transcript-found ensembles without filtering for entrez ID matching?

embltr_matched2 <- remaining3 %>%
  inner_join(ensembl_from_transcript, by=c("Gene_Symbol" = "mgi_symbol")) %>%
  dplyr::select(-c("entrezgene_id")) %>%
  dplyr::distinct()

# Woohoo, 5 more probes...
output_manifest4 <- embltr_matched2 %>%
  dplyr::select("Probe_ID",
                "Probe_Name",
                "Gene_Symbol",
                "ensembl_gene_id",
                "Entrez_ID",
                "Probe_Sequence",
                "Transcripts_Targeted") %>%
  bind_rows(output_manifest3)

check_data_unique(output_manifest4)

###############################################################################################################

remaining4 <- manifest %>%
  dplyr::filter(! Probe_ID %in% output_manifest4$Probe_ID)

# Use org.Mm.eg.db to find ensembl IDs from transcript IDs

orgmmegdb_ensembl_results <- select(
  org.Mm.eg.db,
  keytype = 'ENSEMBLTRANS',
  columns = c('ALIAS','ENSEMBL','ENSEMBLTRANS','ENTREZID','SYMBOL'),
  keys = remaining4$Transcripts_Targeted) %>%
  dplyr::filter(! is.na(ENSEMBL))

orgmatches <- remaining4 %>%
  inner_join(orgmmegdb_ensembl_results, by=c("Transcripts_Targeted" = "ENSEMBLTRANS")) %>%
#  dplyr::filter(Entrez_ID == ENTREZID) %>%
  dplyr::filter(Gene_Symbol == ALIAS) %>%
  dplyr::rename('ensembl_gene_id' = 'ENSEMBL')

# this db has aliases, which is great because genes can have multiple names

check_data_unique(orgmatches)
# No duplicates :)


output_manifest5 <- orgmatches %>%
  dplyr::select("Probe_ID",
                "Probe_Name",
                "Gene_Symbol",
                "ensembl_gene_id",
                "Entrez_ID",
                "Probe_Sequence",
                "Transcripts_Targeted") %>%
  bind_rows(output_manifest4)

check_data_unique(output_manifest5)


###############################################################################################################

remaining5 <- manifest %>%
  dplyr::filter(! Probe_ID %in% output_manifest5$Probe_ID)



# Use org.Mm.eg.db to find ensembl IDs from gene symbols (matching to aliases)

orgmmegdb_alias_results <- select(
  org.Mm.eg.db,
  keytype = 'ALIAS',
  columns = c('ALIAS','ENSEMBL','ENTREZID','SYMBOL'),
  keys = remaining5$Gene_Symbol) %>%
  dplyr::filter(! is.na(ENSEMBL))

# Check for duplicates
check_n_matches5 <- orgmmegdb_alias_results %>% group_by(ALIAS) %>% mutate(n_copies = n())
duplicate_rows5 <- check_n_matches5 %>% filter(n_copies > 1)

# Remove duplicate rows

orgmmegdb_alias_nodup <- orgmmegdb_alias_results %>% 
  dplyr::filter(! ALIAS %in% duplicate_rows5$ALIAS)

orgmmegdb_alias_joined1 <- remaining5 %>%
  inner_join(orgmmegdb_alias_nodup, by=c("Gene_Symbol" = "ALIAS")) %>%
  dplyr::rename('ensembl_gene_id' = 'ENSEMBL')


# There are duplicates, same alias and ensemblID but different Entrez IDs...
# Go with the EntrezID that matches the manifest

orgmatches_alias_dedup <- remaining5 %>%
  inner_join(duplicate_rows5, by=c("Gene_Symbol" = "ALIAS")) %>%
  dplyr::filter(Entrez_ID == ENTREZID) %>%
  dplyr::distinct() %>%
  dplyr::rename('ensembl_gene_id' = 'ENSEMBL') %>%
  dplyr::select(-c("n_copies"))

orgmatches_alias <- bind_rows(orgmmegdb_alias_joined1, orgmatches_alias_dedup)

output_manifest6 <- orgmatches_alias %>%
  dplyr::select("Probe_ID",
                "Probe_Name",
                "Gene_Symbol",
                "ensembl_gene_id",
                "Entrez_ID",
                "Probe_Sequence",
                "Transcripts_Targeted") %>%
  bind_rows(output_manifest5)

check_data_unique(output_manifest6)

###############################################################################################################
# Creating final output manifest
# Not all the probes have an ensemblID, which is annoying, but it's very close!

remaining6 <- manifest %>%
  dplyr::filter(! Probe_ID %in% output_manifest6$Probe_ID)
remaining6$ensembl_gene_id <- NA
nrow(remaining6)

output_manifest_final <- remaining6 %>%
  bind_rows(output_manifest6) %>%
  dplyr::rename("Ensembl_Gene_ID" =  "ensembl_gene_id") %>%
  dplyr::arrange(Probe_ID)

check_data_unique(output_manifest_final)

		
output_file <- file.path(output_directory, output_filename)

write.csv(output_manifest_final, output_file, row.names=F)

###############################################################################################################


# Thing to try in the future to finish the remaining 165 probes without EnsemblIDs: NCBI search with Entrez
# However:
## At the moment it always 400 error times out, even with a system sleep built into the function
## I don't trust all the Entrez IDs to be correct!

# Search NCBI using entrez for those that have it?
# Except I don't trust the entrez...

remaining7 <- manifest %>%
  dplyr::filter(! Probe_ID %in% output_manifest6$Probe_ID)

remaining7wentrez <- remaining7 %>%
  dplyr::filter(! is.na(Entrez_ID))

# Still getting the HTTP 400 failure, arrrg
ncbi_results <- bind_rows(lapply(remaining5wentrez$Entrez_ID, get_gene_info_from_ncbi_using_entrez)) %>% filter(ensembl_id!="NULL")


# Gene symbol found by NCBI matches gene symbol from biomart, but not manifest
# What to do with this? Should I just accept the external_gene_name?
# Pretty sure that's what was done in previous manifest standardization
# But then probe name and gene symbol won't match
symbol_mismatch %>% 
  inner_join(ensembl_results10, by = c("Entrez_ID" = "entrez_id")) %>%
  dplyr::select(Probe_Name, Gene_Symbol, external_gene_name, gene_symbol)




