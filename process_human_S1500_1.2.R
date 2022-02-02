# script to fill in blank IDs in TempO-seq manifests (as much as possible) and unify columns across probe designs

# final design will be: 
# Probe ID, Gene Symbol, Gene Ensembl ID, Entrez ID, Ligated sequence, Aligned transcript IDs

###################
# Human S1500 1.2 #
###################

# The main issue with this set is that it is missing gene IDs completely
# It only has gene symbols and refseq IDs
# I'm not thrilled about either of those as lookup terms for identifiers,
# so I'll require that they both point to the same gene, and manually double-check the rest

library(tidyverse)
library(biomaRt)
library(annotate)
library(seqinr)
library(rentrez)
library(stringr)

# prevents some weird biomaRt SSL errors, as per https://github.com/grimbough/biomaRt/issues/31
# may not be necessary all of the time (I went weeks without needing it)
httr::set_config(httr::config(ssl_verifypeer = FALSE))

#####################################################################################################################
# functions

check_data_unique <- function(x){
  n <- x %>% nrow()
  n_unique <- x %>% unique() %>% nrow()
  stopifnot(n==n_unique)
  n_unique_probeid <- x %>% dplyr::select(Probe.ID) %>% unique() %>% nrow()
  stopifnot(n==n_unique_probeid)
}

#####################################################################################################################

# set up ensembl mart
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#####################################################################################################################

# read in manifest file
original_manifests_directory <- "original_manifests"

human_s1500_1.2_file <- file.path(original_manifests_directory, "181019 Human S1500+ Surrogate 1.2 Manifest.csv")

manifest <- read.csv(human_s1500_1.2_file)


# First, do some data conversion
manifest <- manifest %>%
  separate(Probe.name,sep='[_]',into=c(NA,'Probe.ID'), remove=FALSE) %>%
  rename("Probe.name" = "Probe.Name") %>%
  mutate(Probe.ID = as.numeric(Probe.ID))

output_manifest0 <- data.frame(Probe.ID=integer(),
                               Probe.Name=character(), 
                               Gene.Symbol=character(), 
                               ENSEMBL.Gene.ID=character(), 
                               Entrez.ID=integer(), 
                               Probe.Sequence=character(), 
                               stringsAsFactors=FALSE) 


check_data_unique(output_manifest0)

output_directory <- "output_manifests"
output_filename <- "Human_S1500_1.2_standardized.csv"

#####################################################################################################################

# Start by trying to match both the gene symbol and refseq ID

ensembl_and_entrez_for_missing <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','refseq_mrna'),
                             filters = c('external_gene_name','refseq_mrna'),
                             values = list(manifest$Gene.Symbol,manifest$Reference.Transcript),
                             mart = ensembl)

ensembl_and_entrez_matches <- ensembl_and_entrez_for_missing %>%
  inner_join(manifest, by=c("refseq_mrna"="Reference.Transcript", "external_gene_name"="Gene.Symbol"))

ensembl_and_entrez_single_chr <- ensembl_and_entrez_matches %>%
  group_by(Probe.Name) %>%
  arrange(chromosome_name, by_group=TRUE) %>%
  slice_head(n=1) %>%
  ungroup()

output_manifest1 <- ensembl_and_entrez_single_chr %>%
  dplyr::select(Probe.ID, Gene.Symbol = external_gene_name, ENSEMBL.Gene.ID = ensembl_gene_id, Entrez.ID = entrezgene_id, Probe.Sequence) %>%
  bind_rows(output_manifest0)

check_data_unique(output_manifest1)

#####################################################################################################################

remaining <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest1) %>%
  filter(is.na(Entrez.ID) | is.na(ENSEMBL.Gene.ID)) %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)


# if there is an ensembl record for the NM_ transcript, manually investigate why the gene name didn't match

ensembl_and_entrez_for_missing_2 <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','refseq_mrna'),
                                        filters = c('refseq_mrna'),
                                        values = remaining$Reference.Transcript,
                                        mart = ensembl)

single_chr <- ensembl_and_entrez_for_missing_2 %>%
  group_by(external_gene_name) %>%
  arrange(chromosome_name, by_group=TRUE) %>%
  slice_head(n=1) %>%
  ungroup()

remaining_with_ensembl <- single_chr %>%
  inner_join(remaining, by=c('refseq_mrna' = 'Reference.Transcript'))

#####################################################################################################################

# looks like lots of gene names have been updated. Try looking for the gene name in the synonyms, if there's a match
# we'll assume the ensembl name is correct and it's just been updated.

synonyms <- getBM(attributes = c('ensembl_gene_id','external_synonym'),
                  filters = c('ensembl_gene_id'),
                  values = remaining_with_ensembl$ensembl_gene_id,
                  mart = ensembl)

gene_name_in_synonyms <- synonyms %>% inner_join(remaining_with_ensembl, by=c('ensembl_gene_id', 'external_synonym'='Gene.Symbol'))

output_manifest2 <- gene_name_in_synonyms %>%
  dplyr::select(Probe.ID, Gene.Symbol = external_gene_name, ENSEMBL.Gene.ID = ensembl_gene_id, Entrez.ID = entrezgene_id, Probe.Sequence) %>%
  bind_rows(output_manifest1)

check_data_unique(output_manifest2)



remaining2 <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest2) %>%
  filter(is.na(Entrez.ID) | is.na(ENSEMBL.Gene.ID)) %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)


ensembl_and_entrez_for_missing_3 <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','refseq_mrna'),
                                          filters = c('refseq_mrna'),
                                          values = remaining2$Reference.Transcript,
                                          mart = ensembl)

single_chr_3 <- ensembl_and_entrez_for_missing_3 %>%
  group_by(external_gene_name) %>%
  arrange(chromosome_name, by_group=TRUE) %>%
  slice_head(n=1) %>%
  ungroup()

remaining_with_ensembl_3 <- single_chr_3 %>%
  inner_join(remaining2, by=c('refseq_mrna' = 'Reference.Transcript'))


# load manual corrections
manual_corrections <- read.csv('S1500_1.2_corrections.csv') %>% filter(!is.na(ncbi_id))

output_manifest3 <- manual_corrections %>% 
  left_join(remaining_with_ensembl_3, by="Probe.Name") %>% 
  dplyr::select(Probe.ID, Gene.Symbol = correct_gene_symbol, ENSEMBL.Gene.ID = ensembl_id, Entrez.ID = ncbi_id, Probe.Sequence) %>%
  bind_rows(output_manifest2)

check_data_unique(output_manifest3)


#####################################################################################################################


remaining3 <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest3, by="Probe.ID") %>%
  filter(is.na(Entrez.ID) & is.na(ENSEMBL.Gene.ID)) %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)


# ok, some of the transcripts are ensembl transcript IDs (for some reason...), so check that

ensembl_and_entrez_by_enst <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','refseq_mrna','ensembl_transcript_id'),
                                          filters = c('ensembl_transcript_id'),
                                          values = remaining3$Reference.Transcript,
                                          mart = ensembl)


output_manifest4 <- ensembl_and_entrez_by_enst %>%
  left_join(remaining3, by=c('ensembl_transcript_id' = 'Reference.Transcript')) %>%
  dplyr::select(Probe.ID, Gene.Symbol = external_gene_name, ENSEMBL.Gene.ID = ensembl_gene_id, Entrez.ID = entrezgene_id, Probe.Sequence) %>%
  bind_rows(output_manifest3)

check_data_unique(output_manifest4)

#####################################################################################################################
# start searching by gene name

remaining4 <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest4, by="Probe.ID") %>%
  filter(is.na(Entrez.ID) & is.na(ENSEMBL.Gene.ID)) %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)


ensembl_and_entrez_by_gene_name <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','refseq_mrna','ensembl_transcript_id'),
                                    filters = c('external_gene_name'),
                                    values = remaining4$Gene.Symbol,
                                    mart = ensembl)

single_chr_4 <- ensembl_and_entrez_by_gene_name %>%
  group_by(external_gene_name) %>%
  arrange(chromosome_name, by_group=TRUE) %>%
  slice_head(n=1) %>%
  ungroup()


output_manifest5 <- single_chr_4 %>%
  left_join(remaining4, by=c('external_gene_name' = 'Gene.Symbol')) %>%
  dplyr::select(Probe.ID, Gene.Symbol = external_gene_name, ENSEMBL.Gene.ID = ensembl_gene_id, Entrez.ID = entrezgene_id, Probe.Sequence) %>%
  distinct() %>%
  bind_rows(output_manifest4)

check_data_unique(output_manifest5)

#####################################################################################################################

# I suspect a lot of these are older gene names, so check by synonyms again. This time we'll manually double-check the results

remaining5 <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest5, by="Probe.ID") %>%
  filter(is.na(Entrez.ID) & is.na(ENSEMBL.Gene.ID)) %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)

ensembl_and_entrez_by_synonym <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','refseq_mrna','ensembl_transcript_id','external_synonym'),
                                         filters = c('external_synonym'),
                                         values = remaining5$Gene.Symbol,
                                         mart = ensembl)
single_chr_5 <- ensembl_and_entrez_by_synonym %>%
  group_by(external_gene_name) %>%
  arrange(chromosome_name, by_group=TRUE) %>%
  slice_head(n=1) %>%
  ungroup()


output_manifest6 <- single_chr_5 %>%
  left_join(remaining5, by=c('external_synonym' = 'Gene.Symbol')) %>%
  filter(!(Probe.ID %in% c("HIST1H2AK_20840", "HIST2H2AA4_2963", "HIST2H2AA4_2964", "H2AFV_28424"))) %>%
  dplyr::select(Probe.ID, Gene.Symbol = external_gene_name, ENSEMBL.Gene.ID = ensembl_gene_id, Entrez.ID = entrezgene_id, Probe.Sequence) %>%
  distinct() %>%
  bind_rows(output_manifest5)

check_data_unique(output_manifest6)


#####################################################################################################################

output_manifest_final <- manifest %>%
  left_join(output_manifest6, by=c("Probe.ID"), suffix=(c('','.o'))) %>%
  filter(is.na(Gene.Symbol.o )) %>%
  dplyr::select(Probe.ID) %>%
  left_join(manifest) %>%
  dplyr::select(Probe.ID, Gene.Symbol, Probe.Sequence) %>%
  bind_rows(output_manifest6) %>%
  mutate(Probe.Name = paste0(Gene.Symbol, "_", Probe.ID)) %>%
  inner_join(manifest, by="Probe.ID") 

check_data_unique(output_manifest_final)

output_manifest_formatted <- output_manifest_final %>%
  dplyr::select(Probe_ID = Probe.ID, Probe_Name = Probe.Name.x, Gene_Symbol = Gene.Symbol.x, Ensembl_Gene_ID = ENSEMBL.Gene.ID, Entrez_ID = Entrez.ID, Probe_Sequence = Probe.Sequence.x, Transcripts_Targeted = Transcripts.Targeted) %>%
  arrange(Probe_ID)

output_file <- file.path(output_directory, output_filename)

write.csv(output_manifest_formatted, output_file, row.names=F)




