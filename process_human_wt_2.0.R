# script to fill in blank IDs in TempO-seq manifests (as much as possible) and unify columns across probe designs

# final design will be: 
# Probe ID, Gene Symbol, Gene Ensembl ID, Entrez ID, Ligated sequence, Aligned transcript IDs

#################################
# Human Whole Transcriptome 2.0 #
#################################

# The main issue with this set is missing Ensembl IDs
# A subset have Entrez IDs that we can use to map to Ensembl IDs
# Or vice versa
# There are also 

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

is_chr <- function(data, chr_col){
  chr_str <- data[[chr_col]]
  numeric <- suppressWarnings(as.numeric(chr_str))
  return(case_when(
    numeric<25~TRUE,
    chr_str=="X"~TRUE,
    chr_str=="Y"~TRUE,
    chr_str=="MT"~TRUE,
    TRUE~FALSE))
}

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

#human_wt_2.0_file <- file.path(original_manifests_directory, "191004_Human_Whole_Transcriptome_2.0_Manifest.csv")
human_wt_2.0_file <- file.path(original_manifests_directory, "190620 Human Whole Transcriptome 2.0 Manifest.csv")

manifest <- read.csv(human_wt_2.0_file)


# First, do some data conversion

manifest[manifest == 'NULL'] <- NA
manifest$Entrez.ID <- as.integer(manifest$Entrez.ID)
manifest <- manifest %>% separate(Probe.Genomic.Coordinates,sep='[:-]',into=c('probe_chrom','probe_start','probe_end', NA))

manifest <- manifest %>% separate(Probe.Name,sep='[_]',into=c(NA,'Probe.ID'), remove=FALSE)

output_manifest0 <- manifest %>%
  filter(!is.na(Entrez.ID) & !is.na(ENSEMBL.Gene.ID)) %>%
  dplyr::select(Probe.ID, Probe.Name, Gene.Symbol, ENSEMBL.Gene.ID, Entrez.ID, Probe.Sequence)

check_data_unique(output_manifest0)

output_directory <- "output_manifests"
output_filename <- "Human_Whole_Transcriptome_2.0_standardized.csv"

#####################################################################################################################

# Start with rows missing Ensembl IDs, with Entrez IDs

missing_ensembl_with_entrez <- manifest %>%
  filter(is.na(ENSEMBL.Gene.ID) & !is.na(Entrez.ID)) 

ensembl_for_missing <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'),
                             filters = 'entrezgene_id',
                             values = missing_ensembl_with_entrez$Entrez.ID,
                             mart = ensembl)

entrez_matched <- missing_ensembl_with_entrez %>%
  inner_join(ensembl_for_missing, by=c("Entrez.ID" = "entrezgene_id")) %>%
  filter(probe_chrom == chromosome_name)

# manually fix some problems 
entrez_matched$Gene.Symbol[entrez_matched$Probe.ID == "2955"] <- "H2BC12" # updated gene symbol
entrez_matched$Probe.Name[entrez_matched$Probe.ID == "2955"] <- "H2BC12" # updated gene symbol
entrez_matched$Gene.Symbol[entrez_matched$Probe.ID == "10706"] <- "MARCHF4" # updated gene symbol
entrez_matched$Probe.Name[entrez_matched$Probe.ID == "10706"] <- "MARCHF4" # updated gene symbol

# error checking

stopifnot(all(entrez_matched$Gene.Symbol == entrez_matched$external_gene_name))
check_n_matches <- entrez_matched %>% group_by(Probe.ID) %>% mutate(n_copies = n())
stopifnot(check_n_matches %>% filter(n_copies > 1) %>% nrow() == 0)

output_manifest1 <- entrez_matched %>%
  dplyr::select(Probe.ID, Gene.Symbol, ENSEMBL.Gene.ID = ensembl_gene_id, Entrez.ID, Probe.Sequence) %>%
  bind_rows(output_manifest0)

check_data_unique(output_manifest1)

#####################################################################################################################

# now probes  with Ensembl IDs but missing Entrez IDs

missing_entrez_with_ensembl <- manifest %>%
  filter(is.na(Entrez.ID) & !is.na(ENSEMBL.Gene.ID)) 

entrez_for_missing <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'),
                            filters = 'ensembl_gene_id',
                            values = missing_entrez_with_ensembl %>% dplyr::select(ENSEMBL.Gene.ID),
                            mart = ensembl) %>% filter(!is.na(entrezgene_id))


ensembl_matched <- missing_entrez_with_ensembl %>%
  inner_join(entrez_for_missing, by=c("ENSEMBL.Gene.ID" = "ensembl_gene_id"))

# there seem to be some ensembl genes that match to multiple entrez IDs, but spot checking a few 
# suggests that the larger number is a read-through gene in most cases, so assume that
# the smaller entrez ID is the one we want
# ensembl_matched <- ensembl_matched %>% group_by(ENSEMBL.Gene.ID) %>% slice_max(order_by = entrezgene_id, n=1) %>% ungroup()

# manually fix some problems (all "ACNNNNN" accessions for symbols--one has been assigned a gene symbol so
# we'll use that, the rest don't have ensembl symbols so we'll force it to use the accession numbers)
ensembl_matched$Gene.Symbol[ensembl_matched$Probe.ID == "5295"] <- "PP7080"
ensembl_matched$Gene.Symbol[ensembl_matched$Probe.ID == "12290"] <- "FOXL3"
ensembl_matched$external_gene_name[ensembl_matched$Probe.ID == "18935"] <- "AL449403.2"
ensembl_matched$external_gene_name[ensembl_matched$Probe.ID == "28481"] <- "AC007666.1"
ensembl_matched$external_gene_name[ensembl_matched$Probe.ID == "29034"] <- "AC132008.2"
ensembl_matched$external_gene_name[ensembl_matched$Probe.ID == "88136"] <- "AC008687.4"
ensembl_matched$external_gene_name[ensembl_matched$Probe.ID == "89032"] <- "AC008397.1"
ensembl_matched$external_gene_name[ensembl_matched$Probe.ID == "89097"] <- "AL022312.1"

# error checking

stopifnot(all(ensembl_matched$Gene.Symbol == ensembl_matched$external_gene_name))
stopifnot(all(ensembl_matched$probe_chrom == ensembl_matched$chromosome_name))
check_n_matches <- ensembl_matched %>% group_by(Probe.ID) %>% mutate(n_copies = n())
stopifnot(check_n_matches %>% filter(n_copies > 1) %>% nrow() == 0)

output_manifest2 <- ensembl_matched %>%
  dplyr::select(Probe.ID, Gene.Symbol, ENSEMBL.Gene.ID, Entrez.ID = entrezgene_id, Probe.Sequence) %>%
  bind_rows(output_manifest1)

check_data_unique(output_manifest2)

#####################################################################################################################

# Next we are going to search NCBI for the genes that are still missing either an Ensembl or Entrez ID 

remaining <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest2) %>%
  filter(is.na(Entrez.ID) | is.na(ENSEMBL.Gene.ID)) %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)

remaining_missing_entrez <- remaining %>% filter(is.na(Entrez.ID) & !is.na(ENSEMBL.Gene.ID))

get_gene_info_from_ncbi_using_ensembl <- function(ensembl_id){
  search_result <- entrez_search(db="gene", term=ensembl_id)
  Entrez.ID <- search_result[[1]]
  if(length(Entrez.ID) == 1){
    # extract gene info
    summary <- entrez_summary(id=Entrez.ID, db="gene")
    Gene.Symbol <- summary$name
    entrez_chromosome <- summary$chromosome
    # create data frame of information
    gene.df <- data.frame(ensembl_id, Entrez.ID, Gene.Symbol, entrez_chromosome) 
    gene.df$ensembl_id <- ensembl_id
    gene.df$Entrez.ID <- ifelse(!is.null(Entrez.ID), Entrez.ID, NA)
    gene.df$Gene.Symbol <- ifelse(!is.null(Gene.Symbol), Gene.Symbol, NA)
    gene.df$entrez_chromosome <- ifelse(!is.null(entrez_chromosome), entrez_chromosome, NA)
    return(gene.df)
  }
}

entrez_results <- bind_rows(lapply(remaining_missing_entrez$ENSEMBL.Gene.ID, get_gene_info_from_ncbi_using_ensembl)) %>% unique()

ensembl_matched <- remaining_missing_entrez %>%
  inner_join(entrez_results, by=c("ENSEMBL.Gene.ID" = "ensembl_id"), suffix=c(".remaining", ".entrez"))

# manual corrections

ensembl_matched$Gene.Symbol.remaining[ensembl_matched$Probe.ID == "3844"] <- "LOC728554"
ensembl_matched$Gene.Symbol.remaining[ensembl_matched$Probe.ID == "15246"] <- "LOC728554"
ensembl_matched$Gene.Symbol.remaining[ensembl_matched$Probe.ID == "26933"] <- "SNHG28"
ensembl_matched$Gene.Symbol.remaining[ensembl_matched$Probe.ID == "87748"] <- "FAM90A10"
ensembl_matched$Gene.Symbol.remaining[ensembl_matched$Probe.ID == "87777"] <- "FAM90A22"
ensembl_matched$Gene.Symbol.remaining[ensembl_matched$Probe.ID == "88043"] <- "OR1R1P"
ensembl_matched$Gene.Symbol.remaining[ensembl_matched$Probe.ID == "88151"] <- "PTP4A1"
ensembl_matched$Gene.Symbol.remaining[ensembl_matched$Probe.ID == "88172"] <- "LOC93622"


# error checking

stopifnot(all(ensembl_matched$Gene.Symbol.remaining == ensembl_matched$Gene.Symbol.entrez))
stopifnot(all(ensembl_matched$probe_chrom == ensembl_matched$entrez_chromosome))

ensembl_matched$Entrez.ID <- as.integer(ensembl_matched$Entrez.ID.entrez)

output_manifest3 <- ensembl_matched %>%
  dplyr::select(Probe.ID, Gene.Symbol=Gene.Symbol.remaining, ENSEMBL.Gene.ID, Entrez.ID, Probe.Sequence) %>%
  bind_rows(output_manifest2)

check_data_unique(output_manifest3)

#####################################################################################################################

remaining_missing_ensembl <- remaining %>% filter(!is.na(Entrez.ID) & is.na(ENSEMBL.Gene.ID))

get_gene_info_from_ncbi_using_entrez <- function(Entrez.ID){
  fetch_result <- entrez_fetch(db="gene",id=Entrez.ID, rettype="xml", parsed=TRUE)
  fetch_list <- xmlToList(fetch_result)
  dbtags <- fetch_list$Entrezgene$Entrezgene_gene[["Gene-ref"]][["Gene-ref_db"]]
  external_array <- as.vector(as.matrix(as.data.frame(dbtags)))
  if(length(external_array)>2){
    external_dbs <- external_array[seq(1, length(external_array), 2)]
    external_ids <- external_array[seq(2, length(external_array), 2)]
    external_data <- data.frame(db=external_dbs, id=external_ids)
    ensembl_id <- external_data %>% filter(db=="Ensembl") %>% dplyr::select(id)
    if(nrow(ensembl_id) == 0)
    ensembl_id = NA
  } else {
    ensembl_id = NA
  }
  entrez_chromosome <- fetch_list$Entrezgene$Entrezgene_source$BioSource$BioSource_subtype$SubSource$SubSource_name
  Gene.Symbol <- fetch_list$Entrezgene$Entrezgene_gene[["Gene-ref"]][["Gene-ref_locus"]]
  
  gene.df <- data.frame(Entrez.ID="",ensembl_id="",Gene.Symbol="",entrez_chromosome="") 
  gene.df$Entrez.ID <- Entrez.ID
  gene.df$ensembl_id <- ifelse(length(ensembl_id)>0, ensembl_id, NA)
  gene.df$Gene.Symbol <- Gene.Symbol
  gene.df$entrez_chromosome <- entrez_chromosome
  return(gene.df)
}

ensembl_results <- bind_rows(lapply(remaining_missing_ensembl$Entrez.ID, get_gene_info_from_ncbi_using_entrez)) %>% filter(ensembl_id!="NULL")

# two measly genes, argh
entrez_matched <- remaining_missing_ensembl %>%
  inner_join(ensembl_results, by=c("Entrez.ID"))
# error checking

stopifnot(all(entrez_matched$Gene.Symbol.x == entrez_matched$Gene.Symbol.y))
stopifnot(all(entrez_matched$probe_chrom == entrez_matched$entrez_chromosome))

entrez_matched$ensembl_id <- as.character(entrez_matched$ensembl_id)

output_manifest4 <- entrez_matched %>%
  dplyr::select(Probe.ID, Gene.Symbol=Gene.Symbol.x, ENSEMBL.Gene.ID=ensembl_id, Entrez.ID, Probe.Sequence) %>%
  bind_rows(output_manifest3)

check_data_unique(output_manifest4)

#####################################################################################################################

remaining2 <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest4) %>%
  filter(is.na(Entrez.ID) | is.na(ENSEMBL.Gene.ID)) %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)

# OK, before we start blasting things let's try with the gene names

by_gene_name <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'),
                      filters = 'external_gene_name',
                      values = remaining2$Gene.Symbol,
                      mart = ensembl)

by_gene_name_singles <- by_gene_name %>%
  group_by(external_gene_name) %>%
  mutate(num_results = n())
by_gene_name_singles <- by_gene_name_singles %>% filter(num_results == 1)

matched_by_gene_name <- remaining2 %>%
  inner_join(by_gene_name_singles, by=c('Gene.Symbol' = 'external_gene_name'))

matched_by_gene_name$Entrez.ID <- pmin(matched_by_gene_name$Entrez.ID, matched_by_gene_name$entrezgene_id, na.rm = TRUE)
matched_by_gene_name$ENSEMBL.Gene.ID <- pmin(matched_by_gene_name$ENSEMBL.Gene.ID, matched_by_gene_name$ensembl_gene_id, na.rm = TRUE)
#####################################################################################################################

# Next we are going to take the remaining probes and use BLAT against a set of human mRNAs
# Start the server outside of R using
#  ~/bin/gfServer start 127.0.0.1 1234 -stepSize=5 /home/katecook/Documents/HC_EHSRB/data/genomes/GRCh38_latest_rna.2bit

fasta_file <- "/tmp/probe_sequences.fa"

write.fasta(as.list(remaining2$Probe.Sequence), remaining2$Probe.Name, fasta_file, open = "w", nbchar = 60, as.string = TRUE)

alignment_file <- "/tmp/out.psl"

system(paste("~/bin/gfClient -minScore=20 -minIdentity=0 127.0.0.1 1234 ~/Documents/HC_EHSRB/data/genomes/",fasta_file,alignment_file))

alignments <- read.table(alignment_file, header=FALSE, sep='\t', skip=5)

names(alignments) <- c('match','mismatch','rep.match','Ns','Q_gap_count','Q_gap_bases','T_gap_count','T_gap_bases','strand','Q_name','Q_size','Q_start','Q_end','T_name','T_size','T_start','T_end','block_count','blockSizes','qStarts','tStarts')

alignments <- alignments %>%
  mutate(score = match - mismatch - Q_gap_count - T_gap_count)

best_align_per_probe <- alignments %>%
  group_by(Q_name) %>%
  slice_max(score, n=1) %>%
  slice_head(n=1) %>%
  ungroup()

refseq_matches <- best_align_per_probe %>%
  separate(T_name,into="refseq",sep='[.]') %>%
  dplyr::select(Q_name,refseq)

ensembl_for_refseq <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','refseq_mrna'),
                            filters = 'refseq_mrna',
                            values = refseq_matches$refseq, 
                            mart = ensembl)

# if there are multiple genes, sort by chromosome and pickthe top one (prioritizes standard chromosomes)
ensembl_for_refseq <- ensembl_for_refseq %>%
  group_by(external_gene_name) %>%
  arrange(chromosome_name, by_group=TRUE) %>%
  slice_head(n=1) %>%
  ungroup()

matched_by_sequence <- refseq_matches %>%
  left_join(ensembl_for_refseq, by=c("refseq" = "refseq_mrna")) %>%
  inner_join(remaining2, by=c("Q_name" = "Probe.Name")) %>% 
  dplyr::rename("Probe.Name" = "Q_name") %>%
  filter(is.na(external_gene_name) | external_gene_name == Gene.Symbol) %>%
  filter(is.na(chromosome_name) | is.na(probe_chrom) | chromosome_name == probe_chrom) %>%
  filter(!is.na(ensembl_gene_id))

matched_by_sequence$Entrez.ID <- pmin(matched_by_sequence$Entrez.ID, matched_by_sequence$entrezgene_id, na.rm = TRUE)
matched_by_sequence$ENSEMBL.Gene.ID <- pmin(matched_by_sequence$ENSEMBL.Gene.ID, matched_by_sequence$ensembl_gene_id, na.rm = TRUE)

# ok, all of the matched by sequence results are also in matched by gene name, so to be
# careful let's require that they match to the same ensembl id

# require that the matched by gene name and matched by sequence IDs are identical if they both exist

good_results <- matched_by_gene_name %>%
  left_join(matched_by_sequence, by = c("Probe.Name", "Probe.ID", "Gene.Symbol", "Probe.Sequence", "ENSEMBL.Gene.ID", "Entrez.ID"), suffix=c(".seq", ".symbol")) %>%
  filter(!is.na(Entrez.ID))

output_manifest5 <- good_results %>%
  dplyr::select(Probe.ID, Gene.Symbol, ENSEMBL.Gene.ID, Entrez.ID, Probe.Sequence) %>%
  unique() %>%
  bind_rows(output_manifest4)

check_data_unique(output_manifest5)


#####################################################################################################################

# fill in the rest manually

remaining3 <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest5) %>%
  filter(is.na(Entrez.ID) | is.na(ENSEMBL.Gene.ID)) %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)

# update entrez IDs
remaining3$Entrez.ID[remaining3$ENSEMBL.Gene.ID == "ENSG00000282801"] <- 28299

# update ensembl IDs
remaining3$Gene.Symbol[remaining3$Probe.Name == "LINC01279_28403"] <- "CCDC80"
remaining3$Entrez.ID[remaining3$Probe.Name == "LINC01279_28403"] <- 151887
remaining3$ENSEMBL.Gene.ID[remaining3$Probe.Name == "LINC01279_28403"] <- "ENSG00000091986"

# both are missing, search using gene name or sequence...whatever works

remaining3$ENSEMBL.Gene.ID[remaining3$Gene.Symbol == "ACTBP9"] <- "ENSG00000266920"
remaining3$Entrez.ID[remaining3$Gene.Symbol == "ACTBP9"] <- 69
remaining3$Entrez.ID[remaining3$Gene.Symbol == "ANKRD20A13P"] <- 100132733
remaining3$Entrez.ID[remaining3$Gene.Symbol == "BAGE5"] <- 85316
remaining3$ENSEMBL.Gene.ID[remaining3$Gene.Symbol == "CTAGE8"] <- "ENSG00000289604"
remaining3$Entrez.ID[remaining3$Gene.Symbol == "CTAGE8"] <- 100142659
remaining3$ENSEMBL.Gene.ID[remaining3$Gene.Symbol == "DDTL"] <- "ENSG00000099974"
remaining3$Entrez.ID[remaining3$Gene.Symbol == "DDTL"] <- 100037417
remaining3$Gene.Symbol[remaining3$Probe.Name == "HIST1H3D_17823"] <- "H3C4"
remaining3$ENSEMBL.Gene.ID[remaining3$Probe.Name == "HIST1H3D_17823"] <- "ENSG00000197409"
remaining3$Entrez.ID[remaining3$Probe.Name == "HIST1H3D_17823"] <- 8351
remaining3$ENSEMBL.Gene.ID[remaining3$Gene.Symbol == "HLA-V"] <- "ENSG00000181126"
remaining3$Entrez.ID[remaining3$Gene.Symbol == "HLA-V"] <- 352962
remaining3$ENSEMBL.Gene.ID[remaining3$Gene.Symbol == "HNRNPCL2"] <- "ENSG00000275774"
remaining3$Entrez.ID[remaining3$Gene.Symbol == "HNRNPCL2"] <- 440563
remaining3$ENSEMBL.Gene.ID[remaining3$Gene.Symbol == "HSP90AA2P"] <- "ENSG00000224411"
remaining3$Entrez.ID[remaining3$Gene.Symbol == "HSP90AA2P"] <- 3324
remaining3$Gene.Symbol[remaining3$Probe.Name == "IGHV3_28168"] <- "IGHV3-48"
remaining3$ENSEMBL.Gene.ID[remaining3$Probe.Name == "IGHV3_28168"] <- "ENSG00000211964"
remaining3$Entrez.ID[remaining3$Probe.Name == "IGHV3_28168"] <- 28424
remaining3$Gene.Symbol[remaining3$Probe.Name == "IGHV4_28203"] <- "IGHV4-59"
remaining3$ENSEMBL.Gene.ID[remaining3$Probe.Name == "IGHV4_28203"] <- "ENSG00000224373"
remaining3$Entrez.ID[remaining3$Probe.Name == "IGHV4_28203"] <- 28392
remaining3$Gene.Symbol[remaining3$Probe.Name == "IGKV1_28187"] <- "IGKV1-5"
remaining3$Entrez.ID[remaining3$Probe.Name == "IGKV1_28187"] <- 28299
remaining3$Gene.Symbol[remaining3$Probe.Name == "IGKV3D_28174"] <- "IGKV3D-11"
remaining3$ENSEMBL.Gene.ID[remaining3$Probe.Name == "IGKV3D_28174"] <- "ENSG00000211632"
remaining3$Entrez.ID[remaining3$Probe.Name == "IGKV3D_28174"] <- 28876
remaining3$Gene.Symbol[remaining3$Probe.Name == "IGLV1_28171"] <- "IGLV1-44"
remaining3$ENSEMBL.Gene.ID[remaining3$Probe.Name == "IGLV1_28171"] <- "ENSG00000211651"
remaining3$Entrez.ID[remaining3$Probe.Name == "IGLV1_28171"] <- 28823
remaining3$Gene.Symbol[remaining3$Probe.Name == "KJ904049_34049"] <- "TXLNGY"
remaining3$ENSEMBL.Gene.ID[remaining3$Probe.Name == "KJ904049_34049"] <- "ENSG00000131002"
remaining3$Entrez.ID[remaining3$Probe.Name == "KJ904049_34049"] <- 246126
#####################################################################
# Not sure what to do with this one......bicistronic gene
#####################################################################
remaining3$Gene.Symbol[remaining3$Probe.Name == "LUZP6_19140"] <- "MTPN"
remaining3$ENSEMBL.Gene.ID[remaining3$Probe.Name == "LUZP6_19140"] <- "ENSG00000105887"
remaining3$Entrez.ID[remaining3$Probe.Name == "LUZP6_19140"] <- 136319
#####################################################################
remaining3$Entrez.ID[remaining3$Probe.Name == "MGC70870_28188"] <- 403340
remaining3$Entrez.ID[remaining3$Probe.Name == "PCDHGCT_34000"] <- 56118
remaining3$ENSEMBL.Gene.ID[remaining3$Probe.Name == "PRAMEF33_14999"] <- "ENSG00000237700"
remaining3$Entrez.ID[remaining3$Probe.Name == "PRAMEF33_14999"] <- 645382
remaining3$Entrez.ID[remaining3$Probe.Name == "PRR21_17743"] <- 643905
remaining3$Gene.Symbol[remaining3$Probe.Name == "SSX8_34034"] <- "SSX8P"
remaining3$ENSEMBL.Gene.ID[remaining3$Probe.Name == "SSX8_34034"] <- "ENSG00000157965"
remaining3$Entrez.ID[remaining3$Probe.Name == "SSX8_34034"] <- 280659
remaining3$ENSEMBL.Gene.ID[remaining3$Gene.Symbol == "XR"] <- "ENSG00000197617"
remaining3$Entrez.ID[remaining3$Gene.Symbol == "XR"] <- 317705

output_manifest_final <- remaining3 %>%
  dplyr::select(Probe.ID, Gene.Symbol, ENSEMBL.Gene.ID, Entrez.ID, Probe.Sequence) %>%
  bind_rows(output_manifest5) %>%
  left_join(manifest %>% dplyr::select(Probe.ID, Transcripts.Covered)) 

check_data_unique(output_manifest_final)

output_manifest_formatted <- output_manifest_final %>%
  mutate(Probe_Name = paste0(Gene.Symbol, "_", Probe.ID), Probe_ID = as.integer(Probe.ID)) %>%
  dplyr::select(Probe_ID, Probe_Name, Gene_Symbol = Gene.Symbol, Ensembl_Gene_ID = ENSEMBL.Gene.ID, Entrez_ID = Entrez.ID, Probe_Sequence = Probe.Sequence, Transcripts_Targeted = Transcripts.Covered) %>%
  arrange(Probe_ID)
  
  
  
output_file <- file.path(output_directory, output_filename)

write.csv(output_manifest_formatted, output_file, row.names=F)


