# script to fill in blank IDs in TempO-seq manifests (as much as possible) and unify columns across probe designs

# final design will be: 
# Probe ID, Gene Symbol, Gene Ensembl ID, Entrez ID, Ligated sequence, Aligned transcript IDs

#################################
# Human Whole Transcriptome 2.1 #
#################################


library(tidyverse)
library(biomaRt)
library(annotate)
library(seqinr)
library(rentrez)
library(stringr)
library(XML)

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

human_wt_2.1_file <- file.path(original_manifests_directory, "100739_homo_sapiens_wt_21_full_manifest_revB.csv")

manifest <- read.csv(human_wt_2.1_file)
manifest_2.1_original <- read.csv(human_wt_2.1_file)

# First, do some data conversion

manifest[manifest == 'NULL'] <- NA

# Change column names to match hwt2.0
manifest <- manifest %>%
  dplyr::rename(Probe.ID = PROBE_ID,
                Probe.Name = PROBE_NAME,
                Gene.Symbol = GENE_SYMBOL,
                ENSEMBL.Gene.ID = ENSEMBL_GENE_ID,
                Entrez.ID = ENTREZ_ID,
                Probe.Sequence = PROBE_SEQUENCE,
                Probe.Genomic.Coordinates = PROBE_COORDINATE)

# Remove index number column
manifest <- manifest %>%
  dplyr::select(-X)

manifest$Entrez.ID <- as.integer(manifest$Entrez.ID)
manifest <- manifest %>% separate(Probe.Genomic.Coordinates,sep='[:-]',into=c('probe_chrom','probe_start','probe_end', NA))

manifest <- manifest %>% separate(Probe.Name,sep='[_]',into=c(NA,'Probe.ID'), remove=FALSE)

output_manifest0 <- manifest %>%
  filter(!is.na(Entrez.ID) & ENSEMBL.Gene.ID != "") %>%
  dplyr::select(Probe.ID, Probe.Name, Gene.Symbol, ENSEMBL.Gene.ID, Entrez.ID, Probe.Sequence)

check_data_unique(output_manifest0)

output_directory <- "output_manifests"
output_filename <- "Human_Whole_Transcriptome_2.1_standardized.csv"

#####################################################################################################################
# Count missing info in original manifest

# Missing enseml
manifest %>% filter(ENSEMBL.Gene.ID == "" | ENSEMBL.Gene.ID == "NA") %>% 
  summarise(n = n())

# Missing Entrez
manifest %>% filter(is.na(Entrez.ID)) %>% 
  summarise(n = n())

# Missing both Ensembl IDs and Entrez IDs?
manifest %>% filter((ENSEMBL.Gene.ID == "") & is.na(Entrez.ID)) %>% 
  summarise(n = n())


#####################################################################################################################
# Probes with entrez IDs but not ensembl

missing_ensembl_with_entrez <- manifest %>%
  filter((ENSEMBL.Gene.ID == "") & !is.na(Entrez.ID))


ensembl_for_missing <- getBM(attributes = 
                               c('entrezgene_id',
                                 'ensembl_gene_id',
                                 'external_gene_name',
                                 'chromosome_name',
                                 'start_position',
                                 'end_position'),
                             filters = 'entrezgene_id',
                             values = missing_ensembl_with_entrez$Entrez.ID,
                             mart = ensembl)

entrez_matched <- missing_ensembl_with_entrez %>%
  inner_join(ensembl_for_missing, by=c("Entrez.ID" = "entrezgene_id")) %>%
  filter(probe_chrom == chromosome_name)

# find rows that need manual(ish) fixing
gene_name_mismatched1 <- entrez_matched[entrez_matched$Gene.Symbol != entrez_matched$external_gene_name,] %>% 
  dplyr::select(c(Probe.ID,Gene.Symbol, external_gene_name))

# For all row names with non-matching Gene Symbol and external gene name,
# If there's an external gene name, replace Gene Symbol with that
# If not, replace blank external gene name with existing Gene symbol
# This replaces the rows of manual fixing that Kate was doing

entrez_matched <- entrez_matched %>% 
  mutate(Gene.Symbol = ifelse(Gene.Symbol != external_gene_name & 
                                !is.na(external_gene_name) & 
                                external_gene_name != "", 
                              external_gene_name, Gene.Symbol), 
         external_gene_name = ifelse(external_gene_name == "" & 
                                       !is.na(Gene.Symbol) & 
                                       Gene.Symbol != "", 
                                     Gene.Symbol, external_gene_name)
  )


# Check for multiple ensemble_gene_ids per probe
check_n_matches <- entrez_matched %>% group_by(Probe.ID) %>% mutate(n_copies = n())

check_n_matches %>% dplyr::filter(n_copies > 1) %>% 
  dplyr::select(c(Gene.Symbol,
           Probe.Name,
           probe_chrom,
           Entrez.ID,
           ensembl_gene_id,
           external_gene_name,
           start_position, end_position))

## Manually fix multiple ensemble_gene_ids
## In this case, we get duplicates like this:

# Probe.ID Gene.Symbol Probe.Name   probe_chrom Entrez.ID ensembl_gene_id external_gene_name start_position end_position
#1 13971    ARL17A      ARL17A_13971 17              51326 ENSG00000185829 ARL17A                   46499780     46579695
#2 13971    ARL17B      ARL17A_13971 17              51326 ENSG00000228696 ARL17B                   46274784     4636176

# And I'm just going to go with the ones where the Gene.Symbol/external_gene_name matches the Probe.Name

# Also get duplicates for PDE4C_27660. 
#Following Kate's thought that the larger number is a read-through gene in most cases, 
# I'll take the smaller one (based on start and end positions)

ensembles_to_remove <- c("ENSG00000228696", "ENSG00000285188")

entrez_matched <- entrez_matched %>%
  filter(!ensembl_gene_id %in% ensembles_to_remove)

check_n_matches_postfix <- entrez_matched %>% group_by(Probe.ID) %>% mutate(n_copies = n())

# Add ensembl-fixed rows to output manifest
output_manifest1 <- entrez_matched %>%
  dplyr::select(Probe.ID, Probe.Name, Gene.Symbol, ENSEMBL.Gene.ID = ensembl_gene_id, Entrez.ID, Probe.Sequence) %>%
  bind_rows(output_manifest0)



# error checking
stopifnot(all(entrez_matched$Gene.Symbol == entrez_matched$external_gene_name))

stopifnot(check_n_matches_postfix %>% filter(n_copies > 1) %>% nrow() == 0)
check_data_unique(output_manifest1)

# Code useful for error checking
# Find non-unique probe ids
#output_manifest1 %>% group_by(Probe.ID) %>% filter(n()>1) %>% summarize(n=n()) %>% select(c(Probe.ID,n))

#####################################################################################################################

# now probes  with Ensembl IDs but missing Entrez IDs

missing_entrez_with_ensembl <- manifest %>%
  filter(is.na(Entrez.ID) & !is.na(ENSEMBL.Gene.ID) & ENSEMBL.Gene.ID != "") 

entrez_for_missing <- getBM(attributes = c('entrezgene_id',
                                           'ensembl_gene_id',
                                           'external_gene_name',
                                           'chromosome_name',
                                           'start_position',
                                           'end_position'),
                            filters = 'ensembl_gene_id',
                            values = missing_entrez_with_ensembl %>% 
                              dplyr::select(ENSEMBL.Gene.ID),
                            mart = ensembl) %>% 
                            filter(!is.na(entrezgene_id))


ensembl_matched <- missing_entrez_with_ensembl %>%
  inner_join(entrez_for_missing, by=c("ENSEMBL.Gene.ID" = "ensembl_gene_id"))

# Comment from Kate's wt2.0 code:
# there seem to be some ensembl genes that match to multiple entrez IDs, but spot checking a few 
# suggests that the larger number is a read-through gene in most cases, so assume that
# the smaller entrez ID is the one we want
ensembl_matched <- ensembl_matched %>% 
  group_by(ENSEMBL.Gene.ID) %>% 
  slice_max(order_by = entrezgene_id, n=1) %>% 
  ungroup()

# Find rows that need to be manually (ish) fixed
ensembl_matched_to_fix_manually <- ensembl_matched[ensembl_matched$Gene.Symbol != ensembl_matched$external_gene_name,] %>% 
  dplyr::select(c(Probe.ID,Gene.Symbol, external_gene_name))

# For all row names with non-matching Gene Symbol and external gene name,
# If there's an external gene name (ensembl symbol), replace Gene Symbol with that
# If not, replace blank external gene name (ensembl symbol) with existing Gene symbol (accession number)
# This replaces the rows of manual fixing that Kate was doing

ensembl_matched <- ensembl_matched %>% 
  mutate(Gene.Symbol = ifelse(Gene.Symbol != external_gene_name & 
                                !is.na(external_gene_name) & 
                                external_gene_name != "", 
                              external_gene_name, Gene.Symbol), 
         external_gene_name = ifelse(external_gene_name == "" & 
                                       !is.na(Gene.Symbol) & 
                                       Gene.Symbol != "", 
                                     Gene.Symbol, external_gene_name)
         )


# Find rows where probe chrom and chromosome name don't match
ensembl_matched[ensembl_matched$probe_chrom != ensembl_matched$chromosome_name,] %>% 
  dplyr::select(c(Probe.ID,probe_chrom, chromosome_name))

# All rows with long chromosome_names, the prob_chrom is the same except with CHR_ prepended. Ex.
#Probe.ID probe_chrom                 chromosome_name              
#<chr>    <chr>                       <chr>                        
#  1 93284    CHR_HSCHR6_MHC_QBL_CTG1     HSCHR6_MHC_QBL_CTG1  
# so it should be fine to overwrite probe_chrom with chromosome_name. No loss of info
#Otherwise, unmatched rows are just because probe_chrome is empty. Can overwrite that with info from chromosome_name

ensembl_matched <- ensembl_matched %>% 
  mutate(probe_chrom = ifelse(probe_chrom != chromosome_name &
                                str_starts(probe_chrom, "CHR")  &
                                !is.na(chromosome_name),
                              chromosome_name, probe_chrom))

ensembl_matched <- ensembl_matched %>% 
  mutate(probe_chrom = ifelse(probe_chrom != chromosome_name &
                              probe_chrom == "" &
                              !is.na(chromosome_name),
         chromosome_name, probe_chrom))


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
  filter(is.na(Entrez.ID) | ENSEMBL.Gene.ID == "") %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)

remaining_missing_entrez <- remaining %>% filter(is.na(Entrez.ID) & ENSEMBL.Gene.ID != "")

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

ensembl_matched_remaining <- remaining_missing_entrez %>%
  inner_join(entrez_results, by=c("ENSEMBL.Gene.ID" = "ensembl_id"), suffix=c(".remaining", ".entrez"))

# Correcting mismatched Gene.Symbol.remaining vs Gene.Symbol.entrez
# When mismatched, Gene.Symbol.entrez look like real gene symbols while Gene.Symbol.remaining are accession number
# (or in one weird case, a timestamp...?)
# Overwriting Gene.Symbol.remaining with Gene.Symbol.entrez with 

ensembl_matched_remaining <- ensembl_matched_remaining %>% 
  mutate(Gene.Symbol.remaining = ifelse(Gene.Symbol.remaining != Gene.Symbol.entrez,
                                        Gene.Symbol.entrez, Gene.Symbol.remaining))

# Correcting mismatched prob_chrom and entrez_chromosome
# (there's only 3, and probe_chrom is blank. Overwriting with entrez_chromosome)

ensembl_matched_remaining <- ensembl_matched_remaining %>% 
  mutate(probe_chrom = ifelse(probe_chrom != entrez_chromosome,
                              entrez_chromosome, probe_chrom))

# error checking

stopifnot(all(ensembl_matched_remaining$Gene.Symbol.remaining == ensembl_matched_remaining$Gene.Symbol.entrez))
stopifnot(all(ensembl_matched_remaining$probe_chrom == ensembl_matched_remaining$entrez_chromosome))

ensembl_matched_remaining$Entrez.ID <- as.integer(ensembl_matched_remaining$Entrez.ID.entrez)

output_manifest3 <- ensembl_matched_remaining %>%
  dplyr::select(Probe.ID, Gene.Symbol=Gene.Symbol.remaining, ENSEMBL.Gene.ID, Entrez.ID, Probe.Sequence) %>%
  bind_rows(output_manifest2)

check_data_unique(output_manifest3)

# Code useful for error checking
#ensembl_matched_remaining[ensembl_matched_remaining$Gene.Symbol.remaining != ensembl_matched_remaining$Gene.Symbol.entrez,] %>% select(c(Probe.ID, Gene.Symbol.remaining, Gene.Symbol.entrez))
#ensembl_matched_remaining[ensembl_matched_remaining$probe_chrom != ensembl_matched_remaining$entrez_chromosome,]

#####################################################################################################################

remaining_missing_ensembl <- remaining %>% filter(!is.na(Entrez.ID) & ENSEMBL.Gene.ID == "")

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

entrez_matched <- remaining_missing_ensembl %>%
  inner_join(ensembl_results, by=c("Entrez.ID"))


# error checking

# Two of the probe rows are duplicated, and I can't see any difference between them.
# collapse duplicates

entrez_matched <- entrez_matched %>% dplyr::distinct()

#Fix non-matching Gene.Symbol.X and Gene.Symbol.y
entrez_matched$Gene.Symbol.x[entrez_matched$Probe.ID == "28402"] <- "DPH5-DT"

#Fix non-matching probe_chrom vs entrez_chromosome
#entrez_chromosome are numbers, probe_chrom are scaffolds (ex. CHR_HSCHR6_MHC_COX_CTG1)
#going with the numbers

entrez_matched <- entrez_matched %>% 
  mutate(probe_chrom = ifelse(probe_chrom != entrez_chromosome,
                              entrez_chromosome, probe_chrom))

# error checking
stopifnot(all(entrez_matched$Gene.Symbol.x == entrez_matched$Gene.Symbol.y))
stopifnot(all(entrez_matched$probe_chrom == entrez_matched$entrez_chromosome))

entrez_matched$ensembl_id <- as.character(entrez_matched$ensembl_id)

output_manifest4 <- entrez_matched %>%
  dplyr::select(Probe.ID, Gene.Symbol=Gene.Symbol.x, ENSEMBL.Gene.ID=ensembl_id, Entrez.ID, Probe.Sequence) %>%
  bind_rows(output_manifest3)

check_data_unique(output_manifest4)

#####################################################################################################################

# Remaining to fix: have either ensembl ID OR entrez ID, but not both

remaining2 <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest4) %>%
  filter(is.na(Entrez.ID) | is.na(ENSEMBL.Gene.ID) | ENSEMBL.Gene.ID == "" ) %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)

#####################################################################################################################
# Count missing info now, before we get too fancy

# Missing enseml
remaining2 %>% filter(ENSEMBL.Gene.ID == "" | ENSEMBL.Gene.ID == "NA") %>% 
  summarise(n = n())

# Missing Entrez
remaining2 %>% filter(is.na(Entrez.ID)) %>% 
  summarise(n = n())

# Missing both Ensembl IDs and Entrez IDs?
remaining2 %>% filter((ENSEMBL.Gene.ID == "") & is.na(Entrez.ID)) %>% 
  summarise(n = n())

#####################################################################################################################

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
#/home/bradford/storage/ucsc-utils/blat/gfServer start 127.0.0.1 1234 -stepSize=5 -canStop /home/bradford/storage/dbs/GRCh38/GRCh38_latest_rna.2bit

fasta_file <- "/tmp/probe_sequences.fa"

write.fasta(as.list(remaining2$Probe.Sequence), remaining2$Probe.Name, fasta_file, open = "w", nbchar = 60, as.string = TRUE)

alignment_file <- "/tmp/out.psl"

#system(paste("~/bin/gfClient -minScore=20 -minIdentity=0 127.0.0.1 1234 ~/Documents/HC_EHSRB/data/genomes/",fasta_file,alignment_file))
system(paste("/home/bradford/storage/ucsc-utils/blat/gfClient -minScore=20 -minIdentity=0 127.0.0.1 1234 /home/bradford/storage/dbs/GRCh38/",fasta_file,alignment_file))

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
  dplyr::rename("Probe.Name" = "Q_name") 

# Rows with ensemble_gene_ids have non-matching external_gene_name and Gene.Symbol
#Some of the gene.symbols are just accessions, or uncharacterized loci. 
#When I search the ensembls, the results match the external_gene_name, 
# so let's overwrite the Gene.Symbol with that.

matched_by_sequence <- matched_by_sequence %>% 
  mutate(Gene.Symbol = ifelse(Gene.Symbol != external_gene_name & !is.na(external_gene_name),
                              external_gene_name, Gene.Symbol)) %>%
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

good_results <- good_results %>% 
  mutate(ENSEMBL.Gene.ID = ensembl_gene_id.seq)

output_manifest5 <- good_results %>%
  dplyr::select(Probe.ID, Gene.Symbol, ENSEMBL.Gene.ID, Entrez.ID, Probe.Sequence) %>%
  unique() %>%
  bind_rows(output_manifest4)

check_data_unique(output_manifest5)


#####################################################################################################################

# Remaining missing
remaining3 <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest5) %>%
  filter(is.na(Entrez.ID) | is.na(ENSEMBL.Gene.ID) | ENSEMBL.Gene.ID == "") %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)


# Check standardized hwt2.0 - are same probes missing info there?

human_wt_2_0 <- read.csv("/home/bradford/storage/dbs/biospyder/Human_whole_transcriptome_2_0/Human_Whole_Transcriptome_2.0_standardized.csv")

hwt_2_0_matches <- human_wt_2_0 %>%
                dplyr::select(Probe_Name,Ensembl_Gene_ID,Entrez_ID) %>%
                dplyr::filter(Probe_Name %in% remaining3$Probe.Name)

rem3_2.0_matches <- remaining3 %>%
  inner_join(hwt_2_0_matches, by=c("Probe.Name" = "Probe_Name"))

# Take ensembl info from 2.0 if Entrez IDs match,
# and Entrez info from 2.0 if ensembls match
rem3_fixed <- rem3_2.0_matches %>% 
  mutate(Entrez.ID = ifelse(ENSEMBL.Gene.ID == Ensembl_Gene_ID,
                            Entrez_ID, Entrez.ID)) %>%
  mutate(ENSEMBL.Gene.ID = ifelse(Entrez.ID == Entrez_ID,
                                  Ensembl_Gene_ID, ENSEMBL.Gene.ID)) %>%
  filter(!is.na(Entrez.ID) & !is.na(ENSEMBL.Gene.ID) & ENSEMBL.Gene.ID != "") %>%
  dplyr::select(-c(Ensembl_Gene_ID, Entrez_ID))

output_manifest6 <- rem3_fixed %>%
  dplyr::select(Probe.ID, Gene.Symbol, ENSEMBL.Gene.ID, Entrez.ID, Probe.Sequence) %>%
  unique() %>%
  bind_rows(output_manifest5)


# error checking
check_data_unique(output_manifest6)

##########################################################################################################
# Count unfixable probes

remaining4 <- manifest %>%
  dplyr::select(Probe.ID) %>%
  left_join(output_manifest6) %>%
  filter(is.na(Entrez.ID) | is.na(ENSEMBL.Gene.ID) | ENSEMBL.Gene.ID == "") %>%
  dplyr::select(Probe.ID) %>%
  inner_join(manifest)

# Missing enseml only
remaining4 %>% filter(ENSEMBL.Gene.ID == "" | ENSEMBL.Gene.ID == "NA") %>% 
  filter(!is.na(Entrez.ID)) %>%
  summarise(n = n())

# Missing Entrez only
remaining4 %>% filter(is.na(Entrez.ID)) %>% 
  filter(ENSEMBL.Gene.ID != "" | is.na(ENSEMBL.Gene.ID)) %>%
  summarise(n = n())

# Missing both Ensembl IDs and Entrez IDs?
remaining4 %>% filter((ENSEMBL.Gene.ID == "") & is.na(Entrez.ID)) %>% 
  summarise(n = n())

##########################################################################################################
#Final manifest

output_manifest_final <- remaining4 %>%
  dplyr::select(Probe.ID, Gene.Symbol, ENSEMBL.Gene.ID, Entrez.ID, Probe.Sequence) %>%
  bind_rows(output_manifest6 %>% dplyr::select(Probe.ID, Gene.Symbol, ENSEMBL.Gene.ID, Entrez.ID, Probe.Sequence)) %>%
  left_join(manifest %>% dplyr::select(Probe.ID, ALIGNED_ENSEMBL_TRANSCRIPTS, Probe.Name))

check_data_unique(output_manifest_final)

output_manifest_formatted <- output_manifest_final %>%
  mutate(Probe_ID = as.integer(Probe.ID)) %>%
  dplyr::select(Probe_ID, Probe_Name = Probe.Name, Gene_Symbol = Gene.Symbol, Ensembl_Gene_ID = ENSEMBL.Gene.ID, Entrez_ID = Entrez.ID, Probe_Sequence = Probe.Sequence, Transcripts_Targeted = ALIGNED_ENSEMBL_TRANSCRIPTS) %>%
  arrange(Probe_ID)
  
stopifnot(all(output_manifest_formatted$Probe_ID %in% manifest$Probe.ID))
stopifnot(all(output_manifest_formatted$Probe_Name %in% manifest$Probe.Name))
stopifnot(all(manifest$Probe.ID %in% output_manifest_formatted$Probe_ID))
stopifnot(all(manifest$Probe.Name %in% output_manifest_formatted$Probe_Name))

  
output_file <- file.path(output_directory, output_filename)

write.csv(output_manifest_formatted, output_file, row.names=F)


