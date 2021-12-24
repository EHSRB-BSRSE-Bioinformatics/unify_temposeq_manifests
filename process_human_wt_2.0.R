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
  n_unique_probeid <- x %>% dplyr::select(PROBE_ID) %>% unique() %>% nrow()
  stopifnot(n==n_unique_probeid)
}

#####################################################################################################################

# set up ensembl mart
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# read in manifest file
original_manifests_directory <- "original_manifests"

human_wt_2.0_file <- file.path(original_manifests_directory, "191004_Human_Whole_Transcriptome_2.0_Manifest.csv")

manifest <- read.csv(human_wt_2.0_file)


# First, do some data conversion

manifest[manifest == 'NULL'] <- NA
manifest$ENTREZ_ID <- as.integer(manifest$ENTREZ_ID)
manifest <- manifest %>% separate(LIGATED_COORDINATE,sep='[:-]',into=c('probe_chrom','probe_start','probe_end'))

output_manifest0 <- manifest %>%
  filter(!is.na(ENTREZ_ID) & !is.na(ENSEMBL_GENE_ID)) %>%
  dplyr::select(PROBE_ID, GENE_SYMBOL, ENSEMBL_GENE_ID, ENTREZ_ID, LIGATED_SEQUENCE, ALIGNED_REFSEQ_TRANSCRIPTS)

check_data_unique(output_manifest0)

output_directory <- "output_manifests"
output_filename <- "Human_Whole_Transcriptome_2.0_standardized.csv"

#####################################################################################################################

# Start with rows missing Ensembl IDs, with Entrez IDs

missing_ensembl_with_entrez <- manifest %>%
  filter(is.na(ENSEMBL_GENE_ID) & !is.na(ENTREZ_ID)) 

ensembl_for_missing <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'),
                             filters = 'entrezgene_id',
                             values = missing_ensembl_with_entrez$ENTREZ_ID,
                             mart = ensembl)

entrez_matched <- missing_ensembl_with_entrez %>%
  inner_join(ensembl_for_missing, by=c("ENTREZ_ID" = "entrezgene_id")) %>%
  filter(probe_chrom == chromosome_name)

# manually fix some problems 
entrez_matched <- entrez_matched %>% filter(!(PROBE_ID == 27178 & external_gene_name == ""))
entrez_matched <- entrez_matched %>% filter(!(PROBE_ID == 13971 & external_gene_name == "ARL17B"))
entrez_matched <- entrez_matched %>% filter(!(PROBE_ID == 26852 & external_gene_name == "ARL17B"))
entrez_matched <- entrez_matched %>% filter(!(PROBE_ID == 26853 & external_gene_name == "ARL17B"))
entrez_matched <- entrez_matched %>% filter(!(PROBE_ID == 26854 & external_gene_name == "ARL17B"))
entrez_matched$GENE_SYMBOL[entrez_matched$PROBE_ID == "23267"] <- "MTRFR"
entrez_matched$GENE_SYMBOL[entrez_matched$PROBE_ID == "13931"] <- "TMEM51-AS2"
entrez_matched$external_gene_name[entrez_matched$PROBE_ID == "25035"] <- "LOC100652758"
entrez_matched$external_gene_name[entrez_matched$PROBE_ID == "28535"] <- "LOC101928093"
entrez_matched$external_gene_name[entrez_matched$PROBE_ID == "15219"] <- "LOC389895"
entrez_matched$GENE_SYMBOL[entrez_matched$PROBE_ID == "19140"] <- "LUZP6/MTPN"
entrez_matched$external_gene_name[entrez_matched$PROBE_ID == "19140"] <- "LUZP6/MTPN"
entrez_matched$external_gene_name[entrez_matched$PROBE_ID == "28400"] <- "LOC100505501"

# error checking

stopifnot(all(entrez_matched$GENE_SYMBOL == entrez_matched$external_gene_name))
check_n_matches <- entrez_matched %>% group_by(PROBE_ID) %>% mutate(n_copies = n())
stopifnot(check_n_matches %>% filter(n_copies > 1) %>% nrow() == 0)

output_manifest1 <- entrez_matched %>%
  dplyr::select(PROBE_ID, GENE_SYMBOL, ENSEMBL_GENE_ID = ensembl_gene_id, ENTREZ_ID, LIGATED_SEQUENCE, ALIGNED_REFSEQ_TRANSCRIPTS) %>%
  bind_rows(output_manifest0)

check_data_unique(output_manifest1)

#####################################################################################################################

# now probes  with Ensembl IDs but missing Entrez IDs

missing_entrez_with_ensembl <- manifest %>%
  filter(is.na(ENTREZ_ID) & !is.na(ENSEMBL_GENE_ID)) 

entrez_for_missing <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'),
                            filters = 'ensembl_gene_id',
                            values = missing_entrez_with_ensembl %>% dplyr::select(ENSEMBL_GENE_ID),
                            mart = ensembl) %>% filter(!is.na(entrezgene_id))


ensembl_matched <- missing_entrez_with_ensembl %>%
  inner_join(entrez_for_missing, by=c("ENSEMBL_GENE_ID" = "ensembl_gene_id"))

# there seem to be some ensembl genes that match to multiple entrez IDs, but spot checking a few 
# suggests that the larger number is a read-through gene in most cases, so assume that
# the smaller entrez ID is the one we want
ensembl_matched <- ensembl_matched %>% group_by(ENSEMBL_GENE_ID) %>% slice_max(order_by = entrezgene_id, n=1) %>% ungroup()

# manually fix some problems (all "ACNNNNN" accessions for symbols--one has been assigned a gene symbol so
# we'll use that, the rest don't have ensembl symbols so we'll force it to use the accession numbers)
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "19698"] <- "TRIM67-AS1"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "5295"] <- "PP7080"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "19161"] <- "POLR1G"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "87532"] <- "MSANTD5"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "19732"] <- "MALT1-AS1"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "16581"] <- "LINC02694"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "27401"] <- "KIR3DS1"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "12290"] <- "FOXL3"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "10536"] <- "AC027097.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "14126"] <- "AC011043.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "14857"] <- "AC012236.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "15923"] <- "AC125421.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "17572"] <- "AC055839.2"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "17958"] <- "AC092118.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "17969"] <- "AC007342.3"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "18935"] <- "AL449403.2"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "20202"] <- "AC011525.2"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "23780"] <- "AC020928.2"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "24116"] <- "AC093249.6"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "24596"] <- "AC010255.3"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "25662"] <- "AL353753.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "27956"] <- "AC011525.2"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "28481"] <- "AC007666.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "29034"] <- "AC132008.2"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "34026"] <- "AC095032.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "34027"] <- "AC095032.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "88136"] <- "AC008687.4"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "88210"] <- "AP002495.2"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "88212"] <- "AL845331.2"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "88744"] <- "AL445989.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "88983"] <- "AL354761.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "89032"] <- "AC008397.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "89097"] <- "AL022312.1"
ensembl_matched$external_gene_name[ensembl_matched$PROBE_ID == "89608"] <- "BX276092.9"
ensembl_matched$chromosome_name[ensembl_matched$PROBE_ID == "27401"] <- "CHR_HSCHR19KIR_0019"

# error checking

stopifnot(all(ensembl_matched$GENE_SYMBOL == ensembl_matched$external_gene_name))
stopifnot(all(ensembl_matched$probe_chrom == ensembl_matched$chromosome_name))
check_n_matches <- ensembl_matched %>% group_by(PROBE_ID) %>% mutate(n_copies = n())
stopifnot(check_n_matches %>% filter(n_copies > 1) %>% nrow() == 0)

output_manifest2 <- ensembl_matched %>%
  dplyr::select(PROBE_ID, GENE_SYMBOL, ENSEMBL_GENE_ID, ENTREZ_ID = entrezgene_id, LIGATED_SEQUENCE, ALIGNED_REFSEQ_TRANSCRIPTS) %>%
  bind_rows(output_manifest1)

check_data_unique(output_manifest2)

#####################################################################################################################

# Next we are going to search NCBI for the genes that are still missing either an Ensembl or Entrez ID 

remaining <- manifest %>%
  dplyr::select(PROBE_ID) %>%
  left_join(output_manifest2) %>%
  filter(is.na(ENTREZ_ID) | is.na(ENSEMBL_GENE_ID)) %>%
  dplyr::select(PROBE_ID) %>%
  inner_join(manifest)

remaining_missing_entrez <- remaining %>% filter(is.na(ENTREZ_ID) & !is.na(ENSEMBL_GENE_ID))

get_gene_info_from_ncbi_using_ensembl <- function(ensembl_id){
  search_result <- entrez_search(db="gene", term=ensembl_id)
  entrez_id <- search_result[[1]]
  if(length(entrez_id) == 1){
    # extract gene info
    summary <- entrez_summary(id=entrez_id, db="gene")
    gene_symbol <- summary$name
    entrez_chromosome <- summary$chromosome
    # create data frame of information
    gene.df <- data.frame(ensembl_id, entrez_id, gene_symbol, entrez_chromosome) 
    gene.df$ensembl_id <- ensembl_id
    gene.df$entrez_id <- ifelse(!is.null(entrez_id), entrez_id, NA)
    gene.df$gene_symbol <- ifelse(!is.null(gene_symbol), gene_symbol, NA)
    gene.df$entrez_chromosome <- ifelse(!is.null(entrez_chromosome), entrez_chromosome, NA)
    return(gene.df)
  }
}

entrez_results <- bind_rows(lapply(remaining_missing_entrez$ENSEMBL_GENE_ID, get_gene_info_from_ncbi_using_ensembl)) %>% unique()

ensembl_matched <- remaining_missing_entrez %>%
  inner_join(entrez_results, by=c("ENSEMBL_GENE_ID" = "ensembl_id"))

# manual corrections

ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "3829"] <- "LOC102606465"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "3844"] <- "LOC728554"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "10457"] <- "DELEC1" # excel strikes again
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "12887"] <- "LOC100130370"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "15246"] <- "LOC728554"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "22102"] <- "OR7E85BP"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "22573"] <- "RPH3AL-AS1"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "26692"] <- "LOC389602"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "26933"] <- "SNHG28"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "87489"] <- "TRIM75"
ensembl_matched$gene_symbol[ensembl_matched$PROBE_ID == "87654"] <- "TRBV6-2" # according to hgnc, TRBV6-2 is on the reference chr7 and TRBV6-3 is on the alternate
ensembl_matched$entrez_id[ensembl_matched$PROBE_ID == "87654"] <- "28605" 
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "87691"] <- "TRIM51G"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "87748"] <- "FAM90A10"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "87777"] <- "FAM90A22"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "88043"] <- "OR1R1P"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "88151"] <- "PTP4A1"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "88172"] <- "LOC93622"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "88257"] <- "HNRNPA1L3"
ensembl_matched$GENE_SYMBOL[ensembl_matched$PROBE_ID == "88998"] <- "MCTS2"
ensembl_matched$gene_symbol[ensembl_matched$PROBE_ID == "90891"] <- "IGHV4-31"
ensembl_matched$entrez_id[ensembl_matched$PROBE_ID == "90891"] <- "28396"


# error checking

stopifnot(all(ensembl_matched$GENE_SYMBOL == ensembl_matched$gene_symbol))
stopifnot(all(ensembl_matched$probe_chrom == ensembl_matched$entrez_chromosome))

ensembl_matched$entrez_id <- as.integer(ensembl_matched$entrez_id)

output_manifest3 <- ensembl_matched %>%
  dplyr::select(PROBE_ID, GENE_SYMBOL, ENSEMBL_GENE_ID, ENTREZ_ID = entrez_id, LIGATED_SEQUENCE, ALIGNED_REFSEQ_TRANSCRIPTS) %>%
  bind_rows(output_manifest2)

check_data_unique(output_manifest3)

###

remaining_missing_ensembl <- remaining %>% filter(!is.na(ENTREZ_ID) & is.na(ENSEMBL_GENE_ID))

get_gene_info_from_ncbi_using_entrez <- function(entrez_id){
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

ensembl_results <- bind_rows(lapply(remaining_missing_ensembl$ENTREZ_ID, get_gene_info_from_ncbi_using_entrez)) %>% filter(ensembl_id!="NULL")

# two measly genes, argh
entrez_matched <- remaining_missing_ensembl %>%
  inner_join(ensembl_results, by=c("ENTREZ_ID" = "entrez_id"))
# error checking

stopifnot(all(entrez_matched$GENE_SYMBOL == entrez_matched$gene_symbol))
stopifnot(all(entrez_matched$probe_chrom == entrez_matched$entrez_chromosome))

entrez_matched$ensembl_id <- as.character(entrez_matched$ensembl_id)

output_manifest4 <- entrez_matched %>%
  dplyr::select(PROBE_ID, GENE_SYMBOL, ENSEMBL_GENE_ID=ensembl_id, ENTREZ_ID, LIGATED_SEQUENCE, ALIGNED_REFSEQ_TRANSCRIPTS) %>%
  bind_rows(output_manifest3)

check_data_unique(output_manifest4)

#####################################################################################################################

remaining2 <- manifest %>%
  dplyr::select(PROBE_ID) %>%
  left_join(output_manifest4) %>%
  filter(is.na(ENTREZ_ID) | is.na(ENSEMBL_GENE_ID)) %>%
  dplyr::select(PROBE_ID) %>%
  inner_join(manifest)

# OK, before we start blasting things let's try with the gene names

by_gene_name <- getBM(attributes = c('entrezgene_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'),
                      filters = 'external_gene_name',
                      values = remaining2$GENE_SYMBOL,
                      mart = ensembl)

by_gene_name_singles <- by_gene_name %>%
  group_by(external_gene_name) %>%
  mutate(num_results = n())
by_gene_name_singles <- by_gene_name_singles %>% filter(num_results == 1)

matched_by_gene_name <- remaining2 %>%
  inner_join(by_gene_name_singles, by=c('GENE_SYMBOL' = 'external_gene_name'))
#####################################################################################################################

# Next we are going to take the remaining probes and use BLAT against a set of human mRNAs
# Start the server outside of R using
#  ~/bin/gfServer start 127.0.0.1 1234 -stepSize=5 /home/katecook/Documents/HC_EHSRB/data/genomes/GRCh38_latest_rna.2bit

fasta_file <- "/tmp/probe_sequences.fa"

write.fasta(as.list(remaining2$LIGATED_SEQUENCE), remaining2$PROBE_NAME, fasta_file, open = "w", nbchar = 60, as.string = TRUE)

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
  inner_join(remaining2, by=c("Q_name" = "PROBE_NAME")) %>% 
  dplyr::rename("PROBE_NAME" = "Q_name") %>%
  filter(is.na(external_gene_name) | external_gene_name == GENE_SYMBOL) %>%
  filter(is.na(chromosome_name) | is.na(probe_chrom) | chromosome_name == probe_chrom) %>%
  filter(!is.na(ensembl_gene_id))


# require that the matched by gene name and matched by sequence IDs are identical if they both exist

good_results <- matched_by_sequence %>%
  full_join(matched_by_gene_name, by = c("PROBE_NAME", "PROBE_ID", "GENE_SYMBOL", "LIGATED_SEQUENCE", "ENSEMBL_GENE_ID", "ALIGNED_ENSEMBL_TRANSCRIPTS", "ENTREZ_ID", "ALIGNED_REFSEQ_TRANSCRIPTS"), suffix=c(".seq", ".symbol")) %>%
  filter(is.na(ensembl_gene_id.seq) | is.na(ensembl_gene_id.symbol) | ensembl_gene_id.seq == ensembl_gene_id.symbol)

good_results_symbol <- good_results %>%
  filter(!is.na(ensembl_gene_id.symbol)) %>%
  inner_join(remaining2, by=c("PROBE_NAME")) 

good_results_symbol$entrezgene_id.symbol[!is.na(good_results_symbol$ENTREZ_ID.x)] <- good_results_symbol$ENTREZ_ID.x[!is.na(good_results_symbol$ENTREZ_ID.x)]

# fill in some entrez results manually
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000211592"] <- 3514
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000262619"] <- 100996930
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000224411"] <- 3324
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000211673"] <- 28809
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000224078"] <- 104472715
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000266920"] <- 69
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000197617"] <- 317705
good_results_symbol$entrezgene_id.symbol[good_results_symbol$PROBE_ID.x == 22770] <- 317705
good_results_symbol$entrezgene_id.symbol[good_results_symbol$PROBE_ID.x == 25711] <- 729648
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000284283"] <- 391746
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000224032"] <- 114915
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000198225"] <- 642489
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000236737"] <- 729428
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000283586"] <- 100132413
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000259490"] <- 28318
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000233732"] <- 28306
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000270467"] <- 28304
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000271178"] <- 100287372
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000282651"] <- 28386
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000231292"] <- 28862
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000211676"] <- 28832
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000211678"] <- 28831
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000211643"] <- 28779
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000181355"] <- 266553
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000232535"] <- 79289
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000182257"] <- 55267
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000253649"] <- 346702
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000197364"] <- 645922
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000283740"] <- 112488746
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000249156"] <- 112488739
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000283776"] <- 112488747
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000250782"] <- 112488740
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000284439"] <- 646103
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000284283"] <- 391746
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000284042"] <- 391747
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000283967"] <- 112488737
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000283528"] <- 100533814
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000211808"] <- 28679
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000211829"] <- 28526
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000211694"] <- 6984
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000211693"] <- 6985
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000250913"] <- 101241878
good_results_symbol$entrezgene_id.symbol[good_results_symbol$ensembl_gene_id.symbol == "ENSG00000250091"] <- 642797

output_manifest5 <- good_results_symbol %>%
  dplyr::select(PROBE_ID = PROBE_ID.x, GENE_SYMBOL = GENE_SYMBOL.x, ENSEMBL_GENE_ID=ensembl_gene_id.symbol, ENTREZ_ID = entrezgene_id.symbol, LIGATED_SEQUENCE = LIGATED_SEQUENCE.x, ALIGNED_REFSEQ_TRANSCRIPTS = ALIGNED_REFSEQ_TRANSCRIPTS.x) %>%
  unique() %>%
  bind_rows(output_manifest4)

check_data_unique(output_manifest5)


good_results_seq <- good_results %>%
  filter(!(PROBE_NAME %in% good_results_symbol$PROBE_NAME)) %>%
  inner_join(remaining2, by=c("LIGATED_SEQUENCE"))

output_manifest6 <- good_results_seq %>%
  dplyr::select(PROBE_ID = PROBE_ID.x, GENE_SYMBOL = GENE_SYMBOL.x, ENSEMBL_GENE_ID=ensembl_gene_id.seq, ENTREZ_ID = entrezgene_id.seq, LIGATED_SEQUENCE, ALIGNED_REFSEQ_TRANSCRIPTS = ALIGNED_REFSEQ_TRANSCRIPTS.x) %>%
  bind_rows(output_manifest5)

check_data_unique(output_manifest6)

#####################################################################################################################

# fill in the rest manually

remaining3 <- manifest %>%
  dplyr::select(PROBE_ID) %>%
  left_join(output_manifest6) %>%
  filter(is.na(ENTREZ_ID) | is.na(ENSEMBL_GENE_ID)) %>%
  dplyr::select(PROBE_ID) %>%
  inner_join(manifest)

tofix <-good_results_symbol %>% filter(is.na(entrezgene_id.symbol)) %>% dplyr::select(PROBE_NAME) %>% unique() %>% inner_join(remaining3, by=c("PROBE_NAME"))

# update entrez IDs
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000282627"] <- 28424
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000282691"] <- 28392
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000282758"] <- 28941
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000282801"] <- 28299
remaining3$ENTREZ_ID[remaining3$PROBE_ID == 17995] <- 5414 # not sure what went wrong here but there was a retired ensembl ID. Not sure I should update the gene symbol?
remaining3$ENSEMBL_GENE_ID[remaining3$PROBE_ID == 17995] <- "ENSG00000108387"
remaining3$ENTREZ_ID[remaining3$PROBE_ID == 26910] <- 5414 # ditto. HGNC symbol is also withdrawn
remaining3$ENSEMBL_GENE_ID[remaining3$PROBE_ID == 26910] <- "ENSG00000078114"
remaining3$ENTREZ_ID[remaining3$PROBE_ID == 26911] <- 5414 # ditto. HGNC symbol is also withdrawn, update that?
remaining3$ENSEMBL_GENE_ID[remaining3$PROBE_ID == 26911] <- "ENSG00000078114"
remaining3$ENTREZ_ID[remaining3$PROBE_ID == 26818] <- 445815 # gene symbol should be updated maybe?
remaining3$ENSEMBL_GENE_ID[remaining3$PROBE_ID == 26818] <- "ENSG00000157654"
remaining3$ENTREZ_ID[remaining3$PROBE_ID == 26819] <- 445815 
remaining3$ENSEMBL_GENE_ID[remaining3$PROBE_ID == 26819] <- "ENSG00000157654"
remaining3$ENTREZ_ID[remaining3$PROBE_ID == 17558] <- 445815 
remaining3$ENSEMBL_GENE_ID[remaining3$PROBE_ID == 17558] <- "ENSG00000157654"
remaining3$ENTREZ_ID[remaining3$PROBE_ID == 26023] <- 162699 
remaining3$ENSEMBL_GENE_ID[remaining3$PROBE_ID == 26023] <- "ENSG00000288631"
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000226644"] <- 388780
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000254658"] <- 81169
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000255552"] <- 79136
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000261253"] <- 100287036
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000266920"] <- 69
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000267698"] <- 101927572
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000283967"] <- 112488737
remaining3$ENTREZ_ID[remaining3$ENSEMBL_GENE_ID == "ENSG00000284042"] <- 391747

# update ensembl IDs
remaining3$ENSEMBL_GENE_ID[remaining3$ENTREZ_ID == 11042] <- "ENSG00000196302"
remaining3$ENSEMBL_GENE_ID[remaining3$ENTREZ_ID == 100506621] <- "ENSG00000091986"
remaining3$ENTREZ_ID[remaining3$ENTREZ_ID == 100506621] <- 151887 # also update the entrez ID
remaining3$ENSEMBL_GENE_ID[remaining3$ENTREZ_ID == 102606465] <- "ENSG00000233184"
remaining3$ENSEMBL_GENE_ID[remaining3$ENTREZ_ID == 339862] <- "ENSG00000274840"
remaining3$ENSEMBL_GENE_ID[remaining3$ENTREZ_ID == 100129083] <- "ENSG00000255441"
remaining3$ENTREZ_ID[remaining3$PROBE_ID == 28403] <- 151887
remaining3$ENSEMBL_GENE_ID[remaining3$PROBE_ID == 28403] <- "ENSG00000091986"

# both are missing, search using gene name
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "ANKRD20A9P"] <- "ENSG00000206192"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "ANKRD20A9P"] <- 284232
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "IGHV3-23"] <- "ENSG00000211949"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "IGHV3-23"] <- 28442
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "IGKV1-5"] <- "ENSG00000243466"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "IGKV1-5"] <- 28299
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "IGHG2"] <- "ENSG00000211893"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "IGHG2"] <- 3501
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "TRBC2"] <- "ENSG00000211772"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "TRBC2"] <- 28638
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "BAGE5"] <- 85316
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "FGF11"] <- "ENSG00000161958"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "FGF11"] <- 2256
remaining3$ENSEMBL_GENE_ID[remaining3$PROBE_ID == "2867"] <- "ENSG00000206192"
remaining3$ENTREZ_ID[remaining3$PROBE_ID == "2867"] <- 3021
remaining3$GENE_SYMBOL[remaining3$PROBE_ID == "2867"] <- "H3-3B"
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "LOC100289333"] <- "ENSG00000234773"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "LOC100289333"] <- 100289333
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "LOC100507630"] <- 100507630
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "PRR21"] <- "ENSG00000260610"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "PRR21"] <- 643905
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "PRR33"] <- "ENSG00000283787"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "PRR33"] <- 102724536
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "TBC1D3F"] <- "ENSG00000275954"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "TBC1D3F"] <- 84218
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "TRBC2"] <- "ENSG00000211772"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "TRBC2"] <- 28638
remaining3$ENSEMBL_GENE_ID[remaining3$GENE_SYMBOL == "XR"] <- "ENSG00000197617"
remaining3$ENTREZ_ID[remaining3$GENE_SYMBOL == "XR"] <- 317705

output_manifest_final <- remaining3 %>%
  dplyr::select(PROBE_ID, GENE_SYMBOL, ENSEMBL_GENE_ID, ENTREZ_ID, LIGATED_SEQUENCE, ALIGNED_REFSEQ_TRANSCRIPTS) %>%
  bind_rows(output_manifest6)

check_data_unique(output_manifest_final)


output_file <- file.path(output_directory, output_filename)

write.csv(output_manifest_final, output_file, row.names=F)


