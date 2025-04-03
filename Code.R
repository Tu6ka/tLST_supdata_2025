# To study tLSTs in bacterial populations all complete genomes available on NCBI datasets were downloaded https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/ NCBI Datasets. 
# That was done on 05.02.2024.

# This command was used to download all complete bacterial genomes
$ datasets download genome taxon bacteria --assembly-level complete

# To download all related metadata for complete bacterial genomes
$ datasets summary genome taxon bacteria --assembly-level complete --as-json-lines | dataformat tsv genome > full_metadata.tsv

# To generate core genome allignments:
  # Bakta was used to generate genome annotations https://github.com/oschwengers/bakta
  $ for f in *.fasta ; do bakta --db /data/dmitli/tools/baktaDB/db/ --verbose --output ../bakta/${f//.fasta} --threads 20 $f ; done

  # Panaroo https://github.com/gtonkinhill/panaroo was used for core genome allignment, on gff3 files from bakta
  $ panaroo -i *.gff3 -o ../panaroo_results --clean-mode strict -a core --aligner mafft --core_threshold 0.98 --remove-invalid-genes

  # Phylogenetic trees were generated using IQtree2 http://www.iqtree.org/
  $ iqtree -s core_gene_alignment_filtered.aln -m MFP -nt AUTO

  # SNP distances measured using snp-dists https://github.com/tseemann/snp-dists
  $ snp-dists core_gene_alignment_filtered.aln >snp_dists.tsv


# Following code was used for processing of metadata####
library(dplyr)
library(tidyr)
library(stringr)
library(reactable)
library(tidyverse)
library(taxizedb)
library(pheatmap)

# read tsv file of metadata, full metadata file can be provi
metadata<-read.table(file = "full_metadata.tsv",
                     sep = '\t', quote = "\"", fill = T, header = T)

# put rownames in a column
metadata <- rowid_to_column(metadata, 'rowid')

# find problematic rows
problem_lines <- which(!grepl("GCA|GCF",metadata$`Assembly.Accession`))
problem_lines <- sort(union(problem_lines,problem_lines-1))

# separate them and manualy fix few problematic cases
dataframe <- metadata[problem_lines,]
for(line in c(seq(1,88, by = 2),seq(109,148,by = 2))){
  dataframe[line,59:148]<-dataframe[line+1,3:89]
}
for(line in seq(89,108, by = 2)){
  dataframe[line, 72] <- paste(dataframe[line,72],
                               dataframe[line,73],
                               dataframe[line,115])
  dataframe[line,73:103]<-dataframe[line,118:148]
  dataframe[line,104:148]<-dataframe[line+1,2:46]
}
for(line in seq(149,344, by=7)){
  dataframe[line,57]<-paste(dataframe[line,57],
                            dataframe[line+1,1],
                            dataframe[line+2,1],
                            dataframe[line+3,1],
                            dataframe[line+4,1],
                            dataframe[line+5,1],
                            dataframe[line+6,1])
  dataframe[line,58:144]<-dataframe[line+6,2:88]
}
dataframe <- dataframe[c(seq(1, 148, by = 2), seq(149, 344, by = 7)),] 

#remove problem lines
metadata <- metadata %>% 
  filter(!c(rowid %in% problem_lines)) %>% 
  rbind(dataframe)

# remove all ; for future
metadata <- metadata %>% mutate(across(c(Assembly.BioSample.Attribute.Name,Assembly.BioSample.Attribute.Value), ~ str_replace_all(., ";", ".")))


# Define a function to collapse unique values into a single string
collapse_unique <- function(x) {
  if (length(unique(x)) > 1) {
    paste(as.character(unique(x)), collapse = ';')
  } else {
    as.character(unique(x))
  }
}

# Group by Assembly.Accession's ID and collapse the unique values for each column
collapsed_metadata <- metadata %>% 
  select(-rowid) %>% 
  mutate(Assembly.Accession.bkp = Assembly.Accession) %>% 
  separate(Assembly.Accession.bkp, into = c('DB', 'ID', 'Version'), sep = "[._]") %>% 
  mutate(Assembly.BioSample.Attribute.Name.Value = paste(Assembly.BioSample.Attribute.Name,Assembly.BioSample.Attribute.Value, 
                                                         sep = '$')) %>% 
  select(-c(Assembly.BioSample.Attribute.Name,Assembly.BioSample.Attribute.Value)) %>% 
  group_by(ID) %>% 
  summarize(across(everything(), collapse_unique)) %>%
  ungroup() %>% mutate(Assembly.Accession = 
                         unlist(ifelse(grepl(';', Assembly.Accession),
                                       sapply(strsplit(Assembly.Accession, ';'), 
                                              function(x) {unlist(x)[grep("GCF", unlist(x))]}),
                                       Assembly.Accession))) %>% 
  select( -DB, -ID, -Version)

metadata <- collapsed_metadata

# Write a tsv
write_tsv(collapsed_metadata, "cleaned_metadata_210624.tsv", quote = "needed")

meta <- metadata[c("Assembly.Accession", "Assembly.BioSample.Attribute.Name.Value")]
names(meta) <- c("Assembly", "Attribute.Value")

# Separate data
expanded_data <- meta %>%
  separate_rows(`Attribute.Value`, sep = ";") 
expanded_data <- expanded_data %>%
  separate(`Attribute.Value`, into = c("Attribute", "Value"), sep = "\\$")

# fill empty with 'Empty' so there will be no complanes in future
expanded_data$Attribute[expanded_data$Attribute == ''] <- 'Empty'

# Pivot to wider format
wide_data <- expanded_data %>%
  pivot_wider(names_from = Attribute, 
              values_from = Value, 
              values_fn = list
  ) %>%
  select(Assembly, everything()) %>% 
  mutate(across(everything(), ~replace(., .=="NULL", ''))) %>% 
  select(-Empty)

# deal with some list export problems
wide_data<-as.data.frame(apply(wide_data,2,as.character))

# save as tsv which could be read later to save time
readr::write_tsv(wide_data, "simple_metadata_210624.tsv") # content of SupData 1A



# Code to process taxonomy ####
# Data frame with taxonomic IDs and assembly accession number

metadata <- read.table(file = "cleaned_metadata_210624.tsv", quote = "\"", sep = '\t', header = T)

# with few cases where ANI best hit is messed up or absent "fixed" by using BioSample Description it is not perfect but better than nothing
isolates <- data.frame(
  isolate = metadata$Assembly.Accession,
  taxid = metadata$ANI.Best.ANI.match.Organism,
  fix = metadata$Assembly.BioSample.Description.Organism.Name
) %>% mutate(taxid = case_when(
  taxid == "Bacteroides hominis (ex Liu et al. 2022)" ~ "Bacteroides hominis",
  taxid == "Rhodobaca barguzinensis" ~ "Roseinatronobacter bogoriensis",
  taxid == "[Clostridium] dakarense" ~ "Faecalimicrobium dakarense",
  taxid == "Oscillatoria nigro_viridis" ~ "Oscillatoria nigro-viridis",
  taxid == ";Bacillus velezensis" ~ "Bacillus velezensis",
  taxid == "[Muricauda] yonaguniensis" ~ "Muricauda yonaguniensis",
  taxid == ";Lactococcus cremoris" ~ "Lactococcus cremoris",
  taxid == "Komarekiella delphini_convector" ~ "Komarekiella delphini-convector",
  taxid == "" ~ fix,
  TRUE ~ taxid)) %>% 
  mutate(taxid = case_when(
    taxid == "Candidatus Moranbacteria bacterium" ~ "2045217",
    taxid == "3038980" ~ "Candidatus Lucifugimonas marina",
    TRUE ~ taxid)) %>% 
  select(-fix)

# Download the NCBI taxonomy database (need to be done just once)
db_download_ncbi()

# Function to retrieve taxonomic ranks for a given taxid
get_taxonomy <- function(taxid) {
  taxizedb::classification(taxid, db = 'ncbi')
}

# Get full taxonomy for every unique species name
genera <- as.data.frame(unique(isolates$taxid)) %>% 
  rename("taxid" = "unique(isolates$taxid)") %>%
  rowwise() %>%
  mutate(
    taxonomy = list(get_taxonomy(taxid)[[1]])
  ) %>% unnest(cols = taxonomy)

# Pivot to wider format and ensure unique column names
taxonomic_data_wide <- genera %>%
  select(-id) %>%
  filter(rank != "no rank", rank !="clade") %>% #two problematic ranks
  pivot_wider(
    names_from = rank,
    values_from = name,
    names_repair = "unique"
  )

# Select and rename columns of interest for clarity
taxonomic_data_wide <- taxonomic_data_wide %>%
  select(taxid, phylum, class, order, family, genus, species)

# and get it back to isolats
taxonomy <- isolates %>% inner_join(taxonomic_data_wide, by = "taxid")

# save as tsv for future uses
readr::write_tsv(taxonomy, "taxonomy_240624.tsv") #content of SupData 1B

  
  
  
# Code to choose members of tLST positive species used for tree in figure 2 ####
tLST <- read.delim(file = 'tLST.tab', sep = "\t") #tLST.tab contain same information as SupData 1C

# attach taxonomy
tLST <- tLST %>% inner_join(taxonomy, by =c('name' = 'isolate')) 

# set seed for reproducibility 
set.seed(123)

# extract assembly names of coverage >=80 and id >=90 to a separate list. this list was used to make the tree
unique_assemblies_species <- tLST %>% filter(coverage >= 80, id >=90) %>%
  group_by(species) %>%
  sample_n(1) %>%  # Sample one row per species
  ungroup() %>% # Ungroup the dataframe
  column_to_rownames("name")



# General version of code to generate maps for SupData 1E ####
tLST1 <- tLST %>% filter(GENE == "tLST1_CP025739.1")

# find only non single assembly hits and prepare them
tLST1_2 <- tLST1[duplicated(tLST1$X.FILE)|duplicated(tLST1$X.FILE,fromLast = T),] %>% 
  separate(COVERAGE, into = c('start', 'end', 'total'), sep = "[-/]") %>% 
  mutate(across(c(start,end,total), as.numeric))

# make empty dataframe
df1 <- data.frame(matrix(0, nrow = length(unique(tLST1_2$name)), ncol = max(tLST1_2$total)), 
                 row.names = unique(tLST1_2$name))
names(df1) <- 1:max(tLST1_2$total)

# go through every entry and +1 for every "bp", so if some region appears twice it'll give 2
for(f in 1:nrow(tLST1_2)){
df1[tLST1_2[f,"name"],as.numeric(tLST1_2[f,"start"]):as.numeric(
 tLST1_2[f,"end"])] <- df1[tLST1_2[f,"name"],as.numeric(tLST1_2[f,"start"]):as.numeric(tLST1_2[f,"end"])] + 1
}

tLST1 <- df1 %>%
  select(rev(names(df1)))


pheatmap(tLST1, 
         cluster_cols =  F, 
         cluster_rows = T,
         show_colnames = F,
         fontsize = 4,
         labels_row = tLST1_2$X.FILE[!duplicated(tLST1_2$X.FILE)]
)
# Code for ChatGPT to fix metadata classification ####
# Required packages
library(readr)
library(dplyr)
library(httr)
library(purrr)

# read metadata file from previous sections
metadata <- read.delim(file = "simple_metadata_210624.tsv", 
                       quote = "\"", 
                       sep = '\t', 
                       header = T)

# Load metadata and get rid of all isolates which do not have anything useful
metadata2process_all <- metadata %>% select(Assembly, isolation_source, host) %>% 
  mutate(host = case_when(grepl("missing",host)~'', grepl("Missing",host)~'',
                          grepl("Unknown",host)~'',
                          grepl("unknown",host)~'',
                          grepl("not provided",host)~'',
                          grepl("Not collected",host)~'',
                          grepl("not collected",host)~'',
                          grepl("not available: not collected",host)~'',
                          grepl("Not available",host)~'',
                          grepl("not available",host)~'',
                          grepl("Not Applicable",host)~'',
                          grepl("Not applicable",host)~'',
                          grepl("not applicable",host)~'',
                          grepl("N/A",host)~'',
                          grepl("Missing",host)~'',
                          grepl("missing",host)~'',
                          .default = host
  )) %>% 
  mutate(isolation_source = case_when(grepl("Unknow",isolation_source)~'',
                                      grepl("unknown",isolation_source)~'',
                                      grepl("Unknown",isolation_source)~'',
                                      grepl("UNknown",isolation_source)~'',
                                      grepl("unkown",isolation_source)~'',
                                      grepl("Unkown",isolation_source)~'',
                                      grepl("Not defined",isolation_source)~'',
                                      grepl("not determined",isolation_source)~'',
                                      grepl("not informed",isolation_source)~'',
                                      grepl("not isolated",isolation_source)~'',
                                      grepl("not known",isolation_source)~'',
                                      grepl("not provided",isolation_source)~'',
                                      grepl("not recorded",isolation_source)~'',
                                      grepl("Not Collected",isolation_source)~'',
                                      grepl("Not collected",isolation_source)~'',
                                      grepl("not collected",isolation_source)~'',
                                      grepl("not available: not collected",isolation_source)~'',
                                      grepl("Not available",isolation_source)~'',
                                      grepl("not available",isolation_source)~'',
                                      grepl("Not Applicable",isolation_source)~'',
                                      grepl("Not applicable",isolation_source)~'',
                                      grepl("not applicable",isolation_source)~'',
                                      grepl("not available: to be reported later",isolation_source)~'',
                                      grepl("N/A",isolation_source)~'',
                                      grepl("na",isolation_source)~'',
                                      grepl("none",isolation_source)~'',
                                      grepl("None",isolation_source)~'',
                                      grepl("na",isolation_source)~'',
                                      grepl("MISSING",isolation_source)~'',
                                      grepl("Missing",isolation_source)~'',
                                      grepl("missing",isolation_source)~'',
                                      .default = isolation_source
  )) %>% filter(isolation_source != '' | host != '')

# Classify everything that could be easily classified as human
metadata2process_1 <- metadata2process_all %>%  mutate(source = case_when(
  grepl("Human*",host) ~ "Human",
  grepl("Homo sapiens*",host) ~ "Human",
  grepl("patient",host) ~ "Human",
  grepl("Homosapiens",host) ~ "Human",
  grepl("Homo-sapiens",host) ~ "Human",
  grepl("Homo",host) ~ "Human"
))

# and filter them out as they are processed
metadata2process_2 <- metadata2process_1 %>% filter(is.na(source)) 

# Set your OpenAI API key
api_key <- "Here will go private OpenAI API key"

# Define a function to call OpenAI's GPT model
classify_sample_with_gpt <- function(description) {
  
  # this message will be send to chatgpt
  messages <- list(
    list(
      role = "user",
      content = paste(
        "Classify the following sample isolation source into one of the categories: human, animal, environment, wastewater ",
        "If the description is vague or unclear, and you cannot confidently classify it, respond with 'unknown'. ",
        "Answer with just one term. Sample:", description
      )
    )
  )
  
  # send POST request
  response <- POST(
    url = "https://api.openai.com/v1/chat/completions",
    add_headers(Authorization = paste("Bearer", api_key)),
    body = list(
      model = "gpt-4o-mini",   
      messages = messages,
      max_tokens = 10,
      temperature = 0
    ),
    encode = "json"
  )
  
  # retrieve only message from response
  content <- content(response, "parsed")
  category <- content[["choices"]][[1]][["message"]][["content"]]
  
  return(trimws(category))  # Remove any extra spaces
}



# Split the data into chunks of 100 rows
# this mostly done as a safety measure, but it is not strictly necesary
chunks <- split(metadata2process_2, ceiling(seq_along(metadata2process_2$Assembly) / 100))

# initiate results object
classified_results <- list()

#for loop for chunks
for (i in seq_along(chunks)) {
  message("Processing chunk: ", i)
  
  # Use tryCatch to handle errors
  result <- tryCatch({
    chunks[[i]] %>%
      mutate(category_gpt = sapply(paste0(isolation_source, "; ", host), classify_sample_with_gpt))
  }, error = function(e) {
    # Print an error message and return NULL for this chunk
    message("Error in chunk ", i, ": ", e$message)
    NULL
  })
  
  # Save the result (even if NULL to keep track of processed chunks)
  classified_results[[i]] <- result
  
  # Optionally save intermediate progress to disk
  saveRDS(classified_results, file = "classified_results_progress.rds")
}

# retrieve from disk classification results and 'unchunk' them
clas_results <- do.call(rbind,
                        readRDS('classified_results_progress.rds'))

# write to file for manual curation
write.table(clas_results %>% select(-source), file = "gpt2source.csv", sep = ';', row.names = F)

# everything questionable got manualy curated and reloaded
gpt2source <- read.csv("gpt2source.csv") %>% select(Assembly,category_gpt)

# remove unnecessary columns
metadataProcessed2 <- gpt2source %>% mutate(category_gpt  =  category_gpt4 == '' ~ 'Unknown',
                                                                               category_gpt4 == 'human' ~ 'Human',
                                                                               category_gpt4 == 'animal' ~ 'Animal',
                                                                               category_gpt4 == 'wastewater' ~ 'Wastewater',
                                                                               category_gpt4 == 'environment' ~ 'Environment',
                                                                               category_gpt4 == 'unknown' ~ 'Unknown',
                                                                               .default = category_gpt4))

# Read preprocessed metadata, which was used for big analysis
metadata <- read.table(file = "cleaned_metadata_210624.tsv", quote = "\"", sep = '\t', header = T)


# bind it with processed data
metadata2 <- metadata %>% 
  left_join(metadataProcessed2, by = c('Assembly.Accession' = 'Assembly')) %>% 
  mutate(Source_category = case_when(is.na(Source_category) ~ "Unknown",
                                     .default = Source_category))
