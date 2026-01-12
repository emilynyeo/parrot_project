# Tutorial executable script 
library(tidyverse)
library(broom)
library(ggtext)
setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")
set.seed(19760620)

# Changable files 

# NANOPORE AND MISEQ from original runs
taxa_file <- read.csv("./00_sequence_data/JVB3218_16S/JVB3218-16S-esv-data.csv") %>% 
             rename_all(tolower) %>% 
             mutate(across(kingdom:species, ~ na_if(., ""))) %>%   # treat "" as NA
             mutate(across(kingdom:species, ~ replace_na(., "unclassified"))) %>%
             unite(taxonomy, kingdom:species, sep = ";", remove = FALSE) %>% 
             mutate(taxonomy = str_replace_all(taxonomy, " ", "_"),
                    taxonomy = str_replace_all(taxonomy, "\\[|\\]", ""), 
                    taxonomy = str_replace(taxonomy, ";unclassified", "_unclassified"),
                    taxonomy = str_replace_all(taxonomy, ";unclassified", ""),
                    taxonomy = str_replace_all(taxonomy, ";$", ""),
                    taxonomy = str_replace_all(taxonomy, ".*;", "")) %>% 
             select(esvid, taxonomy)

otu_file <- read.csv("./00_sequence_data/JVB3218_16S/JVB3218-16S-read-data.csv") %>% 
            rename_with(~ gsub("^s0", "S0", .)) %>% 
            select(-starts_with(c("S082895", "S082806","TestId","sequence",
                                  "Kingdom","Phylum","Class","Order","Family",
                                  "Genus","Species","X..match","X..species"))) %>% 
            pivot_longer(-ESVId, names_to = "samples", values_to = "count") %>% 
            rename_all(tolower)
            
  
samples_16 <- read.csv("./00_sequence_data/JVB3218_16S/JVB3218-samples.csv") 

# Original dataset
meta1 <- read.csv("./00_metadata/ParrotFecalSampleMetadata_SerialIDs_McKenzie.csv")
addmeta <- read.csv("./00_metadata/additional_info_seizures.csv")

meta_file  <- left_join(meta1, addmeta) %>%
  unite(serial.prefix, serial.number, sep="", col="SampleId", remove=FALSE) %>%
  select(SampleId, everything()) %>% 
  rename_all(tolower) %>%
  rename(samples = sampleid) %>%
  rename_with(~ str_replace_all(.x, "\\.", "_")) %>% 
  mutate(across(where(is.character), ~ str_replace_all(.x, " ", "_"))) %>% 
  filter(samples != "S082895") %>% 
  filter(samples != "S082806")


composite <- inner_join(otu_file, taxa_file, by ="esvid") %>% 
  group_by(samples, taxonomy) %>% 
  summarise(count = sum(count), .groups = "drop") %>% 
  group_by(samples) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  ungroup() %>%
  select(-count) %>%
  inner_join(., meta_file, by="samples")
  



metadata <- read_tsv("raw_data/baxter.metadata.tsv", 
                     col_types=cols(sample = col_character())) %>%
  rename_all(tolower) %>% 
  rename(group = sample) %>%
  mutate(srn = dx_bin == "Adv Adenoma" | dx_bin == "Cancer",
         lesion = dx_bin == "Adv Adenoma" | dx_bin == "Cancer" | dx_bin == "Adenoma")

composite <- inner_join(shared, taxonomy, by="otu") %>%
  group_by(group, taxonomy) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(group) %>%
  mutate(rel_abund = count / sum(count)) %>% 
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="group")