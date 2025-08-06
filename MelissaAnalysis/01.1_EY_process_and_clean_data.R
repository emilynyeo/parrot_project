#!bin/bash

#### Load libraries ####
library(vegan)
library(MASS)
library(indicspecies)
library(tidyverse)
library(phyloseq)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(scales)
library(calecopal)

setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")

# dir.create("01_process_and_clean_data")
if (!dir.exists("01_process_and_clean_data")) {
  dir.create("01_process_and_clean_data")
}

#### Load data ####
# Original dataset
meta1 <- read.csv("./00_metadata/ParrotFecalSampleMetadata_SerialIDs_McKenzie.csv")
addmeta <- read.csv("./00_metadata/additional_info_seizures.csv")
meta1 <- left_join(meta1, addmeta) %>%
  unite(serial.prefix, serial.number, sep="", col="SampleId", remove=FALSE) %>%
  select(SampleId, everything())

# Renee's extra samples, unknown origin
meta2 <- read.csv("./00_metadata/ReneEb_ParrotSamples_Feb2025_250508_JD.csv") %>%
  rename(SampleId = JonahBarcodeID) %>%
  rename(Date.collected = Sampling.date) %>%
  mutate(Captive.Wild = paste0("Crate ", Shipment, ": Unknown origin")) %>%
  mutate(Tube.type = "5ml RNAlater factory") %>%
  mutate(Sample.type = "Feces", Scientific.Name = "Psittacus erithacus", Common.Name = "African Grey parrot") %>%
  unite(Owner.name,"Shipment",Shipment, col="Location") %>%
  rename(Site.Name = ShipmentCrate) %>%
  mutate(Country = "Unknown")

# Combine the metadata
meta <- meta1 %>%
  full_join(meta2) %>%
  mutate(Given.probiotic= ifelse(is.na(Given.probiotic), "", 
                                 ifelse(as.character(Given.probiotic)=="Yes", "Given probiotics", "")))

# NANOPORE AND MISEQ from original runs
esv_16 <- read.csv("./00_sequence_data/JVB3218_16S/JVB3218-16S-esv-data.csv")
reads_16 <- read.csv("./00_sequence_data/JVB3218_16S/JVB3218-16S-read-data.csv")
samples_16 <- read.csv("./00_sequence_data/JVB3218_16S/JVB3218-samples.csv")

esv_18 <- read.csv("./00_sequence_data/JVB3218_18S_616/JVB3218-18S_616-esv-data.csv")
reads_18 <- read.csv("./00_sequence_data/JVB3218_18S_616/JVB3218-18S_616-read-data.csv")
samples_18 <- read.csv("./00_sequence_data/JVB3218_18S_616/JVB3218-samples.csv")

# NANOPORE- RENEE
esv_16ren <- read.csv("./00_sequence_data/JVB4814_16S/JVB4814-16S-esv-data.csv")
reads_16ren <- read.csv("./00_sequence_data/JVB4814_16S/JVB4814-16S-read-data.csv")
samples_16ren <- read.csv("./00_sequence_data/JVB4814_16S/JVB4814-samples.csv")

esv_18ren <- read.csv("./00_sequence_data/JVB4814_18S_616/JVB4814-18S_616-esv-data.csv")
reads_18ren <- read.csv("./00_sequence_data/JVB4814_18S_616/JVB4814-18S_616-read-data.csv")
samples_18ren <- read.csv("./00_sequence_data/JVB4814_18S_616/JVB4814-samples.csv")


# ### NOTE: The ESV taxonomies kinda suck no matter what databse you use
# because the database is often wrong (e.g. rice is listed as the species under Pythium???)
# Even BLAST is "wrong" because it lacks a lot of info (e.g. no Parrot DNA)
# Going to use Jonah's assignments preliminarily,

###### SampleID names: adjust and edit #########
# Make sampleID column
meta_adj <- meta %>% select(SampleId, everything())

# Make reads into a matrix
# So it turns out both 16S and 18S column names have a weird thing where the first sample ID has a lower case S? Let's fix that.
# Joe note: "On AWS the new Nanopore data is all rep 2, with the exception of 12 samples that already had 2 reps for 18S_616, which I made rep 3."
# Suffixes for 16S are either 0, 1, or 2. 
# Suffixes for 18S are either 0, 1, 2, or 3

# Are set 0's sums 0? ~  Yes; remove these.
reads_16 %>% select(ends_with(".0")) %>% colSums()
reads_18 %>% select(ends_with(".0")) %>% colSums()

colnames(reads_16) <- gsub("^s0","S0", colnames(reads_16)) # 703 cols, 11628 rows 

# Split nano and miseq and renee
# JOE: "AWS new Nanopore data is all rep 2, except 12 samples that already had 2 reps for 18S_616, which I made rep 3."
otu16m_df <- reads_16 %>% select(ESVId, starts_with("S0") & ends_with(".1")) # Run 1 are miseq
otu16n_df <- reads_16 %>% select(ESVId, starts_with("S0") & ends_with(".2")) # Run 2 is Nano

colnames(reads_18) <- gsub("^s0","S0", colnames(reads_18)) # 713 cols
r118 <- reads_18 %>% select(starts_with("S0")) %>% select(ends_with(".1")) %>% colnames() # 344 
r218 <- reads_18 %>% select(starts_with("S0")) %>% select(ends_with(".2")) %>% colnames() # 340 
r318 <- reads_18 %>% select(starts_with("S0")) %>% select(ends_with(".3")) %>% colnames() # 12 
r318_woSuffix <- gsub(".3$","",r318)

otu18m_df <- reads_18 %>%
  select(ESVId, (starts_with("S0") & (ends_with(".1") | one_of(paste0(r318_woSuffix,".2"))))) %>%
  select(-one_of(paste0(r318_woSuffix,".1")))

otu18n_df <- reads_18 %>%
  select(ESVId, (starts_with("S0") & (ends_with(".2") | one_of(paste0(r318_woSuffix,".3"))))) %>%
  select(-one_of(paste0(r318_woSuffix,".2")))

colnames(reads_16ren) <- gsub("^s0","S0", colnames(reads_16ren))
otu16ren_df <- reads_16ren %>%
  select(ESVId, starts_with("S1") & ends_with(".1"))

colnames(reads_18ren) <- gsub("^s0","S0", colnames(reads_18ren))
otu18ren_df <- reads_18ren %>%
  select(ESVId, starts_with("S1") & ends_with(".1"))

# Removing rows with zero sum 
table(rowSums(otu16m_df[,-1])>0)
otu16m_df <- otu16m_df[rowSums(otu16m_df[,-1])>0,] # 16 miseq 3638 removed

table(rowSums(otu16n_df[,-1])>0)
otu16n_df <- otu16n_df[rowSums(otu16n_df[,-1])>0,] # 16 nano 2474 removed

table(rowSums(otu18m_df[,-1])>0) # 1586 removed
otu18m_df <- otu18m_df[rowSums(otu18m_df[,-1])>0,] 

table(rowSums(otu18n_df[,-1])>0) # 10292 removed 
otu18n_df <- otu18n_df[rowSums(otu18n_df[,-1])>0,] 

table(rowSums(otu16ren_df[,-1])>0)
otu16ren_df <- otu16ren_df[rowSums(otu16ren_df[,-1])>0,] # No empty rows 

table(rowSums(otu18ren_df[,-1])>0)
otu18ren_df <- otu18ren_df[rowSums(otu18ren_df[,-1])>0,] # Remove empty rows

### Remove suffixes from samplenames 
colnames(otu16m_df) <- gsub("[.]1","", colnames(otu16m_df))
colnames(otu18m_df) <- gsub("[.]2", "", gsub("[.]1","", colnames(otu18m_df)))

colnames(otu16n_df) <- gsub("[.]2","", colnames(otu16n_df))
colnames(otu18n_df) <- gsub("[.]2", "", gsub("[.]3","", colnames(otu18n_df)))

colnames(otu16ren_df) <- gsub("[.]1","", colnames(otu16ren_df))
colnames(otu18ren_df) <- gsub("[.]1", "",colnames(otu18ren_df))

### Make into a matrix
otu16m_mat <- as.matrix(otu16m_df[,-1])
rownames(otu16m_mat) <- otu16m_df[,1]
# dim(otu16m_mat): 7990  343

otu18m_mat <- as.matrix(otu18m_df[,-1])
rownames(otu18m_mat) <- otu18m_df[,1]
# dim(otu18m_mat): 11868   344

otu16n_mat <- as.matrix(otu16n_df[,-1])
rownames(otu16n_mat) <- otu16n_df[,1]
# dim(otu16n_mat): 9154  344

otu18n_mat <- as.matrix(otu18n_df[,-1])
rownames(otu18n_mat) <- otu18n_df[,1]
# dim(otu18n_mat): 11868   340

otu16ren_mat <- as.matrix(otu16ren_df[,-1])
rownames(otu16ren_mat) <- otu16ren_df[,1]
# dim(otu16ren_mat): 684  23

otu18ren_mat <- as.matrix(otu18ren_df[,-1])
rownames(otu18ren_mat) <- otu18ren_df[,1]
# dim(otu18ren_mat): 200  23

#### Cleaning up any metadata ####

# Add read depth to metadata info

ReadDepth16m <- colSums(otu16m_mat) %>% as.data.frame() %>%
                rownames_to_column(var="SampleId") %>%
                rename(ReadDepth16m=".")
ReadDepth18m <- colSums(otu18m_mat) %>% as.data.frame() %>%
                rownames_to_column(var="SampleId") %>%
                rename(ReadDepth18m=".")
ReadDepth16n <- colSums(otu16n_mat) %>% as.data.frame() %>%
                rownames_to_column(var="SampleId") %>%
                rename(ReadDepth16n=".")
ReadDepth18n <- colSums(otu18n_mat) %>% as.data.frame() %>%
                rownames_to_column(var="SampleId") %>%
                rename(ReadDepth18n=".")
ReadDepth16ren <- colSums(otu16ren_mat) %>% as.data.frame() %>%
                rownames_to_column(var="SampleId") %>%
                rename(ReadDepth16n=".")
ReadDepth16ren$ReadDepth16ren <- as.numeric(ReadDepth16ren$ReadDepth16n)

ReadDepth18ren <- colSums(otu18ren_mat) %>% as.data.frame() %>%
                rownames_to_column(var="SampleId") %>%
                rename(ReadDepth18ren=".")
ReadDepth18ren$ReadDepth18n <- ReadDepth18ren$ReadDepth18ren

meta_adj <- meta_adj %>% 
            left_join(full_join(full_join(full_join(ReadDepth16m,ReadDepth18m), 
                                          full_join(ReadDepth16n, ReadDepth18n)), 
                                full_join(ReadDepth16ren, ReadDepth18ren))) %>%
                      mutate(Captive.Wild = factor(Captive.Wild, 
                            levels=c("Captive","Wild, free ranging", 
                                     "Wild, seized from traffickers", 
                                     "Crate 1: Unknown origin", 
                                     "Crate 2: Unknown origin")))

# Look at diff in read depths 
meta_adj %>% ggplot() + geom_point(aes(x=ReadDepth16m, y=ReadDepth16n))
meta_adj %>% ggplot() + geom_point(aes(x=ReadDepth18m, y=ReadDepth18n))

# What's up with those outliers?
meta_adj %>% filter(ReadDepth18n>100000)

# They are seemingly three random samples? Let's ignore and keep for now.
meta_adj %>% ggplot() + geom_point(aes(x=ReadDepth18n, y=ReadDepth16n, col=Captive.Wild))

# Standardize date format
meta_adj$Date.collected
meta_adj %>% select(Country, Date.collected) %>% table()

# Dashes look like SA dates (MM-DD-YYYY) whereas Nigerian = DD/MM/YY: Let's fix those...
reformatDates <- as.data.frame(matrix(nrow=length(meta_adj$Date.collected), 
                               ncol=5, 
                               dimnames = list(NULL, c("Date.collected","Day","Month","Year","DateAdj"))))

# function to extract dates from strings and re-standardize them
for ( i in 1:length(meta_adj$Date.collected)) {
  val = meta_adj$Date.collected[i]
  if ( i == "no data") {
    newentry <- c(val, NA, NA,NA)
  } else if (length(grep("/", val))>0) { 
    newdate <- ymd(as.Date(val, format = "%d/%m/%Y"))
    newentry <- c(val, day(newdate), month(newdate), year(newdate),as.character(newdate))
  } else if (length(grep("[aA-zZ]",val))>0 ) {
    sepDate <- unlist(strsplit(val,split = "-"))
    newdate <- ymd(as.Date(paste(sepDate[1],match(sepDate[2],month.abb), sepDate[3], sep = "-"), format = "%d-%m-%Y"))
  } else { 
    newdate <- ymd(as.Date(val, format = "%m-%d-%Y"))
    newentry <- c(val, day(newdate), month(newdate), year(newdate),as.character(newdate))
  }
  reformatDates[i, ] <- newentry
}

unique(reformatDates)
meta_adj <- meta_adj %>% left_join(distinct(reformatDates)) 

### Time since captive ###
meta_adj$Seized.date

meta_adj$Seized.date <- ymd(as.Date(meta_adj$Seized.date, format = "%d-%m-%Y"))
meta_adj <- meta_adj %>% mutate(TimeSinceSeizure = as.Date(DateAdj)-Seized.date)

# Time since seizure conflated with probiotic treatment
meta_adj %>% select(Given.probiotic, TimeSinceSeizure) %>%
  table()

## Read Depth by sampling time ##
meta_adj %>%
  ggplot() + geom_point(aes(x=DateAdj, y=ReadDepth16n, col = Captive.Wild, shape=Country))

meta_adj %>%
  ggplot() + geom_point(aes(x=DateAdj, y=ReadDepth18n, col = Captive.Wild, shape=Country))

# Are tube storage options correlated with samples?
meta_adj %>% select(Captive.Wild, Tube.type) %>% table()

# the 10mL conicals were ONLY used in captive. Check later that this isn't an issue
meta_adj %>%
  ggplot() +
  geom_jitter(aes(x=Captive.Wild, y=ReadDepth16n, col=Tube.type))

meta_adj %>%
  ggplot() +
  geom_jitter(aes(x=Captive.Wild, y=ReadDepth18n, col=Tube.type))

# STEP 1: 
#### Get taxonomic assignments original (from Jonah) ####
tax16 <- reads_16 %>% select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) 
tax16 <- tax16 %>%
  mutate(ScientificName = ifelse(Species!="", paste0(Genus, " ", Species),
                          ifelse(Genus !="", paste0(Genus, " sp"),
                          ifelse(Family!="", paste0("Family:", Family),
                          ifelse(Order!="", paste0("Order:", Order),
                          ifelse(Class!="", paste0("Class:", Class),
                          ifelse(Phylum!="", paste0("Phylum:",Phylum),
                          ifelse(Kingdom!="", paste0("Kingdom:", Kingdom),"Unassigned"))))))))
# Split into nanopore and miseq
tax16m <- tax16 %>% filter(ESVId %in% otu16m_df$ESVId)
tax16n <- tax16 %>% filter(ESVId %in% otu16n_df$ESVId)

# Get taxonomic assignments Renees
tax16ren <- reads_16ren %>% select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(ScientificName = ifelse(Species!="", paste0(Genus, " ", Species),
                          ifelse(Genus !="", paste0(Genus, " sp"),
                          ifelse(Family!="", paste0("Family:", Family),
                          ifelse(Order!="", paste0("Order:", Order),
                          ifelse(Class!="", paste0("Class:", Class),
                          ifelse(Phylum!="", paste0("Phylum:",Phylum),
                          ifelse(Kingdom!="", paste0("Kingdom:", Kingdom),"Unassigned"))))))))

## Now 18S original
tax18 <- reads_18 %>% select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(ScientificName = ifelse(Species!="", gsub("_"," ", Species),
                          ifelse((Genus !="") & (Genus!=Family), paste0(Genus, " sp"),
                          ifelse((Family!="")& (Family!=Order), paste0("Family:", Family),
                          ifelse((Order!="")&(Order!=Class), paste0("Order:", Order),
                          ifelse((Class!="")&(Class!=Phylum), paste0("Class:", Class),
                          ifelse((Phylum!="") & (Phylum!=Kingdom), paste0("Phylum:",Phylum),
                          ifelse(Kingdom!="", paste0("Kingdom:", Kingdom),"Unassigned"))))))))

# Split into nanopore and miseq
tax18m <- tax18 %>% filter(ESVId %in% otu18m_df$ESVId)
tax18n <- tax18 %>% filter(ESVId %in% otu18n_df$ESVId)

# 18s Ren
tax18ren <- reads_18ren %>% select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(ScientificName = ifelse(Species!="", gsub("_"," ", Species),
                          ifelse((Genus !="") & (Genus!=Family), paste0(Genus, " sp"),
                          ifelse((Family!="")& (Family!=Order), paste0("Family:", Family),
                          ifelse((Order!="")&(Order!=Class), paste0("Order:", Order),
                          ifelse((Class!="")&(Class!=Phylum), paste0("Class:", Class),
                          ifelse((Phylum!="") & (Phylum!=Kingdom), paste0("Phylum:",Phylum),
                          ifelse(Kingdom!="", paste0("Kingdom:", Kingdom), "Unassigned"))))))))

########## JONAS ASSIGNMENTS ##################
# Processes tax data to classify each ESV as Bacteria, Archaea, Eukaryota, Mitochondria, or Chloroplast, and reshapes the data into a binary presence/absence matrix, where each row represents a taxonomic type and each column an ESV

## 16S misiq original
otu16m_bactornot <- tax16m %>% 
  mutate(Type = ifelse(Kingdom=="Eukaryota",Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
              # ifelse(is.na(Phylum), "No phylum assigned",
                ifelse(Family=="Mitochondria", Family, 
                ifelse(Order == "Chloroplast", "Chloroplast", "Bacteria")))), pres = 1) %>%
    mutate(Type = ifelse(is.na(Type), "Bacteria", Type)) %>%
    select(ESVId, Type, pres) %>%
    filter(ESVId %in% c(rownames(otu16n_mat), rownames(otu16m_mat))) %>%
    pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
    column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

## 16S nanopore original
otu16n_bactornot <- tax16n %>% 
  mutate(Type = ifelse(Kingdom=="Eukaryota", Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
              # ifelse(is.na(Phylum), "No phylum assigned",
                ifelse(Family=="Mitochondria", Family, 
                ifelse(Order == "Chloroplast", "Chloroplast", "Bacteria")))), pres = 1) %>%
  mutate(Type = ifelse(is.na(Type), "Bacteria", Type)) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% c(rownames(otu16n_mat), rownames(otu16m_mat))) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

## 16S ren
otu16ren_bactornot <- tax16ren %>% 
  mutate(Type = ifelse(Kingdom=="Eukaryota",Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
              # ifelse(is.na(Phylum), "No phylum assigned",
                ifelse(Family=="Mitochondria", Family, 
                ifelse(Order == "Chloroplast", "Chloroplast", "Bacteria")))), pres = 1) %>%
  mutate(Type = ifelse(is.na(Type), "Bacteria", Type)) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% c(rownames(otu16ren_mat))) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

# aggregate OTU counts by taxonomic types
otu16n_mat_bactornot <- otu16n_bactornot %*% otu16n_mat[match(colnames(otu16n_bactornot),rownames(otu16n_mat)),]

# calculate relative abundance 
otu16n_mat_bactornot_RA <- apply(otu16n_mat_bactornot, MARGIN=2, FUN=function(x) x/sum(x))


### 18S
# Classify ESV as Bacteria, Archaea, Eukaryota, Mitochondria, or Chloroplast, and 
# Reshapes the data into a binary presence/absence matrix, where each row represents a taxonomic type and each column an ESV

#### 18s miseq Bacteria, Archaea, Eukaryota, Mitochondria, or Chloroplast
otu18m_eukornot <- tax18m %>% 
  mutate(Type = ifelse(Kingdom=="Bacteria",Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
                ifelse(Kingdom=="", "No kingdom assigned",
                ifelse(is.na(Phylum)|(Phylum==""), 
                       "No phylum assigned", 
                       "Eukaryota (wphylum)")))), pres = 1) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% rownames(otu18m_mat)) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

#### 18s nano Bacteria, Archaea, Eukaryota, Mitochondria, or Chloroplast
otu18n_eukornot <- tax18n %>% 
  mutate(Type =ifelse(Kingdom=="Bacteria",Kingdom,
               ifelse(Kingdom=="Archaea", Kingdom,
               ifelse(Kingdom=="", "No kingdom assigned",
               ifelse(is.na(Phylum)|(Phylum==""), 
                      "No phylum assigned", 
                      "Eukaryota (wphylum)")))), pres = 1) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% rownames(otu18n_mat)) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

#### 18s Rene 
otu18ren_eukornot <- tax18ren %>% 
  mutate(Type = ifelse(Kingdom=="Bacteria",Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
                ifelse(Kingdom=="", "No kingdom assigned",
                ifelse(is.na(Phylum)|(Phylum==""), 
                       "No phylum assigned", 
                       "Eukaryota (wphylum)")))), pres = 1) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% c(rownames(otu18ren_mat))) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

# plot 18s miseq original
otu18m_mat_eukornot <- otu18m_eukornot %*% otu18m_mat[match(colnames(otu18m_eukornot),rownames(otu18m_mat)),]
otu18m_mat_eukornot_RA <- apply(otu18m_mat_eukornot, MARGIN=2, FUN=function(x) x/sum(x))

# plot 18s nano original
otu18n_mat_eukornot <- otu18n_eukornot %*% otu18n_mat[match(colnames(otu18n_eukornot),rownames(otu18n_mat)),]
otu18n_mat_eukornot_RA <- apply(otu18n_mat_eukornot, MARGIN=2, FUN=function(x) x/sum(x))

# plot 18s Rene
otu18ren_mat_eukornot <- otu18ren_eukornot %*% otu18ren_mat[match(colnames(otu18ren_eukornot),rownames(otu18ren_mat)),]
otu18ren_mat_eukornot_RA <- apply(otu18ren_mat_eukornot, MARGIN=2, FUN=function(x) x/sum(x))


### 18S-- more detailed
keepPhylum <- c("Cnidaria","Mollusca","Porifera","Vertebrata","Chlorophyta", "Arthropoda","Phragmoplastophyta")

# Miseq
otu18m_eukdetailed <- tax18m %>% 
  mutate(Type = ifelse(Kingdom=="Bacteria",Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
                ifelse(Kingdom=="","No kingdom assigned",
                ifelse(is.na(Phylum)|Phylum=="", "No phylum assigned", 
                ifelse(Phylum %in% keepPhylum, Phylum, "MicroEuk"))))), pres = 1) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% rownames(otu18m_mat)) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

otu18m_mat_eukdetailed <- otu18m_eukdetailed %*% otu18m_mat[match(colnames(otu18m_eukdetailed),rownames(otu18m_mat)),]
otu18m_mat_eukdetailed_RA <- apply(otu18m_mat_eukdetailed, MARGIN=2, FUN=function(x) x/sum(x))

# Nanopore
otu18n_eukdetailed <- tax18n %>% 
  mutate(Type = ifelse(Kingdom=="Bacteria",Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
                ifelse(Kingdom=="","No kingdom assigned",
                ifelse(is.na(Phylum)|Phylum=="", "No phylum assigned", 
                ifelse(Phylum %in% keepPhylum, Phylum, "MicroEuk"))))), pres = 1) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% rownames(otu18n_mat)) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

otu18n_mat_eukdetailed <- otu18n_eukdetailed %*% otu18n_mat[match(colnames(otu18n_eukdetailed),rownames(otu18n_mat)),]
otu18n_mat_eukdetailed_RA <- apply(otu18n_mat_eukdetailed, MARGIN=2, FUN=function(x) x/sum(x))

# Renee
otu18ren_eukdetailed <- tax18ren %>% 
  mutate(Type = ifelse(Kingdom=="Bacteria",Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
                ifelse(Kingdom=="","No kingdom assigned",
                ifelse(is.na(Phylum)|Phylum=="", "No phylum assigned", 
                ifelse(Phylum %in% keepPhylum, Phylum, "MicroEuk"))))), pres = 1) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% rownames(otu18ren_mat)) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

otu18ren_mat_eukdetailed <- otu18ren_eukdetailed %*% otu18ren_mat[match(colnames(otu18ren_eukdetailed),rownames(otu18ren_mat)),]
otu18ren_mat_eukdetailed_RA <- apply(otu18ren_mat_eukdetailed, MARGIN=2, FUN=function(x) x/sum(x))

### 18S-- more detailed 2
tax18 %>% filter(Class =="Embryophyta")
keepClass <- c("Aves","Gastropoda","Insecta","Mammalia","Embryophyta")

# Miseq
otu18m_eukdetailed2 <- tax18m %>% 
  mutate(Type = ifelse(Kingdom=="Bacteria",Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
                ifelse(Kingdom=="", "No kingdom assigned",
                ifelse(is.na(Phylum)|Phylum=="", "No phylum assigned", 
                ifelse(!Phylum %in% keepPhylum, "MicroEuk",
                ifelse(Class %in% keepClass, Class, "Other macroeukaryote")))))), pres = 1) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% rownames(otu18m_mat)) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

otu18m_mat_eukdetailed2 <- otu18m_eukdetailed2 %*% otu18m_mat[match(colnames(otu18m_eukdetailed2),rownames(otu18m_mat)),]
otu18m_mat_eukdetailed2_RA <- apply(otu18m_mat_eukdetailed2, MARGIN=2, FUN=function(x) x/sum(x))

# Nanopore
otu18n_eukdetailed2 <- tax18n %>% 
  mutate(Type = ifelse(Kingdom=="Bacteria",Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
                ifelse(Kingdom=="", "No kingdom assigned",
                ifelse(is.na(Phylum)|Phylum=="", "No phylum assigned", 
                ifelse(!Phylum %in% keepPhylum, "MicroEuk",
                ifelse(Class %in% keepClass, Class, "Other macroeukaryote")))))), pres = 1) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% rownames(otu18n_mat)) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

otu18n_mat_eukdetailed2 <- otu18n_eukdetailed2 %*% otu18n_mat[match(colnames(otu18n_eukdetailed2),rownames(otu18n_mat)),]
otu18n_mat_eukdetailed2_RA <- apply(otu18n_mat_eukdetailed2, MARGIN=2, FUN=function(x) x/sum(x))

# Renee
otu18ren_eukdetailed2 <- tax18ren %>% 
  mutate(Type = ifelse(Kingdom=="Bacteria",Kingdom, 
                ifelse(Kingdom=="Archaea", Kingdom,
                ifelse(Kingdom=="", "No kingdom assigned",
                ifelse(is.na(Phylum)|Phylum=="", "No phylum assigned", 
                ifelse(!Phylum %in% keepPhylum, "MicroEuk",
                ifelse(Class %in% keepClass, Class, "Other macroeukaryote")))))), pres = 1) %>%
  select(ESVId, Type, pres) %>%
  filter(ESVId %in% rownames(otu18ren_mat)) %>%
  pivot_wider(names_from="Type", values_from="pres", values_fill = 0) %>%
  column_to_rownames(var="ESVId") %>% as.matrix() %>% t()

otu18ren_mat_eukdetailed2 <- otu18ren_eukdetailed2 %*% otu18ren_mat[match(colnames(otu18ren_eukdetailed2),rownames(otu18ren_mat)),]
otu18ren_mat_eukdetailed2_RA <- apply(otu18ren_mat_eukdetailed2, MARGIN=2, FUN=function(x) x/sum(x))


#### Remove chloroplasts and mitochondria for 16S ####
# Chloroplasts are Order level; Mitochondria are family level
esv16m_toremove <- tax16m %>% filter(Family == "Mitochondria"| Order =="Chloroplast") %>% pull(ESVId)
esv16n_toremove <- tax16n %>% filter(Family == "Mitochondria"| Order =="Chloroplast") %>% pull(ESVId)
esv16ren_toremove <- tax16ren %>% filter(Family == "Mitochondria"| Order =="Chloroplast") %>% pull(ESVId)

esv_16m_filt <- tax16m %>% filter(!ESVId %in% esv16m_toremove)
esv_16n_filt <- tax16n %>% filter(!ESVId %in% esv16n_toremove)
esv_16ren_filt <- tax16ren %>% filter(!ESVId %in% esv16ren_toremove)

#### Remove host DNA and others for 18S ####
esv18m_toremove <- tax18m %>% filter(Kingdom=="Bacteria" | Kingdom=="Archaea"| Kingdom=="") %>% pull(ESVId)
esv18n_toremove <- tax18n %>% filter(Kingdom=="Bacteria" | Kingdom=="Archaea"| Kingdom=="") %>% pull(ESVId)
esv18ren_toremove <- tax18ren %>% filter(Kingdom=="Bacteria" | Kingdom=="Archaea"| Kingdom=="") %>% pull(ESVId)

esv_18m_filt <- tax18m %>% filter(!ESVId %in% esv18m_toremove)
esv_18n_filt <- tax18n %>% filter(!ESVId %in% esv18n_toremove)
esv_18ren_filt <- tax18ren %>% filter(!ESVId %in% esv18ren_toremove)

# Filter out ESVs from tables based on this
otu16m_mat_filt <- otu16m_mat[which(rownames(otu16m_mat) %in% esv_16m_filt$ESVId),]
otu16n_mat_filt <- otu16n_mat[which(rownames(otu16n_mat) %in% esv_16n_filt$ESVId),]
otu16ren_mat_filt <- otu16ren_mat[which(rownames(otu16ren_mat) %in% esv_16ren_filt$ESVId),]
otu18m_mat_filt <- otu18m_mat[which(rownames(otu18m_mat) %in% esv_18m_filt$ESVId),]
otu18n_mat_filt <- otu18n_mat[which(rownames(otu18n_mat) %in% esv_18n_filt$ESVId),]
otu18ren_mat_filt <- otu18ren_mat[which(rownames(otu18ren_mat) %in% esv_18ren_filt$ESVId),]

# Version with no host DNA and no 'mammal' DNA. Not sure what that is.. if it's actually people or not
esv_18m_nohost_filt <- esv_18m_filt %>% filter(!is.na(Class), Class != "Actinopterygii", Class !="Aves", Class !="Mammalia")
otu18m_nohost_mat_filt <- otu18m_mat[which(rownames(otu18m_mat) %in% esv_18m_nohost_filt$ESVId),]
esv_18n_nohost_filt <- esv_18n_filt %>% filter(!is.na(Class), Class != "Actinopterygii", Class !="Aves", Class !="Mammalia")
otu18n_nohost_mat_filt <- otu18n_mat[which(rownames(otu18n_mat) %in% esv_18n_nohost_filt$ESVId),]
esv_18ren_nohost_filt <- esv_18ren_filt %>% filter(!is.na(Class), Class != "Actinopterygii", Class !="Aves", Class !="Mammalia")
otu18ren_nohost_mat_filt <- otu18ren_mat[which(rownames(otu18ren_mat) %in% esv_18ren_nohost_filt$ESVId),]


# # # #  STEP 2: # # # # 

# Recalculate read depth after host DNA filtering
ReadDepth16m_filt <- colSums(otu16m_mat_filt) %>%
  as.data.frame() %>% rownames_to_column(var = "SampleId") %>%
  rename(ReadDepth16m_filt = ".")

ReadDepth16n_filt <- colSums(otu16n_mat_filt) %>%
  as.data.frame() %>% rownames_to_column(var = "SampleId") %>%
  rename(ReadDepth16n_filt = ".")

ReadDepth16ren_filt <- colSums(otu16ren_mat_filt) %>%
  as.data.frame() %>% rownames_to_column(var = "SampleId") %>%
  rename(ReadDepth16n_filt = ".")
ReadDepth16ren_filt$ReadDepth16ren_filt <- ReadDepth16ren_filt$ReadDepth16n_filt

ReadDepth18m_filt <- colSums(otu18m_nohost_mat_filt) %>%
  as.data.frame() %>% rownames_to_column(var = "SampleId") %>%
  rename(ReadDepth18m_filt = ".")

ReadDepth18n_filt <- colSums(otu18n_nohost_mat_filt) %>%
  as.data.frame() %>% rownames_to_column(var = "SampleId") %>%
  rename(ReadDepth18n_filt = ".")

ReadDepth18ren_filt <- colSums(otu18ren_nohost_mat_filt) %>%
  as.data.frame() %>% rownames_to_column(var = "SampleId") %>%
  rename(ReadDepth18n_filt = ".")
ReadDepth18ren_filt$ReadDepth18ren_filt <- ReadDepth18ren_filt$ReadDepth18n_filt

meta_adj_filt <- meta_adj %>% 
  left_join(full_join(full_join(full_join(ReadDepth16m_filt,ReadDepth18m_filt), 
                                full_join(ReadDepth16n_filt, ReadDepth18n_filt)), 
                      full_join(ReadDepth16ren_filt, ReadDepth18ren_filt)))

# Plot 1: 16m
p1 <- ggplot(meta_adj_filt, aes(x = ReadDepth16m, y = ReadDepth16m_filt)) +
  geom_point(alpha = 0.6) +
  labs(title = "16m", x = "Original Read Depth (16m)", y = "Filtered Read Depth (16m)") +
  theme_minimal()

# Plot 2: 16n
p2 <- ggplot(meta_adj_filt, aes(x = ReadDepth16n, y = ReadDepth16n_filt)) +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = label_number()) +
  scale_y_continuous(labels = label_number()) +
  labs(title = "16n", x = "Original Read Depth (16n)", y = "Filtered Read Depth (16n)") +
  theme_minimal()

# Plot 3: 16ren
p3 <- ggplot(meta_adj_filt, aes(x = ReadDepth16ren, y = ReadDepth16ren_filt)) +
  geom_point(alpha = 0.6) +
  labs(title = "16ren", x = "Original Read Depth (16ren)", y = "Filtered Read Depth (16ren)") +
  theme_minimal()

# Plot 4: 18m
p4 <- ggplot(meta_adj_filt, aes(x = ReadDepth18m, y = ReadDepth18m_filt)) +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = label_number()) +
  scale_y_continuous(labels = label_number()) +
  labs(title = "18m", x = "Original Read Depth (18m)", y = "Filtered Read Depth (18m)") +
  theme_minimal()

# Plot 5: 18n
p5 <- ggplot(meta_adj_filt, aes(x = ReadDepth18n, y = ReadDepth18n_filt)) +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = label_number()) +
  scale_y_continuous(labels = label_number()) +
  labs(title = "18n", x = "Original Read Depth (18n)", y = "Filtered Read Depth (18n)") +
  theme_minimal()

# Plot 6: 18ren
p6 <- ggplot(meta_adj_filt, aes(x = ReadDepth18ren, y = ReadDepth18ren_filt)) +
  geom_point(alpha = 0.6) +
  labs(title = "18ren", x = "Original Read Depth (18ren)", y = "Filtered Read Depth (18ren)") +
  theme_minimal()

post_filt_depths <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
ggsave("01.1_qc_checks/read_depth_pre_post_filtering_host.png", post_filt_depths, height=7, width=6)

# Average read depth 
filt_cols <- meta_adj_filt %>%
  select(c("ReadDepth16m","ReadDepth16m_filt", 
           "ReadDepth18m","ReadDepth18m_filt", 
           "ReadDepth16n", "ReadDepth16n_filt", 
           "ReadDepth18n","ReadDepth18n_filt",
           "ReadDepth16ren", "ReadDepth16ren_filt",
           "ReadDepth18ren", "ReadDepth18ren_filt"))
filt_means <- colMeans(filt_cols, na.rm = TRUE)

# data frame for plot
filt_df <- data.frame(FilteredMetric = names(filt_means),
                      MeanReadDepth = filt_means)
# Bar plot of mean read depth 
bigsur = c("#E4DECE","#F5DECE", "#ECBD95", "#ECBD82", "#8CA1BB", "#8CA1CB", 
           "#79ACBD", "#89BBBD", "#346575", "#277575", "#1F5529", "#3A5139")
mean_depth <- ggplot(filt_df, aes(x = FilteredMetric, y = MeanReadDepth, 
                    fill = FilteredMetric)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(MeanReadDepth, 1)), 
            vjust = -0.3, size = 3.5) +                      
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Average Read Depth Pre & Post Host DNA Filtering",
       x = "Filtered Column",
       y = "Mean Read Depth") +
  scale_fill_discrete(name = "Sequence Type") +
  scale_fill_manual(values = bigsur) 
ggsave("01.1_qc_checks/average_depth_pre_post_filtering_bar.png", mean_depth, height=5, width=8)

# Summary table of Read Depths :

# Create summary table of read depths 
summary_table <- data.frame(
  FilteredMetric = character(),
  Mean = numeric(),
  Median = numeric(),
  Min = numeric(),
  Max = numeric(),
  stringsAsFactors = FALSE)

# Loop through each column and compute stats
for (col in names(filt_cols)) {
  col_data <- filt_cols[[col]]
  summary_table <- rbind(summary_table, data.frame(
    FilteredMetric = col,
    Mean = mean(col_data, na.rm = TRUE),
    Median = median(col_data, na.rm = TRUE),
    Min = min(col_data, na.rm = TRUE),
    Max = max(col_data, na.rm = TRUE)))
}

# View the table
print(summary_table)


### Read Coverage function
analyze_coverage <- function(otu_mat, sample_name = "Sample") {
  sample_counts <- colSums(otu_mat)
  sample_counts <- sample_counts[sample_counts > 0]  # Remove zero-count samples
  
  if (length(sample_counts) == 0) {
    warning(paste("All samples have zero counts in", sample_name))
    return()
  }
  
  log_counts <- log10(sample_counts)
  
  # Summary stats
  Q1 <- quantile(log_counts, 0.25)
  Q3 <- quantile(log_counts, 0.75)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  outliers <- names(sample_counts[log_counts < lower_bound])
  
  # Thresholds
  read_depth_thresh_20 <- 20
  read_depth_thresh_50 <- 50
  read_depth_thresh_500 <- 500
  threshold_20 <- log10(read_depth_thresh_20)
  threshold_50 <- log10(read_depth_thresh_50)
  threshold_500 <- log10(read_depth_thresh_500)
  
  # Counts below thresholds
  low_counts_20 <- sample_counts[sample_counts < read_depth_thresh_20]
  low_counts_50 <- sample_counts[sample_counts < read_depth_thresh_50]
  low_counts_500 <- sample_counts[sample_counts < read_depth_thresh_500]
  
  # Create layout: 2 rows (histogram + ECDF) + 1 for text
  layout(matrix(c(1, 2,
                  3, 0), 
                nrow = 2, byrow = TRUE),
         heights = c(0.5, 0.5))  # Adjust heights if needed
  
  ## Plot 1: Histogram
  hist(log_counts, breaks = 100, main = paste("Coverage Histogram:", sample_name),
       xlab = "log10(Sample Sum)", col = "lightgray", border = "white")
  abline(v = threshold_20, col = "blue", lty = 3, lwd = 2)
  abline(v = threshold_50, col = "#4052D6", lty = 3, lwd = 2)
  abline(v = threshold_500, col = "#1591EA", lty = 3, lwd = 2)
  
  ## Plot 2: ECDF
  plot(ecdf(log_counts), main = paste("ECDF:", sample_name),
       xlab = "log10(Sample Sum)", ylab = "Proportion of Samples")
  abline(v = Q1, col = "red", lty = 2, lwd = 2)
  abline(v = Q3, col = "darkgreen", lty = 2, lwd = 2)
  abline(v = lower_bound, col = "purple", lty = 1, lwd = 2)
  abline(v = threshold_20, col = "blue", lty = 3, lwd = 2)
  abline(v = threshold_50, col = "#4052D6", lty = 3, lwd = 2)
  abline(v = threshold_500, col = "#1591EA", lty = 3, lwd = 2)
  
  legend("bottomright",
         legend = c("Q1", "Q3", "Outlier Cutoff",  
                    "<20 reads", "<50 reads", "<500 reads"),
         col = c("red", "darkgreen", "purple", "blue", "#4052D6", "#1591EA"),
         lty = c(2, 2, 1, 3, 3, 3), lwd = 2, bty = "n")
  
  
  ## Plot 3: Text Panel
  plot.new()
  text(x = 0, y = 1, adj = c(0, 1), cex = 0.9, labels = paste0(
    "--- ", sample_name, " ---\n",
    "Mean reads / samples: ", round(mean(sample_counts), 2), "\n",
    "Median reads / samples: ", round(median(sample_counts, na.rm = TRUE), 2), "\n",
    "Min reads / samples: ", round(min(sample_counts, na.rm = TRUE), 2), "\n",
    "Max reads / samples: ", round(max(sample_counts, na.rm = TRUE), 2), "\n",
    "Samples with < 20 reads: ", length(low_counts_20), "\n",
    "Samples with < 50 reads: ", length(low_counts_50), "\n",
    "Samples with < 500 reads: ", length(low_counts_500), "\n",
    "IQR Outliers (log10 < ", round(lower_bound, 2), "): ", length(outliers)
  ))
}


otu_list <- list(
  "otu16m" = otu16m_mat,
  "otu16m_filt" = otu16m_mat_filt,
  "otu16n" = otu16n_mat,
  "otu16n_filt" = otu16n_mat_filt,
  "otu18m" = otu18m_mat,
  "otu18m_filt" = otu18m_nohost_mat_filt,
  "otu18n" = otu18n_mat,
  "otu18n_filt" = otu18n_nohost_mat_filt,
  "otu16ren" = otu16ren_mat,
  "otu16ren_filt" = otu16ren_mat_filt,
  "otu18ren" = otu18ren_mat, 
  "otu18ren_filt" = otu18ren_nohost_mat_filt)

pdf("01.1_qc_checks/coverage_plots.pdf", width = 8, height = 7)  # Adjust size as needed
#par(mfrow = c(13, 2))
for (name in names(otu_list)) {
  analyze_coverage(otu_list[[name]], sample_name = name)
}
dev.off()

## Ok now remove low read counts 

for (o in c("otu16m_mat_filt", "otu16n_mat_filt", "otu16ren_mat_filt",
            "otu18m_nohost_mat_filt", "otu18n_nohost_mat_filt", "otu18ren_nohost_mat_filt")) {
  
  df <- get(o)  # Retrieve the actual dataframe from the string name
  sample_counts <- colSums(df)
  Q1 <- quantile(sample_counts, 0.25)
  Q3 <- quantile(sample_counts, 0.75)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  outliers <- names(sample_counts[sample_counts < lower_bound])
  
  # Manual Thresholds
  read_depth_thresh_20 <- 20
  read_depth_thresh_50 <- 50
  read_depth_thresh_500 <- 500
  
  # 24 dataframes will be made below
  assign(paste0(o, "_20"),  df[, colnames(df)[which(colSums(df) >= read_depth_thresh_20)]])
  assign(paste0(o, "_50"),  df[, colnames(df)[which(colSums(df) >= read_depth_thresh_50)]])
  assign(paste0(o, "_500"), df[, colnames(df)[which(colSums(df) >= read_depth_thresh_500)]])
  assign(paste0(o, "_outlier"), df[, colnames(df)[which(colSums(df) >= lower_bound)]])
}

# # # # Step 3: # # # # 

# Remove ESVs not seen more than 2 times in at least 10% or 30% of the samples 

base_names <- c("otu16m_mat_filt", "otu1n_mat_filt", "otu16ren_mat_filt",
                "otu18m_nohost_mat_filt", "otu18n_nohost_mat_filt", "otu18ren_nohost_mat_filt")

suffixes <- c("_20", "_50", "_500", "_outlier")
all_names <- as.vector(outer(base_names, suffixes, paste0))

for (name in all_names) {
  cat("Processing:", name, "\n")
  temp <- get(name)
  
  sample_sums <- colSums(temp)
  df_plot <- data.frame(sample = names(sample_sums), reads = sample_sums)
  
  min_count <- 2
  min_sample_fraction_10 <- 0.10 # for 10% 
  min_sample_fraction_20 <- 0.30 # for 30 % 
  min_sample_count_10 <- ceiling(ncol(temp) * min_sample_fraction_10)  
  min_sample_count_30 <- ceiling(ncol(temp) * min_sample_fraction_30)
  
  p <- ggplot(df_plot, aes(x = reads)) +
    geom_histogram(binwidth = 1000, fill = "steelblue", color = "black") +
    geom_vline(xintercept = min_sample_count_10, linetype = "dashed", color = "red") +
    geom_vline(xintercept = min_sample_count_30, linetype = "dashed", color = "purple") +
    theme_minimal() +
    labs(title = paste("OTU count distribution per sample -", o),
         x = "Total OTUs per sample", y = "Frequency")
  
  print(p)
  assign(paste0(name, "_otu_cnt"), p)
  
  presence_matrix <- temp > min_count
  samples_per_taxa <- rowSums(presence_matrix) # samples per taxa meet that criterion
  otu_filtered_10 <- temp[samples_per_taxa >= min_sample_count_10, ]  # Keep taxa meeting threshold
  otu_filtered_30 <- temp[samples_per_taxa >= min_sample_count_30, ]
  dim(otu_filtered_10)
  dim(otu_filtered_30)
  
  # Print filtering summary 10 %
  cat("Filtering", name, " 10% ", "\n")
  cat("Original taxa:", nrow(temp), "\n")
  cat("Removed taxa with 10% thresh:", nrow(temp) - nrow(otu_filtered_10), "\n")
  cat("Remaining taxa:", nrow(otu_filtered_10), "\n\n")
  
  # Print filtering summary 30 %
  cat("Filtering", name, " 30% ", "\n")
  cat("Original taxa:", nrow(temp), "\n")
  cat("Removed taxa with 30% thresh:", nrow(temp) - nrow(otu_filtered_30), "\n")
  cat("Remaining taxa:", nrow(otu_filtered_30), "\n\n")
  
  # Optionally save the filtered matrix back with a new name
  assign(paste0(name, "_count_thresh_10"), otu_filtered_10)
  assign(paste0(name, "_count_thresh_30"), otu_filtered_30)
}

### Final loop 

# Reconstruct the name grid with base names and suffixes
base_names <- c("otu16m_mat_filt", "otu16n_mat_filt", "otu16ren_mat_filt",
                "otu18m_nohost_mat_filt", "otu18n_nohost_mat_filt", "otu18ren_nohost_mat_filt")
esv_list <- c(esv_16m_filt, esv_16n_filt, esv_16ren_filt,
              esv_18m_nohost_filt, esv_18m_nohost_filt, esv_18ren_nohost_filt)
suffixes <- c("_20", "_50", "_500", "_outlier")
reads_list <- c(reads_16, reads_16ren, reads_18, reads_18ren)

# Map each base_name to the correct reads_ object
reads_lookup <- list(
  "otu16m_mat_filt" = reads_16,
  "otu16n_mat_filt" = reads_16, 
  "otu16ren_mat_filt" = reads_16ren,
  "otu18m_nohost_mat_filt" = reads_18, 
  "otu18n_nohost_mat_filt" = reads_18,
  "otu18ren_nohost_mat_filt" = reads_18ren)

esv_list <- list(
  otu16m_mat_filt = esv_16m_filt,
  otu16n_mat_filt = esv_16n_filt,
  otu16ren_mat_filt = esv_16ren_filt,
  otu18m_nohost_mat_filt = esv_18m_nohost_filt,
  otu18n_nohost_mat_filt = esv_18n_nohost_filt,
  otu18ren_nohost_mat_filt = esv_18ren_nohost_filt)

df_info <- expand.grid(
  base = base_names,
  suffix = suffixes,
  stringsAsFactors = FALSE)

# Add the threshold versions
df_info$name_10 <- paste0(df_info$base, df_info$suffix) 
df_info$name_30 <- paste0(df_info$base, df_info$suffix)

head(df_info)# Combine and add base_name info

all_df_info <- rbind(
  data.frame(name = df_info$name_10, base = df_info$base, threshold = 10),
  data.frame(name = df_info$name_30, base = df_info$base, threshold = 30))

# Loop over each filtered name and use correct esv_list entry
for (i in seq_len(nrow(all_df_info))) {
  name <- all_df_info$name[i]
  base_name <- all_df_info$base[i]
  
  suffix <- sub(base_name, "", name, fixed = TRUE)  # extract suffix from the name
  threshold <- all_df_info$threshold[i]
  output_name <- paste0(base_name, suffix, "_thr", threshold)
  
  cat("Processing:", name, "\n")
  df <- get(name)
  toKeepSamples <- colnames(df)
  meta_temp <- meta_adj_filt %>% filter(SampleId %in% toKeepSamples)
  cat("Metadata entries retained:", name, ":", nrow(meta_temp), "\n")
  
  # Match ESV input using the base_name index
  base_index <- which(base_names == base_name)
  esv_filt_temp <- esv_list[[base_name]]

  cat("Initial taxonomy entries (ESV):", nrow(esv_filt_temp), "\n")
  listESV_temp <- rownames(df)
  esv_final_temp <- esv_filt_temp %>% filter(ESVId %in% listESV_temp)
  cat("Retained taxonomy entries (ESV):", nrow(esv_final_temp), "\n")
  cat("Taxonomy entries removed:", nrow(esv_filt_temp) - nrow(esv_final_temp), "\n")
  
  # Set up phyloseq objects
  meta_phyloseq_temp <- meta_temp %>%
    column_to_rownames("SampleId") 
  
  esv_phyloseq_temp <- esv_final_temp %>%
    mutate(Species = ScientificName) %>%
    select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    column_to_rownames(var="ESVId") %>% as.matrix()
  
  # Set up phyloseq objects
  meta_phyloseq_temp <- meta_temp %>%
    column_to_rownames("SampleId") 
  
  esv_phyloseq_temp <- esv_final_temp %>%
    mutate(Species = ScientificName) %>%
    select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    column_to_rownames(var="ESVId") %>% as.matrix()
  
  # Make phyloseq objects
  phyloseq_temp <- phyloseq(otu_table(df, taxa_are_rows=TRUE), 
                            tax_table(esv_phyloseq_temp), 
                            sample_data(meta_phyloseq_temp))
  
  #### Get sequences ####
  reads_df <- reads_lookup[[base_name]]
  seq_temp <- reads_df %>%
    filter(ESVId %in% rownames(df)) %>%
    select(ESVId, sequence)
  cat("Number of sequences retained:", nrow(seq_temp), "\n")
  
  #### Save items ####
  #write_csv(meta_temp, paste0("01.1_qc_checks/meta_",name,".csv"))
  write_csv(meta_temp, paste0("01.1_qc_checks/meta_", output_name, ".csv"))
  
  # Make otu table into a df
  otu_df_final <- df  %>%
    as.data.frame() %>% rownames_to_column(var="SampleId")
  #write_csv(otu_df_final, paste0("01.1_qc_checks/otu_",name,".csv"))
  write_csv(otu_df_final, paste0("01.1_qc_checks/otu_", output_name, ".csv"))
  
  # Save ESV info
  #write_csv(esv_final_temp, paste0("01.1_qc_checks/esv_", name,".csv"))
  write_csv(esv_final_temp, paste0("01.1_qc_checks/esv_", output_name, ".csv"))
  
  ####### R versions #############
  assign(paste0("meta_", output_name), value = meta_temp)
  save(list = paste0("meta_", output_name), file = paste0("01.1_qc_checks/meta_", output_name, ".rds"))
  
  assign(paste0("otu_", output_name), value = df)
  save(list = paste0("otu_", output_name), file = paste0("01.1_qc_checks/otu_", output_name, ".rds"))
  
  assign(paste0("esv_", output_name), value = esv_final_temp)
  save(list = paste0("esv_", output_name), file = paste0("01.1_qc_checks/esv_", output_name, ".rds"))
  
  assign(paste0("phyloseq_", output_name), value = phyloseq_temp)
  save(list = paste0("phyloseq_", output_name), file = paste0("01.1_qc_checks/phyloseq_", output_name, ".rds"))
  
  assign(paste0("seq_", output_name), value = seq_temp)
  save(list = paste0("seq_", output_name), file = paste0("01.1_qc_checks/seq_", output_name, ".rds"))
  
  #assign(paste0("meta_", name), value = meta_temp)
  #save(list=paste0("meta_",name), file = paste0("01.1_qc_checks/meta_",name,".rds"))
  #assign(paste0("otu_",name), value = df)
  #save(list=paste0("otu_",name), file = paste0("01.1_qc_checks/otu_",name,".rds"))
  #assign(paste0("esv_",name), value = esv_final_temp)
  #save(list=paste0("esv_",name), file = paste0("01.1_qc_checks/esv_",name,".rds"))
  #assign(paste0("phyloseq_",name), value = phyloseq_temp)
  #save(list=paste0("phyloseq_",name), file = paste0("01.1_qc_checks/phyloseq_",name,".rds"))
  #assign(paste0("seq_", name), value = seq_temp)
  #save(list=paste0("seq_", name), file = paste0("01.1_qc_checks/seq_",name,".rds"))
  
}

###### Make lists of eukaryotic classes, which we can partition into different types #######
simple18 <- full_join(esv18n, esv18ren) %>%
  select(Kingdom, Phylum, Class) %>%
  distinct()
#write.csv(simple18, file="01_process_and_clean_data/simple18.txt", row.names = FALSE, quote = FALSE)

