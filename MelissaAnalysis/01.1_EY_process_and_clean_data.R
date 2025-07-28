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


# STEP 2: 

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
  log_counts <- log10(sample_counts)
  
  # Calculate summary stats
  Q1 <- quantile(log_counts, 0.25)
  Q3 <- quantile(log_counts, 0.75)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  outliers <- names(sample_counts[log_counts < lower_bound])
  # Define log threshold for read depth of 20
  read_depth_thresh_20 <- 20
  read_depth_thresh_50 <- 50
  read_depth_thresh_500 <- 500
  threshold_20 <- log10(read_depth_thresh_20)
  threshold_50 <- log10(read_depth_thresh_40)
  threshold_500 <- log10(read_depth_thresh_500)
  
  # Histogram
  hist(log_counts, breaks = 100, main = paste("Coverage Histogram:", sample_name),
       xlab = "log10(Sample Sum)", col = "lightgray", border = "white")
  abline(v = threshold_20, col = "blue", lty = 3, lwd = 2)
  abline(v = threshold_50, col = "#4052D6", lty = 3, lwd = 2)
  abline(v = threshold_500, col = "#1591EA", lty = 3, lwd = 2)
  
  # ECDF plot with IQR lines
  plot(ecdf(log_counts), main = paste("ECDF:", sample_name),
       xlab = "log10(Sample Sum)", ylab = "Proportion of Samples")
  abline(v = Q1, col = "red", lty = 2, lwd = 2)
  abline(v = Q3, col = "darkgreen", lty = 2, lwd = 2)
  abline(v = lower_bound, col = "purple", lty = 1, lwd = 2)
  abline(v = threshold_20, col = "blue", lty = 3, lwd = 2)
  abline(v = threshold_50, col = "#4052D6", lty = 3, lwd = 2)
  abline(v = threshold_500, col = "#1591EA", lty = 3, lwd = 2)
  
  legend("bottomright",
         legend = c("Q1", "Q3", "Lower Bound (Outlier Cutoff)", 
                    paste("Read Depth <", threshold_20),
                    paste("Read Depth <", threshold_50),
                    paste("Read Depth <", threshold_500)),
         col = c("red", "darkgreen", "purple", "blue", "#4052D6", "#1591EA"),
         lty = c(2, 2, 1, 3), lwd = 2, bty = "n")
  
  # Print stats
  low_counts_20 <- sample_counts[sample_counts < read_depth_threshold]
  low_counts_1000 <- sample_counts[sample_counts < 1000]
  low_counts_2000 <- sample_counts[sample_counts < 2000]
  
  cat("\n---", sample_name, "---\n")
  cat("Samples with <", threshold_20, "reads:", length(low_counts_20), "\n")
  cat("Samples with <1000 reads:", length(low_counts_1000), "\n")
  cat("Samples with <2000 reads:", length(low_counts_2000), "\n")
  cat("Names of <", threshold_20, " samples:\n"); print(names(low_counts_20))
  cat("Names of <1000 samples:\n"); print(names(low_counts_1000))
  cat("Names of <2000 samples:\n"); print(names(low_counts_2000))
  cat("Outliers by IQR method (log10 < ", round(lower_bound, 2), "):\n", sep = "")
  print(outliers)
}

otu_list <- list(
  "otu16m" = otu16m_mat,
  "otu16m_filt" = otu16m_mat_filt,
  "otu16n" = otu16n_mat,
  "otu16n_filt" = otu16n_mat_filt,
  "otu18m" = otu18m_mat,
  "otu18m_filt" = otu18m_mat_filt,
  "otu18n" = otu18n_mat,
  "otu18n_filt" = otu18n_mat_filt,
  "otu16ren" = otu16ren_mat,
  "otu16ren_filt" = otu16ren_mat_filt,
  "otu18ren" = otu18ren_mat, 
  "otu18ren_filt" = otu18ren_nohost_mat_filt)

pdf("01.1_qc_checks/coverage_plots.pdf", width = 8, height = 30)  # Adjust size as needed
par(mfrow = c(12, 2))
for (name in names(otu_list)) {
  analyze_coverage(otu_list[[name]], sample_name = name)
}
dev.off()


# Step 3: 

#### Looking at and removing (or not) rare sequences ####
# NOTE: phyloseq recommends NOT removing rare sequences... so I'm keeping them here for now, but we can fix this later

# 16S (old)
# Prevalence of OTUs
for ( o in c("otu16n_mat_filt", "otu16m_mat_filt", "otu16ren_mat_filt")) {
  # Get current working file
  temp <- get(o)
  ### NOTE: I removed all filters because Phyloseq recommends keeping all low-count asvs
  # Prevalence 
  prevGrThan2 <- which(rowSums(temp>0) > 0) # Used to be 2; keeping 0 here
  temp2 <- temp[prevGrThan2,]
  # Turn counts of <5 to zero
  lessThan5perSample <- temp2 <0 # Used to be 5, keeping 9 here
  temp2[lessThan5perSample] <- 0
  # How many have total reads less than 10?
  temp2 <- temp2[rowSums(temp2)>0,colSums(temp2)>0]
  assign(paste0(o,"2"), temp2)
}

# 16S
# taxa not seen more than 2 times in at least 10% of the samples can be removed? 
# Prevalence of OTUs
for (o in c("otu16n_mat_filt", "otu16m_mat_filt", "otu16ren_mat_filt")) {
  # Get current working file
  temp <- get(o)# Plot histogram of total reads per sample (column sums)
  sample_sums <- colSums(temp)
  df_plot <- data.frame(sample = names(sample_sums), reads = sample_sums)
  
  # Remove counts 
  min_count <- 2
  min_sample_fraction_10 <- 0.10
  min_sample_fraction_20 <- 0.20
  min_sample_count_10 <- ceiling(ncol(temp) * min_sample_fraction_10)  # of samples (cols) needed
  min_sample_count_20 <- ceiling(ncol(temp) * min_sample_fraction_20)
  
  p <- ggplot(df_plot, aes(x = reads)) +
    geom_histogram(binwidth = 1000, fill = "steelblue", color = "black") +
    geom_vline(xintercept = min_sample_count_10, linetype = "dashed", color = "red") +
    geom_vline(xintercept = min_sample_count_20, linetype = "dashed", color = "purple") +
    theme_minimal() +
    labs(title = paste("OTU count distribution per sample -", o),
         x = "Total OTUs per sample", y = "Frequency")
  print(p)
  
  presence_matrix <- temp > min_count
  samples_per_taxa <- rowSums(presence_matrix) # samples per taxa meet that criterion
  otu_filtered_10 <- temp[samples_per_taxa >= min_sample_count_10, ]  # Keep taxa meeting threshold
  otu_filtered_20 <- temp[samples_per_taxa >= min_sample_count_20, ]
  dim(otu_filtered_10)
  dim(otu_filtered_20)
  
  # Print filtering summary 10 %
  cat("Filtering", o, " 10% ", "\n")
  cat("Original taxa:", nrow(temp), "\n")
  cat("Removed taxa with 10% thresh:", nrow(temp) - nrow(otu_filtered_10), "\n")
  cat("Remaining taxa:", nrow(otu_filtered_10), "\n\n")
  
  # Print filtering summary 20 %
  cat("Filtering", o, " 20% ", "\n")
  cat("Original taxa:", nrow(temp), "\n")
  cat("Removed taxa with 20% thresh:", nrow(temp) - nrow(otu_filtered_20), "\n")
  cat("Remaining taxa:", nrow(otu_filtered_20), "\n\n")
  
  # Optionally save the filtered matrix back with a new name
  assign(paste0(o, "_count_thresh_10"), otu_filtered_10)
  assign(paste0(o, "_count_thresh_20"), otu_filtered_10)
}

# 18S
for ( o in c("otu18n_mat_filt", "otu18m_mat_filt", 
             "otu18ren_mat_filt", "otu18n_nohost_mat_filt", 
             "otu18m_nohost_mat_filt", "otu18ren_nohost_mat_filt")) {
  # Get current working file
  temp <- get(o)
  ### NOTE: I removed all filters because Phyloseq recommends keeping all low-count asvs
  # Prevalence of OTUs
  prev18GrThan2 <- which(rowSums(temp>0) > 0) # Used to be 2, keeping 0 here
  temp2 <- temp[prev18GrThan2,]
  # Turn counts of <5 to zero
  lessThan5perSampl18 <- temp2<0 # Used to be 5, keeping 0 here
  temp2[lessThan5perSampl18] <- 0
  # How many have total reads less than 10?
  which(rowSums(temp2)<0) # None, but let's remove those anyway #10?
  temp2 <- temp2[rowSums(temp2)>0,colSums(temp2)>0]
  assign(paste0(o,"2"), temp2)
}

#### Looking at coverage depth to remove low-coverage samples ####
hist(log10(colSums(otu16m_mat_filt2)), breaks = 100) # Handful of bad samples
hist(log10(colSums(otu16n_mat_filt2)), breaks = 100) # Handful of bad samples
hist(log10(colSums(otu16ren_mat_filt2)), breaks = 100) # Handful of bad samples
hist(log10(colSums(otu18m_mat_filt2)), breaks = 100) # For some reason, 18S distribution is very zero-skewed
hist(log10(colSums(otu18n_mat_filt2)), breaks = 100) # For some reason, 18S distribution is very zero-skewed
hist(log10(colSums(otu18ren_mat_filt2)), breaks = 100) # For some reason, 18S distribution is very zero-skewed
hist(log10(colSums(otu18m_nohost_mat_filt2)), breaks = 100) # 
hist(log10(colSums(otu18n_nohost_mat_filt2)), breaks = 100) # 
hist(log10(colSums(otu18ren_nohost_mat_filt2)), breaks = 100) # 


####  Look at decay curves for sample cutoff ####
# Miseq
plot(ecdf(x =colSums(otu16m_mat_filt) ))
abline(v=4000)
# Nanopore
plot(ecdf(x =colSums(otu16n_mat_filt2) ))
abline(v=17000)
# Ren
plot(ecdf(x =colSums(otu16ren_mat_filt2) ))
abline(v=17000)
# Make cutoff of 17000 for 16S nanopore and 4000 for miseq
cutoff16Smiseq= 4000
cutoff16Snano = 17000


# Look at decay curves for 18s
plot(ecdf(x = log10(colSums(otu18m_mat_filt2) )))
abline(v=log10(100))
plot(ecdf(x = log10(colSums(otu18n_mat_filt2) )))
abline(v=log10(100))
plot(ecdf(x = log10(colSums(otu18ren_mat_filt2) )))
abline(v=log10(100))

# Are the two in miseq and nano the same ones?
which(log10(colSums(otu18n_mat_filt2))<2)
which(log10(colSums(otu18m_mat_filt2))<2)

# Nope. Huh. Super sample dependent.
colSums(otu18m_mat_filt2 [,names(which(log10(colSums(otu18n_mat_filt2))<2))])

# There doesn't seem to be a "good" cutoff for 18S? Very zero-skewed.
# Are those random 2 bad samples the same as the 16S ones?
meta_adj %>% 
  filter(SampleId %in% colnames(otu16n_mat_filt2)[which(colSums(otu16n_mat_filt2)<cutoff16Snano)]) %>%
  pull(SampleId)
which(log10(colSums(otu18n_mat_filt2))<2)
# No, those are completely different samples! Hm. Let's set a cutoff of 100 reads for 18S for now
cutoff18S <- 100

# Filter OTU tables
# which(rownames(otu18n_nohost_mat_filt2) == "ESV_006911") # CHECK THEY ARE STILL GONE

### Final filtering of samples ####
# 16S
for ( o in c("m","n","ren")) {
  otumattemp <- get(paste0("otu16",o,"_mat_filt2"))
  cat("\n-----", o, "-----\n")
  cat("Initial number of samples (columns):", ncol(otumattemp), "\n")
  cat("Initial number of ESVs (rows):", nrow(otumattemp), "\n")
  cutoff <- ifelse(o %in% c("otu16m"), cutoff16Smiseq, cutoff16Snano)
  
  # Filter out samples
  toKeepSamples <- colnames(otumattemp)[which(colSums(otumattemp)>=cutoff)]
  otu_mat_final_temp <- otumattemp[, toKeepSamples]
  cat("Number of samples retained:", length(toKeepSamples), "\n")
  cat("Number of samples removed:", ncol(otumattemp) - length(toKeepSamples), "\n")
  cat("Number of ESVs retained after sample filtering:", nrow(otu_mat_final_temp), "\n")
  
  # Filter out metadata to match OTU table
  meta_temp <- meta_adj %>% filter(SampleId %in% toKeepSamples)
  cat("Metadata entries retained:", nrow(meta_temp), "\n")
  
  # Get list of ESVs to keep
  listESV_temp <- rownames(otu_mat_final_temp)
  esv_filt_temp <- get(paste0("esv_16",o,"_filt"))
  cat("Initial taxonomy entries (ESV):", nrow(esv_filt_temp), "\n")
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
  
  # Make phyloseq objects
  phyloseq_temp <- phyloseq(otu_table(otu_mat_final_temp, taxa_are_rows=TRUE), tax_table(esv_phyloseq_temp), sample_data(meta_phyloseq_temp))
  
  #### Get sequences ####
  seq_temp <- reads_16 %>%
    filter(ESVId %in% rownames(otu_mat_final_temp)) %>%
    select(ESVId, sequence)
  cat("Number of sequences retained:", nrow(seq_temp), "\n")
  #### Save items ####
  #write_csv(meta_temp, paste0("01_process_and_clean_data/meta16",o,".csv"))
  
  # Make otu table into a df
  otu_df_final <- otu_mat_final_temp  %>%
    as.data.frame() %>% rownames_to_column(var="SampleId")
  #write_csv(otu_df_final, paste0("01_process_and_clean_data/otu16",o,".csv"))
  
  # Save ESV info
  #write_csv(esv_final_temp, paste0("01_process_and_clean_data/esv16",o,".csv"))
  
  ####### R versions #############
  assign(paste0("meta16",o), value = meta_temp)
  #save(list=paste0("meta16",o), file = paste0("01_process_and_clean_data/meta16",o,".rds"))
  assign(paste0("otu16",o), value = otu_mat_final_temp)
  #save(list=paste0("otu16",o), file = paste0("01_process_and_clean_data/otu16",o,".rds"))
  assign(paste0("esv16",o), value = esv_final_temp)
  #save(list=paste0("esv16",o), file = paste0("01_process_and_clean_data/esv16",o,".rds"))
  assign(paste0("phyloseq16",o), value = phyloseq_temp)
  #save(list=paste0("phyloseq16",o), file = paste0("01_process_and_clean_data/phyloseq16",o,".rds"))
  assign(paste0("seq16",o), value = seq_temp)
  #save(list=paste0("seq16",o), file = paste0("01_process_and_clean_data/seq16",o,".rds"))
  
}

# for ( o in c("m","n","ren", "m_nohost","n_nohost","ren_nohost")) {
#   
#   cat("\n### Processing 18S:", o, "###\n")
#   
#   # Load OTU matrix
#   otumattemp <- get(paste0("otu18",o,"_mat_filt2"))
#   cat("Initial number of samples (columns):", ncol(otumattemp), "\n")
#   cat("Initial number of ESVs (rows):", nrow(otumattemp), "\n")
#   
#   # Set cutoff
#   cutoff <- cutoff18S
#   
#   # Filter samples
#   toKeepSamples <- colnames(otumattemp)[which(colSums(otumattemp) >= cutoff)]
#   cat("Number of samples retained:", length(toKeepSamples), "\n")
#   cat("Number of samples removed:", ncol(otumattemp) - length(toKeepSamples), "\n")
#   
#   # Subset OTU table
#   otu_mat_final_temp <- otumattemp[, toKeepSamples]
#   cat("Number of ESVs retained after sample filtering:", nrow(otu_mat_final_temp), "\n")
#   cat("Number of ESVs removed:", nrow(otumattemp) - nrow(otu_mat_final_temp), "\n")
#   
#   # Filter metadata
#   meta_temp <- meta_adj %>% filter(SampleId %in% toKeepSamples)
#   cat("Metadata entries retained:", nrow(meta_temp), "\n")
#   
#   # Get list of ESVs to keep
#   listESV_temp <- rownames(otu_mat_final_temp)
#   esv_filt_temp <- get(paste0("esv_18",o,"_filt"))
#   cat("Initial taxonomy entries (ESVs):", nrow(esv_filt_temp), "\n")
#   
#   esv_final_temp <- esv_filt_temp %>% filter(ESVId %in% listESV_temp)
#   cat("Taxonomy entries retained:", nrow(esv_final_temp), "\n")
#   cat("Taxonomy entries removed:", nrow(esv_filt_temp) - nrow(esv_final_temp), "\n")
#   
#   # Set up phyloseq objects
#   meta_phyloseq_temp <- meta_temp %>% column_to_rownames("SampleId")
#   
#   esv_phyloseq_temp <- esv_final_temp %>%
#     mutate(Species = ScientificName) %>%
#     select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
#     column_to_rownames(var = "ESVId") %>% as.matrix()
#   
#   # Make phyloseq object
#   phyloseq_temp <- phyloseq(
#     otu_table(otu_mat_final_temp, taxa_are_rows = TRUE),
#     tax_table(esv_phyloseq_temp),
#     sample_data(meta_phyloseq_temp)
#   )
#   
#   #### Get sequences ####
#   seq_temp <- reads_18 %>%
#     filter(ESVId %in% rownames(otu_mat_final_temp)) %>%
#     select(ESVId, sequence)
#   cat("Number of sequences retained:", nrow(seq_temp), "\n")
#   
#   #### Save items (assign only) ####
#   assign(paste0("meta18", o), value = meta_temp)
#   assign(paste0("otu18", o), value = otu_mat_final_temp)
#   assign(paste0("esv18", o), value = esv_final_temp)
#   assign(paste0("phyloseq18", o), value = phyloseq_temp)
#   assign(paste0("seq18", o), value = seq_temp)
# }


# 18S
# "otu18m","otu18n","otu18ren","otu18m_nohost","otu18n_nohost","otu18ren_nohost"
for ( o in c("m","n","ren", "m_nohost","n_nohost","ren_nohost")) {
  otumattemp <- get(paste0("otu18",o,"_mat_filt2"))
  cutoff <- cutoff18S
  
  # Filter out samples
  toKeepSamples <- colnames(otumattemp)[which(colSums(otumattemp)>=cutoff)]
  otu_mat_final_temp <- otumattemp[, toKeepSamples]
  
  # Filter out metadata to match OTU table
  meta_temp <- meta_adj %>% filter(SampleId %in% toKeepSamples)
  
  # Get list of ESVs to keep
  listESV_temp <- rownames(otu_mat_final_temp)
  esv_filt_temp <- get(paste0("esv_18",o,"_filt"))
  esv_final_temp <- esv_filt_temp %>% filter(ESVId %in% listESV_temp)
  
  # Set up phyloseq objects
  meta_phyloseq_temp <- meta_temp %>%
    column_to_rownames("SampleId") 
  
  esv_phyloseq_temp <- esv_final_temp %>%
    mutate(Species = ScientificName) %>%
    select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    column_to_rownames(var="ESVId") %>% as.matrix()
  
  # Make phyloseq objects
  phyloseq_temp <- phyloseq(otu_table(otu_mat_final_temp, taxa_are_rows=TRUE), tax_table(esv_phyloseq_temp), sample_data(meta_phyloseq_temp))
  
  #### Get sequences ####
  seq_temp <- reads_18 %>%
    filter(ESVId %in% rownames(otu_mat_final_temp)) %>%
    select(ESVId, sequence)
  
  #### Save items ####
  #write_csv(meta_temp, paste0("01_process_and_clean_data/meta18",o,".csv"))

  # Make otu table into a df
  otu_df_final <- otu_mat_final_temp  %>%
    as.data.frame() %>% rownames_to_column(var="SampleId")
  #write_csv(otu_df_final, paste0("01_process_and_clean_data/otu18",o,".csv"))

  # Save ESV info
  #write_csv(esv_final_temp, paste0("01_process_and_clean_data/esv18",o,".csv"))
  
  ####### R versions #############
  assign(paste0("meta18",o), value = meta_temp)
  #save(list=paste0("meta18",o), file = paste0("01_process_and_clean_data/meta18",o,".rds"))
  assign(paste0("otu18",o), value = otu_mat_final_temp)
  #save(list=paste0("otu18",o), file = paste0("01_process_and_clean_data/otu18",o,".rds"))
  assign(paste0("esv18",o), value = esv_final_temp)
  #save(list=paste0("esv18",o), file = paste0("01_process_and_clean_data/esv18",o,".rds"))
  assign(paste0("phyloseq18",o), value = phyloseq_temp)
  #save(list=paste0("phyloseq18",o), file = paste0("01_process_and_clean_data/phyloseq18",o,".rds"))
  assign(paste0("seq18",o), value = seq_temp)
  #save(list=paste0("seq18",o), file = paste0("01_process_and_clean_data/seq18",o,".rds"))

}

###### Make lists of eukaryotic classes, which we can partition into different types #######
simple18 <- full_join(esv18n, esv18ren) %>%
  select(Kingdom, Phylum, Class) %>%
  distinct()
#write.csv(simple18, file="01_process_and_clean_data/simple18.txt", row.names = FALSE, quote = FALSE)

