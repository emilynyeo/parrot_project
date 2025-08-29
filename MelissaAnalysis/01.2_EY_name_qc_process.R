#!bin/bash

#### Load libraries ####
pacman::p_load(tidyverse, vegan, MASS, phyloseq, tibble, ANCOMBC, ggplot2, coin,
               MicrobiotaProcess, patchwork, reshape2, ggnewscale, VennDiagram,
               UpSetR, gridExtra, grid, WGCNA, indicspecies, lubridate, scales,
               calecopal)

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
  select(SampleId, everything()) %>% 
  filter(SampleId != "S082895") %>% 
  filter(SampleId != "S082806")

# NANOPORE AND MISEQ from original runs
esv_16 <- read.csv("./00_sequence_data/JVB3218_16S/JVB3218-16S-esv-data.csv")
reads_16 <- read.csv("./00_sequence_data/JVB3218_16S/JVB3218-16S-read-data.csv")
samples_16 <- read.csv("./00_sequence_data/JVB3218_16S/JVB3218-samples.csv")

esv_18 <- read.csv("./00_sequence_data/JVB3218_18S_616/JVB3218-18S_616-esv-data.csv")
reads_18 <- read.csv("./00_sequence_data/JVB3218_18S_616/JVB3218-18S_616-read-data.csv")
samples_18 <- read.csv("./00_sequence_data/JVB3218_18S_616/JVB3218-samples.csv")

# Remove odd looking samples S082895 and S082806
colnames(reads_16) <- gsub("^s0","S0", colnames(reads_16)) # Make it start with capital S
colnames(reads_18) <- gsub("^s0","S0", colnames(reads_18)) 

reads16 <- reads_16 %>% dplyr::select(!starts_with("S082895")) %>% select(!(starts_with("S082806")))
reads18 <- reads_18 %>% dplyr::select(!starts_with("S082895")) %>% select(!(starts_with("S082806")))

samples_16 <- samples_16 %>% 
                   dplyr::filter(SampleId != "S082895") %>% 
                          filter(SampleId != "S082806")

samples_18 <- samples_18 %>% 
                   dplyr::filter(SampleId != "S082895") %>% 
                          filter(SampleId != "S082806")

# Remove zero sum read rows  (the ones ending in 0)
reads_16 %>% select(ends_with(".0")) %>% colSums()
reads_18 %>% select(ends_with(".0")) %>% colSums()

# Split reaads into nano and miseq
otu16m_df <- reads_16 %>% select(ESVId, starts_with("S0") & ends_with(".1")) # Run 1 are miseq
otu16n_df <- reads_16 %>% select(ESVId, starts_with("S0") & ends_with(".2")) # Run 2 is Nano

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

# Remove zero sum read rows  (the ones ending in 0)
reads_16 %>% select(ends_with(".0")) %>% colSums()
reads_18 %>% select(ends_with(".0")) %>% colSums()

table(rowSums(otu16m_df[,-1])>0)
otu16m_df <- otu16m_df[rowSums(otu16m_df[,-1])>0,] # 16 miseq 3638 removed

table(rowSums(otu16n_df[,-1])>0)
otu16n_df <- otu16n_df[rowSums(otu16n_df[,-1])>0,] # 16 nano 2474 removed

table(rowSums(otu18m_df[,-1])>0) # 1586 removed
otu18m_df <- otu18m_df[rowSums(otu18m_df[,-1])>0,] 

table(rowSums(otu18n_df[,-1])>0) # 10292 removed 
otu18n_df <- otu18n_df[rowSums(otu18n_df[,-1])>0,] 

### Remove suffixes from samplenames 
colnames(otu16m_df) <- gsub("[.]1","", colnames(otu16m_df))
colnames(otu18m_df) <- gsub("[.]2", "", gsub("[.]1","", colnames(otu18m_df)))
colnames(otu16n_df) <- gsub("[.]2","", colnames(otu16n_df))
colnames(otu18n_df) <- gsub("[.]2", "", gsub("[.]3","", colnames(otu18n_df)))

### Make into an OTU matrix
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

### Clean up Meta data 

# Incorporate Metadata Read Data 
ReadDepth16m <- colSums(otu16m_mat) %>% as.data.frame() %>%
                rownames_to_column(var="SampleId") %>% rename(ReadDepth16m=".")

ReadDepth18m <- colSums(otu18m_mat) %>% as.data.frame() %>%
                rownames_to_column(var="SampleId") %>% rename(ReadDepth18m=".")

ReadDepth16n <- colSums(otu16n_mat) %>% as.data.frame() %>%
                rownames_to_column(var="SampleId") %>% rename(ReadDepth16n=".")
  
ReadDepth18n <- colSums(otu18n_mat) %>% as.data.frame() %>%
                rownames_to_column(var="SampleId") %>% rename(ReadDepth18n=".")

meta_adj <- meta1 %>% 
  left_join(full_join(full_join(ReadDepth16m,ReadDepth18m), 
                                full_join(ReadDepth16n, ReadDepth18n)))%>% 
  dplyr::mutate(Captive.Wild = factor(Captive.Wild, 
                               levels=c("Captive", "Wild, free ranging", 
                                        "Wild, seized from traffickers")))

# Clean up dates 

# Standardize date format
meta_adj$Date.collected
meta_adj %>% select(Country, Date.collected) %>% table()

# Dashes look like SA dates (MM-DD-YYYY) whereas Nigerian = DD/MM/YY: Let's fix those...
reformatDates <- as.data.frame(matrix(nrow=length(meta_adj$Date.collected), 
                                      ncol=5, 
                                      dimnames = list(NULL, c("Date.collected","Day","Month","Year","DateAdj"))))

# extract dates from strings and re-standardize them
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

# Time sincw captive 
meta_adj$Seized.date
meta_adj$Seized.date <- ymd(as.Date(meta_adj$Seized.date, format = "%d-%m-%Y"))
meta_adj <- meta_adj %>% mutate(TimeSinceSeizure = as.Date(DateAdj)-Seized.date)

### Taxonomic Assignments ### 
#tax16 <- reads_16 %>% select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) 
tax16 <- reads_16 %>% 
  select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(ScientificName = ifelse(Species!="", paste0(Genus, " ", Species),
                                 ifelse(Genus !="", paste0("Genus:", Genus),
                                        ifelse(Family!="", paste0("Family:", Family),
                                               ifelse(Order!="", paste0("Order:", Order),
                                                      ifelse(Class!="", paste0("Class:", Class),
                                                             ifelse(Phylum!="", paste0("Phylum:",Phylum),
                                                                    ifelse(Kingdom!="", paste0("Kingdom:", Kingdom),"Unassigned")))))))) %>% 
  mutate(across(Kingdom:Species, ~na_if(., ""))) %>% # Replace "" with NA 
  mutate(Kingdom = coalesce(Kingdom, "Unassigned"), # Fill from highest to lowest rank
    Phylum  = coalesce(Phylum, Kingdom),
    Class   = coalesce(Class, Phylum),
    Order   = coalesce(Order, Class),
    Family  = coalesce(Family, Order),
    Genus   = coalesce(Genus, Family),
    Species = coalesce(Species, Genus))

tax16m <- tax16 %>% filter(ESVId %in% otu16m_df$ESVId)
tax16n <- tax16 %>% filter(ESVId %in% otu16n_df$ESVId)

#tax18 <- reads_18 %>% select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) 
tax18 <- reads_18 %>% select(ESVId, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  mutate(ScientificName = ifelse(Species!="", paste0(Genus, " ", gsub("_"," ", Species)),
                                 ifelse((Genus !="") & (Genus!=Family), paste0("Genus:", Genus),
                                        ifelse((Family!="")& (Family!=Order), paste0("Family:", Family),
                                               ifelse((Order!="")&(Order!=Class), paste0("Order:", Order),
                                                      ifelse((Class!="")&(Class!=Phylum), paste0("Class:", Class),
                                                             ifelse((Phylum!="") & (Phylum!=Kingdom), paste0("Phylum:",Phylum),
                                                                    ifelse(Kingdom!="", paste0("Kingdom:", Kingdom),"Unassigned")))))))) %>% 
  mutate(across(Kingdom:Species, ~na_if(., ""))) %>% # Replace "" with NA 
  mutate(Kingdom = coalesce(Kingdom, "Unassigned"), # Fill from highest to lowest rank
         Phylum  = coalesce(Phylum, Kingdom),
         Class   = coalesce(Class, Phylum),
         Order   = coalesce(Order, Class),
         Family  = coalesce(Family, Order),
         Genus   = coalesce(Genus, Family),
         Species = coalesce(Species, Genus))

tax18m <- tax18 %>% filter(ESVId %in% otu18m_df$ESVId)
tax18n <- tax18 %>% filter(ESVId %in% otu18n_df$ESVId)

# Remove 16s chloroplasts & mitochondria : Chloroplasts are Order ; Mitochondria are family level
esv16_toremove <- tax16 %>% filter(Family == "Mitochondria"| Order =="Chloroplast") %>% pull(ESVId)
esv_16_filt <- tax16 %>% filter(!ESVId %in% esv16_toremove)

# Remove host DNA and others for 18S 
esv18_toremove <- tax18 %>% filter(Kingdom=="Bacteria" | Kingdom=="Archaea"| Kingdom=="") %>% pull(ESVId)
esv_18_filt <- tax18 %>% filter(!ESVId %in% esv18_toremove)

# remove birds 
allBirdESVs <- read_csv("00_metadata/allBirdESVs.csv")
bird_esv <- allBirdESVs$ESVId
esv_16_filt <- esv_16_filt %>% filter(!ESVId %in% bird_esv)
esv_18_filt <- esv_18_filt %>% filter(!ESVId %in% bird_esv)

# Filter out ESVs from tables based on this
otu16m_mat_filt <- otu16m_mat[which(rownames(otu16m_mat) %in% esv_16_filt$ESVId),]
otu16n_mat_filt <- otu16n_mat[which(rownames(otu16n_mat) %in% esv_16_filt$ESVId),]

otu18m_mat_filt <- otu18m_mat[which(rownames(otu18m_mat) %in% esv_18_filt$ESVId),]
otu18n_mat_filt <- otu18n_mat[which(rownames(otu18n_mat) %in% esv_18_filt$ESVId),]

# Version with no host DNA and no 'mammal' DNA. Not sure what that is.. if it's actually people or not
esv_18_nohost_filt <- esv_18_filt %>% filter(!is.na(Class), Class != "Actinopterygii", Class !="Aves", Class !="Mammalia")
otu18m_nohost_mat_filt <- otu18m_mat[which(rownames(otu18m_mat) %in% esv_18_nohost_filt$ESVId),]
otu18n_nohost_mat_filt <- otu18n_mat[which(rownames(otu18n_mat) %in% esv_18_nohost_filt$ESVId),]

# Recalculate read depth after host DNA filtering
ReadDepth16m_filt <- colSums(otu16m_mat_filt) %>%
  as.data.frame() %>% rownames_to_column(var = "SampleId") %>%
  rename(ReadDepth16m_filt = ".")

ReadDepth16n_filt <- colSums(otu16n_mat_filt) %>%
  as.data.frame() %>% rownames_to_column(var = "SampleId") %>%
  rename(ReadDepth16n_filt = ".")

ReadDepth18m_filt <- colSums(otu18m_nohost_mat_filt) %>%
  as.data.frame() %>% rownames_to_column(var = "SampleId") %>%
  rename(ReadDepth18m_filt = ".")

ReadDepth18n_filt <- colSums(otu18n_nohost_mat_filt) %>%
  as.data.frame() %>% rownames_to_column(var = "SampleId") %>%
  rename(ReadDepth18n_filt = ".")

meta_adj_filt <- meta_adj %>% 
  left_join(full_join(full_join(ReadDepth16m_filt,ReadDepth18m_filt), 
                                full_join(ReadDepth16n_filt, ReadDepth18n_filt)))

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
  read_depth_thresh_500 <- 500
  threshold_20 <- log10(read_depth_thresh_20)
  threshold_500 <- log10(read_depth_thresh_500)
  # Counts below thresholds
  low_counts_20 <- sample_counts[sample_counts < read_depth_thresh_20]
  low_counts_500 <- sample_counts[sample_counts < read_depth_thresh_500]
  # Create layout: 2 rows (histogram + ECDF) + 1 for text
  layout(matrix(c(1, 2, 3, 0), nrow = 2, byrow = TRUE),
         heights = c(0.5, 0.5))  # Adjust heights if needed
  ## Plot 1: Histogram
  hist(log_counts, breaks = 100, main = paste("Coverage Histogram:", sample_name),
       xlab = "log10(Sample Sum)", col = "lightgray", border = "white")
  abline(v = threshold_20, col = "#4052D6", lty = 3, lwd = 2)
  abline(v = threshold_500, col = "#1591EA", lty = 3, lwd = 2)
  ## Plot 2: ECDF
  plot(ecdf(log_counts), main = paste("ECDF:", sample_name),
       xlab = "log10(Sample Sum)", ylab = "Proportion of Samples")
  abline(v = Q1, col = "red", lty = 2, lwd = 2)
  abline(v = Q3, col = "darkgreen", lty = 2, lwd = 2)
  abline(v = lower_bound, col = "purple", lty = 1, lwd = 2)
  abline(v = threshold_20, col = "blue", lty = 3, lwd = 2)
  abline(v = threshold_500, col = "#1591EA", lty = 3, lwd = 2)
  legend("bottomright", legend = c("Q1", "Q3", "Outlier Cutoff",  
                    "<20 reads", "<500 reads"),
         col = c("red", "darkgreen", "purple", "blue", "#1591EA"),
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
  "otu18n_filt" = otu18n_nohost_mat_filt)

pdf("01.1_qc_checks/coverage_plots_no_parrot.pdf", width = 8, height = 7)  # Adjust size as needed
#par(mfrow = c(13, 2))
for (name in names(otu_list)) {
  analyze_coverage(otu_list[[name]], sample_name = name)
}
dev.off()

## Ok now remove low read counts 
for (o in c("otu16m_mat_filt", "otu16n_mat_filt",
            "otu18m_nohost_mat_filt", "otu18n_nohost_mat_filt")) {
  df <- get(o)  # Retrieve the actual dataframe from the string name
  sample_counts <- colSums(df)
  Q1 <- quantile(sample_counts, 0.25)
  Q3 <- quantile(sample_counts, 0.75)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  outliers <- names(sample_counts[sample_counts < lower_bound])
  # Manual Thresholds
  read_depth_thresh_20 <- 20
  read_depth_thresh_500 <- 500
  read_depth_max_9500 <- 9500
  # dataframes will be made below
  assign(paste0(o, "_20"),  df[, colnames(df)[which(colSums(df) >= read_depth_thresh_20)]])
  assign(paste0(o, "_500"), df[, colnames(df)[which(colSums(df) >= read_depth_thresh_500)]])
  assign(paste0(o, "_outlier"), df[, colnames(df)[which(colSums(df) >= lower_bound)]])
  assign(paste0(o, "_max_9500"), df[, colnames(df)[which(colSums(df) <= read_depth_max_9500)]])
  assign(paste0(o, "data.frame(name = df_info$name_30, base = df_info$base, threshold = 30))"), 
         df[, colnames(df)[which(colSums(df) <= read_depth_max_9500)]])
}


# Remove ESVs not seen more than 2 0r 10 times in at least 10% or 30% of the samples 
base_names <- c("otu16m_mat_filt", "otu16n_mat_filt",
                "otu18m_nohost_mat_filt", "otu18n_nohost_mat_filt")
suffixes <- c("_20", "_500", "_outlier", "_max_9500")
all_names <- as.vector(outer(base_names, suffixes, paste0))

for (name in all_names) {
  cat("Processing:", name, "\n")
  temp <- get(name)
  sample_sums <- colSums(temp)
  df_plot <- data.frame(sample = names(sample_sums), reads = sample_sums)
  min_count <- 2
  min_sample_fraction_10 <- 0.10 # for 10% 
  min_sample_fraction_30 <- 0.30  # for 30%
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
  otu_filtered_2cnt_10pct <- temp[samples_per_taxa >= min_sample_count_10, ]  
  otu_filtered_2cnt_30pct <- temp[samples_per_taxa >= min_sample_count_30, ]
  # Print filtering summary 10 %
  cat("Filtering", name, " 10% ", "\n")
  cat("Original taxa:", nrow(temp), "\n")
  cat("Removed 2 count taxa with 10% thresh:", nrow(temp) - nrow(otu_filtered_2cnt_10pct), "\n")
  cat("Remaining 2count taxa 10%:", nrow(otu_filtered_2cnt_10pct), "\n\n")
  # Print filtering summary 30%
  cat("Filtering", name, " 30% ", "\n")
  cat("Original taxa:", nrow(temp), "\n")
  cat("Removed 2 count taxa with 30% thresh:", nrow(temp) - nrow(otu_filtered_2cnt_30pct), "\n")
  cat("Remaining 2 count taxa with 30% thresh: ", nrow(otu_filtered_2cnt_30pct), "\n\n")
  # Optionally save the filtered matrix back with a new name
  assign(paste0(name, "_2count_thresh_10"), otu_filtered_2cnt_10pct)
  assign(paste0(name, "_2count_thresh_30"), otu_filtered_2cnt_30pct)
}

### Final loop 

# Reconstruct the name grid with base names and suffixes
base_names <- c("otu16m_mat_filt", "otu16n_mat_filt",
                "otu18m_nohost_mat_filt", "otu18n_nohost_mat_filt")
esv_list <- c(esv_16_filt, esv_18_nohost_filt)
suffixes <- c("_20", "_500", "_outlier", "_max_9500")
reads_list <- c(reads_16, reads_18)
all_phyloseq_names <- c() 

# Map each base_name to the correct reads_ object
reads_lookup <- list(
  "otu16m_mat_filt" = reads_16,
  "otu16n_mat_filt" = reads_16, 
  "otu18m_nohost_mat_filt" = reads_18, 
  "otu18n_nohost_mat_filt" = reads_18)

esv_list <- list(
  otu16m_mat_filt = esv_16_filt,
  otu16n_mat_filt = esv_16_filt,
  otu18m_nohost_mat_filt = esv_18_nohost_filt,
  otu18n_nohost_mat_filt = esv_18_nohost_filt)

df_info <- expand.grid(
  base = base_names,
  suffix = suffixes,
  stringsAsFactors = FALSE)

# Add the threshold versions just doing 10% 
df_info$name_10 <- paste0(df_info$base, df_info$suffix) 
#df_info$name_30 <- paste0(df_info$base, df_info$suffix)

head(df_info)# Combine and add base_name info

all_df_info <- rbind(data.frame(name = df_info$name_10, 
                                base = df_info$base, threshold = 10))
  #data.frame(name = df_info$name_30, base = df_info$base, threshold = 30))

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
  
  tax_tab <- phyloseq::tax_table(as.matrix(esv_phyloseq_temp))
  otu_tab <- phyloseq::otu_table(df, taxa_are_rows=TRUE)
  samp_tab <- phyloseq::sample_data(meta_phyloseq_temp)
  
  # Make phyloseq objects
  phyloseq_temp <- phyloseq(otu_tab, tax_tab,
                            samp_tab)
  
  #### Get sequences ####
  reads_df <- reads_lookup[[base_name]]
  seq_temp <- reads_df %>%
    filter(ESVId %in% rownames(df)) %>%
    select(ESVId, sequence)
  cat("Number of sequences retained:", nrow(seq_temp), "\n")
  
  ####### R versions #############
  assign(paste0("meta_", output_name), value = meta_temp)
  save(list = paste0("meta_", output_name), file = paste0("01.1_qc_checks/meta_", output_name, ".rds"))
  
  assign(paste0("otu_", output_name), value = df)
  save(list = paste0("otu_", output_name), file = paste0("01.1_qc_checks/otu_", output_name, ".rds"))
  
  assign(paste0("esv_", output_name), value = esv_final_temp)
  save(list = paste0("esv_", output_name), file = paste0("01.1_qc_checks/esv_", output_name, ".rds"))
  
  assign(paste0("phyloseq_", output_name), value = phyloseq_temp)
  save(list = paste0("phyloseq_", output_name), file = paste0("01.1_qc_checks/phyloseq_", output_name, ".rds"))
  phyloseq_name <- paste0("phyloseq_", output_name)
  all_phyloseq_names <- c(all_phyloseq_names, phyloseq_name)
  
  assign(paste0("seq_", output_name), value = seq_temp)
  save(list = paste0("seq_", output_name), file = paste0("01.1_qc_checks/seq_", output_name, ".rds"))
  
}

saveRDS(all_phyloseq_names, file = "01.1_qc_checks/all_phyloseq_object_names.rds")



