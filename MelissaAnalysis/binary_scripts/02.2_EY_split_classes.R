# Split 16S and 18S data into classes/categories
# Accepts a list of phyloseq objects and processes them flexibly

library(phyloseq)
library(tidyverse)

setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")

# Create output directory
if (!dir.exists("02_split_18_classes")) dir.create("02_split_18_classes")

#### Helper Functions ####

# Add Captive.Wild.probiotic factor to sample data
add_probiotic_factor <- function(ps) {
  newsampdata <- sample_data(ps) %>%
    data.frame() %>%
    mutate(Captive.Wild.probiotic = ifelse(
      Captive.Wild == "Wild, seized from traffickers",
      ifelse(Given.probiotic == "Given probiotics",
             "Wild, seized from traffickers\nGiven probiotics",
             as.character(Captive.Wild)),
      as.character(Captive.Wild))) %>%
    mutate(Captive.Wild.probiotic = factor(
      Captive.Wild.probiotic,
      levels = c("Captive", "Wild, free ranging",
                "Wild, seized from traffickers",
                "Wild, seized from traffickers\nGiven probiotics",
                "Crate 1: Unknown origin", "Crate 2: Unknown origin")))
  sample_data(ps) <- sample_data(newsampdata)
  return(ps)
}

# Process 18S data with categorization
process_18s_data <- function(ps, name, catdat, listPlantSpeciesFILLED = NULL) {
  cat("Processing 18S data:", name, "\n")
  
  # Add probiotic factor
  ps <- add_probiotic_factor(ps)
  
  # Add categorization
  # Ensure catdat has ESVId column
  if (!"ESVId" %in% colnames(catdat)) {
    if ("ESVId" %in% rownames(catdat)) {
      catdat <- catdat %>% rownames_to_column(var = "ESVId")
    } else {
      stop("catdat must have ESVId column or row names")
    }
  }
  
  allTaxWithCat <- tax_table(ps) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ESVId") %>%
    left_join(catdat, by = "ESVId")
  
  catTax <- allTaxWithCat %>%
    select(ESVId, Phylum, BroadCat2, BroadCat) %>%
    column_to_rownames(var = "ESVId") %>%
    as.matrix()
  
  ps_broadcat <- phyloseq(
    otu_table(ps, taxa_are_rows = TRUE),
    tax_table(catTax),
    sample_data(ps))
  
  ps_broadcat_glom <- tax_glom(ps_broadcat, taxrank = "BroadCat2")
  ps_broadcat_glom_nohost <- subset_taxa(
    ps_broadcat_glom,
    BroadCat2 != "Presumed_host" & BroadCat2 != "Human")
  
  ps_broadcat_glom_RA <- transform_sample_counts(
    ps_broadcat_glom,
    fun = function(x) x / sum(x))
  ps_broadcat_glom_RA_nohost <- transform_sample_counts(
    ps_broadcat_glom_nohost,
    fun = function(x) x / sum(x))
  
  # Create plots
  color_leg_broadcat2 <- c(
    Arachnid = "black", Presumed_host = "darkgrey",
    Mammal = "beige", Insect = "brown",
    MicroEukaryote = "lightblue", Mollusk = "darkblue",
    Plant = "darkolivegreen", Worm = "salmon")
  
  gg_broadcat <- plot_bar(ps_broadcat_glom_RA, fill = "BroadCat2") +
    geom_bar(aes(col = BroadCat2, fill = BroadCat2), stat = "identity", position = "stack") +
    facet_grid(~Captive.Wild.probiotic, scales = "free", drop = TRUE) +
    scale_fill_manual(values = color_leg_broadcat2) +
    scale_color_manual(values = color_leg_broadcat2) +
    labs(col = "Broad category", fill = "Broad category") +
    ylab("Relative Abundance") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  ggsave(
    filename = paste0("02_split_18_classes/gg_broadcat_", name, ".png"),
    gg_broadcat, height = 5, width = 10)
  
  gg_broadcat_nohost <- plot_bar(ps_broadcat_glom_RA_nohost, fill = "BroadCat2") +
    geom_bar(aes(col = BroadCat2, fill = BroadCat2), stat = "identity", position = "stack") +
    facet_grid(~Captive.Wild.probiotic, scales = "free", drop = TRUE) +
    scale_fill_manual(values = color_leg_broadcat2) +
    scale_color_manual(values = color_leg_broadcat2) +
    labs(col = "Broad category", fill = "Broad category") +
    ylab("Relative Abundance (excluding host and human DNA)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  ggsave(
    filename = paste0("02_split_18_classes/gg_broadcat_", name, "_nohost.png"),
    gg_broadcat_nohost, height = 5, width = 10)
  
  # Split into plants
  plant_ESVs <- allTaxWithCat %>%
    filter(BroadCat2 == "Plant") %>%
    pull(ESVId)
  ps_plants <- prune_taxa(taxa = plant_ESVs, ps)
  
  if (!is.null(listPlantSpeciesFILLED) && nrow(listPlantSpeciesFILLED) > 0) {
    taxtablenewplants <- tax_table(ps_plants) %>%
      as.data.frame() %>%
      rownames_to_column(var = "ESVId") %>%
      left_join(listPlantSpeciesFILLED %>% select(-Colour), by = "ESVId") %>%
      column_to_rownames(var = "ESVId") %>%
      as.matrix()
    tax_table(ps_plants) <- tax_table(taxtablenewplants)
    
    ps_plants_glom <- tax_glom(ps_plants, taxrank = "CommonName")
    ps_plants_glom_RA <- transform_sample_counts(
      ps_plants_glom,
      fun = function(x) x / sum(x))
    
    color_leg_plants <- listPlantSpeciesFILLED$Colour
    names(color_leg_plants) <- listPlantSpeciesFILLED$CommonName
    
    gg_plantsonly <- plot_bar(ps_plants_glom_RA, fill = "CommonName") +
      geom_bar(aes(col = CommonName, fill = CommonName), stat = "identity", position = "stack") +
      facet_grid(~Captive.Wild.probiotic, scales = "free", drop = TRUE) +
      scale_fill_manual(values = color_leg_plants) +
      scale_color_manual(values = color_leg_plants) +
      labs(col = "Plant common name", fill = "Plant common name") +
      ylab("Relative Abundance") +
      theme(axis.text.x = element_blank())
    
    ggsave(
      filename = paste0("02_split_18_classes/gg_plantsonly_", name, ".png"),
      gg_plantsonly, height = 5, width = 14)
    
    # Grey out rare plants
    toKeepColours <- names(which(
      rowSums(otu_table(ps_plants_glom_RA) > 0.2, na.rm = TRUE) > 0))
    color_leg_plants2 <- tax_table(ps_plants_glom_RA) %>%
      as.data.frame() %>%
      rownames_to_column(var = "ESVId") %>%
      mutate(Colour = ifelse(
        ESVId %in% toKeepColours,
        color_leg_plants[CommonName],
        "grey")) %>%
      pull(Colour)
    names(color_leg_plants2) <- tax_table(ps_plants_glom_RA) %>%
      as.data.frame() %>%
      pull(CommonName)
    
    gg_plantsonly_wgrey <- plot_bar(ps_plants_glom_RA, fill = "CommonName") +
      geom_bar(aes(col = CommonName, fill = CommonName), stat = "identity", position = "stack") +
      facet_grid(~Captive.Wild.probiotic, scales = "free", drop = TRUE) +
      scale_fill_manual(values = color_leg_plants2) +
      scale_color_manual(values = color_leg_plants2) +
      labs(col = "Plant common name", fill = "Plant common name") +
      ylab("Relative Abundance") +
      theme(axis.text.x = element_blank())
    
    ggsave(
      filename = paste0("02_split_18_classes/gg_plantsonly_wgrey_", name, ".png"),
      gg_plantsonly_wgrey, height = 5, width = 14)
  }
  
  # Split into microeukaryotes
  listMicroEuk <- tax_table(ps_broadcat) %>%
    as.data.frame() %>%
    filter(BroadCat2 == "MicroEukaryote") %>%
    rownames()
  
  ps_microeuk <- prune_taxa(x = ps, taxa = listMicroEuk)
  ps_microeuk_glom_RA <- transform_sample_counts(
    tax_glom(ps_microeuk, taxrank = "Phylum"),
    fun = function(x) x / sum(x))
  
  set.seed(0923)
  colorlegmicroeuk <- sample(
    colors()[-grep("white|grey|gray|beige", colors())],
    size = length(taxa_names(ps_microeuk_glom_RA)))
  
  gg_microeuk_phylum <- plot_bar(ps_microeuk_glom_RA, fill = "Phylum") +
    geom_bar(aes(col = Phylum, fill = Phylum), stat = "identity", position = "stack") +
    facet_grid(~Captive.Wild.probiotic, scales = "free", drop = TRUE) +
    scale_fill_manual(values = colorlegmicroeuk) +
    scale_color_manual(values = colorlegmicroeuk) +
    labs(col = "Phylum\n(Microeukaryotes only)", fill = "Phylum\n(Microeukaryotes only)") +
    ylab("Relative Abundance") +
    theme(axis.text.x = element_blank())
  
  ggsave(filename = paste0("02_split_18_classes/gg_microeuk_phylum_", name, ".png"),
    gg_microeuk_phylum, height = 5, width = 8)
  
  # Create nohost version
  nohost_ESVs <- allTaxWithCat %>%
    filter(!BroadCat2 %in% c("Presumed_host", "Mammal")) %>%
    pull(ESVId)
  ps_nohost <- prune_taxa(taxa = nohost_ESVs, ps)
  
  return(list(
    plants = ps_plants,
    microeuk = ps_microeuk,
    nohost = ps_nohost))
}

# Process 16S data (simpler - no plant/microeuk splits)
process_16s_data <- function(ps, name) {
  cat("Processing 16S data:", name, "\n")
  
  # Add probiotic factor
  ps <- add_probiotic_factor(ps)
  
  # For 16S, we might just want a basic nohost version
  # (16S typically doesn't have the same categorization needs)
  # Return the processed object
  return(list(nohost = ps))
}

# Detect if phyloseq is 16S or 18S based on name or taxonomy
detect_region <- function(name, ps) {
  # Check name first
  if (grepl("16", name, ignore.case = TRUE)) return("16S")
  if (grepl("18", name, ignore.case = TRUE)) return("18S")
  
  # Check taxonomy table for clues
  tax_tab <- tax_table(ps)
  if (ncol(tax_tab) > 0) {
    # 18S often has more diverse kingdoms
    kingdoms <- unique(tax_tab[, "Kingdom"])
    if (any(grepl("Plant|Animal|Fungi", kingdoms, ignore.case = TRUE))) {
      return("18S")
    }
  }
  
  # Default to 16S if uncertain
  return("16S")
}

#### Main Processing Function ####

# Process a list of phyloseq objects
process_phyloseq_list <- function(ps_list, 
                                   catdat_18s = NULL,
                                   listPlantSpeciesFILLED = NULL,
                                   merge_list = NULL) {
  # merge_list: named list of phyloseq objects to merge before processing
  # e.g., list(nano = list(ps1, ps2)) will merge ps1 and ps2 into "nano"
  
  results <- list()
  
  # Handle merging if specified
  # merge_list format: list(merged_name = list(ps1, ps2, ...))
  if (!is.null(merge_list)) {
    for (merged_name in names(merge_list)) {
      ps_to_merge <- merge_list[[merged_name]]
      if (is.list(ps_to_merge) && length(ps_to_merge) > 0) {
        merged_ps <- ps_to_merge[[1]]
        if (length(ps_to_merge) > 1) {
          for (i in 2:length(ps_to_merge)) {
            merged_ps <- merge_phyloseq(merged_ps, ps_to_merge[[i]])
          }
        }
        ps_list[[merged_name]] <- merged_ps
        cat("Merged", length(ps_to_merge), "objects into:", merged_name, "\n")
      }
    }
  }
  
  # Process each phyloseq object
  for (name in names(ps_list)) {
    ps <- ps_list[[name]]
    
    if (is.null(ps)) {
      warning(paste("Skipping", name, "- phyloseq object is NULL"))
      next
    }
    
    region <- detect_region(name, ps)
    cat("\n=== Processing", name, "(", region, ") ===\n")
    
    if (region == "18S") {
      if (is.null(catdat_18s)) {
        warning(paste("No categorization data provided for 18S, skipping", name))
        next
      }
      processed <- process_18s_data(ps, name, catdat_18s, listPlantSpeciesFILLED)
    } else {
      processed <- process_16s_data(ps, name)
    }
    
    # Save outputs
    for (subset_name in names(processed)) {
      output_name <- paste0(name, "_", subset_name)
      output_file <- paste0("02_split_18_classes/phyloseq_", output_name, ".rds")
      
      assign(output_name, processed[[subset_name]])
      save(list = output_name, file = output_file)
      cat("Saved:", output_file, "\n")
      
      results[[output_name]] <- processed[[subset_name]]
    }
  }
  
  return(results)
}

#### Load Data and Run ####

# Load categorization data (if exists)
catdat_18s <- NULL
if (file.exists("01_process_and_clean_data/simple18_labelled2.csv")) {
  catdat_18s <- read.csv("01_process_and_clean_data/simple18_labelled2.csv")
}

# Load plant species list (if exists)
listPlantSpeciesFILLED <- NULL
if (file.exists("02_split_18_classes/listPlantSpeciesFILLED.txt")) {
  listPlantSpeciesFILLED <- read.csv("02_split_18_classes/listPlantSpeciesFILLED.txt") %>%
    as.data.frame()
}

#### Example Usage ####

# Example 1: Load phyloseq objects from 01.1_qc_checks
# Uncomment and modify paths as needed

# ps_list <- list(
#   otu16m_mat_filt_20_thr10 = readRDS("01.1_qc_checks/phyloseq_otu16m_mat_filt_20_thr10.rds"),
#   otu16n_mat_filt_20_thr10 = readRDS("01.1_qc_checks/phyloseq_otu16n_mat_filt_20_thr10.rds"),
#   otu18m_nohost_mat_filt_20_thr10 = readRDS("01.1_qc_checks/phyloseq_otu18m_nohost_mat_filt_20_thr10.rds"),
#   otu18n_nohost_mat_filt_20_thr10 = readRDS("01.1_qc_checks/phyloseq_otu18n_nohost_mat_filt_20_thr10.rds")
# )
# 
# results <- process_phyloseq_list(ps_list, catdat_18s, listPlantSpeciesFILLED)

# Example 2: With merging (e.g., merge Renee's data with nanopore)
# merge_list <- list(
#   otu18nano = list(
#     readRDS("01.1_qc_checks/phyloseq_otu18n_nohost_mat_filt_20_thr10.rds"),
#     readRDS("01.1_qc_checks/phyloseq_otu18ren_nohost_mat_filt_20_thr10.rds")
#   )
# )
# 
# ps_list <- list(
#   otu18miseq = readRDS("01.1_qc_checks/phyloseq_otu18m_nohost_mat_filt_20_thr10.rds")
# )
# 
# results <- process_phyloseq_list(ps_list, catdat_18s, listPlantSpeciesFILLED, merge_list)

cat("\n=== Script ready! ===\n")
cat("To use, create a named list of phyloseq objects and call:\n")
cat("  results <- process_phyloseq_list(ps_list, catdat_18s, listPlantSpeciesFILLED)\n")
