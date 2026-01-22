# Combine LDA and ANCOMBC Binary Results
# Works with binary comparison results (Captive vs Wild, free ranging)
# Filters out "Wild, seized from traffickers" samples

#### Load libraries ####
pacman::p_load(tidyverse, vegan, MASS, phyloseq, tibble, ANCOMBC, ggplot2, coin,
               MicrobiotaProcess, patchwork, reshape2, ggnewscale, VennDiagram,
               UpSetR, gridExtra, grid, WGCNA)

# Set working directory to MelissaAnalysis (parent directory)
# This allows the script to work when run from binary_scripts/ or MelissaAnalysis/
if (basename(getwd()) == "binary_scripts") {
  setwd("..")
}
setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")

# Source utils.R for helper functions
source("utils.R")

# Create output directory
if (!dir.exists("05_taxa_driving_groups/combined_lda_ancombc/binary")) {
  dir.create("05_taxa_driving_groups/combined_lda_ancombc/binary", recursive = TRUE)
}

# Read in ANCOMBC and LDA binary results 
cat("Loading binary LDA results...\n")
load("05_taxa_driving_groups/LDA/binary/lda_binary_results.RData")

cat("Loading binary ANCOMBC results...\n")
load("05_taxa_driving_groups/ANCOMBC/binary/ancombc_binary_results.RData")

# Binary color scheme (only Captive and Wild, no Seized)
colorlist_binary <- c("Captive" = "#00AED7", "Wild" = "#FD9347")

### Format LDA Results ###
cat("\n=== Formatting LDA Results ===\n")
LDA_plot_dfs <- list()
for (taxrank in names(lda_diff_results)) { # (Species, Genus, Family)
  LDA_plot_dfs[[taxrank]] <- list()
  for (dataset in names(lda_diff_results[[taxrank]])) {
    cat("Processing LDA:", taxrank, "-", dataset, "\n")
    deres <- lda_diff_results[[taxrank]][[dataset]]
    
    # Select top features by LDAmean
    top_feat <- deres@result %>%
      arrange(desc(LDAmean)) %>%
      pull(f)
    
    # Create a new object containing only the top 30 features
    deres_top <- deres
    deres_top@result <- deres@result %>% filter(f %in% top_feat)
    deres_top@result$label <- gsub("__", ": ", deres_top@result$f)
    LDA_plot_dfs[[taxrank]][[dataset]] <- deres_top
  }
}

### Format ANCOMBC Results ###
cat("\n=== Formatting ANCOMBC Results ===\n")

# Loop through datasets and taxonomic levels
anc_plot_dfs <- list()
for (dataset_name in names(anc_results_trimmed)) {
  anc_plot_dfs[[dataset_name]] <- list()
  for (tax_level in names(anc_results_trimmed[[dataset_name]])) {
    
    cat("Fetching ANCOMBC results for:", dataset_name, "at", tax_level, "\n")
    res_obj <- anc_results_trimmed[[dataset_name]][[tax_level]]
    
    # Check if result object contains `res` element
    if (!"res" %in% names(res_obj)) {
      warning(paste("No 'res' in", dataset_name, tax_level))
      next
    }
    res_df <- res_obj$res
    
    # Check for valid results (don't filter yet - filter after reshaping)
    if (nrow(res_df) == 0) {
      cat("  No data in res_df\n")
      next
    }
    cat("  res_df size:", format(object.size(res_df), units = "MB"), 
        "with", nrow(res_df), "rows\n")
    
    # Reshape the lfc, q, and diff values separately
    lfc_df <- res_df %>%
      dplyr::select(taxon, starts_with("lfc_")) %>%
      pivot_longer(cols = -taxon, names_to = "group", values_to = "lfc") %>%
      mutate(group = sub("lfc_", "", group))
    
    q_df <- res_df %>%
      dplyr::select(taxon, starts_with("q_")) %>%
      pivot_longer(cols = -taxon, names_to = "group", values_to = "q") %>%
      mutate(group = sub("q_", "", group))
    
    diff_df <- res_df %>%
      dplyr::select(taxon, starts_with("diff_")) %>%
      pivot_longer(cols = -taxon, names_to = "group", values_to = "diff") %>%
      mutate(group = sub("diff_", "", group))
    
    # Join all into a single long-form dataframe
    plot_df <- reduce(list(lfc_df, q_df, diff_df), left_join, by = c("taxon", "group"))
    
    cat("  Before filtering: ", nrow(plot_df), "rows\n")
    
    # Filter significant results (diff == TRUE and q < 0.05)
    plot_df <- plot_df %>% filter(diff == TRUE, q < 0.05)
    
    cat("  After filtering (diff==TRUE & q<0.05): ", nrow(plot_df), "rows\n")
    
    if (nrow(plot_df) == 0) {
      cat("  Skipping - no significant taxa after filtering\n")
      next  # Skip if no significant taxa
    }
    
    # Clean and format for plotting (binary version - only Captive and Wild)
    plot_df$group <- str_replace_all(plot_df$group, "captive_wild_binary", "")
    plot_df$taxon <- gsub("_", " ", plot_df$taxon)
    plot_df <- plot_df %>%
      mutate(group = ifelse(group == "(Intercept)", "Captive", group))
    plot_df$group <- sub("^Captive\\.", "", plot_df$group)
    plot_df$group <- sub("WildWild", "Wild", plot_df$group)
    
    cat("  plot_df size:", format(object.size(plot_df), units = "MB"), "\n")
    anc_plot_dfs[[dataset_name]][[tax_level]] <- plot_df
  }
}

### Combine LDA and ANCOMBC Results ###
cat("\n=== Combining LDA and ANCOMBC Results ===\n")

# Show what datasets and tax levels are available
cat("Available ANCOMBC datasets:", paste(names(anc_plot_dfs), collapse = ", "), "\n")
if (length(anc_plot_dfs) > 0) {
  for (ds in names(anc_plot_dfs)) {
    cat("  ", ds, ":", paste(names(anc_plot_dfs[[ds]]), collapse = ", "), "\n")
  }
}

cat("Available LDA tax ranks:", paste(names(LDA_plot_dfs), collapse = ", "), "\n")
if (length(LDA_plot_dfs) > 0) {
  for (tr in names(LDA_plot_dfs)) {
    cat("  ", tr, ":", paste(names(LDA_plot_dfs[[tr]]), collapse = ", "), "\n")
  }
}

# Check if data exists
if (length(anc_plot_dfs) == 0) {
  stop("No ANCOMBC data available. All results were filtered out or empty.")
}

# Process all available dataset and taxonomic level combinations
# Get all unique combinations
all_combinations <- list()
for (dataset_name in names(anc_plot_dfs)) {
  for (tax_level in names(anc_plot_dfs[[dataset_name]])) {
    # Capitalize first letter to match stored names
    tax_level_cap <- paste0(toupper(substring(tax_level, 1, 1)), 
                            tolower(substring(tax_level, 2)))
    
    # Check if LDA data exists for this combination
    if (tax_level_cap %in% names(LDA_plot_dfs) &&
        dataset_name %in% names(LDA_plot_dfs[[tax_level_cap]])) {
      all_combinations[[length(all_combinations) + 1]] <- list(
        dataset = dataset_name,
        tax_level = tax_level,
        tax_level_cap = tax_level_cap
      )
    }
  }
}

cat("\nFound", length(all_combinations), "valid combinations to process\n")
cat("Processing all combinations...\n\n")

# Loop through all combinations
for (i in seq_along(all_combinations)) {
  combo <- all_combinations[[i]]
  dataset_name <- combo$dataset
  tax_level <- combo$tax_level
  tax_level_cap <- combo$tax_level_cap
  
  cat("\n", "=", rep("=", 60), "\n", sep = "")
  cat("Processing:", dataset_name, "-", tax_level, "\n")
  cat("=", rep("=", 60), "\n", sep = "")
  
  # Create descriptive base filename with clear dataset and tax level info
  # Extract dataset type (16m, 16n, 18m, 18n) from dataset name
  dataset_short <- gsub("phyloseq", "", dataset_name)
  base_filename <- paste0("05_taxa_driving_groups/combined_lda_ancombc/binary/",
                          "LDA_ANCOMBC_", dataset_short, "_", tax_level)
  
  # From ANCOMBC (use capitalized tax level name)
  ancom <- anc_plot_dfs[[dataset_name]][[tax_level_cap]]
ancom$taxon <- as.character(ancom$taxon)

# From LDA
  lda_obj <- LDA_plot_dfs[[tax_level_cap]][[dataset_name]]
lda <- lda_obj@result
lda$taxon <- gsub("^[a-z]:\\s*", "", lda$label)
lda$taxon <- sub("^[a-z]__+", "", lda$taxon)

# ANCOMBC dataframe
ancom_df <- ancom %>%
    dplyr::select(taxon, group, value = lfc) %>%
  mutate(stat = "LFC")

# LDA dataframe (binary version - only Captive and Wild)
lda_df <- lda %>%
    dplyr::select(taxon, group = captive_wild_binary, value = LDAmean) %>%
  mutate(stat = "LDA")

# Ensure group names match
lda_df <- lda_df %>%
  mutate(group = as.character(group),
         group = ifelse(group == "Wild, free ranging", "Wild", group))

combined_df <- bind_rows(ancom_df, lda_df)

# Filter to taxa that appear in both LFC and LDA
combined_df <- combined_df %>%
  dplyr::group_by(taxon) %>%
  filter(n() > 1) %>%  # Appears in both LFC and LDA
  dplyr::ungroup()

cat("Combined dataframe has", nrow(combined_df), "rows\n")
cat("Unique taxa in combined:", length(unique(combined_df$taxon)), "\n")

# Create combined plot
p_combined <- ggplot(combined_df, aes(x = reorder(taxon, value), y = value, fill = group)) +
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), aes(alpha = stat)) +
  coord_flip() +
  facet_wrap(~stat, scales = "free_x") +  # Separate panels for LFC and LDA
  scale_fill_manual(values = colorlist_binary) +
  scale_alpha_manual(values = c("LFC" = 1, "LDA" = 0.6)) +
  labs(title = "Differential Abundance (Binary): ANCOMBC LFC and LDA Scores",
         subtitle = paste(dataset_name, "-", tax_level),
       x = "Taxon", y = "Score (LFC or LDA)", fill = "Group", alpha = "Stat Type") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave(
    filename = paste0(base_filename, "_combined_plot.png"),
  plot = p_combined,
  width = 16,
  height = 10,
  dpi = 300
)

  cat("Combined plot saved:", paste0(base_filename, "_combined_plot.png"), "\n")
  
  ### Save Input Data as CSV Files ###
  cat("Saving CSV files...\n")
  
  # 1. Save LDA results (from the LDA object used for plotting)
  lda_results_to_save <- lda_obj@result
  write.csv(lda_results_to_save, 
            file = paste0(base_filename, "_LDA_input_data.csv"),
            row.names = FALSE)
  cat("  Saved LDA input data:", paste0(base_filename, "_LDA_input_data.csv"), "\n")
  
  # 2. Save ANCOMBC results (from anc_plot_dfs)
  ancom_results_to_save <- anc_plot_dfs[[dataset_name]][[tax_level_cap]]
  write.csv(ancom_results_to_save,
            file = paste0(base_filename, "_ANCOMBC_input_data.csv"),
            row.names = FALSE)
  cat("  Saved ANCOMBC input data:", paste0(base_filename, "_ANCOMBC_input_data.csv"), "\n")
  
  # 3. Save combined dataframe
  write.csv(combined_df,
            file = paste0(base_filename, "_combined_LDA_ANCOMBC.csv"),
            row.names = FALSE)
  cat("  Saved combined dataframe:", paste0(base_filename, "_combined_LDA_ANCOMBC.csv"), "\n")

### Create Separate Combined Plots (p1, p2, p3) ###
  cat("  Creating separate combined plots (p1, p2, p3)...\n")

  # Prepare LDA object for plotting (use the same lda_obj from above)
  lda_for_plots <- lda_obj

# Clean feature names in LDA object
  if (!is.null(lda_for_plots@originalD)) {
    original_names <- colnames(lda_for_plots@originalD)
  cleaned_names <- sub("^[a-z]__+", "", original_names)
    colnames(lda_for_plots@originalD) <- cleaned_names
}

# Clean taxon names in result
  lda_for_plots@result$taxon <- gsub("^[a-z]:\\s*", "", lda_for_plots@result$f)
  lda_for_plots@result$taxon <- sub("^[a-z]__+", "", lda_for_plots@result$f)
  lda_for_plots@result$f <- gsub("^[a-z]:\\s*", "", lda_for_plots@result$f)
  lda_for_plots@result$f <- sub("^[a-z]__+", "", lda_for_plots@result$f)

# Get feature lists
  # LDA features - use the cleaned names from @result$f (these match @originalD column names after cleaning)
  lda_feature_list <- unique(lda_for_plots@result$f)
  
  # ANCOMBC features - convert spaces back to underscores to match LDA format
  anc_feature_list_raw <- unique(anc_plot_dfs[[dataset_name]][[tax_level_cap]]$taxon)
  anc_feature_list <- gsub(" ", "_", anc_feature_list_raw)  # Convert spaces to underscores

  # Find overlap (comparing cleaned names)
  overlap_lda <- lda_feature_list[lda_feature_list %in% anc_feature_list]
  overlap_anc <- anc_feature_list[anc_feature_list %in% lda_feature_list]
  
  unique_to_ancombc <- setdiff(anc_feature_list, lda_feature_list)
  unique_to_lda <- setdiff(lda_feature_list, anc_feature_list)
  
  # Use LDA feature names for the ordered list (these match @originalD column names)
  ordered_features_all <- c(overlap_lda, unique_to_ancombc, unique_to_lda)
  
  # Verify that all features exist in originalD
  available_features <- colnames(lda_for_plots@originalD)
  missing_features <- setdiff(ordered_features_all, available_features)
  if (length(missing_features) > 0) {
    cat("    Warning: Some features not found in originalD:", paste(missing_features, collapse = ", "), "\n")
    ordered_features_all <- intersect(ordered_features_all, available_features)
  }
  
  cat("    Overlap:", length(overlap_lda), "taxa\n")
  cat("    Unique to ANCOMBC:", length(unique_to_ancombc), "taxa\n")
  cat("    Unique to LDA:", length(unique_to_lda), "taxa\n")
  cat("    Total features for plotting:", length(ordered_features_all), "\n")
  
  # Save feature overlap information (only overlapping features)
  feature_overlap_info <- data.frame(
    feature = overlap_lda,
    in_LDA = TRUE,
    in_ANCOMBC = TRUE,
    category = "Overlap"
  )
  write.csv(feature_overlap_info,
            file = paste0(base_filename, "_feature_overlap.csv"),
            row.names = FALSE)
  cat("    Saved feature overlap info (", nrow(feature_overlap_info), " overlapping features):", 
      paste0(base_filename, "_feature_overlap.csv"), "\n")

# Plot p1 - Boxplot (LDA)
p1 <- make_p1_boxplot(
    obj = lda_for_plots,
  featurelist = ordered_features_all,
  colorlist = colorlist_binary,
  box_notch = FALSE,
  xlabtext = "Relative Abundance"
)

# Plot p1_ordered - Ordered boxplot
p1_ord <- make_p1_boxplot_ordered(
    obj = lda_for_plots,
  featurelist = ordered_features_all,
  colorlist = colorlist_binary,
  box_notch = FALSE,
  xlabtext = "Relative Abundance"
)

# Plot p2 - Effect size (LDA)
  nodedfres_all <- lda_for_plots@result
missing_in_p2 <- setdiff(ordered_features_all, nodedfres_all$f)
if (length(missing_in_p2) > 0) {
  na_rows <- data.frame(f = missing_in_p2)
  na_rows[, setdiff(colnames(nodedfres_all), "f")] <- NA
  nodedfres_all <- rbind(nodedfres_all, na_rows)
}
nodedfres_all$f <- factor(nodedfres_all$f, levels = rev(ordered_features_all))

p2 <- make_p2_effectsize(
    obj = lda_for_plots,
  nodedfres = nodedfres_all, 
  colorlist = colorlist_binary, 
  xlim_range = c(2.5, 5.5)
)

p2 <- p2 + theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank())

# Plot p3 - ANCOMBC
  # Convert ordered_features_all (underscores) back to spaces for ANCOMBC matching
  ordered_features_all_spaces <- gsub("_", " ", ordered_features_all)
  anc_df <- anc_plot_dfs[[dataset_name]][[tax_level_cap]]
missing_in_p3 <- setdiff(ordered_features_all_spaces, anc_df$taxon)
if (length(missing_in_p3) > 0) {
  na_rows <- data.frame(taxon = missing_in_p3)
  na_rows[, setdiff(colnames(anc_df), "taxon")] <- NA
  anc_df <- rbind(anc_df, na_rows)
}

anc_df$taxon <- factor(anc_df$taxon, levels = rev(ordered_features_all_spaces))

p3 <- make_p3_ancombc(
  ancombc_df = anc_df, 
    dataset_name = dataset_short,
    tax_level = tax_level,
  featurelist = ordered_features_all_spaces,  # Use space-separated names for ANCOMBC
  colorlist = colorlist_binary
)

p3 <- p3 + theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 plot.title = element_blank())

# Combine p1, p2, p3
p_combined_separate <- p1 | p2 | p3

ggsave(
    filename = paste0(base_filename, "_separate_plots_p1_p2_p3.png"),
  plot = p_combined_separate,
  width = 20,
  height = 10,
  dpi = 300
)

  cat("    Separate combined plot saved:", paste0(base_filename, "_separate_plots_p1_p2_p3.png"), "\n")
  
  cat("\nCompleted processing:", dataset_name, "-", tax_level, "\n")
}

cat("\n", "=", rep("=", 60), "\n", sep = "")
cat("=== Analysis Complete ===\n")
cat("Processed", length(all_combinations), "combinations\n")
cat("Output files saved to: 05_taxa_driving_groups/combined_lda_ancombc/binary/\n")
cat("=", rep("=", 60), "\n", sep = "")