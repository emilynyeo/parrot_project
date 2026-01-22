# ANCOMBC2 Differential Abundance Analysis - Binary Comparison (Captive vs Wild, free ranging)
# Filters out "Wild, seized from traffickers" for binary comparison
# Accepts a list of phyloseq objects

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

# Create output directory
if (!dir.exists("05_taxa_driving_groups/ANCOMBC/binary")) {
  dir.create("05_taxa_driving_groups/ANCOMBC/binary", recursive = TRUE)
}

#### Helper Functions ####

# Filter phyloseq object to remove "Wild, seized from traffickers"
filter_binary_groups <- function(ps, group_var = "Captive.Wild") {
  # Get sample data - convert to regular data.frame
  sd <- as.data.frame(sample_data(ps))
  
  # Check if group variable exists
  if (!group_var %in% colnames(sd)) {
    warning(paste("Group variable", group_var, "not found in sample data"))
    return(ps)
  }
  
  # Filter out "Wild, seized from traffickers" - use base R subsetting to avoid issues
  keep_rows <- sd[[group_var]] != "Wild, seized from traffickers" & !is.na(sd[[group_var]])
  sd_filtered <- sd[keep_rows, , drop = FALSE]
  
  # Check if we have samples left
  if (nrow(sd_filtered) == 0) {
    warning("No samples remaining after filtering")
    return(ps)
  }
  
  # Create binary group variable
  sd_filtered$captive_wild_binary <- factor(
    sd_filtered[[group_var]],
    levels = c("Captive", "Wild, free ranging"),
    labels = c("Captive", "Wild")
  )
  
  # Get sample names to keep
  samples_to_keep <- rownames(sd_filtered)
  
  # Subset phyloseq object to only include filtered samples
  ps_filtered <- prune_samples(
    sample_names(ps) %in% samples_to_keep,
    ps
  )
  
  # Update sample data with filtered data (convert back to sample_data)
  sample_data(ps_filtered) <- sample_data(sd_filtered)
  
  # Remove samples with zero counts
  ps_filtered <- prune_samples(sample_sums(ps_filtered) > 0, ps_filtered)
  
  return(ps_filtered)
}

# Run ANCOMBC2 analysis on a single phyloseq object at a specific taxonomic level
run_ancombc_analysis <- function(ps,
                                 tax_level = "Genus",
                                 group_var = "captive_wild_binary",
                                 p_adj_method = "holm",
                                 prv_cut = 0.01,
                                 lib_cut = 50,
                                 s0_perc = 0.05) {
  
  cat("Running", tax_level, "level ANCOMBC2...\n")
  
  # Run ANCOMBC2
  result <- ancombc2(
    data = ps,
    tax_level = tax_level,
    fix_formula = group_var,
    rand_formula = NULL,
    p_adj_method = p_adj_method,
    pseudo_sens = TRUE,
    prv_cut = prv_cut,
    lib_cut = lib_cut,
    s0_perc = s0_perc,
    group = group_var,
    struc_zero = TRUE,
    neg_lb = TRUE
  )
  
  return(result)
}

# Create ANCOMBC plot
create_ancombc_plot <- function(res_df,
                                dataset_name,
                                tax_level,
                                output_name = NULL,
                                save_plot = TRUE,
                                colorlist = c("#00AED7", "#FD9347")) {
  
  # Filter significant results (q < 0.05 and diff == TRUE)
  res_df_filtered <- res_df %>%
    filter(`q_(Intercept)` < 0.05)
  
  if (nrow(res_df_filtered) == 0) {
    warning(paste("No significant taxa for", dataset_name, "at", tax_level))
    return(NULL)
  }
  
  # Reshape the lfc, q, and diff values separately
  lfc_df <- res_df_filtered %>%
    select(taxon, starts_with("lfc_")) %>%
    pivot_longer(cols = -taxon, names_to = "group", values_to = "lfc") %>%
    mutate(group = sub("lfc_", "", group))
  
  q_df <- res_df_filtered %>%
    select(taxon, starts_with("q_")) %>%
    pivot_longer(cols = -taxon, names_to = "group", values_to = "q") %>%
    mutate(group = sub("q_", "", group))
  
  diff_df <- res_df_filtered %>%
    select(taxon, starts_with("diff_")) %>%
    pivot_longer(cols = -taxon, names_to = "group", values_to = "diff") %>%
    mutate(group = sub("diff_", "", group))
  
  # Join all into a single long-form dataframe
  plot_df <- reduce(list(lfc_df, q_df, diff_df), left_join, by = c("taxon", "group"))
  
  # Filter significant results
  plot_df <- plot_df %>% filter(diff == TRUE, q < 0.05)
  
  if (nrow(plot_df) == 0) {
    warning(paste("No significant differential taxa for", dataset_name, "at", tax_level))
    return(NULL)
  }
  
  # Clean and format for plotting
  plot_df$group <- str_replace_all(plot_df$group, "captive_wild_binary", "")
  plot_df$taxon <- gsub("_", " ", plot_df$taxon)
  plot_df <- plot_df %>%
    mutate(group = ifelse(group == "(Intercept)", "Captive", group))
  plot_df$group <- sub("^Captive\\.", "", plot_df$group)
  plot_df$group <- sub("WildWild", "Wild", plot_df$group)
  
  # Create plot
  p <- ggplot(plot_df, aes(x = reorder(taxon, lfc), y = lfc, fill = group)) +
    geom_col(position = position_dodge(width = 0.9)) +
    coord_flip() +
    labs(title = paste("Significant DA (ANCOMBC2):", dataset_name, "-", tax_level),
         x = "Taxon", y = "Log Fold Change", fill = "Group") +
    scale_fill_manual(values = colorlist) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14))
  
  # Save plot if requested
  if (save_plot && !is.null(output_name)) {
    ggsave(
      filename = paste0("05_taxa_driving_groups/ANCOMBC/binary/", output_name, ".png"),
      plot = p,
      height = 15,
      width = 18,
      dpi = 300
    )
  }
  
  return(list(plot = p, plot_df = plot_df))
}

#### Main Processing Function ####

# Process a list of phyloseq objects
process_ancombc_binary <- function(ps_list,
                                    tax_levels = c("Family", "Genus", "Species"),
                                    group_var = "Captive.Wild",
                                    p_adj_method = "holm",
                                    prv_cut = 0.01,
                                    lib_cut = 50,
                                    s0_perc = 0.05,
                                    prv_cut_species = 0.1,  # Different cutoff for species
                                    colorlist = c("#00AED7", "#FD9347"),
                                    save_results = TRUE) {
  
  # Initialize nested list to store results
  anc_results <- list()
  anc_plots <- list()
  
  # Filter each phyloseq object
  ps_list_filtered <- list()
  for (name in names(ps_list)) {
    cat("Filtering", name, "...\n")
    ps <- ps_list[[name]]
    
    if (is.null(ps)) {
      warning(paste("Skipping", name, "- phyloseq object is NULL"))
      next
    }
    
    # Filter to binary groups
    ps_filtered <- filter_binary_groups(ps, group_var = group_var)
    
    # Check if we have at least 2 groups
    sd_check <- as.data.frame(sample_data(ps_filtered))
    if (!"captive_wild_binary" %in% colnames(sd_check)) {
      warning(paste("Skipping", name, "- could not create binary group variable"))
      next
    }
    
    groups_present <- unique(sd_check$captive_wild_binary)
    if (length(groups_present) < 2) {
      warning(paste("Skipping", name, "- need at least 2 groups, found:", 
                    paste(groups_present, collapse = ", ")))
      next
    }
    
    cat("  Groups present:", paste(groups_present, collapse = ", "), "\n")
    cat("  Samples:", nsamples(ps_filtered), "\n")
    
    ps_list_filtered[[name]] <- ps_filtered
  }
  
  if (length(ps_list_filtered) == 0) {
    stop("No phyloseq objects could be processed after filtering")
  }
  
  # Run ANCOMBC analysis for each taxonomic level and dataset
  for (name in names(ps_list_filtered)) {
    cat("\n=== Processing", name, "===\n")
    ps_filtered <- ps_list_filtered[[name]]
    anc_results[[name]] <- list()
    anc_plots[[name]] <- list()
    
    for (tax_level in tax_levels) {
      cat("  Running", tax_level, "level...\n")
      
      tryCatch({
        # Use different prv_cut for species level
        current_prv_cut <- ifelse(tax_level == "Species", prv_cut_species, prv_cut)
        
        # Run ANCOMBC analysis
        result <- run_ancombc_analysis(
          ps = ps_filtered,
          tax_level = tax_level,
          group_var = "captive_wild_binary",
          p_adj_method = p_adj_method,
          prv_cut = current_prv_cut,
          lib_cut = lib_cut,
          s0_perc = s0_perc
        )
        
        anc_results[[name]][[tax_level]] <- result
        
        # Extract results dataframe
        res_df <- result$res
        
        # Create plot
        plot_name <- paste0("ancombc_", tax_level, "_", name)
        plot_result <- create_ancombc_plot(
          res_df = res_df,
          dataset_name = name,
          tax_level = tax_level,
          output_name = plot_name,
          save_plot = save_results,
          colorlist = colorlist
        )
        
        if (!is.null(plot_result)) {
          anc_plots[[name]][[tax_level]] <- plot_result$plot
          cat("    Found", nrow(plot_result$plot_df), "significant taxa\n")
        } else {
          cat("    No significant taxa found\n")
        }
        
      }, error = function(e) {
        warning(paste("Error processing", name, "at", tax_level, ":", e$message))
      })
    }
  }
  
  # Trim results (keep only res component to save space)
  anc_results_trimmed <- list()
  for (dataset_name in names(anc_results)) {
    anc_results_trimmed[[dataset_name]] <- list()
    for (tax_level in names(anc_results[[dataset_name]])) {
      res_obj <- anc_results[[dataset_name]][[tax_level]]
      anc_results_trimmed[[dataset_name]][[tax_level]] <- list(
        res = res_obj$res
      )
    }
  }
  
  # Save results if requested
  if (save_results) {
    save(anc_results_trimmed,
         file = "05_taxa_driving_groups/ANCOMBC/binary/ancombc_binary_results.RData")
    cat("\nResults saved to: 05_taxa_driving_groups/ANCOMBC/binary/ancombc_binary_results.RData\n")
  }
  
  return(list(
    results = anc_results_trimmed,
    plots = anc_plots,
    filtered_phyloseq = ps_list_filtered
  ))
}

#### Example Usage ####

# Example 1: Load phyloseq objects and process
# Uncomment and modify as needed:

# ps_list <- list(
#   phyloseq16m = readRDS("01.1_qc_checks/phyloseq_otu16m_mat_filt_20_thr10.rds"),
#   phyloseq16n = readRDS("01.1_qc_checks/phyloseq_otu16n_mat_filt_20_thr10.rds"),
#   phyloseq18m = readRDS("01.1_qc_checks/phyloseq_otu18m_nohost_mat_filt_20_thr10.rds"),
#   phyloseq18n = readRDS("01.1_qc_checks/phyloseq_otu18n_nohost_mat_filt_20_thr10.rds")
# )
# 
# results <- process_ancombc_binary(
#   ps_list = ps_list,
#   tax_levels = c("Family", "Genus", "Species"),
#   save_results = TRUE
# )

cat("\n=== Script ready! ===\n")
cat("To use, create a named list of phyloseq objects and call:\n")
cat("  results <- process_ancombc_binary(ps_list)\n")
cat("\nThis will:\n")
cat("  1. Filter out 'Wild, seized from traffickers' samples\n")
cat("  2. Create binary comparison (Captive vs Wild, free ranging)\n")
cat("  3. Run ANCOMBC2 analysis at Family, Genus, and Species levels\n")
cat("  4. Generate plots and save results\n")
