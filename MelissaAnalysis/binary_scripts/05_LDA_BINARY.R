# LDA Differential Analysis - Binary Comparison (Captive vs Wild, free ranging)
# Filters out "Wild, seized from traffickers" for binary comparison
# Accepts a list of phyloseq objects

#### Load libraries ####
pacman::p_load(tidyverse, vegan, MASS, phyloseq, tibble, ANCOMBC, ggplot2, coin,
               MicrobiotaProcess, patchwork, reshape2, ggnewscale, VennDiagram,
               UpSetR, gridExtra, grid, WGCNA)
setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")

# Create output directory
if (!dir.exists("05_taxa_driving_groups/LDA/binary")) {
  dir.create("05_taxa_driving_groups/LDA/binary", recursive = TRUE)
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

# Run LDA analysis on a single phyloseq object
run_lda_analysis <- function(ps, 
                             taxrank = "Genus",
                             group_var = "captive_wild_binary",
                             firstalpha = 0.05,
                             secondalpha = 0.001,
                             ldascore = 4,
                             seed = 1024) {
  set.seed(seed)
  
  # Agglomerate to taxonomic level
  ps_agglom <- tax_glom(ps, taxrank = taxrank)
  
  # Run diff_analysis
  result <- diff_analysis(
    obj = ps_agglom,
    classgroup = group_var,
    mlfun = "lda",
    filtermod = "pvalue",
    firstcomfun = "kruskal_test",
    firstalpha = firstalpha,
    strictmod = TRUE,
    secondcomfun = "wilcox_test",
    subclmin = 3,
    subclwilc = TRUE,
    secondalpha = secondalpha,
    ldascore = ldascore
  )
  
  return(result)
}

# Create LDA plot
create_lda_plot <- function(deres,
                            top_n = 30,
                            colorlist = c("#00AED7", "#FD9347"),
                            output_name = NULL,
                            save_plot = TRUE) {
  # Select top features by LDAmean
  top_feat <- deres@result %>%
    arrange(desc(LDAmean)) %>%
    slice_head(n = top_n) %>%
    pull(f)
  
  # Create a new object containing only the top features
  deres_top <- deres
  deres_top@result <- deres@result %>%
    filter(f %in% top_feat)
  deres_top@result$label <- gsub("__", ": ", deres_top@result$f)
  
  # Create plot
  diffbox <- ggdiffbox(
    obj = deres_top,
    box_notch = FALSE,
    lineheight = 0.1,
    linewidth = 0.3,
    colorlist = colorlist,
    l_xlabtext = "Relative Abundance"
  )
  
  # Relabel facets
  diffbox <- diffbox + scale_y_discrete(labels = deres_top@result$label)
  
  # Save plot if requested
  if (save_plot && !is.null(output_name)) {
    ggsave(
      filename = paste0("05_taxa_driving_groups/LDA/binary/", output_name, ".png"),
      plot = diffbox,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  return(list(plot = diffbox, result = deres_top))
}

#### Main Processing Function ####

# Process a list of phyloseq objects
process_lda_binary <- function(ps_list,
                                tax_levels = c("Species", "Genus", "Family"),
                                group_var = "Captive.Wild",
                                firstalpha = 0.05,
                                secondalpha = 0.001,
                                ldascore = 4,
                                top_n_plots = 30,
                                colorlist = c("#00AED7", "#FD9347"),
                                seed = 1024,
                                save_results = TRUE) {
  
  # Initialize nested list to store results
  lda_diff_results <- list()
  diffbox_plots <- list()
  
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
  
  # Run LDA analysis for each taxonomic level and dataset
  set.seed(seed)
  for (taxrank in tax_levels) {
    cat("\n=== Processing", taxrank, "level ===\n")
    lda_diff_results[[taxrank]] <- list()
    diffbox_plots[[taxrank]] <- list()
    
    for (name in names(ps_list_filtered)) {
      cat("  Processing", name, "...\n")
      ps_filtered <- ps_list_filtered[[name]]
      
      tryCatch({
        # Run LDA analysis
        result <- run_lda_analysis(
          ps = ps_filtered,
          taxrank = taxrank,
          group_var = "captive_wild_binary",
          firstalpha = firstalpha,
          secondalpha = secondalpha,
          ldascore = ldascore,
          seed = seed
        )
        
        lda_diff_results[[taxrank]][[name]] <- result
        
        # Create plot
        plot_name <- paste0("diffbox_", taxrank, "_", name)
        plot_result <- create_lda_plot(
          deres = result,
          top_n = top_n_plots,
          colorlist = colorlist,
          output_name = plot_name,
          save_plot = save_results
        )
        
        diffbox_plots[[taxrank]][[name]] <- plot_result$plot
        
        cat("    Found", nrow(result@result), "significant taxa\n")
        
      }, error = function(e) {
        warning(paste("Error processing", name, "at", taxrank, ":", e$message))
      })
    }
  }
  
  # Save results if requested
  if (save_results) {
    save(lda_diff_results,
         file = "05_taxa_driving_groups/LDA/binary/lda_binary_results.RData")
    cat("\nResults saved to: 05_taxa_driving_groups/LDA/binary/lda_binary_results.RData\n")
  }
  
  return(list(
    results = lda_diff_results,
    plots = diffbox_plots,
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
# results <- process_lda_binary(
#   ps_list = ps_list,
#   tax_levels = c("Species", "Genus", "Family"),
#   top_n_plots = 30,
#   save_results = TRUE
# )

cat("\n=== Script ready! ===\n")
cat("To use, create a named list of phyloseq objects and call:\n")
cat("  results <- process_lda_binary(ps_list)\n")
cat("\nThis will:\n")
cat("  1. Filter out 'Wild, seized from traffickers' samples\n")
cat("  2. Create binary comparison (Captive vs Wild, free ranging)\n")
cat("  3. Run LDA analysis at Family, Genus, and Species levels\n")
cat("  4. Generate plots and save results\n")
