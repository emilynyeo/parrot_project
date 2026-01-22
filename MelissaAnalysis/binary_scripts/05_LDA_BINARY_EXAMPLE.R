# Example usage of 05_LDA_BINARY.R
# This shows how to use the binary LDA analysis script

source("binary_scripts/05_LDA_BINARY.R")

#### Load your phyloseq objects ####
# Files are saved with save(), so we need to use load() not readRDS()
# Load and rename directly
load("01.1_qc_checks/phyloseq_otu16m_mat_filt_20_thr10.rds")
load("01.1_qc_checks/phyloseq_otu16n_mat_filt_20_thr10.rds")
load("01.1_qc_checks/phyloseq_otu18m_nohost_mat_filt_20_thr10.rds")
load("01.1_qc_checks/phyloseq_otu18n_nohost_mat_filt_20_thr10.rds")

# Create list with renamed objects
ps_list <- list(
  phyloseq16m = phyloseq_otu16m_mat_filt_20_thr10,
  phyloseq16n = phyloseq_otu16n_mat_filt_20_thr10,
  phyloseq18m = phyloseq_otu18m_nohost_mat_filt_20_thr10,
  phyloseq18n = phyloseq_otu18n_nohost_mat_filt_20_thr10
)

# Clean up - remove the original objects if desired
rm(phyloseq_otu16m_mat_filt_20_thr10, 
   phyloseq_otu16n_mat_filt_20_thr10,
   phyloseq_otu18m_nohost_mat_filt_20_thr10,
   phyloseq_otu18n_nohost_mat_filt_20_thr10)

#### Run binary LDA analysis ####
results <- process_lda_binary(
  ps_list = ps_list,
  tax_levels = c("Species", "Genus", "Family"),
  top_n_plots = 30,
  colorlist = c("#00AED7", "#FD9347"),  # Colors for Captive and Wild
  save_results = TRUE
)

#### Access results ####
# LDA results
# results$results$Genus$phyloseq16m  # Genus-level results for 16S Miseq
# results$results$Species$phyloseq18n  # Species-level results for 18S Nanopore

# Plots
# results$plots$Genus$phyloseq16m  # Plot for Genus-level 16S Miseq

# Filtered phyloseq objects (with binary groups)
# results$filtered_phyloseq$phyloseq16m

#### Customize analysis ####
# You can customize various parameters:
# results_custom <- process_lda_binary(
#   ps_list = ps_list,
#   tax_levels = c("Genus", "Family"),  # Only analyze these levels
#   firstalpha = 0.01,  # Stricter first filter
#   secondalpha = 0.0001,  # Stricter second filter
#   ldascore = 3,  # Lower LDA score threshold
#   top_n_plots = 20,  # Show top 20 taxa in plots
#   colorlist = c("blue", "orange"),  # Custom colors
#   seed = 1234,  # Different random seed
#   save_results = TRUE
# )

cat("\n=== Analysis complete! ===\n")
cat("Output files saved to: 05_taxa_driving_groups/LDA/binary/\n")
cat("Access results via: results$results$[tax_level]$[dataset_name]\n")
cat("Example: results$results$Genus$phyloseq16m\n")
