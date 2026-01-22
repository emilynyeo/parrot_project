# Example usage of 06_ANCOMBC_BINARY.R
# This shows how to use the binary ANCOMBC analysis script
# Can be run from binary_scripts/ or MelissaAnalysis/ directory

# Set working directory to MelissaAnalysis if running from binary_scripts/
if (basename(getwd()) == "binary_scripts") {
  setwd("..")
}

# Source the main script
source("binary_scripts/06_ANCOMBC_BINARY.R")

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

#### Run binary ANCOMBC analysis ####
results <- process_ancombc_binary(
  ps_list = ps_list,
  tax_levels = c("Family", "Genus", "Species"),
  colorlist = c("#00AED7", "#FD9347"),  # Colors for Captive and Wild
  save_results = TRUE
)

#### Access results ####
# ANCOMBC results
# results$results$phyloseq16m$Genus  # Genus-level results for 16S Miseq
# results$results$phyloseq18n$Species  # Species-level results for 18S Nanopore

# Plots
# results$plots$phyloseq16m$Genus  # Plot for Genus-level 16S Miseq

# Filtered phyloseq objects (with binary groups)
# results$filtered_phyloseq$phyloseq16m

#### Customize analysis ####
# You can customize various parameters:
# results_custom <- process_ancombc_binary(
#   ps_list = ps_list,
#   tax_levels = c("Genus", "Family"),  # Only analyze these levels
#   p_adj_method = "fdr",  # Use FDR instead of Holm
#   prv_cut = 0.05,  # Stricter prevalence cutoff
#   lib_cut = 100,  # Higher library size cutoff
#   prv_cut_species = 0.2,  # Different cutoff for species level
#   colorlist = c("blue", "orange"),  # Custom colors
#   save_results = TRUE
# )

cat("\n=== Analysis complete! ===\n")
cat("Output files saved to: 05_taxa_driving_groups/ANCOMBC/binary/\n")
cat("Access results via: results$results$[dataset_name]$[tax_level]\n")
cat("Example: results$results$phyloseq16m$Genus\n")
