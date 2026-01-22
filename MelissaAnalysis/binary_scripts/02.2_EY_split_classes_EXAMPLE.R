# Example usage of 02.2_EY_split_classes.R
# This shows how to use the flexible splitting script

source("binary_scripts/02.2_EY_split_classes.R")

#### Load your phyloseq objects ####
# Files are saved with save(), so we need to use load() not readRDS()
# Load and rename directly
load("01.1_qc_checks/phyloseq_otu16m_mat_filt_20_thr10.rds")
load("01.1_qc_checks/phyloseq_otu16n_mat_filt_20_thr10.rds")
load("01.1_qc_checks/phyloseq_otu18m_nohost_mat_filt_20_thr10.rds")
load("01.1_qc_checks/phyloseq_otu18n_nohost_mat_filt_20_thr10.rds")

# Create list with renamed objects
ps_list <- list(
  otu16m_mat_filt_20_thr10 = phyloseq_otu16m_mat_filt_20_thr10,
  otu16n_mat_filt_20_thr10 = phyloseq_otu16n_mat_filt_20_thr10,
  otu18m_nohost_mat_filt_20_thr10 = phyloseq_otu18m_nohost_mat_filt_20_thr10,
  otu18n_nohost_mat_filt_20_thr10 = phyloseq_otu18n_nohost_mat_filt_20_thr10
)

# Clean up - remove the original objects if desired
rm(phyloseq_otu16m_mat_filt_20_thr10, 
   phyloseq_otu16n_mat_filt_20_thr10,
   phyloseq_otu18m_nohost_mat_filt_20_thr10,
   phyloseq_otu18n_nohost_mat_filt_20_thr10)

#### Load categorization data ####
catdat_18s <- read.csv("01_process_and_clean_data/simple18_labelled2.csv")
listPlantSpeciesFILLED <- read.csv("02_split_18_classes/listPlantSpeciesFILLED.txt") %>%
  as.data.frame()

#### Process all phyloseq objects ####
results <- process_phyloseq_list(
  ps_list = ps_list,
  catdat_18s = catdat_18s,
  listPlantSpeciesFILLED = listPlantSpeciesFILLED
)

#### Optional: Merge datasets before processing ####
# Example: Merge Renee's nanopore data with original nanopore
# merge_list <- list(
#   otu18nano = list(
#     load("01.1_qc_checks/phyloseq_otu18n_nohost_mat_filt_20_thr10.rds"),
#     load("01.1_qc_checks/phyloseq_otu18ren_nohost_mat_filt_20_thr10.rds")
#   )
# )
# 
# ps_list_merged <- list(
#   otu18miseq = load("01.1_qc_checks/phyloseq_otu18m_nohost_mat_filt_20_thr10.rds")
# )
# 
# results_merged <- process_phyloseq_list(
#   ps_list = ps_list_merged,
#   catdat_18s = catdat_18s,
#   listPlantSpeciesFILLED = listPlantSpeciesFILLED,
#   merge_list = merge_list
# )

cat("\n=== Processing complete! ===\n")
cat("Output files saved to: 02_split_18_classes/\n")
cat("Access results via: results$[name]_[subset]\n")
cat("Example: results$otu18m_nohost_mat_filt_20_thr10_plants\n")
