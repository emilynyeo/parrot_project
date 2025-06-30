# Which Taxa drive the differences between groups? 

#### Load libraries ####
pacman::p_load(tidyverse, vegan, MASS, phyloseq, tibble, ANCOMBC)
               
#dir.create("04_beta_diversity")
#### Load data ####
# Also converting to matrices with rownames for phyloseq
load("01_process_and_clean_data/phyloseq16n.rds")
load("01_process_and_clean_data/phyloseq16ren.rds")

load("01_process_and_clean_data/phyloseq18n.rds")
load("01_process_and_clean_data/phyloseq18ren_nohost.rds")
load("01_process_and_clean_data/phyloseq18m.rds")

# Split into classes phylo
load("02_split_18_classes/phyloseq18nano_nohost.rds")
load("02_split_18_classes/phyloseq18nano_plants.rds")
load("02_split_18_classes/phyloseq18nano_microeuk.rds")

load("02_split_18_classes/phyloseq18miseq_nohost.rds")
load("02_split_18_classes/phyloseq18miseq_plants.rds")
load("02_split_18_classes/phyloseq18miseq_microeuk.rds")

# Merge phyloseq objects 
phyloseq16all <- merge_phyloseq(phyloseq16n, phyloseq16ren)
phyloseq18all <- phyloseq18nano_nohost
phyloseq18all <- prune_samples(sample_sums(phyloseq18all) >0, phyloseq18all)
phyloseq18miseq_nohost <- prune_samples(sample_sums(phyloseq18miseq_nohost) >0, phyloseq18miseq_nohost)

# Combine Crate Data into one
sample_df <- as.data.frame(sample_data(phyloseq16all))

# Create new variable with recoded values
sample_df$cap_wild_crate <- as.character(sample_df$Captive.Wild)
sample_df$cap_wild_crate[sample_df$Captive.Wild %in% c("Crate 1: Unknown origin", 
                                                       "Crate 2: Unknown origin")] <- "New crates"
# factor with levels
sample_df$cap_wild_crate <- factor(sample_df$cap_wild_crate,
                                   levels = c("Captive", 
                                              "Wild, free ranging", 
                                              "Wild, seized from traffickers", 
                                              "New crates"))
# Replace with updated data
sample_data(phyloseq16all) <- sample_data(sample_df)
levels(sample_data(phyloseq16all)$cap_wild_crate)

#### ANCOM-BC2 for all 16s nano data ####
?ancombc2

set.seed(123)
p16n_all_anc <-  ancombc2(data = phyloseq16all, tax_level = "Genus",
                  fix_formula = "Captive.Wild", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0, lib_cut = 1000, s0_perc = 0.05,
                  group = "Captive.Wild", struc_zero = TRUE, neg_lb = TRUE)

p16n_all_anc_1crate <-  ancombc2(data = phyloseq16all, tax_level = "Genus",
                          fix_formula = "cap_wild_crate", rand_formula = NULL,
                          p_adj_method = "holm", pseudo_sens = TRUE,
                          prv_cut = 0, lib_cut = 1000, s0_perc = 0.05,
                          group = "cap_wild_crate", struc_zero = TRUE, neg_lb = TRUE)
                  
res <- p16n_all_anc_1crate$res

# https://www.yanh.org/2021/01/01/microbiome-r/#differential-abundance-analysis
# select the bottom 20 with lowest p values
res.or_p <- rownames(res$q_val["body.sitetongue"])[base::order(res$q_val[,"body.sitetongue"])]
taxa_sig <- res.or_p[1:20]
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)

# Only keep gut and tongue samples 
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.taxa.sub)), ps.taxa.rel.sig)

matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Species"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))

# Define the annotation color for columns and rows
annotation_col = data.frame(
  Subject = as.factor(metadata_sub$subject), 
  `Body site` = as.factor(metadata_sub$body.site), 
  check.names = FALSE
)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(
  Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"])
)
rownames(annotation_row) = rownames(matrix)

# ann_color should be named vectors
phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)
ann_colors = list(
  Subject = c(`subject-1` = "red", `subject-2` = "blue"),
  `Body site` = c(`gut` = "purple", tongue = "yellow"),
  Phylum = phylum_col
)
ComplexHeatmap::pheatmap(matrix, 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row, 
                         annotation_colors = ann_colors)
# Plot 
plot_ancombc_pq(
  phyloseq16all,
  p16n_all_anc,
  filter_passed = TRUE,
  filter_diff = TRUE,
  min_abs_lfc = 0,
  tax_col = "Genus",
  tax_label = "Species",
  add_marginal_vioplot = TRUE,
  add_label = TRUE,
  add_hline_cut_lfc = NULL
)

save(p16n_all_anc, file = "05_taxa_driving_groups/ancom_16sn.RData")

# MICROBIOME WORKSHOP
### https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html

library(coin)
library(MicrobiotaProcess)
library(patchwork)