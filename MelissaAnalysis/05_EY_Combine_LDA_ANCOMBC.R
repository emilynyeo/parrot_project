# COmbine LDA and ANCOMBC 

#### Load libraries ####
pacman::p_load(tidyverse, vegan, MASS, phyloseq, tibble, ANCOMBC, ggplot2, coin,
               MicrobiotaProcess, patchwork, reshape2, ggnewscale, VennDiagram,
               UpSetR, gridExtra, grid, WGCNA)
source("utils.R")
setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")

# Read in ANCOMBC and LDA results 
load("05_taxa_driving_groups/LDA/lda_16m_16n_18m_18n.RData")
load("05_taxa_driving_groups/ancom_16m_16n_18m_18n.RData")

### Format LDA PLOTS 
diffbox_plots <- list()
LDA_plot_dfs <- list()
for (taxrank in names(lda_diff_results_2)) { # (Species, Genus, Family)
  for (dataset in names(lda_diff_results_2[[taxrank]])) {
    deres <- lda_diff_results_2[[taxrank]][[dataset]]
    
    # Select top features by LDAmean
    top_feat <- deres@result %>%
      arrange(desc(LDAmean)) %>%
      pull(f)
    
    # Create a new object containing only the top 30 features
    deres_top <- deres
    deres_top@result <- deres@result %>% filter(f %in% top_feat)
    deres_top@result$label <- gsub("__", ": ", deres_top@result$f)
    LDA_plot_dfs[[taxrank]][[dataset]] <- deres_top
    
    # Create LDA plot 
    #diffbox <- ggdiffbox(
    #  obj = deres_top,
    #  box_notch = FALSE,
    #  lineheight = 0.1,
    #  linewidth = 0.3,
    #  colorlist = c("#00AED7", "#FD9347", "seagreen", "purple4"),
    #  l_xlabtext = "Relative Abundance")
    
    # Optional: relabel facets or axes after plot is created
    #diffbox <- diffbox + scale_y_discrete(labels = deres_top@result$label)
    
    # Save or store the plot
    #plot_name <- paste0("diffbox_", taxrank, "_", dataset)
    #diffbox_plots[[taxrank]][[dataset]] <- diffbox
  }
}

### Format ANCOMBC plots 


# Trim ANCOM results after the fact
anc_results_trimmed <- list()
for (dataset_name in names(anc_results)) {
  anc_results_trimmed[[dataset_name]] <- list()
  
  for (tax_level in names(anc_results[[dataset_name]])) {
    res_obj <- anc_results[[dataset_name]][[tax_level]]
    
    anc_results_trimmed[[dataset_name]][[tax_level]] <- list(
      res = res_obj$res
      #res_global = res_obj$res_global,
      #res_dunn = res_obj$res_dunn
    )
  }
}


# Loop through datasets and taxonomic levels
anc_plot_dfs <- list()
for (dataset_name in names(anc_results_trimmed)) {
  for (tax_level in names(anc_results_trimmed[[dataset_name]])) {
    
    cat("Fetching ancom results for:", dataset_name, "at", tax_level, "\n")
    res_obj <- anc_results_trimmed[[dataset_name]][[tax_level]]
    
    # Check if result object contains `res` element
    if (!"res" %in% names(res_obj)) {
      warning(paste("No 'res' in", dataset_name, tax_level))
      next
    }
    res_df <- res_obj$res
    res_df <- res_df %>% 
      dplyr::filter(`q_(Intercept)` < 0.05)
    cat("res_df size:", format(object.size(res_df), units = "MB"), "\n")
    
    # Check for valid results
    if (nrow(res_df) == 0) next
    
    # Reshape the lfc, q, and diff values separately
    lfc_df <- res_df %>%
      select(taxon, starts_with("lfc_")) %>%
      pivot_longer(cols = -taxon, names_to = "group", values_to = "lfc") %>%
      mutate(group = sub("lfc_", "", group))
    
    q_df <- res_df %>%
      select(taxon, starts_with("q_")) %>%
      pivot_longer(cols = -taxon, names_to = "group", values_to = "q") %>%
      mutate(group = sub("q_", "", group))
    
    diff_df <- res_df %>%
      select(taxon, starts_with("diff_")) %>%
      pivot_longer(cols = -taxon, names_to = "group", values_to = "diff") %>%
      mutate(group = sub("diff_", "", group))
    
    # Join all into a single long-form dataframe
    plot_df <- reduce(list(lfc_df, q_df, diff_df), left_join, by = c("taxon", "group"))
    
    # Filter significant results
    plot_df <- plot_df %>% filter(diff == TRUE, q < 0.05)
    
    if (nrow(plot_df) == 0) next  # Skip if no significant taxa
    
    # Clean and format for plotting
    plot_df$group <- str_replace_all(plot_df$group, "captive_wild", "")
    plot_df$taxon <- gsub("_", " ", plot_df$taxon)
    plot_df <- plot_df %>%
      mutate(group = ifelse(group == "(Intercept)", "Captive", group)) 
    plot_df$group <- sub("^Captive\\.", "", plot_df$group)
    plot_df$group <- sub("WildWild", "Wild", plot_df$group)
    
    cat("plot_df size:", format(object.size(plot_df), units = "MB"), "\n")
    anc_plot_dfs[[dataset_name]][[tax_level]] <- plot_df
    
    # Create ANCOMBC plot
    p <- ggplot(plot_df, aes(x = reorder(taxon, lfc), y = lfc, fill = group)) +
      geom_col(position = position_dodge(width = 0.9)) +
      coord_flip() +
      labs(title = paste("Significant DA (ANCOMBC2):", dataset_name, "-", tax_level),
           x = "Taxon", y = "Log Fold Change", fill = "Group") +
      scale_fill_manual(
        values = c("Captive" = "#00AED7", "Wild" = "seagreen",
                   "Seized" = "purple4")) +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14))
  }
}


### Combine the two plot types 

# From ANCOMBC
ancom <- anc_plot_dfs[["phyloseq16n"]]$genus
ancom$taxon <- as.character(ancom$taxon)

# From LDA
lda <- LDA_plot_dfs[["Genus"]]$phyloseq16n@result
lda$taxon <- lda$label  # Or extract from `label` with gsub if needed
lda$taxon <- gsub("^[a-z]:\\s*", "", lda$label)

# ANCOMBC
ancom_df <- ancom %>%
  select(taxon, group, value = lfc) %>%
  mutate(stat = "LFC")

# LDA
lda_df <- lda %>%
  select(taxon, group = captive_wild, value = LDAmean) %>%
  mutate(stat = "LDA")
  
lda_df <- lda_df %>%
  mutate(group = recode(group,
                        "Wild, seized from traffickers" = "Seized",
                        "Wild, free ranging" = "Wild"))

combined_df <- bind_rows(ancom_df, lda_df)

combined_df <- combined_df %>%
  dplyr::group_by(taxon) %>%
  filter(n() > 1) %>%  # Appears in both LFC and LDA
  dplyr::ungroup()


library(ggplot2)

ggplot(combined_df, aes(x = reorder(taxon, value), y = value, fill = group)) +
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), aes(alpha = stat)) +
  coord_flip() +
  facet_wrap(~stat, scales = "free_x") +  # Separate panels for LFC and LDA
  scale_fill_manual(values = c("Captive" = "#00AED7", "Wild" = "seagreen", "Seized" = "purple4")) +
  scale_alpha_manual(values = c("LFC" = 1, "LDA" = 0.6)) +
  labs(title = "Differential Abundance: ANCOMBC LFC and LDA Scores",
       x = "Taxon", y = "Score (LFC or LDA)", fill = "Group", alpha = "Stat Type") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


length(unique(ancom_df$taxon %in% lda_df$taxon))
length(unique(lda_df$taxon %in% ancom_df$taxon))



lda_genus_16n <- lda_diff_results_2[["Genus"]]$phyloseq16n

# Get the current column names
original_names <- colnames(lda_genus_16n@originalD)
cleaned_names <- sub("^[a-z]__+", "", original_names)
colnames(lda_genus_16n@originalD) <- cleaned_names


lda_genus_16n@result$taxon <- gsub("^[a-z]:\\s*", "", lda_genus_16n@result$f)
lda_genus_16n@result$taxon <- sub("^[a-z]__+", "", lda_genus_16n@result$f)
lda_genus_16n@result$f <- gsub("^[a-z]:\\s*", "", lda_genus_16n@result$f)
lda_genus_16n@result$f <- sub("^[a-z]__+", "", lda_genus_16n@result$f)

#lda_genus_16n@result$Group <- lda_genus_16n@result$Captive.Wild %>% 
#                              recode("Wild, seized from traffickers" = "Seized",
#                                     "Wild, free ranging" = "Wild")
#lda_genus_16n@result$Captive.Wild <- lda_genus_16n@result$Captive.Wild %>% 
#  recode("Wild, seized from traffickers" = "Seized",
#         "Wild, free ranging" = "Wild")



### Trying to plot seperately 
lda_genus_16n_feature_list <- unique(lda_genus_16n@result$f) 
anc_genus_16n_feature_list <- unique(anc_plot_dfs[["phyloseq16n"]]$genus$taxon)
colorlist <- c("Captive" = "#00AED7", "Wild" = "seagreen", "Seized" = "purple4")

#### Taxa overlap
overlap <- lda_genus_16n_feature_list[lda_genus_16n_feature_list %in% anc_genus_16n_feature_list]
unique_to_ancombc <- setdiff(anc_genus_16n_feature_list, overlap)
unique_to_lda <- setdiff(lda_genus_16n_feature_list, overlap)

ordered_features_ancombc <- c(overlap, unique_to_ancombc)
ordered_features_all <- c(overlap, unique_to_ancombc, unique_to_lda)


# Plot p1 
p1 <- make_p1_boxplot(
  obj = lda_genus_16n,
  featurelist = ordered_features_all,
  colorlist = colorlist,
  box_notch = FALSE,
  xlabtext = "Relative Abundance"
)
p1

p1_ord <- make_p1_boxplot_ordered(
  obj = lda_genus_16n,
  featurelist = ordered_features_all,
  colorlist = colorlist,
  box_notch = FALSE,
  xlabtext = "Relative Abundance"
)

p1_ord

# Plot effect size (p2)
# Add missing rows to nodedfres_all (for p2):

# Fill missing features for effect size plot with NA
nodedfres_all <- lda_genus_16n@result
missing_in_p2 <- setdiff(ordered_features_all, nodedfres_all$f)
if (length(missing_in_p2) > 0) {
  na_rows <- data.frame(f = missing_in_p2)
  # Fill other columns with NA or suitable values
  na_rows[, setdiff(colnames(nodedfres_all), "f")] <- NA
  nodedfres_all <- rbind(nodedfres_all, na_rows)
}
# Reorder levels to match p1
nodedfres_all$f <- factor(nodedfres_all$f, levels = rev(ordered_features_all))

p2 <- make_p2_effectsize(obj = lda_genus_16n, 
                         nodedfres = nodedfres_all, 
                         colorlist = colorlist, 
                         xlim_range = c(2.5, 5.5))

p2 <- p2 + theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank())
p2
# Plot p3 ancombc 

# Add missing rows to anc_plot_dfs[["phyloseq16n"]]$genus (for p3):

anc_df <- anc_plot_dfs[["phyloseq16n"]]$genus
missing_in_p3 <- setdiff(ordered_features_all, anc_df$taxon)
if (length(missing_in_p3) > 0) {
  na_rows <- data.frame(taxon = missing_in_p3)
  na_rows[, setdiff(colnames(anc_df), "taxon")] <- NA
  anc_df <- rbind(anc_df, na_rows)
}

anc_df$taxon <- factor(anc_df$taxon, levels = rev(ordered_features_all))

p3 <- make_p3_ancombc(ancombc_df = anc_df, 
                      dataset_name = "16N",
                      tax_level = "Genus", 
                      featurelist = ordered_features_all,
                      colorlist = colorlist)

p3 <- p3 + theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 plot.title = element_blank())

p3

p1 | p2 | p3 

