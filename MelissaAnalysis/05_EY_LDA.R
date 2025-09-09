### Use LDA diff taxa analysis 


#### Load libraries ####
pacman::p_load(tidyverse, vegan, MASS, phyloseq, tibble, ANCOMBC, ggplot2, coin,
               MicrobiotaProcess, patchwork, reshape2, ggnewscale, VennDiagram,
               UpSetR, gridExtra, grid, WGCNA)
setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")
#BiocManager::install(c("GO.db", "impute", "preprocessCore"))

#### Load data ####
# Also converting to matrices with rownames for phyloseq
# load("01_process_and_clean_data/phyloseq16m.rds")
# load("01_process_and_clean_data/phyloseq16n.rds")
load("01.1_qc_checks/phyloseq_otu16m_mat_filt_20_thr10.rds")
load("01.1_qc_checks/phyloseq_otu16n_mat_filt_20_thr10.rds")

sample_df <- as.data.frame(sample_data(phyloseq_otu16m_mat_filt_20_thr10))
head(sample_df$Captive.Wild)

# Split into classes phylo
#load("02_split_18_classes/phyloseq18nano_nohost.rds")
#load("02_split_18_classes/phyloseq18miseq_nohost.rds")
load("01.1_qc_checks/phyloseq_otu18m_nohost_mat_filt_20_thr10.rds")
load("01.1_qc_checks/phyloseq_otu18n_nohost_mat_filt_20_thr10.rds")

phyloseq16m <- prune_samples(sample_sums(phyloseq_otu16m_mat_filt_20_thr10) >0, phyloseq_otu16m_mat_filt_20_thr10)
phyloseq16n <- prune_samples(sample_sums(phyloseq_otu16n_mat_filt_20_thr10) >0, phyloseq_otu16n_mat_filt_20_thr10)
phyloseq18m <- prune_samples(sample_sums(phyloseq_otu18m_nohost_mat_filt_20_thr10) >0, phyloseq_otu18m_nohost_mat_filt_20_thr10)
phyloseq18n <- prune_samples(sample_sums(phyloseq_otu18n_nohost_mat_filt_20_thr10) >0, phyloseq_otu18n_nohost_mat_filt_20_thr10)

rm(phyloseq_otu16m_mat_filt_20_thr10, phyloseq_otu16n_mat_filt_20_thr10,
   phyloseq_otu18m_nohost_mat_filt_20_thr10, phyloseq_otu18n_nohost_mat_filt_20_thr10)

phyloseq_list <- list(
  phyloseq16m = phyloseq16m,
  phyloseq16n = phyloseq16n,
  phyloseq18m = phyloseq18m,
  phyloseq18n = phyloseq18n)

# rename wild cap groups 
for (i in seq_along(phyloseq_list)) {
  phy <- phyloseq_list[[i]]
  sd <- as.data.frame(sample_data(phy))
  sd$captive_wild <- factor(sd$Captive.Wild,
                            levels = c("Captive", "Wild, free ranging", "Wild, seized from traffickers"),
                            labels = c("Captive", "Wild", "Seized"))
  sample_data(phy) <- sample_data(sd)
  phyloseq_list[[i]] <- phy
}


# BIOMARKER DISCOVERY USING LDA

tax_levels <- c("Species", "Genus", "Family")

# Initialize nested list to store results
lda_diff_results_2 <- list()
set.seed(1024)
for (taxrank in tax_levels) {
  lda_diff_results_2[[taxrank]] <- list()
  for (name in names(phyloseq_list)) {
    ps_raw <- phyloseq_list[[name]]
    ps_agglom <- tax_glom(ps_raw, taxrank = taxrank)# Agglomerate to taxa level
    
    # Run diff_analysis
    result <- diff_analysis(obj = ps_agglom,
                            classgroup = "captive_wild",
                            mlfun = "lda",
                            filtermod = "pvalue",
                            firstcomfun = "kruskal_test",
                            firstalpha = 0.05,
                            strictmod = TRUE,
                            secondcomfun = "wilcox_test",
                            subclmin = 3,
                            subclwilc = TRUE,
                            secondalpha = 0.001,
                            ldascore = 3)
    lda_diff_results_2[[taxrank]][[name]] <- result
  }
}

save(lda_diff_results_2, file = "05_taxa_driving_groups/LDA/lda3_16m_16n_18m_18n.RData")
dim(lda_diff_results_2$Species$phyloseq16m@result)

# Initialize a list to store the plots (optional)
diffbox_plots <- list()
for (taxrank in names(lda_diff_results_2)) { # (Species, Genus, Family)
  for (dataset in names(lda_diff_results_2[[taxrank]])) {
    deres <- lda_diff_results_2[[taxrank]][[dataset]]
    
    # Select top 30 features by LDAmean
    top30_feat <- deres@result %>%
      arrange(desc(LDAmean)) %>%
      #slice_head(n = 30) %>%
      pull(f)
    
    # Create a new object containing only the top 30 features
    deres_top30 <- deres
    deres_top30@result <- deres@result %>%
      filter(f %in% top30_feat)
    #deres_top30@result$f <- gsub("__", ": ", deres_top30@result$f)
    deres_top30@result$label <- gsub("__", ": ", deres_top30@result$f)
    
    
    diffbox <- ggdiffbox(
      obj = deres_top30,
      box_notch = FALSE,
      lineheight = 0.1,
      linewidth = 0.3,
      colorlist = c("#00AED7", "#FD9347", "seagreen", "purple4"),
      l_xlabtext = "Relative Abundance"
    )
    
    # Optional: relabel facets or axes after plot is created
    diffbox <- diffbox + scale_y_discrete(labels = deres_top30@result$label)
    
    # Save or store the plot
    plot_name <- paste0("diffbox_", taxrank, "_", dataset)
    diffbox_plots[[taxrank]][[dataset]] <- diffbox
    
    # Optionally: Print or save plot
    print(diffbox)
    ggsave(paste0("05_taxa_driving_groups/LDA/", plot_name, ".png"), 
           diffbox, width = 8, height = 6, dpi = 300)
  }
}


## CORRELATION OF TAXONOMY
ps <- phyloseq16m
genustab <- get_taxadf(ps, taxlevel=6)
genustab <- data.frame(t(otu_table(genustab)), check.names=FALSE)
genustab <- data.frame(apply(genustab, 2, function(x)x/sum(x)), check.names=FALSE)

cortest <- WGCNA::corAndPvalue(genustab, method="spearman", alternative="two.sided")
cortest$cor[upper.tri(cortest$cor, diag = TRUE)] <- NA
cortest$p[upper.tri(cortest$p, diag = TRUE)] <- NA
cortab1 <- na.omit(melt(t(cortest$cor))) %>% rename(from=Var1,to=Var2,cor=value)
corptab1 <- na.omit(melt(t(cortest$p))) %>% rename(pvalue=value)
cortab1$fdr <- p.adjust(corptab1$pvalue, method="fdr")

cortab1 <- cortab1 %>% mutate(correlation=case_when(cor>0 ~ "positive",cor < 0 ~ "negative",TRUE ~ "No"))
cortab2 <- cortab1 %>% filter(fdr <= 0.05) %>% filter(cor <= -0.5 | cor >= 0.8)


p <- ggdiffclade(
  obj=deres,
  alpha=0.3,
  linewd=0.25,
  skpointsize=0.2,
  layout="inward_circular",
  taxlevel=7,
  cladetext=0,
  setColors=FALSE,
  xlim=16
) +
  scale_fill_manual(values=c("#00AED7", "#FD9347", "seagreen", "purple4"),
                    guide=guide_legend(keywidth=0.5,
                                       keyheight=0.5,
                                       order=3,
                                       override.aes=list(alpha=1))
  ) +
  scale_size_continuous(range=c(1, 3),
                        guide=guide_legend(keywidth=0.5,keyheight=0.5,order=4,
                                           override.aes=list(shape=21))) +
  scale_colour_manual(values=rep("white", 100),guide="none")

p2 <- p +
  new_scale_color() +
  new_scale("size") +
  geom_tiplab(size=1, hjust=1) +
  geom_taxalink(
    data=cortab2,
    mapping=aes(taxa1=from,
                taxa2=to,
                colour=correlation,
                size=abs(cor)),
    alpha=0.4,
    ncp=10,
    hratio=1,
    offset=1.2) +
  scale_size_continuous(range = c(0.2, 1),
                        guide=guide_legend(keywidth=1, keyheight=0.5,
                                           order=1, override.aes=list(alpha=1))) +
  scale_colour_manual(values=c("chocolate2", "#009E73", "seagreen", "purple4"),
                      guide=guide_legend(keywidth=0.5, keyheight=0.5,
                                         order=2, override.aes=list(alpha=1, size=1)))
p2
cortab2
