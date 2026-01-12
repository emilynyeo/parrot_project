# Which Taxa drive the differences between groups? 

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

anc_results <- list()

sample_df <- as.data.frame(sample_data(phyloseq16m))
head(sample_df$Captive.Wild)

for (name in names(phyloseq_list)) {
  cat("Processing:", name, "\n")
  phy <- phyloseq_list[[name]]
  print(nsamples(phy)) # Sample Count 
  
  # ANCOMBC Genus
  cat("Running family level ANCOMBC2 on: ", name, "\n")
  anc_family <-  ancombc2(data = phy, tax_level = "Family",
                         fix_formula = "Captive.Wild", rand_formula = NULL,
                         p_adj_method = "holm", pseudo_sens = TRUE,
                         prv_cut = 0.01, lib_cut = 50, s0_perc = 0.05,
                         group = "Captive.Wild", struc_zero = TRUE, neg_lb = TRUE)
  
  cat("Running genus level ANCOMBC2 on: ", name, "\n")
  anc_genus <-  ancombc2(data = phy, tax_level = "Genus",
                         fix_formula = "Captive", rand_formula = NULL,
                         p_adj_method = "holm", pseudo_sens = TRUE,
                         prv_cut = 0.01, lib_cut = 50, s0_perc = 0.05,
                         group = "Captive", struc_zero = TRUE, neg_lb = TRUE)
  
  cat("Running species level ANCOMBC2 on: ", name, "\n")
  anc_species <-  ancombc2(data = phy, tax_level = "Species",
                         fix_formula = "Captive", rand_formula = NULL,
                         p_adj_method = "holm", pseudo_sens = TRUE,
                         prv_cut = 0.1, lib_cut = 50, s0_perc = 0.05,
                         group = "Captive", struc_zero = TRUE, neg_lb = TRUE)
  
  anc_results[[name]] <- list(
    family = anc_family,
    genus = anc_genus,
    species = anc_species)
}

save(anc_results, 
     file = "05_taxa_driving_groups/ancom_16m_16n_18m_18n.RData")

load("05_taxa_driving_groups/ancom_16m_16n_18m_18n.RData")

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


# Optional: Save the trimmed version
save(anc_results_trimmed, file = "05_taxa_driving_groups/ancom_16m_16n_18m_trimmed.RData")
rm(anc_results)

# Loop through datasets and taxonomic levels
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
    plot_df$group <- str_replace_all(plot_df$group, "cap_wild_crate", "")
    plot_df$taxon <- gsub("_", " ", plot_df$taxon)
    plot_df <- plot_df %>%
      mutate(group = ifelse(group == "(Intercept)", "Captive", group)) 
    plot_df$group <- sub("^Captive\\.", "", plot_df$group)
    plot_df$group <- sub("WildWild", "Wild", plot_df$group)
    
    cat("plot_df size:", format(object.size(plot_df), units = "MB"), "\n")
    
    # Create plot
    p <- ggplot(plot_df, aes(x = reorder(taxon, lfc), y = lfc, fill = group)) +
      geom_col(position = position_dodge(width = 0.9)) +
      coord_flip() +
      labs(title = paste("Significant DA (ANCOMBC2):", dataset_name, "-", tax_level),
           x = "Taxon", y = "Log Fold Change", fill = "Group") +
      scale_fill_manual(
        values = c("Captive" = "#00AED7",
                   "Wild, free ranging" = "seagreen",
                   "Wild, seized from traffickers" = "purple4")) +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
    
    # Define output file name
    output_dir <- "05_taxa_driving_groups"
    plot_file <- file.path(output_dir, 
                           paste0(dataset_name, "_ancombc_", tolower(tax_level), ".png"))
    ggsave(filename = plot_file, plot = p, height = 15, width = 18)
  }
}



# MICROBIOME WORKSHOP
### https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html

# ALPHA D and READ DEPTH
set.seed(1024)
ps <- phyloseq16all
rareres <- get_rarecurve(obj=ps, chunks=400)

prare1 <- ggrarecurve(obj=rareres, factorNames="cap_wild_crate",
                      indexNames=c("Observe", "Chao1", "ACE")
) +
  scale_fill_manual(values=c("#00AED7", "#FD9347","seagreen","purple4"))+
  scale_color_manual(values=c("#00AED7", "#FD9347","seagreen","purple4"))+
  theme_bw()+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"))          

prare2 <- ggrarecurve(obj=rareres,
                      factorNames="cap_wild_crate",
                      shadow=FALSE,
                      indexNames=c("Observe", "Chao1", "ACE")) +
  scale_color_manual(values=c("#00AED7", "#FD9347","seagreen","purple4"))+
  coord_cartesian(ylim = c(0, 400)) +
  theme_bw()+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"))
ggsave(filename = "05_taxa_driving_groups/16s_alpha_reads.png",prare2, height=5, width=8)
prare1 / prare2


# OVERLAPPING TAXA

vennlist <- get_vennlist(obj=ps, factorNames="cap_wild_crate")
upsetda <- get_upset(obj=ps, factorNames="cap_wild_crate")

vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#00AED7", "#FD9347", "seagreen", "purple4"),
                      cat.col=c("#00AED7", "#FD9347", "seagreen", "purple4"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.25,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
grid::grid.draw(vennp)
ggsave(filename = "05_taxa_driving_groups/16s_taxa_overlap.png",vennp, height=5, width=8)

# PCOA's

# distmethod: "unifrac","wunifrac","manhattan","euclidean","canberra","bray","kulczynski" (vegdist, dist options)
pcoares <- get_pcoa(obj=ps, distmethod="bray", method="hellinger")

# Visualizing the result
pcoaplot1 <- ggordpoint(obj=pcoares, biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("cap_wild_crate"), ellipse=TRUE) +
  scale_color_manual(values=c("#00AED7", "#FD9347", "seagreen", "purple4")) +
  scale_fill_manual(values=c("#00AED7", "#FD9347", "seagreen", "purple4"))
ggsave(filename = "05_taxa_driving_groups/16s_pcoaplot1.png",pcoaplot1, height=5, width=8)

# first and third principal co-ordinates
pcoaplot2 <- ggordpoint(obj=pcoares, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("cap_wild_crate"), ellipse=TRUE) +
  scale_color_manual(values=c("#00AED7", "#FD9347", "seagreen", "purple4")) +
  scale_fill_manual(values=c("#00AED7", "#FD9347", "seagreen", "purple4"))
ggsave(filename = "05_taxa_driving_groups/16s_pcoaplot2.png",pcoaplot2, height=5, width=8)
pcoaplot1 | pcoaplot2

# PERMANOVA
distme <- get_dist(ps, distmethod ="bray", method="hellinger")
sampleda <- data.frame(sample_data(ps), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
sampleda$Group <- factor(sampleda$cap_wild_crate)
set.seed(1024)
adores <- adonis2(distme ~ cap_wild_crate, data=sampleda, permutation=9999)
data.frame(adores$aov.tab)
print(adores)

Table_grob <- tableGrob(adores, rows = NULL, theme = ttheme_default(
  core = list(fg_params = list(hjust = 0.5)),
  colhead = list(fg_params = list(fontface = 2, hjust = 0.5))
))

# Save to PNG
png("05_taxa_driving_groups/permanova_16all.png", width = 7, height = 5)
grid.draw(Table_grob)
dev.off()

# BIOMARKER DISCOVERY 

set.seed(1024)
deres <- diff_analysis(obj = ps, classgroup = "cap_wild_crate",
                       mlfun = "lda",
                       filtermod = "pvalue",
                       firstcomfun = "kruskal_test",
                       firstalpha = 0.05,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                       subclmin = 3,
                       subclwilc = TRUE,
                       secondalpha = 0.01,
                       lda=3)

deres_rf <- diff_analysis(obj = ps, classgroup = "cap_wild_crate",
                       mlfun = "rf",
                       filtermod = "pvalue",
                       firstcomfun = "kruskal_test",
                       firstalpha = 0.05,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                       subclmin = 3,
                       subclwilc = TRUE,
                       secondalpha = 0.01,
                       lda=3)
deres
deres_rf

# Enriched group:
result_df <- deres_rf@result
head(result_df[, c("f", "MDAmean", "cap_wild_crate")])

# Biomarker Abundance (r) & conf. interval of effect size (LDA or MDA) of biomarker. The bigger confident interval shows that the biomarker is more fluctuant, owing to the influence of samples number.

top30_feat <- deres@result %>%
  arrange(desc(LDAmean)) %>%
  slice_head(n = 30) %>% pull(f)
deres_top30 <- deres
deres_top30@result <- deres@result %>% filter(f %in% top30_feat)

diffbox <- ggdiffbox(obj=deres_top30, box_notch=FALSE, 
                     colorlist=c("#00AED7", "#FD9347", "seagreen", "purple4"), 
                     l_xlabtext="relative abundance")
diffbox
ggsave(filename = "05_taxa_driving_groups/16s_biomarker_abundance_LDA.png",diffbox, height=7, width=8)

# Confidence interval of effect

top40_feat <- deres@result %>%
  arrange(desc(LDAmean)) %>%
  slice_head(n = 40) %>% pull(f)
deres_top40 <- deres
deres_top40@result <- deres@result %>% filter(f %in% top40_feat)

es_p <- ggeffectsize(obj=deres_top40, 
                     lineheight=0.1,
                     linewidth=0.3) + 
  scale_color_manual(values=c("#00AED7", "#FD9347", "seagreen", "purple4")) 

es_p
ggsave(filename = "05_taxa_driving_groups/16s_biomarker_top_LDA.png",es_p, height=5, width=8)

## CORRELATION OF TAXONOMY

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
    offset=1.2
  ) +
  scale_size_continuous(range = c(0.2, 1),
                        guide=guide_legend(keywidth=1, keyheight=0.5,
                                           order=1, override.aes=list(alpha=1))
  ) +
  scale_colour_manual(values=c("chocolate2", "#009E73", "seagreen", "purple4"),
                      guide=guide_legend(keywidth=0.5, keyheight=0.5,
                                         order=2, override.aes=list(alpha=1, size=1)))
p2
cortab2
