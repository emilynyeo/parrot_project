# Which features drive the differences between groups for 18s? 

#### Load libraries ####
setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")
pacman::p_load(tidyverse, vegan, MASS, phyloseq, tibble, ANCOMBC, ggplot2, coin,
               MicrobiotaProcess, patchwork, reshape2, ggnewscale, VennDiagram,
               UpSetR, gridExtra, grid, WGCNA, gt)
#BiocManager::install(c("GO.db", "impute", "preprocessCore"))

#### Load data ####
#load("01_process_and_clean_data/phyloseq18n.rds")
# Split into classes phylo
load("02_split_18_classes/phyloseq18nano_nohost.rds")
#load("02_split_18_classes/phyloseq18nano_plants.rds")
#load("02_split_18_classes/phyloseq18nano_microeuk.rds")

head(phyloseq::tax_table(phyloseq18n))
dim(phyloseq::tax_table(phyloseq18n))
head(phyloseq::tax_table(phyloseq18nano_microeuk))
dim(phyloseq::tax_table(phyloseq18nano_microeuk))
head(phyloseq::tax_table(phyloseq18nano_plants))
dim(phyloseq::tax_table(phyloseq18nano_plants))
head(phyloseq::tax_table(phyloseq18nano_nohost))
dim(phyloseq::tax_table(phyloseq18nano_nohost))


# ** Merge phyloseq objects **
phyloseq18all <- phyloseq18nano_nohost
phyloseq18all <- prune_samples(sample_sums(phyloseq18all) >0, phyloseq18all)

# ** Combine Crate Data into one **
sample_df <- as.data.frame(sample_data(phyloseq18all))

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
sample_data(phyloseq18all) <- sample_data(sample_df)
levels(sample_data(phyloseq18all)$cap_wild_crate)

# Create new phyloseq object without "New crates"
phyloseq_no_new_crates <- subset_samples(phyloseq18all, cap_wild_crate != "New crates")
levels(sample_data(phyloseq_no_new_crates)$cap_wild_crate)
taxa_names(phyloseq_no_new_crates)
phyloseq::tax_table(phyloseq_no_new_crates)
dim(phyloseq::tax_table(phyloseq_no_new_crates))

#### ANCOM-BC2 for all 18s nano data ####
?ancombc2

#set.seed(123)
p18n_all_anc_1crate <-  ancombc2(data = phyloseq18all, tax_level = "Genus",
                                 fix_formula = "cap_wild_crate", rand_formula = NULL,
                                 p_adj_method = "holm", pseudo_sens = TRUE,
                                 prv_cut = 0, lib_cut = 1000, s0_perc = 0.05,
                                 group = "cap_wild_crate", struc_zero = TRUE, neg_lb = TRUE)

p18n_all_anc_nocrate <-  ancombc2(data = phyloseq_no_new_crates, tax_level = "Genus",
                                  fix_formula = "cap_wild_crate", rand_formula = NULL,
                                  p_adj_method = "holm", pseudo_sens = TRUE,
                                  prv_cut = 0, lib_cut = 1000, s0_perc = 0.05,
                                  group = "cap_wild_crate", struc_zero = TRUE, neg_lb = TRUE)
save(p18n_all_anc_1crate, 
     p18n_all_anc_nocrate,
     file = "05_taxa_driving_groups/18s/ancom_18sn.RData")
load("05_taxa_driving_groups/18s/ancom_18sn.RData")       

res_df <- p18n_all_anc_1crate$res
res_df_nocrate <- p18n_all_anc_nocrate$res

# Reshape: gather LFC and FDR (q-value) columns for 1 Crate 
plot_df <- res_df %>%
  select(taxon, starts_with("lfc_"), starts_with("q_"), starts_with("diff_")) %>%
  pivot_longer(cols = -taxon, names_to = "metric", values_to = "value") %>%
  mutate(
    type = str_extract(metric, "^[^_]+"),  # Extract everything before first underscore
    group = str_replace(metric, "^[^_]+_", "")  # Remove the prefix
  ) %>%
  pivot_wider(names_from = type, values_from = value)
  
plot_df$group <- str_replace_all(plot_df$group, "cap_wild_crate", "")
plot_df$taxon <- gsub(":", ": ", plot_df$taxon)
plot_df <- plot_df %>%
  mutate(group = ifelse(group == "(Intercept)", "Captive", group))

# Plot: barplot of LFCs
n18_1crate <- ggplot(plot_df, aes(x = reorder(taxon, lfc), y = lfc, fill = group)) +
  geom_col(position = position_dodge(width = 0.9)) +
  coord_flip() +
  labs(title = "Significant Differential Abundance (ANCOMBC2)",
       x = "Taxon",
       y = "Log Fold Change",
       fill = "Group") +
  scale_fill_manual(
    values = c("Captive" = "purple4",
               "New crates" = "#CD2C7A",
               "Wild, free ranging" = "seagreen",
               "Wild, seized from traffickers" = "goldenrod")) +
  theme_bw(base_size = 14) +  # increases base text size and uses white background
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

ggsave(filename = "05_taxa_driving_groups/18s/18n_ancombc_1crate.png",n18_1crate, height=15, width=18)

# Reshape: gather LFC and FDR (q-value) columns for NO Crates 
plot_df_nc <- res_df_nocrate %>%
  dplyr::select(taxon, starts_with("lfc_"), starts_with("q_"), starts_with("diff_")) %>%
  pivot_longer(cols = -taxon, names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("type", "group"), sep = "_(?=cap|\\(Intercept\\))", extra = "merge") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  filter(diff == TRUE, q < 0.05)  # Keep only significant taxa


plot_df_nc$group <- str_replace_all(plot_df_nc$group, "cap_wild_crate", "")
plot_df_nc$taxon <- gsub("_", " ", plot_df_nc$taxon)
plot_df_nc <- plot_df_nc %>%
  mutate(group = ifelse(group == "(Intercept)", "Captive", group))

# Plot: barplot of LFCs
n18_no_crate <- ggplot(plot_df_nc, aes(x = reorder(taxon, lfc), y = lfc, fill = group)) +
  geom_col(position = position_dodge(width = 0.9)) +
  coord_flip() +
  labs(title = "Significant Differential Abundance (ANCOMBC2)",
       x = "Taxon",
       y = "Log Fold Change",
       fill = "Group") +
  scale_fill_manual(
    values = c("Captive" = "purple4",
               "Wild, free ranging" = "seagreen",
               "Wild, seized from traffickers" = "goldenrod")) +
  theme_bw(base_size = 14) +  # increases base text size and uses white background
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

ggsave(filename = "05_taxa_driving_groups/18s/18n_ancombc_no_crate.png",n18_no_crate, height=15, width=18)

# MICROBIOME WORKSHOP
### https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html

# ALPHA D and READ DEPTH
set.seed(1024)
ps <- phyloseq18all
rareres <- MicrobiotaProcess::get_rarecurve(obj=ps, chunks=400)

prare1 <- ggrarecurve(obj=rareres, factorNames="cap_wild_crate",
                      indexNames=c("Observe", "Chao1", "ACE")
) +
  scale_fill_manual(values=c("purple4", "seagreen", "goldenrod", "#CD2C7A"))+
  scale_color_manual(values=c("purple4", "seagreen", "goldenrod", "#CD2C7A"))+
  theme_bw()+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"))          

prare2 <- ggrarecurve(obj=rareres,
                      factorNames="cap_wild_crate",
                      shadow=FALSE,
                      indexNames=c("Observe", "Chao1", "ACE")) +
  scale_color_manual(values=c("purple4", "seagreen", "goldenrod", "#CD2C7A"))+
  coord_cartesian(ylim = c(0, 400)) +
  theme_bw()+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"))
ggsave(filename = "05_taxa_driving_groups/18s/18s_alpha_reads.png",prare2, height=5, width=8)
prare1 / prare2


# OVERLAPPING TAXA

vennlist <- get_vennlist(obj=ps, factorNames="cap_wild_crate")
upsetda <- get_upset(obj=ps, factorNames="cap_wild_crate")

vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("purple4", "seagreen", "goldenrod", "#CD2C7A"),
                      cat.col=c("purple4", "seagreen", "goldenrod", "#CD2C7A"),
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
ggsave(filename = "05_taxa_driving_groups/18s/18s_euk_overlap.png",vennp, height=5, width=8)

# *** PCOA's ***

# special function to get ESVs
ggordpoint_extract <- function(obj, pc = c(1, 2), mapping = NULL, sampleda = NULL, 
                               factorNames = NULL, factorLevels = NULL,
                               poinsize = 2, linesize = 0.3, arrowsize = 1.5,
                               arrowlinecolour = "grey", ellipse = FALSE,
                               showsample = FALSE, ellipse_pro = 0.9,
                               ellipse_alpha = 0.2, ellipse_linewd = 0.5,
                               ellipse_lty = 3, biplot = FALSE, topn = 5,
                               settheme = TRUE, speciesannot = FALSE,
                               fontsize = 2.5, labelfactor = NULL,
                               stroke = 0.1, fontface = "bold.italic",
                               fontfamily = "sans", textlinesize = 0.02, ...) {
  
  # 1. Coordinate extraction
  plotcoordclass <- get_coord(obj, pc)
  plotcoord <- plotcoordclass@coord
  xlab_text <- plotcoordclass@xlab
  ylab_text <- plotcoordclass@ylab
  title_text <- plotcoordclass@title
  
  defaultmapping <- aes_string(x = colnames(plotcoord)[1], y = colnames(plotcoord)[2])
  if (is.null(mapping)) {
    mapping <- defaultmapping
  } else {
    mapping <- modifyList(defaultmapping, mapping)
  }
  
  if (!is.null(sampleda)) {
    plotcoord <- merge(plotcoord, sampleda, by = 0)
    if (!is.null(factorNames)) {
      tmpfactormap <- get_factormap(factorNames)
      ellipsemapping <- get_ellipsemap(factorNames)
    } else {
      tmpfactormap <- get_factormap(colnames(sampleda))
      ellipsemapping <- get_ellipsemap(colnames(sampleda))
    }
    mapping <- modifyList(mapping, tmpfactormap)
    ellipsemapping <- modifyList(mapping, ellipsemapping)
    ellipsemapping <- modifyList(ellipsemapping, aes_string(starshape = NULL, size = NULL))
    if (!is.null(factorLevels)) {
      plotcoord <- setfactorlevels(plotcoord, factorLevels)
    }
  } else {
    ellipsemapping <- NULL
  }
  
  # 2. Biplot handling
  biplotcoord <- NULL
  if (biplot) {
    varcontrib <- get_varct(obj)
    varcontr <- varcontrib$VarContribution[, pc]
    tmpvars <- names(sort(rowSums(varcontr), decreasing = TRUE))
    #varlist <- get_varlist(namevector = tmpvars, n = topn)
    
    if (identical(topn, "all")) {
      varlist <- tmpvars
    } else {
      varlist <- get_varlist(namevector = tmpvars, n = topn)
    }
    
    biplotcoord <- varcontrib$VarCoordinates[match(varlist, rownames(varcontrib$VarCoordinates)), pc, drop = FALSE]
    biplotcoord <- data.frame(biplotcoord, check.names = FALSE)
    if (speciesannot) {
      biplotcoord$tax <- rownames(biplotcoord)
    }
  }
  
  # Return all internal data
  return(list(
    plotcoord = plotcoord,
    mapping = mapping,
    ellipsemapping = ellipsemapping,
    xlab = xlab_text,
    ylab = ylab_text,
    title = title_text,
    biplotcoord = biplotcoord
  ))
}

# function for the function above 
get_varlist <- function(namevector, n){
  if (inherits(n, "integer") || inherits(n, "numeric")){
    if (length(n)==2){
      varnames <- namevector[n[1]:n[2]]
    }
    if (length(n)> 2){
      varnames <- namevector[n]
    }
    if (length(n)==1){
      varnames <- namevector[seq_len(n)]
    }
  }
  if (inherits(n, "character")){
    if (length(n)>1){
      varnames <- n[match(n, namevector)]
    }else{
      stop("the n should be a character vector, integer  or integer vector")
    }
  }
  return(varnames)
}


# distmethod: "unifrac","wunifrac","manhattan","euclidean","canberra","bray","kulczynski" (vegdist, dist options)
pcoares <- get_pcoa(obj=ps, distmethod="bray", method="hellinger")
pca_obj <- pcoares@pca

# Visualizing the result for pca 1 and 2
pcoaplot1 <- ggordpoint(obj=pcoares, biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("cap_wild_crate"), 
                        ellipse=TRUE, topn = 10) +
  scale_color_manual(values=c("purple4", "seagreen", "goldenrod", "#CD2C7A")) +
  scale_fill_manual(values=c("purple4", "seagreen", "goldenrod", "#CD2C7A"))
ggsave(filename = "05_taxa_driving_groups/18s/18s_pcoaplot_1_2.png",pcoaplot1, height=5, width=8)

# first and third principal co-ordinates
pcoaplot2 <- ggordpoint(obj=pcoares, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("cap_wild_crate"), 
                        ellipse=TRUE, topn = 10) +
  scale_color_manual(values=c("purple4", "seagreen", "goldenrod", "#CD2C7A")) +
  scale_fill_manual(values=c("purple4", "seagreen", "goldenrod", "#CD2C7A"))
ggsave(filename = "05_taxa_driving_groups/18s/18s_pcoaplot_1_3.png",pcoaplot2, height=5, width=8)

# second and third principal co-ordinates
pcoaplot3 <- ggordpoint(obj=pcoares, pc=c(2, 3), biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("cap_wild_crate"), 
                        ellipse=TRUE, topn = 10) +
  scale_color_manual(values=c("purple4", "seagreen", "goldenrod", "#CD2C7A")) +
  scale_fill_manual(values=c("purple4", "seagreen", "goldenrod", "#CD2C7A"))
ggsave(filename = "05_taxa_driving_groups/18s/18s_pcoaplot_2_3.png",pcoaplot3, height=5, width=8)

# Get PCA ESVs of PCA 1 and 2
ESV_result_1_2 <- ggordpoint_extract(obj = pca_obj, biplot = TRUE, speciesannot = TRUE,
                                     factorNames = c("cap_wild_crate"), ellipse = TRUE,
                                     topn = 50)
ESV_result_1_3 <- ggordpoint_extract(obj = pca_obj, biplot = TRUE, speciesannot = TRUE,
                                     factorNames = c("cap_wild_crate"), ellipse = TRUE,
                                     topn = 50, pc=c(1, 3))
ESV_result_2_3 <- ggordpoint_extract(obj = pca_obj, biplot = TRUE, speciesannot = TRUE,
                                     factorNames = c("cap_wild_crate"), ellipse = TRUE,
                                     topn = 50, pc=c(2, 3))

ESV_coords_1_2 <- as.data.frame(ESV_result_1_2$biplotcoord) 
ESV_coords_1_3 <- as.data.frame(ESV_result_1_3$biplotcoord)
ESV_coords_2_3 <- as.data.frame(ESV_result_2_3$biplotcoord)

# taxonomy from ps
tax_table_df <- as.data.frame(phyloseq::tax_table(ps))  # Convert to data.frame
tax_table_df$ESV_ID <- rownames(tax_table_df)  # Add ESV IDs as a column

ESV_coords_tax_1_2 <- merge(ESV_coords_1_2, tax_table_df, by.x = "tax", by.y = "ESV_ID")
ESV_coords_tax_1_3 <- merge(ESV_coords_1_3, tax_table_df, by.x = "tax", by.y = "ESV_ID")
ESV_coords_tax_2_3 <- merge(ESV_coords_2_3, tax_table_df, by.x = "tax", by.y = "ESV_ID")

write.csv(ESV_coords_tax_1_2, file = "05_taxa_driving_groups/18s/PCoA_ESV_coords_1_2.csv", row.names = FALSE)
write.csv(ESV_coords_tax_1_3, file = "05_taxa_driving_groups/18s/PCoA_ESV_coords_1_3.csv", row.names = FALSE)
write.csv(ESV_coords_tax_2_3, file = "05_taxa_driving_groups/18s/PCoA_ESV_coords_2_3.csv", row.names = FALSE)

# PERMANOVA
distme <- get_dist(ps, distmethod ="bray", method="hellinger")
sampleda <- data.frame(sample_data(ps), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
sampleda$Group <- factor(sampleda$cap_wild_crate)
set.seed(1024)
adores <- adonis2(distme ~ cap_wild_crate, data=sampleda, permutation=9999)
data.frame(adores$aov.tab)
print(adores)
adores_df <- as.data.frame(adores)
adores_df %>%
  gt(rownames_to_stub = TRUE) %>%
  tab_header(title = "Adonis2 PERMANOVA Results") %>%
  fmt_number(columns = where(is.numeric), decimals = 2) %>%
  gtsave("05_taxa_driving_groups/18s/adonis2_results.html")

############ Stats for beta diversity; all ########
dist16_bc_transformed <- vegdist(wisconsin(sqrt(t(otu_table(ps)))), method="bray")
meta16 <- data.frame(sample_data(ps))
sink("05_taxa_driving_groups/18s/permanova18_captivewildcrate.txt")
adonis2(dist16_bc_transformed ~ cap_wild_crate, data=meta16)
sink()

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
sink("05_taxa_driving_groups/18s/18s_biomarker_discovery_LDA.txt")
deres
sink()
sink("05_taxa_driving_groups/18s/18s_biomarker_discovery_RF.txt")
deres_rf
sink()

# Biomarker Abundance (r) & conf. interval of effect size (LDA or MDA) of biomarker. The bigger confident interval shows that the biomarker is more fluctuant, owing to the influence of samples number.

# Enriched group using LDA :
top20_feat <- deres@result %>%
  arrange(desc(LDAmean)) %>%
  slice_head(n = 20) %>% pull(f)
deres_top20 <- deres
deres_top20@result <- deres@result %>% filter(f %in% top20_feat)

diffbox <- ggdiffbox(obj=deres_top20, box_notch=FALSE, 
                     colorlist=c("purple4", "seagreen", "goldenrod", "#CD2C7A"), 
                     l_xlabtext="relative abundance")
diffbox
ggsave(filename = "05_taxa_driving_groups/18s/18s_biomarker_abundance_top20_LDA.png",diffbox, height=7, width=8)

# Confidence interval of effect LDA :
es_p <- ggeffectsize(obj=deres_top20, 
                     lineheight=0.1,
                     linewidth=0.3) + 
  scale_color_manual(values=c("purple4", "seagreen", "goldenrod", "#CD2C7A")) 
es_p
ggsave(filename = "05_taxa_driving_groups/18s/18s_biomarker_top20_LDA.png",es_p, height=5, width=8)


# Enriched group using Random Forest :
top20_feat_rf <- deres_rf@result %>%
  arrange(desc(MDAmean)) %>%
  slice_head(n = 20) %>% pull(f)
deres_rf_top20 <- deres_rf
deres_rf_top20@result <- deres_rf@result %>% filter(f %in% top20_feat_rf)

diffbox_rf <- ggdiffbox(obj=deres_rf_top20, box_notch=FALSE, 
                        colorlist=c("purple4", "seagreen", "goldenrod", "#CD2C7A"), 
                        l_xlabtext="Relative Abundance")
diffbox_rf
ggsave(filename = "05_taxa_driving_groups/18s/18s_biomarker_abundance_top20_RF.png",diffbox_rf, height=7, width=8)

# Confidence interval of effect RF :
es_p_rf <- ggeffectsize(obj=deres_rf_top20, 
                        lineheight=0.1,
                        linewidth=0.3) + 
  scale_color_manual(values=c("purple4", "#CD2C7A", "goldenrod"))  #"seagreen",
es_p_rf
ggsave(filename = "05_taxa_driving_groups/18s/18s_biomarker_top20_RF.png",es_p_rf, height=5, width=8)


