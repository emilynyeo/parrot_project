# Which Taxa drive the differences between groups? - USING MICROBIOME PROCESS 

### https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html

# Editted functions 
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

generate_pcoa_plot <- function(physeq_obj, dataset_name, 
                               taxrank = "Genus", 
                               color_palette = c("#00AED7", "#FD9347", "seagreen", "purple4")) {
  # 1: Agglomerate at specified taxonomic rank
  physeq_tax <- tax_glom(physeq_obj, taxrank = taxrank)
  # 2: Replace taxon names with chosen rank
  tax_names <- as.character(phyloseq::tax_table(physeq_tax)[, taxrank])
  tax_names <- make.unique(tax_names)
  taxa_names(physeq_tax) <- tax_names
  # 3: Hellinger transform
  physeq_tax_hell <- transform_sample_counts(physeq_tax, function(x) sqrt(x / sum(x)))
  # 4: PCoA
  pcoa_res <- get_pcoa(obj = physeq_tax_hell, distmethod = "bray", method = "hellinger")
  pca_obj <- pcoa_res@pca
  assign(paste0("pca_obj_", tolower(dataset_name), "_", tolower(taxrank)), pca_obj, envir = .GlobalEnv)
  assign(paste0("pca_res_", tolower(dataset_name), "_", tolower(taxrank)), pcoa_res, envir = .GlobalEnv)
  # 5: Plot
  plot <- ggordpoint(obj = pcoa_res, biplot = TRUE, speciesannot = TRUE,
                     factorNames = c("Captive.Wild"), ellipse = TRUE) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    ggtitle(paste("PCoA Plot -", taxrank, "-", dataset_name))
  
  # Create output directory if it doesn't exist
  out_dir <- file.path("05_taxa_driving_groups", paste0("pcoas_", tolower(taxrank)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save plot
  out_path <- file.path(out_dir, paste0("PCoA_", tolower(taxrank), "_", dataset_name, ".png"))
  ggsave(filename = out_path, plot = plot, width = 8, height = 6, dpi = 300)
  return(plot)
}

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

# PCOA's
# Dataset names as character strings
datasets <- c("phyloseq16m", "phyloseq16n", "phyloseq18m", "phyloseq18n")
# Loop through and apply function
pcoa_plots_family <- list()
pcoa_plots_genus <- list()
pcoa_plots_species <- list()
for (ds_name in datasets) {
  physeq_obj <- get(ds_name)  # Get the object by name
  pcoa_plots_family[[ds_name]] <- generate_pcoa_plot(physeq_obj, dataset_name = ds_name, taxrank = "Family")
  pcoa_plots_genus[[ds_name]] <- generate_pcoa_plot(physeq_obj, dataset_name = ds_name, taxrank = "Genus")
  pcoa_plots_species[[ds_name]] <- generate_pcoa_plot(physeq_obj, dataset_name = ds_name, taxrank = "Species")
}

grid_plot_family <- (
  pcoa_plots_family$phyloseq16m + pcoa_plots_family$phyloseq16n +
  pcoa_plots_family$phyloseq18n + pcoa_plots_family$phyloseq18m) + 
  plot_layout(ncol = 2) + plot_annotation(tag_levels = "A")

grid_plot_genus <- (
  pcoa_plots_genus$phyloseq16m + pcoa_plots_genus$phyloseq16n +
  pcoa_plots_genus$phyloseq18n + pcoa_plots_genus$phyloseq18m) + 
  plot_layout(ncol = 2) + plot_annotation(tag_levels = "A")

grid_plot_species <- (
  pcoa_plots_species$phyloseq16m + pcoa_plots_species$phyloseq16n +
  pcoa_plots_species$phyloseq18n + pcoa_plots_species$phyloseq18m) + 
  plot_layout(ncol = 2) + plot_annotation(tag_levels = "A")

ggsave("05_taxa_driving_groups/pcoa_family_grid.png", grid_plot_family, width = 12, height = 10, dpi = 300)
ggsave("05_taxa_driving_groups/pcoa_genus_grid.png", grid_plot_genus, width = 12, height = 10, dpi = 300)
ggsave("05_taxa_driving_groups/pcoa_species_grid.png", grid_plot_species, width = 12, height = 10, dpi = 300)


### Ok now I want all the ESVs for the plots above 

pcoa_objs <- c("pca_obj_family", "pca_obj_genus", "pca_obj_species")


# Loop over all pca_obj_*_* objects in the environment
for (obj_name in ls(pattern = "^pca_obj_phyloseq.*_.*$")) {
  
  # Extract dataset and taxrank from object name
  parts <- strsplit(obj_name, "_")[[1]]
  dataset_name <- parts[3]  # assuming format is pca_obj_dataset_taxrank
  taxrank <- parts[4]
  
  # Get the actual pca object
  pca_obj <- get(obj_name)
  
  # Run ggordpoint_extract for different PC combinations
  ESV_result_1_2 <- ggordpoint_extract(obj = pca_obj, biplot = TRUE, speciesannot = TRUE,
                                       factorNames = c("Captive.Wild"), ellipse = TRUE,
                                       topn = 50)
  ESV_result_1_3 <- ggordpoint_extract(obj = pca_obj, biplot = TRUE, speciesannot = TRUE,
                                       factorNames = c("Captive.Wild"), ellipse = TRUE,
                                       topn = 50, pc = c(1, 3))
  ESV_result_2_3 <- ggordpoint_extract(obj = pca_obj, biplot = TRUE, speciesannot = TRUE,
                                       factorNames = c("Captive.Wild"), ellipse = TRUE,
                                       topn = 50, pc = c(2, 3))
  
  # Convert biplot coordinates to data frames
  ESV_coords_1_2 <- as.data.frame(ESV_result_1_2$biplotcoord) 
  ESV_coords_1_3 <- as.data.frame(ESV_result_1_3$biplotcoord)
  ESV_coords_2_3 <- as.data.frame(ESV_result_2_3$biplotcoord)
  
  # Get taxonomy table and convert to data.frame
  tax_table_df <- as.data.frame(tax_table(ps))
  tax_table_df$ESV_ID <- rownames(tax_table_df)
  
  # Merge coords with taxonomy
  ESV_coords_tax_1_2 <- merge(ESV_coords_1_2, tax_table_df, by.x = "tax", by.y = "ESV_ID")
  ESV_coords_tax_1_3 <- merge(ESV_coords_1_3, tax_table_df, by.x = "tax", by.y = "ESV_ID")
  ESV_coords_tax_2_3 <- merge(ESV_coords_2_3, tax_table_df, by.x = "tax", by.y = "ESV_ID")
  
  # Define output directory (same as used for PNG)
  out_dir <- file.path("05_taxa_driving_groups", paste0("pcoas_", tolower(taxrank)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save CSVs in that directory
  write.csv(ESV_coords_tax_1_2, 
            file = file.path(out_dir, paste0("PCoA_ESV_coords_1_2_", dataset_name, "_", taxrank, ".csv")), row.names = FALSE)
  write.csv(ESV_coords_tax_1_3, 
            file = file.path(out_dir, paste0("PCoA_ESV_coords_1_3_", dataset_name, "_", taxrank, ".csv")), row.names = FALSE)
  write.csv(ESV_coords_tax_2_3, 
            file = file.path(out_dir, paste0("PCoA_ESV_coords_2_3_", dataset_name, "_", taxrank, ".csv")), row.names = FALSE)
}


for (obj_name in ls(pattern = "^pca_obj_phyloseq.*_.*$")) {
  
  # Extract parts from the object name: pca_obj_<phyloseq_obj>_<taxrank>
  parts <- strsplit(obj_name, "_")[[1]]
  dataset_name <- parts[3]     # e.g., "phyloseq16m"
  taxrank <- parts[4]          # e.g., "genus"
  taxrank_cap <- tools::toTitleCase(taxrank)
  # Get the actual pca object
  pca_obj <- get(obj_name)
  
  # Get the corresponding phyloseq object
  if (!exists(dataset_name)) {
    warning(paste("Phyloseq object", dataset_name, "not found. Skipping."))
    next
  }
  physeq_obj <- get(dataset_name)
  
  # Run ggordpoint_extract for different PC combinations
  ESV_result_1_2 <- ggordpoint_extract(obj = pca_obj, biplot = TRUE, speciesannot = TRUE,
                                       factorNames = c("Captive.Wild"), ellipse = TRUE,
                                       topn = 50)
  ESV_result_1_3 <- ggordpoint_extract(obj = pca_obj, biplot = TRUE, speciesannot = TRUE,
                                       factorNames = c("Captive.Wild"), ellipse = TRUE,
                                       topn = 50, pc = c(1, 3))
  ESV_result_2_3 <- ggordpoint_extract(obj = pca_obj, biplot = TRUE, speciesannot = TRUE,
                                       factorNames = c("Captive.Wild"), ellipse = TRUE,
                                       topn = 50, pc = c(2, 3))
  
  # Convert biplot coordinates to data frames
  ESV_coords_1_2 <- as.data.frame(ESV_result_1_2$biplotcoord) 
  ESV_coords_1_3 <- as.data.frame(ESV_result_1_3$biplotcoord)
  ESV_coords_2_3 <- as.data.frame(ESV_result_2_3$biplotcoord)
  
  # Get taxonomy table from correct phyloseq object
  tax_table_df <- as.data.frame(phyloseq::tax_table(physeq_obj))
  tax_table_df$ESV_ID <- rownames(tax_table_df)
  
  # Merge biplot coordinates with taxonomy
  ESV_coords_tax_1_2 <- merge(ESV_coords_1_2, tax_table_df, by.x = "tax", by.y = taxrank_cap)
  ESV_coords_tax_1_3 <- merge(ESV_coords_1_3, tax_table_df, by.x = "tax", by.y = taxrank_cap)
  ESV_coords_tax_2_3 <- merge(ESV_coords_2_3, tax_table_df, by.x = "tax", by.y = taxrank_cap)
  
  # Define output directory
  out_dir <- file.path("05_taxa_driving_groups", paste0("pcoas_", tolower(taxrank)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save merged CSVs
  write.csv(ESV_coords_tax_1_2,
            file = file.path(out_dir, paste0("PCoA_ESV_coords_1_2_", dataset_name, "_", taxrank, ".csv")),
            row.names = FALSE)
  write.csv(ESV_coords_tax_1_3,
            file = file.path(out_dir, paste0("PCoA_ESV_coords_1_3_", dataset_name, "_", taxrank, ".csv")),
            row.names = FALSE)
  write.csv(ESV_coords_tax_2_3,
            file = file.path(out_dir, paste0("PCoA_ESV_coords_2_3_", dataset_name, "_", taxrank, ".csv")),
            row.names = FALSE)
}


# PERMANOVA
ps <- phyloseq16m
distme <- get_dist(phyloseq16m, distmethod ="bray", method="hellinger")
sampleda <- data.frame(sample_data(ps), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
sampleda$Group <- factor(sampleda$Captive.Wild)
set.seed(1024)
adores <- adonis2(distme ~ Captive.Wild, data=sampleda, permutation=9999)
data.frame(adores$aov.tab)
print(adores)

Table_grob <- tableGrob(adores, rows = NULL, theme = ttheme_default(
  core = list(fg_params = list(hjust = 0.5)),
  colhead = list(fg_params = list(fontface = 2, hjust = 0.5))
))

# Save to PNG
#png("05_taxa_driving_groups/permanova_16all.png", width = 7, height = 5)
#grid.draw(Table_grob)
#dev.off()
