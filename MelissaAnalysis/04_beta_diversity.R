#!bin/bash

#### Load libraries ####
library(vegan)
library(MASS)
library(phyloseq)
library(tibble)
library(tidyverse)

dir.create("04_beta_diversity")
#### Load data ####
# Also converting to matrices with rownames for phyloseq

load("01_process_and_clean_data/phyloseq16n.rds")
load("01_process_and_clean_data/phyloseq16ren.rds")

load("02_split_18_classes/phyloseq18nano_nohost.rds")
load("02_split_18_classes/phyloseq18nano_plants.rds")
load("02_split_18_classes/phyloseq18nano_microeuk.rds")


load("02_split_18_classes/phyloseq18miseq_nohost.rds")
load("02_split_18_classes/phyloseq18miseq_plants.rds")
load("02_split_18_classes/phyloseq18miseq_microeuk.rds")

load("01_process_and_clean_data/phyloseq16ren.rds")
load("01_process_and_clean_data/phyloseq18ren_nohost.rds")

# Merge
phyloseq16all <- merge_phyloseq(phyloseq16n, phyloseq16ren)
phyloseq18all <- phyloseq18nano_nohost
phyloseq18all <- prune_samples(sample_sums(phyloseq18all) >0, phyloseq18all)
phyloseq18miseq_nohost <- prune_samples(sample_sums(phyloseq18miseq_nohost) >0, phyloseq18miseq_nohost)
captivewild_legend <- c(Captive="purple4", "Wild, free ranging"="seagreen", "Wild, seized from traffickers"="goldenrod",
                        "Crate 1: Unknown origin" = "magenta3", "Crate 2: Unknown origin" = "firebrick3")
location_legend <- c(Nigeria ="goldenrod4",`South Africa`="darkred")

#### Make ordinations #####
set.seed(423)
phylo16n_ord <- ordinate(phyloseq16n, method = "NMDS", distance = "bray")

gg_bray16n_captivetube <- plot_ordination(phyloseq16n, phylo16n_ord, type = "samples", color = "Captive.Wild", shape = "Tube.type")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Bacteria (Stress: ", round(phylo16n_ord$stress, 2),")")) +
  labs(shape="Tube type", col="Captive\nor wild")
gg_bray16n_captivetube
ggsave("04_beta_diversity/gg_bray16n_captivetube.png", gg_bray16n_captivetube, height=5, width=7)

gg_bray16n_captiveprobiotic <- plot_ordination(phyloseq16n, phylo16n_ord, type = "samples", color = "Captive.Wild", shape = "Given.probiotic")+
  scale_shape_manual(values=c(19,3)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Bacteria (Stress: ", round(phylo16n_ord$stress, 2),")")) +
  labs(shape="Microbiome\ntreatment", col="Captive\nor wild")
gg_bray16n_captiveprobiotic
ggsave("04_beta_diversity/gg_bray16n_captiveprobiotic.png", gg_bray16n_captiveprobiotic, height=5, width=7)


gg_bray16n_captivelocation <- plot_ordination(phyloseq16n, phylo16n_ord, type = "samples", color = "Captive.Wild", shape = "Country")+
  scale_shape_manual(values=c(19,21))+
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Bacteria (Stress: ", round(phylo16n_ord$stress, 2),")")) +
  labs(shape="Collection\nlocation", col="Captive\nor wild")
gg_bray16n_captivelocation
ggsave("04_beta_diversity/gg_bray16n_captivelocation.png", gg_bray16n_captivelocation, height=5, width=7)


gg_bray16n_captivelocation <- plot_ordination(phyloseq16n, phylo16n_ord, type = "samples", color = "Site.name", shape = "Captive.Wild") +
  scale_shape_manual(values=c(19,21,15))+
  # scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Bacteria (Stress: ", round(phylo16n_ord$stress, 2),")")) +
  labs(shape="Captive\nor wild", col="Site Name")
gg_bray16n_captivelocation
ggsave("04_beta_diversity/gg_bray16n_captivelocation.png", gg_bray16n_captivelocation, height=5, width=8)

### WITH RENEES SAMPLES
phylo16all_ord <- ordinate(phyloseq16all, method = "NMDS", distance = "bray")
gg_bray16WITHREN_givenprobiotic <- plot_ordination(phyloseq16all, phylo16all_ord, type = "samples", color = "Captive.Wild", shape = "Given.probiotic")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Bacteria (Stress: ", round(phylo16all_ord$stress, 2),")")) +
  labs(shape="Known or unknown", col="Captive\nor wild")
gg_bray16WITHREN_givenprobiotic
ggsave("04_beta_diversity/gg_bray16WITHREN_givenprobiotic.png", gg_bray16WITHREN_givenprobiotic, height=5, width=7)

gg_bray16WITHREN_givenprobiotic_samples_blanked <- plot_ordination(phyloseq16all, phylo16all_ord, type = "samples", color = "Captive.Wild", shape="Given.probiotic")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=c(Captive="purple4", "Wild, free ranging"="seagreen", "Wild, seized from traffickers"="goldenrod",
                              "Crate 1: Unknown origin" = rgb(0,0,0,0), "Crate 2: Unknown origin" = rgb(0,0,0,0)))+
  labs(subtitle=paste0("Bacteria (Stress: ", round(phylo16all_ord$stress, 2),")")) +
  labs(shape="Known or unknown", col="Captive\nor wild")
gg_bray16WITHREN_givenprobiotic_samples_blanked
ggsave("04_beta_diversity/gg_bray16WITHREN_givenprobiotic_samples_blanked.png", gg_bray16WITHREN_givenprobiotic_samples_blanked, height=5, width=7)

###### Val requests ######
# Captive, wild, seized (probioitic), no mystery
captivewild_legend
gg_bray16WITHREN_givenprobiotics_nomystery <- plot_ordination(phyloseq16all, phylo16all_ord, type = "samples", color = "Captive.Wild", shape = "Given.probiotic")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=c(Captive="purple4", "Wild, free ranging"="seagreen", "Wild, seized from traffickers"="goldenrod",
                              "Crate 1: Unknown origin" = rgb(0,0,0,0), "Crate 2: Unknown origin" = rgb(0,0,0,0)))+
  labs(subtitle=paste0("Bacteria (Stress: ", round(phylo16all_ord$stress, 2),")")) +
  labs(shape="Given probiotic?", col="Captive\nor wild") 
gg_bray16WITHREN_givenprobiotics_nomystery
ggsave("04_beta_diversity/gg_bray16WITHREN_givenprobiotics_nomystery.png", gg_bray16WITHREN_givenprobiotics_nomystery, height=5, width=7)

# Captive and wild only
gg_bray16WITHREN_givenprobiotics_nomysterynoseized <- plot_ordination(phyloseq16all, phylo16all_ord, type = "samples", color = "Captive.Wild", shape = "Given.probiotic")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=c(Captive="purple4", "Wild, free ranging"="seagreen", "Wild, seized from traffickers"=rgb(0,0,0,0),
                              "Crate 1: Unknown origin" = rgb(0,0,0,0), "Crate 2: Unknown origin" = rgb(0,0,0,0)))+
  labs(subtitle=paste0("Bacteria (Stress: ", round(phylo16all_ord$stress, 2),")")) +
  labs(shape="Given probiotic?", col="Captive\nor wild") 
gg_bray16WITHREN_givenprobiotics_nomysterynoseized
ggsave("04_beta_diversity/gg_bray16WITHREN_givenprobiotics_nomysterynoseized.png", gg_bray16WITHREN_givenprobiotics_nomysterynoseized, height=5, width=7)

# Captive wild mystery
gg_bray16WITHREN_givenprobiotics_noseized <- plot_ordination(phyloseq16all, phylo16all_ord, type = "samples", color = "Captive.Wild", shape = "Given.probiotic")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=c(Captive="purple4", "Wild, free ranging"="seagreen", "Wild, seized from traffickers"=rgb(0,0,0,0),
                              "Crate 1: Unknown origin" = "magenta3", "Crate 2: Unknown origin" = "firebrick3"))+
  labs(subtitle=paste0("Bacteria (Stress: ", round(phylo16all_ord$stress, 2),")")) +
  labs(shape="Given probiotic?", col="Captive\nor wild") 
gg_bray16WITHREN_givenprobiotics_noseized
ggsave("04_beta_diversity/gg_bray16WITHREN_givenprobiotics_noseized.png", gg_bray16WITHREN_givenprobiotics_noseized, height=5, width=7)



#### CURRENT PROGRESS##########
# 18s is weird. miseq works but you need to use the 01_ version. 
# I think the nanopore one has less diversity?
# Basically, lots of samples are very similar. If I remove low-read samples, it gets better. suggests that low-read samples are very similar to each other, possible just one read.
##### 18S all samples #########
load("01_process_and_clean_data/phyloseq18m.rds")
phylo18m_ord <- ordinate(phyloseq18m, method = "NMDS", distance = "bray")
plot_ordination(phyloseq18m, phylo18m_ord, type = "samples", color = "Captive.Wild", shape = "Tube.type")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Eukaryote (Stress: ", round(phylo18m_ord$stress, 2),")")) +
  labs(shape="Tube type", col="Captive\nor wild") 
load("01_process_and_clean_data/phyloseq18n.rds")

phyloseq18n_c1000 <- prune_samples(sample_sums(phyloseq18n)>1000, phyloseq18n)
phylo18n_ord <- ordinate(phyloseq18n_c1000, method = "NMDS", distance = "bray")
plot_ordination(phyloseq18n_c1000, phylo18n_ord, type = "samples", color = "Captive.Wild", shape = "Tube.type")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Eukaryote (Stress: ", round(phylo18m_ord$stress, 2),")")) +
  labs(shape="Tube type", col="Captive\nor wild") 



phylo18n_ord <- ordinate(phyloseq18all, method = "NMDS", distance = "bray")
plot_ordination(phyloseq18all, phylo18n_ord, type = "samples", color = "Captive.Wild", shape = "Tube.type")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Eukaryote (Stress: ", round(phylo18n_ord$stress, 2),")")) +
  labs(shape="Tube type", col="Captive\nor wild") 
plot_ordination(phyloseq18all, phylo18n_ord, type = "samples", color = "Captive.Wild", shape = "Tube.type")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Eukaryote (Stress: ", round(phylo18n_ord$stress, 2),")")) +
  labs(shape="Tube type", col="Captive\nor wild") +
  xlim(-100,100) + ylim(-100,100)

# There are a few samples that are REALLY out there? 
weirdsamps <- phylo18_ord$points %>% as.data.frame() %>%
  filter(MDS1< -0.4 | MDS2 > 0.4 | MDS1> 0.4 | MDS2< -0.4) %>% rownames()
onlyweirdsamp <- prune_taxa( taxa_sums(prune_samples(samples=weirdsamps, phyloseq18_nohost))>0, prune_samples(samples=weirdsamps, phyloseq18_nohost))
tax_table(onlyweirdsamp)
otu_table(onlyweirdsamp)
sample_data(onlyweirdsamp)

# NOTE: This sample is really weird because it's extremely different from everything else. going to manually remove it here
phyloseq18_nohost <- subset_samples(phyloseq18_nohost, !sample_names(phyloseq18_nohost) %in% weirdsamps )
# phyloseq18_microeuk <- subset_samples(phyloseq18_microeuk, sample_names(phyloseq18_microeuk) !=weirdsamps )
phylo18_ord <- ordinate(phyloseq18_nohost, method = "NMDS", distance = "bray")
gg_bray18_captivetube <- plot_ordination(phyloseq18_nohost, phylo18_ord, type = "samples", color = "Captive.Wild", shape = "Tube.type")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Eukaryote (Stress: ", round(phylo18_ord$stress, 2),")")) +
  labs(shape="Tube type", col="Captive\nor wild") 
gg_bray18_captivetube
ggsave("04_beta_diversity/gg_bray18_captivetube.png", gg_bray18_captivetube, height=5, width=7)

gg_bray18_captiveprobiotic <- plot_ordination(phyloseq18_nohost, phylo18_ord, type = "samples", color = "Captive.Wild", shape = "Given.probiotic")+
  scale_shape_manual(values=c(19,3)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Eukaryote (Stress: ", round(phylo18_ord$stress, 2),")")) +
  labs(shape="Microbiome\ntreatment", col="Captive\nor wild")
gg_bray18_captiveprobiotic
ggsave("04_beta_diversity/gg_bray18_captiveprobiotic.png", gg_bray18_captiveprobiotic, height=5, width=7)


gg_bray18_captivelocation <- plot_ordination(phyloseq18_nohost, phylo18_ord, type = "samples", color = "Captive.Wild", shape = "Country")+
  scale_shape_manual(values=c(19,21))+
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Eukaryote (Stress: ", round(phylo18_ord$stress, 2),")")) +
  labs(shape="Collection\nlocation", col="Captive\nor wild")
gg_bray18_captivelocation
ggsave("04_beta_diversity/gg_bray18_captivelocation.png", gg_bray18_captivelocation, height=5, width=7)


gg_bray18_captivelocation2 <- plot_ordination(phyloseq18_nohost, phylo18_ord, type = "samples", color = "Site.name", shape = "Captive.Wild") +
  scale_shape_manual(values=c(19,21,15))+
  # scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Eukaryote (Stress: ", round(phylo18_ord$stress, 2),")")) +
  labs(shape="Captive\nor wild", col="Site Name")
gg_bray18_captivelocation2
ggsave("04_beta_diversity/gg_bray18_captivelocation2.png", gg_bray18_captivelocation2, height=5, width=8)



##### 18S microeukrs #########
# remove zeros
phyloseq18_microeuk_nozeros <- prune_taxa(taxa_sums(phyloseq18_microeuk)>0, prune_samples(sample_sums(phyloseq18_microeuk) > 0, phyloseq18_microeuk))
phylo18_ord_microeuk <- ordinate(phyloseq18_microeuk_nozeros, method = "NMDS", distance = "bray")
plot_ordination(phyloseq18_microeuk, phylo18_ord_microeuk, type = "samples", color = "Captive.Wild", shape = "Tube.type")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Microeukaryotes only (Stress: ", round(phylo18_ord_microeuk$stress, 2),")")) +
  labs(shape="Tube type", col="Captive\nor wild") 
# There are some samples that look weird again

# There are a few samples thatare REALLY out there? 
weirdsamps_microeuk <- phylo18_ord_microeuk$points %>% as.data.frame() %>%
  filter(MDS1< -0.5 | MDS2 > 0.5 | MDS1>0.5 | MDS2< -0.5) %>% rownames()
onlyweirdsamp <- prune_taxa( taxa_sums(prune_samples(samples=weirdsamps_microeuk, phyloseq18_microeuk))>0, prune_samples(samples=weirdsamps_microeuk, phyloseq18_microeuk))
tax_table(onlyweirdsamp)
otu_table(onlyweirdsamp)
sample_data(onlyweirdsamp)
# Okay so these are weird because there's only 1 ASV in each sample-- let's remove these samples
phyloseq18_microeuk <- subset_samples(phyloseq18_microeuk, !sample_names(phyloseq18_microeuk)  %in% weirdsamps_microeuk )
phyloseq18_microeuk <- prune_samples(sample_sums(prune_taxa(taxa_sums(phyloseq18_microeuk)>0, phyloseq18_microeuk))>0, prune_taxa(taxa_sums(phyloseq18_microeuk)>0, phyloseq18_microeuk))
phylo18_ord_microeuk <- ordinate(phyloseq18_microeuk, method = "NMDS", distance = "bray")


gg_bray18_captivetube_microeuk <- plot_ordination(phyloseq18_microeuk, phylo18_ord_microeuk, type = "samples", color = "Captive.Wild", shape = "Tube.type")+
  scale_shape_manual(values=c(19,21)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Microeukaryotes only (Stress: ", round(phylo18_ord_microeuk$stress, 2),")")) +
  labs(shape="Tube type", col="Captive\nor wild") 
gg_bray18_captivetube_microeuk
ggsave("04_beta_diversity/gg_bray18_captivetube_microeuk.png", gg_bray18_captivetube_microeuk, height=5, width=7)

gg_bray18_captiveprobiotic_microeuk <- plot_ordination(phyloseq18_microeuk_nozeros, phylo18_ord_microeuk, type = "samples", color = "Captive.Wild", shape = "Given.probiotic")+
  scale_shape_manual(values=c(19,3)) +
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Microeukaryote only (Stress: ", round(phylo18_ord_microeuk$stress, 2),")")) +
  labs(shape="Microbiome\ntreatment", col="Captive\nor wild")
gg_bray18_captiveprobiotic_microeuk
ggsave("04_beta_diversity/gg_bray18_captiveprobiotic_microeuk.png", gg_bray18_captiveprobiotic_microeuk, height=5, width=7)


gg_bray18_captivelocation_microeuk<- plot_ordination(phyloseq18_microeuk_nozeros, phylo18_ord_microeuk, type = "samples", color = "Captive.Wild", shape = "Country")+
  scale_shape_manual(values=c(19,21))+
  scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Eukaryote (Stress: ", round(phylo18_ord_microeuk$stress, 2),")")) +
  labs(shape="Collection\nlocation", col="Captive\nor wild")
gg_bray18_captivelocation_microeuk
ggsave("04_beta_diversity/gg_bray18_captivelocation_microeuk.png", gg_bray18_captivelocation_microeuk, height=5, width=7)


gg_bray18_captivelocation2_microeuk <- plot_ordination(phyloseq18_microeuk_nozeros, phylo18_ord_microeuk, type = "samples", color = "Site.name", shape = "Captive.Wild") +
  scale_shape_manual(values=c(19,21,15))+
  # scale_color_manual(values=captivewild_legend)+
  labs(subtitle=paste0("Eukaryote (Stress: ", round(phylo18_ord$stress, 2),")")) +
  labs(shape="Captive\nor wild", col="Site Name")
gg_bray18_captivelocation2_microeuk
ggsave("04_beta_diversity/gg_bray18_captivelocation2_microeuk.png", gg_bray18_captivelocation2_microeuk, height=5, width=8)


############ Stats for beta diversity; all ########
dist16_bc_transformed <- vegdist(wisconsin(sqrt(t(otu_table(phyloseq16)))), method="bray")
dist18_bc_transformed <- vegdist(wisconsin(sqrt(t(otu_table(phyloseq18_nohost)))), method="bray")
dist18_bc_transformed_microeuk <- vegdist(wisconsin(sqrt(t(otu_table(phyloseq18_microeuk_nozeros)))), method="bray")

### YUP, double checked they look the same!
# test <- isoMDS(dist16_bc_transformed, k=2)
# test$points %>% as.data.frame() %>%
#   rownames_to_column(var="sampleid") %>%
#   left_join(data.frame(sample_data(phyloseq16)) %>% rownames_to_column(var="sampleid")) %>%
#   ggplot() +
#   geom_point(aes(x=V1, y=V2, col=Captive.Wild))
meta16 <- data.frame(sample_data(phyloseq16))
sink("04_beta_diversity/permanova16_captivewild.txt")
adonis2(dist16_bc_transformed ~ Captive.Wild, data=meta16)
sink()

meta18 <- data.frame(sample_data(phyloseq18_nohost))
sink("04_beta_diversity/permanova18_captivewild.txt")
adonis2(dist18_bc_transformed ~ Captive.Wild, data=meta18)
sink()

meta18_microeuk <- data.frame(sample_data(phyloseq18_microeuk_nozeros))
sink("04_beta_diversity/permanova18_captivewild_microeuk.txt")
adonis2(dist18_bc_transformed_microeuk ~ Captive.Wild, data=meta18_microeuk)
sink()

# Filter to only include captive ones to assess effect of location
phyloseq16_captiveonly <- subset_samples(phyloseq16, Captive.Wild=="Captive")
dist16_bc_transformed_captiveonly <- vegdist(wisconsin(sqrt(t(otu_table(phyloseq16_captiveonly)))), method="bray")
adonis2(dist16_bc_transformed_captiveonly ~ Site.name, data=data.frame(sample_data(phyloseq16_captiveonly)))
adonis2(dist16_bc_transformed_captiveonly ~ Tube.type, data=data.frame(sample_data(phyloseq16_captiveonly)))
adonis2(dist16_bc_transformed_captiveonly ~ rand, data=data.frame(sample_data(phyloseq16_captiveonly)) %>% mutate(rand = sample(c("yes","no"), size =  160, replace = TRUE)))

