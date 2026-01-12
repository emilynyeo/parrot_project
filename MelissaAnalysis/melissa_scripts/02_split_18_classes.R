#!bin/bash
library(phyloseq)
library(tidyverse)

setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")
dir.create("02_split_18_classes")

load("01_process_and_clean_data/phyloseq18n.rds")
load("01_process_and_clean_data/phyloseq18n_nohost.rds")

load("01_process_and_clean_data/phyloseq18m.rds")
load("01_process_and_clean_data/phyloseq18m_nohost.rds")
phyloseq18miseq <- phyloseq18m
phyloseq18miseq_nohost <- phyloseq18m_nohost
remove(phyloseq18m)
remove(phyloseq18m_nohost)

load("01_process_and_clean_data/phyloseq18ren.rds")
load("01_process_and_clean_data/phyloseq18ren_nohost.rds")

catdat <- read.csv("01_process_and_clean_data/simple18_labelled2.csv")

## Merge into one
phyloseq18nano <- merge_phyloseq(phyloseq18n, phyloseq18ren)
phyloseq18nano_nohost <- merge_phyloseq(phyloseq18n_nohost, phyloseq18ren_nohost)

# manually make a legend of common names
### NOTE THIS NEEDS TO BE MANUALLY UPDATED
# listPlantSpecies <- tax_table(merge(phyloseq18n_plants)) %>% data.frame() %>%
#   select(Species) %>% distinct() %>% pull(Species)
# write.table(data.frame(Species=listPlantSpecies), file = "02_split_18_classes/listPlantSpecies.txt", quote=FALSE, sep="\t", row.names = FALSE)
listPlantSpeciesFILLED <- read.csv("02_split_18_classes/listPlantSpeciesFILLED.txt") %>%
  as.data.frame()



# Do for each nano and miseq
for ( o in c("nano","miseq")) {
  phyloseq18 <- get(paste0("phyloseq18",o))
  # Add in a probiotics one
  newsampdata <- sample_data(phyloseq18) %>%
    data.frame() %>%
    mutate(Captive.Wild.probiotic = ifelse(Captive.Wild == "Wild, seized from traffickers",
                                           ifelse(Given.probiotic == "Given probiotics", 
                                                  "Wild, seized from traffickers\nGiven probiotics", 
                                                  as.character(Captive.Wild)), as.character(Captive.Wild))) %>%
    mutate(Captive.Wild.probiotic = factor(Captive.Wild.probiotic, 
                                           levels=c("Captive","Wild, free ranging", "Wild, seized from traffickers",
                                                    "Wild, seized from traffickers\nGiven probiotics",
                                                    "Crate 1: Unknown origin", "Crate 2: Unknown origin")))
  sample_data(phyloseq18) <- sample_data(newsampdata)
  
  ###
  # Split?
  allTaxWithCat <- tax_table(phyloseq18) %>%
    as.data.frame() %>% rownames_to_column(var="ESVId") %>%
    left_join(catdat)
  
  catTax <- allTaxWithCat %>% select(ESVId, Phylum, BroadCat2, BroadCat ) %>%
    column_to_rownames(var="ESVId") %>% as.matrix()
  
  phyloseq18_broadcat <- phyloseq(otu_table(phyloseq18, taxa_are_rows = TRUE), tax_table(catTax), sample_data(phyloseq18))
  
  phyloseq18_broadcat_glom <- tax_glom(phyloseq18_broadcat, taxrank="BroadCat2")
  phyloseq18_broadcat_glom_nohost <- subset_taxa(phyloseq18_broadcat_glom, BroadCat2!="Presumed_host" & BroadCat2!="Human")
  
  phyloseq18_broadcat_glom_RA <- transform_sample_counts(phyloseq18_broadcat_glom, fun = function(x) x/sum(x))
  phyloseq18_broadcat_glom_RA_nohost <- transform_sample_counts(phyloseq18_broadcat_glom_nohost, fun = function(x) x/sum(x))
  
  color_leg_broadcat2 <- c(Arachnid="black", Presumed_host="darkgrey",Mammal="beige", Insect="brown", MicroEukaryote="lightblue", Mollusk="darkblue",Plant="darkolivegreen", Worm="salmon" )
  gg_broadcat <- plot_bar(phyloseq18_broadcat_glom_RA, fill="BroadCat2") +
    geom_bar(aes(col=BroadCat2, fill=BroadCat2), stat ="identity", position="stack") +
    facet_grid(~Captive.Wild.probiotic, scales="free", drop=TRUE) +
    scale_fill_manual(values=color_leg_broadcat2)+
    scale_color_manual(values=color_leg_broadcat2)+
    labs(col="Broad category",fill="Broad category") +
    ylab("Relative Abundance") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  gg_broadcat
  ggsave(filename = paste0("02_split_18_classes/gg_broadcat_",o,".png"), gg_broadcat, height=5, width=10)

  gg_broadcat_nohost <- plot_bar(phyloseq18_broadcat_glom_RA_nohost, fill="BroadCat2") +
    geom_bar(aes(col=BroadCat2, fill=BroadCat2), stat ="identity", position="stack") +
    facet_grid(~Captive.Wild.probiotic, scales="free", drop=TRUE) +
    scale_fill_manual(values=color_leg_broadcat2)+
    scale_color_manual(values=color_leg_broadcat2)+
    labs(col="Broad category",fill="Broad category") +
    ylab("Relative Abundance (excluding host and human DNA)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  gg_broadcat_nohost
  ggsave(filename = paste0("02_split_18_classes/gg_broadcat_",o,"_nohost.png"), gg_broadcat_nohost, height=5, width=10)

  #### Look at plants only ####
  plant_ESVs <- allTaxWithCat %>% filter(BroadCat2=="Plant") %>% pull(ESVId)
  phyloseq18_plants <- prune_taxa(taxa=plant_ESVs, phyloseq18)
  
  # Add in plant label
  taxtablenewplants <- tax_table(phyloseq18_plants) %>% as.data.frame()  %>%
    left_join(listPlantSpeciesFILLED %>% select(-Colour)) %>% as.matrix()
  rownames(taxtablenewplants) <- rownames( tax_table(phyloseq18_plants))
  tax_table(phyloseq18_plants) <- tax_table(taxtablenewplants)
  phyloseq18_plants_glom <- tax_glom(phyloseq18_plants, taxrank = "CommonName")
  phyloseq18_plants_glom_RA <- transform_sample_counts(phyloseq18_plants_glom, fun = function(x) x/sum(x))
  
  # Custom legend
  color_leg_plants <- listPlantSpeciesFILLED$Colour
  names(color_leg_plants) <- listPlantSpeciesFILLED$CommonName
  
  gg_plantsonly <- plot_bar(phyloseq18_plants_glom_RA, fill="CommonName") +
    geom_bar(aes(col=CommonName, fill=CommonName), stat ="identity", position="stack") +
    facet_grid(~Captive.Wild.probiotic, scales="free", drop=TRUE) +
    scale_fill_manual(values=color_leg_plants)+
    scale_color_manual(values=color_leg_plants)+
    labs(col="Plant common name",fill="Plant common name") +
    ylab("Relative Abundance") +
    theme(axis.text.x = element_blank())
  gg_plantsonly
  ggsave(filename = paste0("02_split_18_classes/gg_plantsonly_",o,".png"), gg_plantsonly, height=5, width=14)
  
  ######### 
  # Change colours of plants that are VERY rare to grey-- so there's not so many colours.
  toKeepColours <- names(which(rowSums(otu_table(phyloseq18_plants_glom_RA)>0.2, na.rm = TRUE)>0))
  color_leg_plants2 <- tax_table(phyloseq18_plants_glom_RA) %>% as.data.frame() %>%
    rownames_to_column(var="ESVId") %>%
    mutate(Colour = ifelse(ESVId %in% toKeepColours, color_leg_plants[CommonName], "grey")) %>%
    pull(Colour)
  names(color_leg_plants2) <- tax_table(phyloseq18_plants_glom_RA) %>% as.data.frame() %>% pull(CommonName)
  
  gg_plantsonly_wgrey <- plot_bar(phyloseq18_plants_glom_RA, fill="CommonName") +
    geom_bar(aes(col=CommonName, fill=CommonName), stat ="identity", position="stack") +
    facet_grid(~Captive.Wild.probiotic, scales="free", drop=TRUE) +
    scale_fill_manual(values=color_leg_plants2)+
    scale_color_manual(values=color_leg_plants2)+
    labs(col="Plant common name",fill="Plant common name") +
    ylab("Relative Abundance") +
    theme(axis.text.x = element_blank())
  gg_plantsonly_wgrey
  ggsave(filename = paste0("02_split_18_classes/gg_plantsonly_wgrey_",o,".png"), gg_plantsonly_wgrey, height=5, width=14)
  
  #### Look at microeukaryotes ####
  phyloseq18_broadcat 
  listMicroEuk <- tax_table(phyloseq18_broadcat) %>% as.data.frame() %>%
    filter(BroadCat2=="MicroEukaryote") %>%
    rownames()
  
  phyloseq18_microeuk <-prune_taxa(x=phyloseq18, taxa=listMicroEuk)
  phyloseq18_microeuk_glom_RA <- transform_sample_counts( tax_glom(phyloseq18_microeuk, taxrank = "Phylum"), fun=function(x) x/sum(x))
  
  set.seed(0923)
  colorlegmicroeuk <- sample(colors()[-grep("white|grey|gray|beige", colors())], size = length(taxa_names(phyloseq18_microeuk_glom_RA)))
  gg_microeuk_phylum <- plot_bar(phyloseq18_microeuk_glom_RA, fill="Phylum") +
    geom_bar(aes(col=Phylum, fill=Phylum), stat ="identity", position="stack") +
    facet_grid(~Captive.Wild.probiotic, scales="free", drop=TRUE) +
    scale_fill_manual(values=colorlegmicroeuk)+
    scale_color_manual(values=colorlegmicroeuk)+
    labs(col="Phylum\n(Microeukaryotes only)",fill="Phylum\n(Microeukaryotes only)") +
    ylab("Relative Abundance") +
    theme(axis.text.x = element_blank())
  gg_microeuk_phylum
  ggsave(paste0("02_split_18_classes/gg_microeuk_phylum_",o,".png"), gg_microeuk_phylum, height=5, width=8)
  
  # Get a no host version (also no mammal)
  nohost_ESVs <- allTaxWithCat %>% filter(!BroadCat2 %in% c("Presumed_host","Mammal")) %>% pull(ESVId)
  phyloseq18_nohost <- prune_taxa(taxa=nohost_ESVs, phyloseq18)
  
  ##### Save different versions of 18S stuff ####
  assign(paste0("phyloseq18",o,"_plants"),value = phyloseq18_plants)
  assign(paste0("phyloseq18",o,"_microeuk"),value = phyloseq18_microeuk)
  assign(paste0("phyloseq18",o,"_nohost"),value = phyloseq18_nohost)
  
  save(list=c(paste0("phyloseq18",o,"_plants")), file=paste0("02_split_18_classes/phyloseq18",o,"_plants.rds"))
  save(list=c(paste0("phyloseq18",o,"_microeuk")), file=paste0("02_split_18_classes/phyloseq18",o,"_microeuk.rds"))
  save(list=c(paste0("phyloseq18",o,"_nohost")), file=paste0("02_split_18_classes/phyloseq18",o,"_nohost.rds"))
  # 
  # if ( o == "nano") {
  #   phyloseq18n_plants <- phyloseq18_plants
  #   phyloseq18n_microeuk <- phyloseq18_microeuk
  #   phyloseq18n_nohost <- phyloseq18_nohost
  #   save(phyloseq18n_plants, file="02_split_18_classes/phyloseq18n_plants.rds")
  #   save(phyloseq18n_microeuk, file="02_split_18_classes/phyloseq18n_microeuk.rds")
  #   save(phyloseq18n_nohost, file="02_split_18_classes/phyloseq18n_nohost.rds")
  #   
  # } else if (o=="miseq") {
  #   phyloseq18m_plants <- phyloseq18_plants
  #   phyloseq18m_microeuk <- phyloseq18_microeuk
  #   phyloseq18m_nohost <- phyloseq18_nohost
  #   save(phyloseq18m_plants, file="02_split_18_classes/phyloseq18m_plants.rds")
  #   save(phyloseq18m_microeuk, file="02_split_18_classes/phyloseq18m_microeuk.rds")
  #   save(phyloseq18m_nohost, file="02_split_18_classes/phyloseq18m_nohost.rds")
  # }

}

