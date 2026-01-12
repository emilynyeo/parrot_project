#!bin/bash

#### Load libraries ####
library(vegan)
library(MASS)
library(phyloseq)
# library(indicspecies)
library(tibble)
library(tidyverse)
library(ggpubr)

setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")
if (!dir.exists("03_alpha_diversity")) {
  dir.create("03_alpha_diversity")
}
#### Load data ####
# Also converting to matrices with rownames for phyloseq

load("01_process_and_clean_data/phyloseq16n.rds")
load("01_process_and_clean_data/phyloseq16ren.rds")
# load("02_taxonomy_add/phyloseq18_detailtax.rds")
load("02_split_18_classes/phyloseq18nano_nohost.rds")
load("02_split_18_classes/phyloseq18nano_plants.rds")
load("02_split_18_classes/phyloseq18nano_microeuk.rds")

# Merge phyloseq16 and ren
phyloseq16nano <- merge_phyloseq(phyloseq16n, phyloseq16ren)

#### Colour legends ####
captivewild_legend <- c(Captive="purple4", "Wild, free ranging"="seagreen", "Wild, seized from traffickers"="goldenrod"
                        , "Crate 1: Unknown origin" = "magenta3", "Crate 2: Unknown origin" = "firebrick3")
location_legend <- c(Nigeria ="goldenrod4",`South Africa`="darkred", Unknown="grey")


####  Alpha diversity ####
##### 16S #####

# All samples
alpha_16 <- plot_richness(phyloseq16nano, measures = c("Chao1","Shannon","InvSimpson"))$data%>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) 
gg_alpha16_sourceprobiotic <- alpha_16 %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild))+
  geom_boxplot() + 
  geom_point(aes(group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Bacterial Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21))
gg_alpha16_sourceprobiotic
ggsave(filename = "03_alpha_diversity/gg_alpha16_sourceprobiotic.png",gg_alpha16_sourceprobiotic,  height=5, width=8)

# Just original ones
gg_alpha16n_sourceprobiotic <- alpha_16 %>%
  filter(Captive.Wild %in% c("Captive","Wild, free ranging","Wild, seized from traffickers")) %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild))+
  geom_boxplot() + 
  geom_point(aes(group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Bacterial Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21))
gg_alpha16n_sourceprobiotic
ggsave(filename = "03_alpha_diversity/gg_alpha16n_sourceprobiotic.png",gg_alpha16n_sourceprobiotic,  height=5, width=8)

alpha16n_stats <- compare_means(value_adj ~ Captive.Wild, data = alpha_16 %>%   filter(Captive.Wild %in% c("Captive","Wild, free ranging","Wild, seized from traffickers")) 
, p.adjust.method="fdr",
              method="t.test",
              group.by = "variable_adj") %>%
  mutate(p.signif.adj = ifelse(p.adj>=0.05, "ns", ifelse(p.adj>=0.01, "*", ifelse(p.adj>=0.001, "**", ifelse(p.adj>=0.001, "***", "****"))))) %>%
  mutate(y.position = c(350,400,200,5,5.5,4.5,1.75,2,1.5))
write.table(alpha16n_stats, file="03_alpha_diversity/alpha16_stats.txt", quote = FALSE, row.names = FALSE, sep="\t")
gg_alpha16n_sourceprobiotic_wstats <- alpha_16 %>%
  filter(Captive.Wild %in% c("Captive","Wild, free ranging","Wild, seized from traffickers")) %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, col=Captive.Wild))+
  geom_boxplot(aes(pch=Given.probiotic)) +
  stat_pvalue_manual(alpha16n_stats, label="p.signif.adj") +
  geom_point(aes(group=Given.probiotic, pch=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Bacterial Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21))
gg_alpha16n_sourceprobiotic_wstats
ggsave(filename = "03_alpha_diversity/gg_alpha16_sourceprobiotic_wstats.png",gg_alpha16_sourceprobiotic_wstats,  height=5, width=8)

# alpha16_comp <- list(c("Captive", "Wild, free ranging"), c("Wild, free ranging", "Wild, seized from traffickers"), c("Captive", "Wild, seized from traffickers"))
# gg_alpha16_sourceprobiotic_wstats <- ggboxplot(alpha_16, x="Captive.Wild", y="value_adj", pch="Given.probiotic", color = "Captive.Wild")+
#   stat_compare_means(comparisons =alpha16_comp
#                      , method = "wilcox.test"
#                      , p.adjust.method="fdr"
#                      , label = "p.signif"
#                      # , label = "p.adj"
#                      # , ref.group=".all."
#                      ) +
#   # stat_compare_means(method="anova") +
#   geom_point(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild, group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
#   facet_wrap(.~variable_adj, scales="free", nrow=1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
#   labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Bacterial Diversity") + xlab("Captive or wild") +
#   scale_color_manual(values=captivewild_legend) +
#   scale_shape_manual(values=c(19,21))
# gg_alpha16_sourceprobiotic_wstats


##### 18S #####
####### FIX SIMPSON, there are inf values #########
alpha_18_nohost <- plot_richness(phyloseq18nano_nohost, measures = c("Chao1","Shannon","InvSimpson"))$data %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  mutate(Type = "All eukaryotes")
alpha_18_plants <- plot_richness(phyloseq18nano_plants, measures = c("Chao1","Shannon","InvSimpson"))$data %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable)))%>%
  mutate(Type = "Plants only")
alpha_18_microeuk <- plot_richness(phyloseq18nano_microeuk, measures = c("Chao1","Shannon","InvSimpson"))$data %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable)))%>%
  mutate(Type = "Microeukaryotes only")
# 
# alpha_18 <- full_join(alpha_18_nohost, alpha_18_plants ) %>%
#   full_join(alpha_18_microeuk)

gg_alpha18_sourceprobiotic <- alpha_18_nohost %>%
  filter(Type == "All eukaryotes") %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild))+
  geom_boxplot() + 
  geom_point(aes(group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Eukaryote Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21))
gg_alpha18_sourceprobiotic
ggsave(filename = "03_alpha_diversity/gg_alpha18_sourceprobiotic.png", gg_alpha18_sourceprobiotic, height=5, width=8)


gg_alpha18n_sourceprobiotic <- alpha_18_nohost %>%
  filter(Captive.Wild %in% c("Captive","Wild, free ranging","Wild, seized from traffickers")) %>%
  filter(Type == "All eukaryotes") %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild))+
  geom_boxplot() + 
  geom_point(aes(group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Eukaryote Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21))
gg_alpha18n_sourceprobiotic
ggsave(filename = "03_alpha_diversity/gg_alpha18n_sourceprobiotic.png", gg_alpha18n_sourceprobiotic, height=5, width=8)

gg_alpha18_sourceprobiotic_plants <- alpha_18_plants %>%
  filter(Type == "Plants only") %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild))+
  geom_boxplot() + 
  geom_point(aes(group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Plant Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21))
gg_alpha18_sourceprobiotic_plants
ggsave(filename = "03_alpha_diversity/gg_alpha18_sourceprobiotic_plants.png", gg_alpha18_sourceprobiotic_plants, height=5, width=8)


gg_alpha18n_sourceprobiotic_plants <- alpha_18_plants %>%
  filter(Captive.Wild %in% c("Captive","Wild, free ranging","Wild, seized from traffickers")) %>%
  filter(Type == "Plants only") %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild))+
  geom_boxplot() + 
  geom_point(aes(group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Plant Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21))
gg_alpha18n_sourceprobiotic_plants
ggsave(filename = "03_alpha_diversity/gg_alpha18n_sourceprobiotic_plants.png", gg_alpha18n_sourceprobiotic_plants, height=5, width=8)

gg_alpha18_sourceprobiotic_microeuk <- alpha_18_microeuk %>%
  filter(Type == "Microeukaryotes only") %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild))+
  geom_boxplot() + 
  geom_point(aes(group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Microeukaryote Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21))
gg_alpha18_sourceprobiotic_microeuk
ggsave(filename = "03_alpha_diversity/gg_alpha18_sourceprobiotic_microeuk.png", gg_alpha18_sourceprobiotic_microeuk, height=5, width=8)

gg_alpha18n_sourceprobiotic_microeuk <- alpha_18_microeuk %>%
  filter(Captive.Wild %in% c("Captive","Wild, free ranging","Wild, seized from traffickers")) %>%
  filter(Type == "Microeukaryotes only") %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild))+
  geom_boxplot() + 
  geom_point(aes(group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Microeukaryote Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21))
gg_alpha18n_sourceprobiotic_microeuk
ggsave(filename = "03_alpha_diversity/gg_alpha18n_sourceprobiotic_microeuk.png", gg_alpha18n_sourceprobiotic_microeuk, height=5, width=8)


##### FIX STATS lov inv simpson has infinity ####
# 
# alpha18_stats_nohost <- compare_means(value_adj ~ Captive.Wild, data = alpha_18_nohost %>% filter(Captive.Wild %in% c("Captive","Wild, free ranging","Wild, seized from traffickers")) , p.adjust.method="fdr",
#                                method="t.test",
#                                group.by = "variable_adj") %>%
#   mutate(p.signif.adj = ifelse(p.adj>=0.05, "ns", ifelse(p.adj>=0.01, "*", ifelse(p.adj>=0.001, "**", ifelse(p.adj>=0.001, "***", "****"))))) %>%
#   mutate(y.position = c(225,250,200,5,5.5,4.5,1.75,2,1.5))
# write.table(alpha18_stats_nohost, file="03_alpha_diversity/alpha18_stats_nohost.txt", quote = FALSE, row.names = FALSE, sep="\t")
# gg_alpha18_sourceprobiotic_wstats <- alpha_18_nohost %>%
#   ggplot(aes(x=Captive.Wild, y=value_adj, col=Captive.Wild))+
#   geom_boxplot(aes(pch=Given.probiotic)) +
#   stat_pvalue_manual(alpha18_stats_nohost, label="p.signif.adj") +
#   geom_point(aes(group=Given.probiotic, pch=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
#   facet_wrap(.~variable_adj, scales="free", nrow=1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
#   labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Eukaryote Diversity") + xlab("Captive or wild") +
#   scale_color_manual(values=captivewild_legend) +
#   scale_shape_manual(values=c(19,21))
# gg_alpha18_sourceprobiotic_wstats
# ggsave(filename = "03_alpha_diversity/gg_alpha18_sourceprobiotic_wstats.png",gg_alpha18_sourceprobiotic_wstats,  height=5, width=8)
# 
# alpha_18_microeuk <- alpha_18_microeuk %>% filter(!is.na(value_adj), is.finite(value_adj))
# alpha18_stats_microeuk <- compare_means(value_adj ~ Captive.Wild, data = alpha_18_microeuk, p.adjust.method="fdr",
#                                       method="t.test",
#                                       group.by = "variable_adj") %>%
#   mutate(p.signif.adj = ifelse(p.adj>=0.05, "ns", ifelse(p.adj>=0.01, "*", ifelse(p.adj>=0.001, "**", ifelse(p.adj>=0.001, "***", "****"))))) %>%
#   mutate(y.position = c(225,250,200,5,5.5,4.5,1.75,2,1.5))
# write.table(alpha18_stats_microeuk, file="03_alpha_diversity/alpha18_stats_microeuk.txt", quote = FALSE, row.names = FALSE, sep="\t")
# gg_alpha18_sourceprobiotic_microeuk_wstats <- alpha_18_microeuk %>%
#   ggplot(aes(x=Captive.Wild, y=value_adj, col=Captive.Wild))+
#   geom_boxplot(aes(pch=Given.probiotic)) +
#   stat_pvalue_manual(alpha18_stats_microeuk, label="p.signif.adj") +
#   geom_point(aes(group=Given.probiotic, pch=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
#   facet_wrap(.~variable_adj, scales="free", nrow=1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
#   labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Microeukaryote Diversity") + xlab("Captive or wild") +
#   scale_color_manual(values=captivewild_legend) +
#   scale_shape_manual(values=c(19,21))
# gg_alpha18_sourceprobiotic_microeuk_wstats
# ggsave(filename = "03_alpha_diversity/gg_alpha18_sourceprobiotic_microeuk_wstats.png",gg_alpha18_sourceprobiotic_microeuk_wstats,  height=5, width=8)
# 
# alpha18_comp <- list(c("Captive", "Wild, free ranging"), c("Wild, free ranging", "Wild, seized from traffickers"), c("Captive", "Wild, seized from traffickers"))
# gg_alpha18_sourceprobiotic_wstats <- ggboxplot(alpha_18, x="Captive.Wild", y="value_adj", pch="Given.probiotic", color = "Captive.Wild")+
#   stat_compare_means(comparisons =alpha18_comp
#                      , method = "wilcox.test"
#                      , p.adjust.method="fdr"
#                      , label = "p.signif"
#                      # , label = "p.adj"
#                      # , ref.group=".all."
#   ) +
#   # stat_compare_means(method="anova") +
#   geom_point(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild, group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
#   facet_wrap(.~variable_adj, scales="free", nrow=1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
#   labs(col="Captive\nor wild", shape="Microbiome\ntreatment") +ylab("Bacterial Diversity") + xlab("Captive or wild") +
#   scale_color_manual(values=captivewild_legend) +
#   scale_shape_manual(values=c(19,21))
# gg_alpha18_sourceprobiotic_wstats


#### Is read depth correlated with diversity ####
summary(lm(value ~ ReadDepth16n, data=alpha_16 %>% filter(variable == "Chao1")))
summary(lm(value ~ ReadDepth16n, data=alpha_16 %>% filter(variable == "Shannon")))
summary(lm(value ~ ReadDepth16n, data=alpha_16 %>% filter(variable == "InvSimpson")))

gg_readdepth_div_16 <- alpha_16 %>%
  ggplot(aes(x=ReadDepth16m, y=value)) +
  geom_point() + facet_grid(variable~., scales="free") + 
  geom_smooth(method="lm") +
  ylab("16S diversity or eveness") + xlab("Read depth") +
  scale_x_log10()
gg_readdepth_div_16
ggsave(filename="03_alpha_diversity/gg_readdepth_div_16.png", gg_readdepth_div_16, height=8, width=5)

summary(lm(value ~ ReadDepth18n, data=alpha_18 %>%
             filter(Type=="All eukaryotes") %>% filter(variable == "Chao1")))
summary(lm(value ~ ReadDepth18n, data=alpha_18 %>%
             filter(Type=="All eukaryotes") %>% filter(variable == "Shannon")))
summary(lm(value ~ ReadDepth18n, data=alpha_18 %>%
             filter(Type=="All eukaryotes") %>% filter(variable == "InvSimpson")))

gg_readdepth_div_18 <- alpha_18 %>%
  filter(Type=="All eukaryotes")  %>%
  ggplot(aes(x=ReadDepth18n, y=value)) +
  geom_point() + facet_grid(variable~., scales="free")+ 
  geom_smooth(method="lm") +
  ylab("18S diversity or eveness") + xlab("Read depth") +
  scale_x_log10()
gg_readdepth_div_18
ggsave(filename="03_alpha_diversity/gg_readdepth_div_18.png", gg_readdepth_div_18, height=8, width=5)

##### Try rarefying #####
phyloseq16_r17000 <- rarefy_even_depth(phyloseq16nano, sample.size = 17000, rngseed = 42345, trimOTUs = TRUE)

alpha_16_r17000 <- plot_richness(phyloseq16_r17000, measures = c("Chao1","Shannon","InvSimpson"))$data
gg_alpha16_r17000_sourceprobiotic <- alpha_16_r17000 %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild))+
  geom_boxplot() + 
  geom_point(aes(group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Bacterial Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21)) +
  labs(pch="Microbiome\ntreatment", col="Captive\nor wild")
gg_alpha16_r17000_sourceprobiotic
ggsave(filename = "03_alpha_diversity/gg_alpha16_r17000_sourceprobiotic.png", gg_alpha16_r17000_sourceprobiotic, height=5, width=8)

### Doesn't look any different from non-rarefied, 18S threshold is sort of arbitrary here because values go super low
hist(colSums(otu_table(phyloseq18nano_nohost)))
phyloseq18_nohost_r1000 <- rarefy_even_depth(phyloseq18nano_nohost, sample.size = 1000, rngseed = 42345, trimOTUs = TRUE)

alpha_18_r1000 <- plot_richness(phyloseq18nano_nohost, measures = c("Chao1","Shannon","InvSimpson"))$data
gg_alpha18_r1000_sourceprobiotic <- alpha_18_r1000 %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  ggplot(aes(x=Captive.Wild, y=value_adj, pch=Given.probiotic, col=Captive.Wild))+
  geom_boxplot() + 
  geom_point(aes(group=Given.probiotic),position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Eukaryote Diversity") + xlab("Captive or wild") +
  scale_color_manual(values=captivewild_legend) +
  scale_shape_manual(values=c(19,21)) +
  labs(pch="Microbiome\ntreatment", col="Captive\nor wild")
gg_alpha18_r1000_sourceprobiotic
ggsave(filename = "03_alpha_diversity/gg_alpha18_r1000_sourceprobiotic.png", gg_alpha18_r1000_sourceprobiotic, height=5, width=8)


######### FIX COLOURS FOR UNKNOWN TUBES ############
#### Tube preservative? location?
gg_alpha16_tubelocation <- alpha_16 %>%
  mutate(Country = ifelse(is.na(Country), "Unknown", Country)) %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  ggplot(aes(x=Tube.type, y=value_adj, col=Country))+
  geom_boxplot() + 
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Collection\nlocation") +ylab("Bacterial Diversity") + xlab("DNA preservation tube type") +
  scale_color_manual(values=location_legend)
gg_alpha16_tubelocation
ggsave("03_alpha_diversity/gg_alpha16_tubelocation.png", gg_alpha16_tubelocation, height=5, width=8)

tube16_stats <- compare_means(value_adj ~ Tube.type, data = alpha_16,
                              # p.adjust.method="fdr",
                               method="t.test",
                               group.by = "variable_adj") %>%
  mutate(p.signif.adj = ifelse(p.adj>=0.05, "ns", ifelse(p.adj>=0.01, "*", ifelse(p.adj>=0.001, "**", ifelse(p.adj>=0.001, "***", "****"))))) %>%
  mutate(y.position = c(250,4.5,1.75))
write.table(tube16_stats, file="03_alpha_diversity/tube16_stats.txt", quote = FALSE, row.names = FALSE, sep="\t")
# tube16_comp <- list(c("5ml RNAlater factory", "SANBI 10ml conical"))
gg_alpha16_tubelocation_wstats <- alpha_16 %>%
  mutate(Country = ifelse(is.na(Country), "Unknown", Country)) %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  ggplot(aes(x=Tube.type, y=value_adj))+
  geom_boxplot(aes(col=Country)) + 
  stat_pvalue_manual(data=tube16_stats, label="p.signif.adj") +
  geom_point(aes( col=Country), position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Collection\nlocation") +ylab("Bacterial Diversity") + xlab("DNA preservation tube type") +
  scale_color_manual(values=location_legend)
gg_alpha16_tubelocation_wstats
ggsave(filename = "03_alpha_diversity/gg_alpha16_tubelocation_wstats.png",gg_alpha16_tubelocation_wstats,  height=5, width=8)


gg_alpha18_tubelocation_nohost <- alpha_18_nohost %>%
  mutate(Country = ifelse(is.na(Country), "Unknown", Country)) %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  ggplot(aes(x=Tube.type, y=value_adj, col=Country))+
  geom_boxplot() + 
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Collection\nlocation") +ylab("Eukaryote Diversity") + xlab("DNA preservation tube type") +
  scale_color_manual(values=location_legend)
gg_alpha18_tubelocation_nohost
ggsave("03_alpha_diversity/gg_alpha18_tubelocation_nohost.png", gg_alpha18_tubelocation_nohost, height=5, width=8)

gg_alpha18_tubelocation_microeuk <- alpha_18_microeuk %>%
  mutate(Country = ifelse(is.na(Country), "Unknown", Country)) %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  ggplot(aes(x=Tube.type, y=value_adj, col=Country))+
  geom_boxplot() + 
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Collection\nlocation") +ylab("Microeukaryote Diversity") + xlab("DNA preservation tube type") +
  scale_color_manual(values=location_legend)
gg_alpha18_tubelocation_microeuk
ggsave("03_alpha_diversity/gg_alpha18_tubelocation_microeuk.png", gg_alpha18_tubelocation_microeuk, height=5, width=8)


gg_alpha18_tubelocation_plants <- alpha_18_plants %>%
  mutate(Country = ifelse(is.na(Country), "Unknown", Country)) %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  ggplot(aes(x=Tube.type, y=value_adj, col=Country))+
  geom_boxplot() + 
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Collection\nlocation") +ylab("Plant Diversity") + xlab("DNA preservation tube type") +
  scale_color_manual(values=location_legend)
gg_alpha18_tubelocation_plants
ggsave("03_alpha_diversity/gg_alpha18_tubelocation_plants.png", gg_alpha18_tubelocation_plants, height=5, width=8)

######## STATS BELOW ARE BROKEN ########
tube18_nohost_stats <- compare_means(value_adj ~ Tube.type, data = alpha_18_nohost %>% filter(Country=="South Africa"),
                              # p.adjust.method="fdr",
                              method="t.test",
                              group.by = "variable_adj") %>%
  mutate(p.signif.adj = ifelse(p.adj>=0.05, "ns", ifelse(p.adj>=0.01, "*", ifelse(p.adj>=0.001, "**", ifelse(p.adj>=0.001, "***", "****"))))) %>%
  mutate(y.position = c(250,4,1.5))
write.table(tube18_nohost_stats, file="03_alpha_diversity/tube18_nohost_stats.txt", quote = FALSE, row.names = FALSE, sep="\t")
# tube18_comp <- list(c("5ml RNAlater factory", "SANBI 10ml conical"))
gg_alpha18_tubelocation_nohost_wstats <- alpha_18_nohost %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  ggplot(aes(x=Tube.type, y=value_adj))+
  geom_boxplot(aes(col=Country)) + 
  stat_pvalue_manual(data=tube18_nohost_stats, label="p.signif.adj") +
  geom_point(aes( col=Country), position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Collection\nlocation") +ylab("Bacterial Diversity") + xlab("DNA preservation tube type") +
  scale_color_manual(values=location_legend)
gg_alpha18_tubelocation_nohost_wstats
ggsave(filename = "03_alpha_diversity/gg_alpha18_tubelocation_nohost_wstats.png",gg_alpha18_tubelocation_nohost_wstats,  height=5, width=8)


gg_alpha16_datesource <- alpha_16 %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  ggplot(aes(x=as.Date(DateAdj), y=value_adj, col=Captive.Wild))+
  geom_jitter(height=0, width=0.1)+
  # geom_boxplot() + 
  # geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild") +ylab("Bacterial Diversity") + xlab("Date sampled") +
  scale_color_manual(values=captivewild_legend)
gg_alpha16_datesource
ggsave("03_alpha_diversity/gg_alpha16_datesource.png",gg_alpha16_datesource,  height=5, width=8)

gg_alpha18_datesource <- alpha_18_nohost %>%
  mutate(value_adj = ifelse(variable=="InvSimpson", log10(value), as.numeric(value))) %>%
  mutate(variable_adj = ifelse(variable=="InvSimpson", "Log InvSimpson", as.character(variable))) %>%
  ggplot(aes(x=as.Date(DateAdj), y=value_adj, col=Captive.Wild))+
  geom_jitter(height=0, width=0.1)+
  # geom_boxplot() + 
  # geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width = 0.75)) +
  facet_wrap(.~variable_adj, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  labs(col="Captive\nor wild") +ylab("Eukaryote Diversity") + xlab("Date sampled") +
  scale_color_manual(values=captivewild_legend)
gg_alpha18_datesource
ggsave("03_alpha_diversity/gg_alpha18_datesource.png",gg_alpha18_datesource,  height=5, width=8)

