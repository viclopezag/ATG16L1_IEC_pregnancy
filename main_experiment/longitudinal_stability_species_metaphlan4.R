
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(glue)


# Code for - Longitudinal changing species: shotgun metagenomics-------

# Setting Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loading Longitudinal per-feature test results from Maaslin2
path = "outputs/differential_abundance_species_metaphlan4/lmm_maaslin2"

complete_model <- read_tsv(paste0(path,"/species/alldata/all_results.tsv")) %>%
  mutate(group = "complete_model")

long_WT_w3 <- read_tsv(paste0(path,"/species/WT_long/all_results.tsv"))%>%
  filter(value == "w3") %>%
  mutate(group = "longWT_w3")

long_WT_w6 <- read_tsv(paste0(path,"/species/WT_long/all_results.tsv"))%>%
  filter(value == "w6") %>%
  mutate(group = "longWT_w6")

long_KO_w3 <- read_tsv(paste0(path,"/species/KO_long/all_results.tsv"))%>%
  filter(value == "w3") %>%
  mutate(group = "longKO_w3")

long_KO_w6 <- read_tsv(paste0(path,"/species/KO_long/all_results.tsv"))%>%
  filter(value == "w6") %>%
  mutate(group = "longKO_w6")


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Changing of Taxa in WT #

results_WT <- rbind(long_WT_w3,long_WT_w6)

fdr_filtered_features_WT <- rbind(long_WT_w3,long_WT_w6) %>%
  mutate(dif_abundant = ifelse(coef > 1 & qval <= 0.05,"up",
                               ifelse(coef < -1 & qval <= 0.05,"down","stable"))) %>%
  mutate(label = ifelse(dif_abundant == "up",feature,
                        ifelse(dif_abundant == "down",feature,NA))) %>%
  mutate(label = str_remove(label,"s__"))

write.table(fdr_filtered_features_WT,file="outputs/differential_abundance_species_metaphlan4/lmm_maaslin2/species/WT_long/long_WT_fdr_log2_stability.txt",row.names=F,sep="\t")

# Volcano Plot #
myvolcanoplot <- ggplot(data = fdr_filtered_features_WT, 
                        aes(x = coef, y = -log10(qval), 
                            col = dif_abundant,shape = value, label = label)) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed')  +
  geom_point(size = 1.0) +
  scale_color_manual(values = c("down"="darkblue", "stable"="grey", "up"="darkred"), # to set the colours of our variable
                     labels = c("Decrease", "Stable", "Increase")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  scale_shape_manual(values = c("w3"=17,"w6"=18)) +
  coord_cartesian(ylim = c(0, 10), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Stability', shape = "pregnancy/post-pregnancy", #legend_title
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) +
  #scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis<br />
  ggtitle('Stability of microbiota - WT mice') + # Plot title<br />
  geom_text_repel(size = 1.0,
                  box.padding = unit(0.5, "lines"),
                  max.overlaps = Inf) + # To show all labels
  theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=7),
        legend.title = element_text(size=0), plot.title = element_text(size=9, face="bold", hjust = 0.5))


plot(myvolcanoplot)

ggsave("supp_figures/Supp_Fig2C_changing_of_species_long_WT.pdf", width= 12, height=10, units="cm")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Changing of Taxa in KO #
results_KO <- rbind(long_KO_w3,long_KO_w6)

fdr_filtered_features_KO <- rbind(long_KO_w3,long_KO_w6) %>%
  mutate(dif_abundant = ifelse(coef > 1 & qval <= 0.05,"up",
                               ifelse(coef < -1 & qval <= 0.05,"down","stable"))) %>%
  mutate(label = ifelse(dif_abundant == "up",feature,
                        ifelse(dif_abundant == "down",feature,NA))) %>%
  mutate(label = str_remove(label,"s__"))

write.table(fdr_filtered_features_KO,file="outputs/differential_abundance_species_metaphlan4/lmm_maaslin2/species/KO_long/long_KO_fdr_log2_stability.txt",row.names=F,sep="\t")


# Volcano Plot #
myvolcanoplot <- ggplot(data = fdr_filtered_features_KO, 
                        aes(x = coef, y = -log10(qval), 
                            col = dif_abundant,shape = value, label = label)) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed')  +
  geom_point(size = 1) +
  scale_color_manual(values = c("down"="darkblue", "stable"="grey", "up"="darkred"), # to set the colours of our variable
                     labels = c("Decrease", "Stable", "Increase")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  scale_shape_manual(values = c("w3"=17,"w6"=18)) +
  coord_cartesian(ylim = c(0, 10), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Stability', shape = "pregnancy/post-pregnancy", #legend_title
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) +
  ggtitle('Stability of microbiota - KO mice') + # Plot title<br />
  geom_text_repel(size = 1.0,
                  box.padding = unit(0.5, "lines"),
                  max.overlaps = Inf) + # To show all labels
  theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=7),
        legend.title = element_text(size=0), plot.title = element_text(size=9, face="bold", hjust = 0.5))


plot(myvolcanoplot)

ggsave("supp_figures/Supp_Fig2D_changing_of_species_long_KO.pdf", width= 12, height=10, units="cm")


#------------------------------------------------------------------------------
# Longitudinal Plots #

# metadata  
preg_metadata <- read.table("metaInfo/metainformation_16S_and_shotgun_metagenomics.txt", check.names = F, header = T, sep = "\t") %>% 
  filter(!is.na(barcode_shotgun)) %>%
  column_to_rownames(var = "barcode_shotgun") %>%
  mutate(time_num = ifelse(timepoint == "BL",0,
                           ifelse(timepoint == "w3",3,6)))

# Relative Abundance Data #

metaphlan <- read.table("inputs/metaphlan_rel_abundances.txt", sep = "\t", check.names = F, header = T) %>%
  separate(col = clade_name, into = c("kingdom","phylum","class","order","family","genus","species","SGBs"), sep = '[|]')

species <- metaphlan %>%
  filter(is.na(SGBs) & !is.na(species)) %>%
  select(-kingdom,-phylum,-class,-order,-family,-genus,-SGBs,-clade_taxid) %>%
  column_to_rownames(var="species")/100

sample_ids <- substr(colnames(species), 1, 16) # extracting only the first 6 strings.

# change colnames
colnames(species) <- sample_ids

s_met_and_taxa <- merge(preg_metadata, t(species), by = 0) %>%
  mutate(Condition = paste0(genotype," ",timepoint))
colnames(s_met_and_taxa)[1] <- "sample_id"

# Loading Significant Unstable taxa #

unstable_taxa_to_plot <- c("s__Lactobacillus_johnsonii","s__Limosilactobacillus_reuteri",
                           "s__Ligilactobacillus_murinus","s__Muribaculum_gordoncarteri",
                           "s__Muribaculaceae_bacterium_Isolate_110_HZI")

long_s_met_and_taxa <- s_met_and_taxa %>%
  select(time_num,genotype,unstable_taxa_to_plot) %>%
  pivot_longer(., cols = starts_with("s__"), names_to = "species", values_to = "RAs")


dir.create("supp_figures/Supp_Fig2E")

for (i in 1:length(unstable_taxa_to_plot)){
  
  # Longitudinal Plot #
  
  ggline(long_s_met_and_taxa %>% filter(species == glue("{unstable_taxa_to_plot[i]}")),x = "time_num", y = "RAs", add = c("mean_se","jitter"), color = "genotype",
         add.params = list(size=1,alpha=0.4)) +
    stat_compare_means(aes(group = genotype), label = "p.signif", method = "wilcox.test", size = 2) +
    stat_summary(mapping = aes(fill=genotype), fun = "mean", geom = "point", color = "black", pch = 21, size = 2)+
    ylab("Relative Abundances") + xlab("Weeks") +
    scale_fill_manual(values = c("WT"="#5da1d4ff", "KO"="#c88cb2ff")) +
    scale_color_manual(values = c("WT"="#5da1d4ff", "KO"="#c88cb2ff")) +
    theme_bw() + theme(axis.text=element_text(size=7, color = "black"), 
                       axis.title=element_text(size=7), 
                       legend.text = element_text(size=7), 
                       legend.title = element_text(size=0), 
                       plot.title = element_text(size=9, face="bold", hjust = 0.5)) + 
    ggtitle(glue("{unstable_taxa_to_plot[i]}"))

  ggsave(file = glue("supp_figures/Supp_Fig2E/{unstable_taxa_to_plot[i]}.pdf"), width = 6, height = 5, units="cm")
  
}
