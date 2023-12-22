
library(tidyverse)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(ggrepel)
library(vegan)
library(microbiome)
library(ggpubr)
library(dplyr)
library(Maaslin2)
library(glue)

# Code for - Longitudinal changing genus: 16S data-------

# Setting Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ps <- readRDS("inputs/phyloseq_object.rds")

# Cutoff #

(ps_fr <- subset_samples(ps, sample_sums(ps) > 5000))

# Deleting taxa with relative abundance lowest than 1e-5
minTotRelAbun <- 1e-5          
x <- taxa_sums(ps_fr)
keepTaxa <- (x / sum(x)) > minTotRelAbun
ps_f <- prune_taxa(keepTaxa, ps_fr)

# Agglomarating taxa 
ps_genus <- tax_glom(ps_f, taxrank = "Genus", NArm = T)

# Creating a curated name of genus #
taxa_annotation <- data.frame(ps_genus@tax_table) %>%
  cbind(data.frame(ASVs = rownames(.)),.) %>%
  mutate(tax_name = paste0(ASVs," | ",Family," | ",Genus)) %>%
  select(-Species)

#:::::::::::: Linear mixed model with maaslin2 :::::::::::#
# Creating directory 
dir.create("outputs/longitudinal_changes_genus")
dir.create("outputs/longitudinal_changes_genus/lmm_maaslin2")
dir.create("outputs/longitudinal_changes_genus/lmm_maaslin2/genus")
dir.create("outputs/longitudinal_changes_genus/lmm_maaslin2/genus/complete_model")

# longitudinal comparison - complete model (WT and KO)
m2_g <- Maaslin2(
    input_data = data.frame(ps_genus@otu_table),
    input_metadata = data.frame(ps_genus@sam_data),
    output="outputs/longitudinal_changes_genus/lmm_maaslin2/genus/complete_model",
    normalization = "CLR",
    transform = "NONE",
    analysis_method = "LM",
    fixed_effects = c("genotype","timepoint"),
    random_effects = c("mouse_id"),
    reference = c("genotype,WT","timepoint,BL"),
    min_abundance = 1e-8,
    min_prevalence = 0.12, # 80% modified rule (min(samples_in_time)/total_samples)*80
    correction = "BH",
    max_significance = 0.05,
    plot_heatmap = TRUE,
    heatmap_first_n = 50,
    plot_scatter = TRUE
)

# longitudinal comparison - Only WT

m2_g_WT <- Maaslin2(
  input_data = data.frame(subset_samples(ps_genus,genotype == "WT")@otu_table),
  input_metadata = data.frame(subset_samples(ps_genus,genotype == "WT")@sam_data),
  output="outputs/longitudinal_changes_genus/lmm_maaslin2/genus/longitudinal_WT",
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("timepoint"),
  random_effects = c("mouse_id"),
  reference = c("timepoint,BL"),
  min_abundance = 1e-8,
  min_prevalence = 0.12, # 80% modified rule (min(samples_in_time)/total_samples)*80
  correction = "BH",
  max_significance = 0.05,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)

#::::::::::::::::::::::::::::::::::::::::::::#
# Longitudinal Stability  - WT Pregnant mice #
# Stability compared to BL #

results_WT <- m2_g_WT$results %>%
  left_join(.,taxa_annotation,by = c("feature"="ASVs"))

# Filtering significant abundant taxa #
fdr_filtered_features_WT <- results_WT %>%
  #filter(value != "WT") %>%
  subset(., Genus != "unknown") %>%
  mutate(dif_abundant = ifelse(coef > 1 & qval <= 0.05,"up",
                               ifelse(coef < -1 & qval <= 0.05,"down","stable"))) %>%
  mutate(label = ifelse(dif_abundant == "up",Genus,
                        ifelse(dif_abundant == "down",Genus,NA))) %>%
  mutate(group = "WT")

write.table(fdr_filtered_features_WT,file="outputs/longitudinal_changes_genus/lmm_maaslin2/genus/longitudinal_WT/long_WT_fdr_log2_stability.txt",row.names=F,sep="\t")

# Volcano Plot #
myvolcanoplot <- ggplot(data = fdr_filtered_features_WT, 
                        aes(x = coef, y = -log10(qval), 
                            col = dif_abundant,shape = value, label = label)) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed')  +
  geom_point(size = 1) +
  scale_color_manual(values = c("down"="darkblue", "stable"="darkgrey", "up"="darkred"), # to set the colours of our variable
                     labels = c("Decrease", "Stable", "Increase")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  scale_shape_manual(values = c("w3"=17,"w6"=18)) +
  coord_cartesian(ylim = c(0, 7), xlim = c(-5, 5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Stability', shape = "pregnancy/post-pregnancy", #legend_title
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) +
  ggtitle('Stability of microbiota - WT mice') + # Plot title<br />
  geom_text_repel(size = 1.0,
                  max.overlaps = Inf) + # To show all labels
  theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), 
        axis.title=element_text(size=9), 
        legend.text = element_text(size=7),
        legend.title = element_text(size=0), 
        plot.title = element_text(size=9, face="bold", hjust = 0.5))


plot(myvolcanoplot)

ggsave("supp_figures/Supp_Fig2A_changing_of_genus_long_WT.pdf", width= 12, height=10, units="cm")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# longitudinal comparison - Only KO

m2_g_KO <- Maaslin2(
  input_data = data.frame(subset_samples(ps_genus,genotype == "KO")@otu_table),
  input_metadata = data.frame(subset_samples(ps_genus,genotype == "KO")@sam_data),
  output="outputs/longitudinal_changes_genus/lmm_maaslin2/genus/longitudinal_KO",
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("timepoint"),
  random_effects = c("mouse_id"),
  reference = c("timepoint,BL"),
  min_abundance = 1e-8,
  min_prevalence = 0.12, # 80% modified rule (min(samples_in_time)/total_samples)*80
  correction = "BH",
  max_significance = 0.05,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)

#::::::::::::::::::::::::::::::::::::::::::::#
# Longitudinal Stability  - KO Pregnant mice #
# Stability compared to BL #

results_KO <- m2_g_KO$results %>%
  left_join(.,taxa_annotation,by = c("feature"="ASVs"))

# Filtering significant abundant taxa #
fdr_filtered_features_KO <- results_KO %>%
  subset(., Genus != "unknown") %>%
  mutate(dif_abundant = ifelse(coef > 1 & qval <= 0.05,"up",
                               ifelse(coef < -1 & qval <= 0.05,"down","stable"))) %>%
  mutate(label = ifelse(dif_abundant == "up",Genus,
                        ifelse(dif_abundant == "down",Genus,NA))) %>%
  mutate(group = "KO")

write.table(fdr_filtered_features_KO,file="outputs/longitudinal_changes_genus/lmm_maaslin2/genus/longitudinal_KO/long_KO_fdr_log2_stability.txt",row.names=F,sep="\t")

# Volcano Plot #
myvolcanoplot <- ggplot(data = fdr_filtered_features_KO, 
                        aes(x = coef, y = -log10(qval), 
                            col = dif_abundant,shape = value, label = label)) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed')  +
  geom_point(size = 1) +
  scale_color_manual(values = c("down"="darkblue", "stable"="darkgrey", "up"="darkred"), # to set the colours of our variable
                     labels = c("Decrease", "Stable", "Increase")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  scale_shape_manual(values = c("w3"=17,"w6"=18)) +
  coord_cartesian(ylim = c(0, 12), xlim = c(-5,5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Stability', shape = "pregnancy/post-pregnancy", #legend_title
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) +
  #scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis<br />
  ggtitle('Stability of microbiota - KO mice') + # Plot title<br />
  geom_text_repel(size = 1.0,
                  max.overlaps = Inf) + # To show all labels
  theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=7),
        legend.title = element_text(size=0), plot.title = element_text(size=9, face="bold", hjust = 0.5))


plot(myvolcanoplot)

ggsave("supp_figures/Supp_Fig2B_changing_of_genus_long_KO.pdf", width= 12, height=10, units="cm")

#::::::::::::::::::::::::::::::::::::::::::::::#
# Plotting Triangular Dotplot #

# Only including FDR < 0.05 taxa with FC > 2 or FC < -2

combined_results <- rbind(results_WT %>% mutate(group = "WT"),results_KO %>% mutate(group = "KO")) %>%
  subset(., Genus != "unknown") %>%
  mutate(dif_abundant = ifelse(coef > 1 & qval <= 0.05,"up",
                               ifelse(coef < -1 & qval <= 0.05,"down","stable"))) %>%
  mutate(label = ifelse(dif_abundant == "up",Genus,
                        ifelse(dif_abundant == "down",Genus,NA))) %>%
  mutate(group_comparison = paste0(group,"-",value))
  

combined_results$group_comparison <- factor(combined_results$group_comparison, levels = c("WT-w3","KO-w3","WT-w6","KO-w6"))

combined_results %>%
  filter(!is.na(label)) %>%
  ggplot (aes (x = group_comparison, y =reorder(Genus,coef))) + 
  geom_point (aes (fill = coef, shape = as.factor(sign (coef)), size = -log10(qval)), 
              color = "#000000", alpha = 1) + 
  geom_point (aes (fill = coef, shape = as.factor(sign (coef)), size = -log10(qval)), 
              color = "#777777", alpha = 0.5) + 
  theme_bw() + theme(legend.position="right") +
  theme (axis.text.x = element_text (angle = 45, hjust = 1, size = 7, face = "italic"),
         axis.text.y = element_text (size = 9)) + scale_shape_manual (values = c (25, 24)) +
  scale_fill_gradient2 (mid = "white", high = "darkred", low = "darkblue", 
                        midpoint = 0, guide = guide_colorbar (CLRter = F), breaks = c(-1,1,2,-2,3,-3)) + 
  xlab ("") + ylab ("") + labs(fill = "Log2FC", shape = "Effect sense", size = "-log10(FDR)") +
  ggtitle('') +
  theme(axis.text=element_text(size=7, color = "black"),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title=element_text(size=7), 
        legend.text = element_text(size=7), 
        legend.title = element_text(size=7), 
        plot.title = element_text(size=7, face="bold", hjust = 0.5))

ggsave("figures/Fig1E_changing_genus_triangle_plot.pdf", width= 10, height=12, units="cm")  

##########################################

# Ploting high unstable taxa #

# Testing the differences of top genera #

genus_CLR_ps <- microbiome::transform(ps_genus, transform = "clr")
genus_CLR <- data.frame(genus_CLR_ps@otu_table)

genus_CLR_table <- cbind(data.frame(Genus = taxa_annotation[,"Genus"]),t(genus_CLR)) %>%
  pivot_longer(.,cols=2:ncol(.),names_to ="samples",values_to="CLR") %>%
  group_by(samples,Genus) %>%
  summarize(CLR = sum(CLR)) %>% # summing up all the "unknown" Genus
  filter(Genus != "unknown") # Deleting unknown Genus taxa.

genus_CLR_matrix <- genus_CLR_table %>%
  pivot_wider(.,names_from = Genus, values_from = CLR) %>%
  column_to_rownames(var="samples")

# Unique Unstable Taxa #

unstable_taxa <- c("Enterorhabdus","Lachnospiraceae UCG-001","Lachnospiraceae FCS020 group",
                   "Turicibacter","Parabacteroides","Prevotellaceae UCG-001")

genus_RA_to_plot <- genus_CLR_matrix %>%
  cbind(data.frame(genus_CLR_ps@sam_data),.) %>%
  select(genotype,timepoint,unstable_taxa) %>%
  pivot_longer(cols = Enterorhabdus:ncol(.), names_to = "taxa", values_to="CLR") %>%
  mutate(Condition = paste0(genotype,"-",timepoint),
         time_num = ifelse(timepoint == "BL",0,
                           ifelse(timepoint == "w3",3,6)))



genus_RA_to_plot$Condition <- factor(genus_RA_to_plot$Condition, levels= c("WT-BL","KO-BL",
                                                                           "WT-w3","KO-w3",                                                                           "WT T2","KO T2",
                                                                           "WT-w6","KO-w6"), ordered = T)

dir.create("figures/Fig1F_highly_changing_genus")

for (i in 1:length(unstable_taxa)){

# Longitudinal Plot #

ggline(genus_RA_to_plot %>% filter(taxa == glue("{unstable_taxa[i]}")),x = "time_num", y = "CLR", add = c("mean_se","jitter"), color = "genotype",
       add.params = list(size=1,alpha=0.4)) +
  stat_compare_means(aes(group = genotype), label = "p.signif", method = "wilcox.test", size = 2) +
  stat_summary(mapping = aes(fill=genotype), fun = "mean", geom = "point", color = "black", pch = 21, size = 2)+
  ylab("Center-log-ratio") + xlab("Weeks") +
  scale_fill_manual(values = c("WT"="#5da1d4ff", "KO"="#c88cb2ff")) +
  scale_color_manual(values = c("WT"="#5da1d4ff", "KO"="#c88cb2ff")) +
  theme_bw() + theme(axis.text=element_text(size=7, color = "black"), 
                          axis.title=element_text(size=7), 
                          legend.text = element_text(size=7), 
                          legend.title = element_text(size=0), 
                          plot.title = element_text(size=9, face="bold", hjust = 0.5)) + 
  ggtitle(glue("{unstable_taxa[i]}"))
#facet_wrap(~taxa, scales = "free")

ggsave(file = glue("figures/Fig1F_highly_changing_genus/{unstable_taxa[i]}.pdf"), width = 6, height = 5, units="cm")

}
