

library(tidyverse)
library(phyloseq)
library(Biostrings)
library(ggplow2)
library(ggrepel)
library(vegan)
library(microbiome)
library(ggpubr)
library(dplyr)
library(Maaslin2)
library(ComplexHeatmap)
library(circlize)


# Code for - Differential Abundance Analysis: 16S data : validation experiment -------

# Setting working directory #
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ps <- readRDS("inputs/phyloseq_object_validation.rds")


# Cutoff #

(ps_fr <- subset_samples(ps, sample_sums(ps) > 5000))
# (ps_fr <- subset_samples(ps, sample_sums(ps) > 5000 & pregnancy_status == "P"))

# Deleting taxa with relative abundance lowest than 1e-5
minTotRelAbun <- 1e-5          
x <- taxa_sums(ps_fr)
keepTaxa <- (x / sum(x)) > minTotRelAbun
ps_f <- prune_taxa(keepTaxa, ps_fr)


#:::::::::::: Maaslin2 :::::::::::#
# Creating directory 
dir.create("outputs/differential_abundance_analysis_16S")
dir.create("outputs/differential_abundance_analysis_16S/lmm_maaslin2")
dir.create("outputs/differential_abundance_analysis_16S/lmm_maaslin2/ASVs")

# Pair-wise comparison
m2_s <- Maaslin2(
    input_data = data.frame(ps_f@otu_table),
    input_metadata = data.frame(ps_f@sam_data),
    output="outputs/differential_abundance_analysis_16S/lmm_maaslin2/ASVs/complete_model",
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
#::::::::::::::::::::::::::::::::::::::#
# Comparison : BL : WT vs KO #

m2_s <- Maaslin2(
  input_data = data.frame(subset_samples(ps_f, timepoint == "BL")@otu_table),
  input_metadata = data.frame(ps_f@sam_data) %>% filter(timepoint == "BL"),
  output="outputs/differential_abundance_analysis_16S/lmm_maaslin2/ASVs/BL_WTvsKO",
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("genotype"),
  #random_effects = c("mouse_id"),
  reference = c("genotype,WT"),
  min_abundance = 1e-8,
  min_prevalence = 0.12, # 80% modified rule (min(samples_in_time)/total_samples)*80
  correction = "BH",
  max_significance = 0.05,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)

#::::::::::::::::::::::::::::::::::::::#
# Comparison : w3 : WT vs KO #

m2_s <- Maaslin2(
  input_data = data.frame(subset_samples(ps_f, timepoint == "w3")@otu_table),
  input_metadata = data.frame(ps_f@sam_data) %>% filter(timepoint == "w3"),
  output="outputs/differential_abundance_analysis_16S/lmm_maaslin2/ASVs/w3_WTvsKO",
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("genotype"),
  #random_effects = c("mouse_id"),
  reference = c("genotype,WT"),
  min_abundance = 1e-8,
  min_prevalence = 0.12, # 80% modified rule (min(samples_in_time)/total_samples)*80
  correction = "BH",
  max_significance = 0.05,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)
#::::::::::::::::::::::::::::::::::::::#
# Comparison : w1 : WT vs KO #

m2_s <- Maaslin2(
  input_data = data.frame(subset_samples(ps_f, timepoint == "w1")@otu_table),
  input_metadata = data.frame(ps_f@sam_data) %>% filter(timepoint == "w1"),
  output="outputs/differential_abundance_analysis_16S/lmm_maaslin2/ASVs/w1_WTvsKO",
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("genotype"),
  #random_effects = c("mouse_id"),
  reference = c("genotype,WT"),
  min_abundance = 1e-8,
  min_prevalence = 0.12, # 80% modified rule (min(samples_in_time)/total_samples)*80
  correction = "BH",
  max_significance = 0.05,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)

#::::::::::::::::::::::::::::::::::::::#
# Comparison : w2 : WT vs KO #

m2_s <- Maaslin2(
  input_data = data.frame(subset_samples(ps_f, timepoint == "w2")@otu_table),
  input_metadata = data.frame(ps_f@sam_data) %>% filter(timepoint == "w2"),
  output="outputs/differential_abundance_analysis_16S/lmm_maaslin2/ASVs/w2_WTvsKO",
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("genotype"),
  #random_effects = c("mouse_id"),
  reference = c("genotype,WT"),
  min_abundance = 1e-8,
  min_prevalence = 0.12, # 80% modified rule (min(samples_in_time)/total_samples)*80
  correction = "BH",
  max_significance = 0.05,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)

#------------------------------------------------------------------------------
# Heatmap Plot #
ps_clr <- microbiome::transform(ps_f, transform = "clr")

# Loading Longitudinal per-feature test results from Maaslin2
path = "outputs/differential_abundance_analysis_16S/lmm_maaslin2"

complete_model <- read_tsv(paste0(path,"/ASVs/complete_model/all_results.tsv")) %>%
  mutate(group = "complete_model")

WTvsKO_BL <- read_tsv(paste0(path,"/ASVs/BL_WTvsKO/all_results.tsv"))%>%
  mutate(group = "WTvsKOBL")

WTvsKO_w1 <- read_tsv(paste0(path,"/ASVs/w1_WTvsKO/all_results.tsv"))%>%
  mutate(group = "WTvsKOw1")

WTvsKO_w2 <- read_tsv(paste0(path,"/ASVs/w2_WTvsKO/all_results.tsv"))%>%
  mutate(group = "WTvsKOw2")

WTvsKO_w3 <- read_tsv(paste0(path,"/ASVs/w3_WTvsKO/all_results.tsv"))%>%
  mutate(group = "WTvsKOw3")

# Using FDR cut-off qval < 0.1 #

fdr_filtered_features <- rbind(complete_model %>% filter(value=="KO"),WTvsKO_BL,WTvsKO_w1,WTvsKO_w2,WTvsKO_w3) %>%
  subset(., qval <= 0.05) %>%
  select(feature) %>%
  distinct()


# Loading Metadata #
preg_metadata <- data.frame(ps_clr@sam_data) %>%
  mutate(genotype = factor(genotype,levels = c("WT","KO"), ordered = T)) %>%
  arrange(timepoint,genotype)

# Loading Taxonomy Table #
taxonomy_table <- data.frame(ps_clr@tax_table) %>%
  cbind(data.frame(ASVs=rownames(.)),.) %>%
  mutate(taxa_name = paste0(ASVs,"|",Family,"|",Species)) %>%
  left_join(fdr_filtered_features,.,by = c("feature"="ASVs")) %>%
  filter(Family != "unknown") # filtering Unclassified Bacteria

fdr_filtered_features <- taxonomy_table %>%
  select(feature)

# FDR filtered Features #
qval_filtered_features <- rbind(complete_model %>% filter(value=="KO"),WTvsKO_BL,WTvsKO_w1,WTvsKO_w2,WTvsKO_w3) %>%
  as_tibble() %>%
  select(feature,qval,group) %>%
  pivot_wider(names_from=group, values_from = qval) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 1)) %>%
  left_join(fdr_filtered_features,., by = "feature")

colnames(qval_filtered_features) <- c("feature","FDR_WTvsKO","FDR_WTvsKO_BL","FDR_WTvsKO_w1","FDR_WTvsKO_w2","FDR_WTvsKO_w3")


# Prevalence Computation #
prevalence_filtered_features <- rbind(complete_model %>% filter(value=="WT"),WTvsKO_BL,WTvsKO_w1,WTvsKO_w2,WTvsKO_w3) %>%
  as_tibble() %>%
  select(feature,N,N.not.0) %>% 
  group_by(feature) %>%
  mutate(prevalence = sum(N.not.0)/sum(N)) %>%
  select(feature, prevalence) %>%
  distinct() %>%
  left_join(fdr_filtered_features,., by = "feature")


# Mean Abundance table
abundance_table <- data.frame(ps_clr@sam_data,ps_clr@otu_table) %>%
  select(timepoint,genotype,fdr_filtered_features$feature) %>%
  pivot_longer(.,cols=starts_with("ASV"),names_to = "taxa",values_to="CLR") %>%
  group_by(timepoint,genotype,taxa) %>%
  summarize(mean = mean(CLR)) %>%
  left_join(.,taxonomy_table, by = c("taxa"="feature")) %>%
  select(timepoint,genotype,taxa,taxa_name,mean) %>%
  arrange(.,desc(genotype))

mean_metadata <- abundance_table %>%
  left_join(fdr_filtered_features,., by = c("feature"="taxa")) %>%
  select(-feature) %>%
  pivot_wider(.,names_from = "taxa_name",values_from="mean") %>%
  select(timepoint,genotype) %>%
  mutate(time_genotype = paste0(timepoint,"-",genotype))

mean_abundance_table <- abundance_table %>%
  left_join(fdr_filtered_features,., by = c("feature"="taxa")) %>%
  select(-feature) %>%
  pivot_wider(.,names_from = "taxa_name",values_from="mean") %>%
  mutate(time_genotype = paste0(timepoint,"-",genotype)) %>%
  column_to_rownames(var="time_genotype") %>%
  select(starts_with("ASV")) %>%
  t(.)

#::::::::::::::::::::::::::::::#
# Building Heatmap per ASVs # 


# Computing zscore for better visualization #

matrix_zscore <- t(scale(t(mean_abundance_table), center = T, scale = T))

# Adding * to significant features in the comparison

# fdrs
is.sig_fdr <- qval_filtered_features$FDR_WTvsKO < 0.05
pch_fdr = rep("*", length(is.sig_fdr))
pch_fdr[!is.sig_fdr] = NA


is.sig_fdr_BL <- qval_filtered_features$FDR_WTvsKO_BL < 0.05
pch_fdr_BL = rep("*", length(is.sig_fdr_BL))
pch_fdr_BL[!is.sig_fdr_BL] = NA

is.sig_fdr_w1 <- qval_filtered_features$FDR_WTvsKO_w1 < 0.05
pch_fdr_w1 = rep("*", length(is.sig_fdr_w1))
pch_fdr_w1[!is.sig_fdr_w1] = NA

is.sig_fdr_w2 <- qval_filtered_features$FDR_WTvsKO_w2 < 0.05
pch_fdr_w2 = rep("*", length(is.sig_fdr_w2))
pch_fdr_w2[!is.sig_fdr_w2] = NA

is.sig_fdr_w3 <- qval_filtered_features$FDR_WTvsKO_w3 < 0.05
pch_fdr_w3 = rep("*", length(is.sig_fdr_w3))
pch_fdr_w3[!is.sig_fdr_w3] = NA




# FDR color bar function

fdr_col_fun = colorRamp2(c(0,2,3,4), c("#55C667FF","#20A387FF","#287D8EFF","#404788FF"))

# Color Ramp for mouse id #
preg_metadata$animal <- as.numeric(preg_metadata$animal)

col_mouse = colorRamp2(c(1005, 1045, 1081), c("darkgreen", "yellow", "darkorange"))


# Top Annotation
tha = HeatmapAnnotation(Time = mean_metadata$timepoint,
                        Genotype = mean_metadata$genotype,
                        col = list(#Mouse = col_mouse,
                          Time = c("BL"="white","w3"="white","w1"="white","w2"="white"),
                          Genotype = c("WT"="#5da1d4ff","KO"="#c88cb2ff")),
                        annotation_name_gp = gpar(fontsize = 9),
                        annotation_name_side = "right",
                        border = TRUE)

# Phylum annotation curation #

phylum = taxonomy_table %>%
  select(feature,Phylum) %>%
  left_join(fdr_filtered_features,.,by = "feature") %>%
  column_to_rownames(var = "feature")


rha = rowAnnotation( Phylum = phylum$Phylum,
                     Prevalence = anno_barplot(prevalence_filtered_features[,c("prevalence")]),
                     `FDR WTvsKO`=anno_simple(-log10(qval_filtered_features[,c("FDR_WTvsKO")]),col = fdr_col_fun, pch = pch_fdr),
                     `FDR WTvsKO-BL` = anno_simple(-log10(qval_filtered_features[,c("FDR_WTvsKO_BL")]),col = fdr_col_fun, pch = pch_fdr_BL),
                     `FDR WTvsKO-w1` = anno_simple(-log10(qval_filtered_features[,c("FDR_WTvsKO_w1")]),col = fdr_col_fun, pch = pch_fdr_w1),
                     `FDR WTvsKO-w2` = anno_simple(-log10(qval_filtered_features[,c("FDR_WTvsKO_w2")]),col = fdr_col_fun, pch = pch_fdr_w2),
                     `FDR WTvsKO-w3` = anno_simple(-log10(qval_filtered_features[,c("FDR_WTvsKO_w3")]),col = fdr_col_fun, pch = pch_fdr_w3),
                     col = list(Phylum= c("Actinobacteriota"="#FFDC9199",
                                          "Bacteroidota"="#0072B599",
                                          "Desulfobacterota"="#7876B199",
                                          "Firmicutes"="#BC3C2999",
                                          "Proteobacteria"="#20854E99",
                                          "Deferribacteres"="#E1872799",
                                          "Verrucomicrobia"="lightblue",
                                          "Bacteria_unclassified"="grey")),
                     annotation_name_gp = gpar(fontsize = 9),
                     annotation_name_side = "top",
                     border = TRUE)

# now we generate two legends, one for the p-value
# see how we define the legend for p value
lgd_fdr = Legend(title = "FDR value", col_fun = fdr_col_fun, at = c(0, 1, 2), 
                 labels = c("0.05","0.01","0.001"))
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")

zscore_col_fun = circlize::colorRamp2(c(-1.8, 0, 2.3), c("darkblue", "white", "darkred"))


heat_map_1 <- Heatmap(matrix_zscore, name = "z-score",
                      top_annotation = tha,
                      left_annotation = rha,
                      col = zscore_col_fun,
                      na_col = "darkgrey",
                      row_names_side = "right",
                      column_title_gp = gpar(fontsize = 9),
                      column_names_gp = gpar(fontsize = 9),
                      column_names_side = "top",
                      column_names_rot = 90,
                      column_names_centered = TRUE,
                      row_names_gp = gpar(fontsize = 7),
                      cluster_rows = TRUE,
                      cluster_row_slices = FALSE,
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = TRUE,
                      show_row_dend = TRUE,
                      row_dend_side = "left",
                      column_dend_side = "bottom",
                      show_parent_dend_line = FALSE,
                      column_labels = rep("",ncol(mean_abundance_table)),
                      row_title = NULL,
                      split = phylum$Phylum,
                      column_split = factor(mean_metadata$timepoint,levels = c("BL","w1","w2","w3")),
                      border = TRUE,
                      border_gp = gpar(col = "black"),
                      gap = unit(1, "mm"),
                      row_gap = unit(1, "mm"),
                      column_gap = unit(1, "mm"),

)

# PDF heatmap
pdf("supp_figures/Supp_Fig5B_Heatmap_ASVs_mean_CLR_ZSCORE.pdf", width= 7, height=8)

draw(heat_map_1,
     annotation_legend_list = list(lgd_fdr, 
                                   lgd_sig),
     padding = unit(c(2, 2, 2, 25), "mm"),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

dev.off()
