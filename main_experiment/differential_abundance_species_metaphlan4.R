
library(vegan)
library(ggplot2)
library(tidyverse)
library(glue)
library(ggpubr)
library(EnvStats)
library(Maaslin2)

# Code for - Differential Abundance Analysis - from Metaphlan4 - Species level: Shot-gun data-------

# Setting working directory #
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Relative Abundance Data #

metaphlan <- read.table("inputs/metaphlan_rel_abundances.txt", sep = "\t", check.names = F, header = T) %>%
  separate(col = clade_name, into = c("kingdom","phylum","class","order","family","genus","species","SGBs"), sep = '[|]')


# Taxonomy Table #

taxonomy_table <- metaphlan %>%
  filter(!is.na(species) & !is.na(SGBs)) %>%
  select(kingdom,phylum,class,order,family,genus,species,SGBs)

write.table(taxonomy_table, file = "inputs/taxonomy_table_species.txt", sep = "\t", quote = F, row.names = F)

# relative abundance data
species <- metaphlan %>%
  filter(is.na(SGBs) & !is.na(species)) %>%
  select(-kingdom,-phylum,-class,-order,-family,-genus,-SGBs,-clade_taxid) %>%
  #rename(taxa = species) %>%
  column_to_rownames(var="species")

# metadata  
preg_metadata <- read.csv("metaInfo/shotgun_samples.csv", check.names = F, header = T)[,1:8] %>% 
  separate(col= Name, into = c("experiment","type_of_sample","mouse_id","time","genotype","preg_status"), sep ='_') %>%
  column_to_rownames(var = "Barcode") %>%
  mutate(trimester = if_else(time == "week0","BL",
                             if_else(time == "week3","w3","w6")))

sample_ids <- substr(colnames(species), 1, 16) # extracting only the first 6 strings.

# change colnames
colnames(species) <- sample_ids

# Deleting taxa with relative abundance lowest than 1e-5
minTotRelAbun <- 1e-5        
x <- rowSums(species)
keepTaxa <- x > minTotRelAbun
(species_filt <- species[keepTaxa, ])  

# Combine metadata with relative abundance table

s_met_and_taxa <- merge(preg_metadata[,c(3,8,9,10,11,13)], t(species_filt), by = 0)
colnames(s_met_and_taxa)[1] <- "sample_id"

s_metadata <- s_met_and_taxa %>%
  select(sample_id, mouse_id,genotype,trimester) %>%
  column_to_rownames("sample_id")

s_matrix <- s_met_and_taxa %>%
  select(-mouse_id,-time, -genotype,-Status,-preg_status,-trimester) %>%
  column_to_rownames("sample_id")

#:::::::::::: lmm_maaslin2 :::::::::::#
# Pair-wise comparison
dir.create("outputs/differential_abundance_species_metaphlan4")
dir.create("outputs/differential_abundance_species_metaphlan4/lmm_maaslin2")
dir.create("outputs/differential_abundance_species_metaphlan4/lmm_maaslin2/species")


m2_s <- Maaslin2(
    input_data = s_matrix,
    input_metadata = s_metadata,
    output="outputs/differential_abundance_species_metaphlan4/lmm_maaslin2/species/alldata",
    normalization = "NONE",
    transform = "NONE",
    analysis_method = "LM",
    fixed_effects = c("trimester","genotype"),
    random_effects = c("mouse_id"),
    reference = c("trimester,BL","genotype,WT"),
    min_abundance = 0,
    min_prevalence = 0.2,
    correction = "BH",
    max_significance = 0.25,
    plot_heatmap = TRUE,
    heatmap_first_n = 50,
    plot_scatter = TRUE
)
#::::::::::::::::::::::::::::::::::::::#
# Comparison : BL : WT vs KO #

s_metadata_BL_WTvsKO <- subset(s_met_and_taxa, trimester == "BL")%>%
  as_tibble() %>%
  select(sample_id, mouse_id,genotype,trimester) %>%
  column_to_rownames("sample_id")

s_matrix_BL_WTvsKO <- subset(s_met_and_taxa, trimester == "BL") %>%
  as_tibble() %>%
  select(-mouse_id,-time, -genotype,-Status,-preg_status,-trimester) %>%
  column_to_rownames("sample_id")

m2_s_BL_WTvsKO <- Maaslin2(
  input_data = s_matrix_BL_WTvsKO,
  input_metadata = s_metadata_BL_WTvsKO,
  output="outputs/differential_abundance_species_metaphlan4/lmm_maaslin2/species/BL_WTvsKO",
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("genotype"),
  reference = c("genotype,WT"),
  min_abundance = 0.0,
  min_prevalence = 0.2,
  correction = "BH",
  max_significance = 0.25,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)

#::::::::::::::::::::::::::::::::::::::#
# Comparison : w3 : WT vs KO #

s_metadata_w3_WTvsKO <- subset(s_met_and_taxa, trimester == "w3") %>%
  select(sample_id, mouse_id,genotype,trimester) %>% 
  as_tibble() %>%
  column_to_rownames("sample_id")

s_matrix_w3_WTvsKO <- subset(s_met_and_taxa, trimester == "w3") %>%
  select(-mouse_id,-time, -genotype,-Status,-preg_status,-trimester) %>%
  as_tibble() %>%
  column_to_rownames("sample_id")

m2_s_w3_WTvsKO <- Maaslin2(
  input_data = s_matrix_w3_WTvsKO,
  input_metadata = s_metadata_w3_WTvsKO,
  output="outputs/differential_abundance_species_metaphlan4/lmm_maaslin2/species/w3_WTvsKO",
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("genotype"),
  # random_effects = c("Mouse"),
  reference = c("genotype,WT"),
  min_abundance = 0.0,
  min_prevalence = 0.2,
  correction = "BH",
  max_significance = 0.25,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)


# Comparison : w6 : WT vs KO #

s_metadata_w6_WTvsKO <- subset(s_met_and_taxa, trimester == "w6") %>%
  select(sample_id, mouse_id,genotype,trimester) %>% 
  as_tibble() %>%
  column_to_rownames("sample_id")

s_matrix_w6_WTvsKO <- subset(s_met_and_taxa, trimester == "w6") %>%
  select(-mouse_id,-time, -genotype,-Status,-preg_status,-trimester) %>%
  as_tibble() %>%
  column_to_rownames("sample_id")

m2_s_w6_WTvsKO <- Maaslin2(
  input_data = s_matrix_w6_WTvsKO,
  input_metadata = s_metadata_w6_WTvsKO,
  output="outputs/differential_abundance_species_metaphlan4/lmm_maaslin2/species/w6_WTvsKO",
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("genotype"),
  reference = c("genotype,WT"),
  min_abundance = 0.0,
  min_prevalence = 0.2,
  correction = "BH",
  max_significance = 0.25,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)



#::::::::::::::::::::::::::::::::::::::#
# Comparison : WT: longitudinal #

s_metadata_WT_long <- subset(s_met_and_taxa, genotype == "WT") %>%
  select(sample_id, mouse_id,genotype,trimester) %>% 
  as_tibble() %>%
  column_to_rownames("sample_id")

s_matrix_WT_long <- subset(s_met_and_taxa, genotype == "WT") %>%
  select(-mouse_id,-time, -genotype,-Status,-preg_status,-trimester) %>%
  as_tibble() %>%
  column_to_rownames("sample_id")

m2_s_WT_long <- Maaslin2(
  input_data = s_matrix_WT_long,
  input_metadata = s_metadata_WT_long,
  output="outputs/differential_abundance_species_metaphlan4/lmm_maaslin2/species/WT_long",
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("trimester"),
  random_effects = c("mouse_id"),
  reference = c("trimester,BL"),
  min_abundance = 0.0,
  min_prevalence = 0.2,
  correction = "BH",
  max_significance = 0.25,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)

#::::::::::::::::::::::::::::::::::::::#
# Comparison : KO: longitudinal #

s_metadata_KO_long <- subset(s_met_and_taxa, genotype == "KO") %>%
  select(sample_id, mouse_id,genotype,trimester) %>% 
  as_tibble() %>%
  column_to_rownames("sample_id")

s_matrix_KO_long <- subset(s_met_and_taxa, genotype == "KO") %>%
  select(-mouse_id,-time, -genotype,-Status,-preg_status,-trimester) %>%
  as_tibble() %>%
  column_to_rownames("sample_id")

m2_s_KO_long <- Maaslin2(
  input_data = s_matrix_KO_long,
  input_metadata = s_metadata_KO_long,
  output="outputs/differential_abundance_species_metaphlan4/lmm_maaslin2/species/KO_long",
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("trimester"),
  random_effects = c("mouse_id"),
  reference = c("trimester,BL"),
  min_abundance = 0.0,
  min_prevalence = 0.2,
  correction = "BH",
  max_significance = 0.25,
  plot_heatmap = TRUE,
  heatmap_first_n = 50,
  plot_scatter = TRUE
)

#------------------------------------------------------------------------------
# Heatmap Plot #

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

WTvsKO_BL <- read_tsv(paste0(path,"/species/BL_WTvsKO/all_results.tsv"))%>%
  mutate(group = "WTvsKOBL")

WTvsKO_w3 <- read_tsv(paste0(path,"/species/w3_WTvsKO/all_results.tsv"))%>%
  mutate(group = "WTvsKOw3")

WTvsKO_w6 <- read_tsv(paste0(path,"/species/w6_WTvsKO/all_results.tsv"))%>%
  mutate(group = "WTvsKOw6")


# Using FDR cut-off qval < 0.05 #

fdr_filtered_features <- rbind(complete_model %>% filter(value=="KO"),WTvsKO_BL,WTvsKO_w3,WTvsKO_w6) %>%
  subset(., qval <= 0.05) %>%
  select(feature) %>%
  distinct()


# Loading Metadata #
preg_metadata <- read.table("metaInfo/metainformation_16S_and_shotgun_metagenomics.txt", check.names = F, header = T, sep = "\t") %>% 
  filter(!is.na(barcode_shotgun)) %>%
  mutate(genotype = factor(genotype,levels = c("WT","KO"), ordered = T)) %>%
  arrange(timepoint,genotype)

# Loading Taxonomy Table #
taxonomy_table <- read.table("inputs/taxonomy_table_species.txt", sep = "\t", 
                             header = T, check.names = F) %>%
  select(-SGBs) %>%
  distinct() %>%
  left_join(fdr_filtered_features,.,by = c("feature"="species")) %>%
  mutate(taxa_name = case_when(str_detect(feature, 's__GGB') ~ paste0(family," | ",feature),
                               TRUE ~ feature)) %>%
  filter(phylum != "p__Bacteria_unclassified") # filtering Unclassified Bacteria

fdr_filtered_features <- taxonomy_table %>%
  select(feature)

# FDR filtered Features #
qval_filtered_features <- rbind(complete_model %>% filter(value=="KO"),WTvsKO_BL,WTvsKO_w3,WTvsKO_w6) %>%
  as_tibble() %>%
  select(feature,qval,group) %>%
  pivot_wider(names_from=group, values_from = qval) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 1)) %>%
  left_join(fdr_filtered_features,., by = "feature")

colnames(qval_filtered_features) <- c("feature","FDR_WTvsKO","FDR_WTvsKO_BL","FDR_WTvsKO_w3","FDR_WTvsKO_w6")

# Prevalence Computation #
prevalence_filtered_features <- rbind(WTvsKO_BL,WTvsKO_w3,WTvsKO_w6) %>%
  as_tibble() %>%
  select(feature,N,N.not.0) %>% 
  group_by(feature) %>%
  mutate(prevalence = sum(N.not.0)/sum(N)) %>%
  select(feature, prevalence) %>%
  distinct() %>%
  left_join(fdr_filtered_features,., by = "feature")

# Loading Abundance Table #

abundance_table <- read.table("inputs/metaphlan_rel_abundances.txt", sep = "\t", check.names = F, header = T) %>%
  separate(col = clade_name, into = c("kingdom","phylum","class","order","family","genus","species","SGBs"), sep = '[|]') %>%
  filter(is.na(SGBs) & !is.na(species)) %>%
  select(-kingdom,-phylum,-class,-order,-family,-genus,-SGBs,-clade_taxid) %>%
  left_join(fdr_filtered_features,., by = c("feature"="species")) %>%
  column_to_rownames(var="feature")

sample_names <- str_sub(colnames(abundance_table),1,16)
colnames(abundance_table) = sample_names
# Matching order of samples with metadata #
abundance_table <- abundance_table[,preg_metadata$barcode_shotgun]
#rownames(abundance_table) = taxonomy_table$taxa_name

# Mean Relative Abundances #

mean_abundances <- t(abundance_table) %>%
  cbind(preg_metadata,.) %>%
  pivot_longer(.,cols=contains("s__"), names_to = "species",values_to="RAs") %>%
  select(timepoint,genotype,species,RAs) %>%
  mutate(time_genotype = paste0(timepoint,"-",genotype)) %>%
  group_by(time_genotype,species) %>%
  summarize(mean = mean(RAs)) %>%
  pivot_wider(.,names_from = "species",values_from="mean") %>%
  column_to_rownames(var="time_genotype") %>%
  t(.)

mean_abundances <- mean_abundances[fdr_filtered_features$feature,c("BL-WT","BL-KO","w3-WT","w3-KO","w6-WT","w6-KO")]

mean_metadata <- data.frame(time_genotype = colnames(mean_abundances)) %>%
  separate(., col = "time_genotype", into = c("timepoint","genotype"), sep = "-") %>%
  mutate(time_genotype = paste0(timepoint,"-",genotype))

#::::::::::::::::::::::::::::::#
# Building Heatmap per species # 

# Computing zscore for better visualization #
matrix_zscore <- t(scale(t(mean_abundances), center = T, scale = T))

# Adding * to significant features in the comparison

# fdrs
is.sig_fdr <- qval_filtered_features$FDR_WTvsKO < 0.05
pch_fdr = rep("*", length(is.sig_fdr))
pch_fdr[!is.sig_fdr] = NA


is.sig_fdr_BL <- qval_filtered_features$FDR_WTvsKO_BL < 0.05
pch_fdr_BL = rep("*", length(is.sig_fdr_BL))
pch_fdr_BL[!is.sig_fdr_BL] = NA

is.sig_fdr_w3 <- qval_filtered_features$FDR_WTvsKO_w3 < 0.05
pch_fdr_w3 = rep("*", length(is.sig_fdr_w3))
pch_fdr_w3[!is.sig_fdr_w3] = NA

is.sig_fdr_w6 <- qval_filtered_features$FDR_WTvsKO_w6 < 0.05
pch_fdr_w6 = rep("*", length(is.sig_fdr_w6))
pch_fdr_w6[!is.sig_fdr_w6] = NA


# FDR color bar function
fdr_col_fun = colorRamp2(c(0,2,3,4), c("#55C667FF","#20A387FF","#287D8EFF","#404788FF"))

# Color Ramp for mouse id #
preg_metadata$mouse_id <- as.numeric(preg_metadata$mouse_id)

col_mouse = colorRamp2(c(2104, 2183, 2258), c("darkgreen", "yellow", "darkorange"))


# Top Annotation
tha = HeatmapAnnotation(Time = mean_metadata$timepoint,
                        Genotype = mean_metadata$genotype,
                        col = list(
                          Time = c("BL"="white","w3"="white","w6"="white"),
                          Genotype = c("WT"="#5da1d4ff","KO"="#c88cb2ff")),
                        annotation_name_gp = gpar(fontsize = 7),
                        annotation_name_side = "right",
                        border = TRUE)

# Phylum annotation curation #

phylum = taxonomy_table %>%
  select(taxa_name,phylum) %>%
  mutate_at("phylum", str_replace, "p__", "") %>%
  column_to_rownames(var = "taxa_name")

rownames(matrix_zscore) <- rownames(phylum)

rha = rowAnnotation( Phylum = phylum$phylum,
                     Prevalence = anno_barplot(prevalence_filtered_features[,c("prevalence")]),
                     `FDR WTvsKO`=anno_simple(-log10(qval_filtered_features[,c("FDR_WTvsKO")]),col = fdr_col_fun, pch = pch_fdr),
                     `FDR WTvsKO-BL` = anno_simple(-log10(qval_filtered_features[,c("FDR_WTvsKO_BL")]),col = fdr_col_fun, pch = pch_fdr_BL),
                     `FDR WTvsKO-w3` = anno_simple(-log10(qval_filtered_features[,c("FDR_WTvsKO_w3")]),col = fdr_col_fun, pch = pch_fdr_w3),
                     `FDR WTvsKO-w6` = anno_simple(-log10(qval_filtered_features[,c("FDR_WTvsKO_w6")]),col = fdr_col_fun, pch = pch_fdr_w6),
                     col = list(Phylum= c("Actinobacteria"="#FFDC9199",
                                          "Bacteroidetes"="#0072B599",
                                          "Desulfobacterota"="#7876B199",
                                          "Firmicutes"="#BC3C2999",
                                          "Proteobacteria"="#20854E99",
                                          "Deferribacteres"="#E1872799",
                                          "Verrucomicrobia"="lightblue",
                                          "Bacteria_unclassified"="grey")),
                     annotation_name_gp = gpar(fontsize = 7),
                     annotation_name_side = "top",
                     border = TRUE)

# now we generate two legends, one for the p-value
# see how we define the legend for p value
lgd_fdr = Legend(title = "FDR value", col_fun = fdr_col_fun, at = c(0, 1, 2), 
                 labels = c("0.05","0.01","0.001"))
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")

zscore_col_fun = circlize::colorRamp2(c(-2, 0, 2), c("darkblue", "white", "darkred"))


heat_map_1 <- Heatmap(matrix_zscore, name = "z-score",
                      top_annotation = tha,
                      left_annotation = rha,
                      col = zscore_col_fun,
                      na_col = "darkgrey",
                      row_names_side = "right",
                      column_title_gp = gpar(fontsize = 7),
                      column_names_gp = gpar(fontsize = 7),
                      column_names_side = "top",
                      column_names_rot = 90,
                      column_names_centered = TRUE,
                      row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                      cluster_rows = TRUE,
                      cluster_row_slices = FALSE,
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = TRUE,
                      show_row_dend = TRUE,
                      row_dend_side = "left",
                      column_dend_side = "bottom",
                      show_parent_dend_line = FALSE,
                      column_labels = rep("",6),
                      row_title = NULL,
                      split = phylum$phylum,
                      column_split = factor(mean_metadata$timepoint,levels = c("BL","w3","w6")),
                      border = TRUE,
                      border_gp = gpar(col = "black"),
                      gap = unit(1, "mm"),
                      row_gap = unit(1, "mm"),
                      column_gap = unit(1, "mm"),
                      # column_order = metadata$sample_id,
                      
)

# PDF heatmap
pdf("figures/Fig2C_Heatmap_metaphlan4_meanRA_ZSCORE.pdf", width= 6, height=5)

draw(heat_map_1,
     annotation_legend_list = list(lgd_fdr, 
                                   lgd_sig),
     padding = unit(c(2, 2, 2, 25), "mm"),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

dev.off()

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Individual Plots #

# Relative Abundance Data #

metaphlan <- read.table("inputs/metaphlan_rel_abundances.txt", sep = "\t", check.names = F, header = T) %>%
  separate(col = clade_name, into = c("kingdom","phylum","class","order","family","genus","species","SGBs"), sep = '[|]')

species <- metaphlan %>%
  filter(is.na(SGBs) & !is.na(species)) %>%
  select(-kingdom,-phylum,-class,-order,-family,-genus,-SGBs,-clade_taxid) %>%
  rename(taxa = species) %>%
  column_to_rownames(var="taxa")/100

sample_ids <- substr(colnames(species), 1, 16) # extracting only the first 6 strings.

# change colnames
colnames(species) <- sample_ids

s_met_and_taxa <- merge(preg_metadata %>% column_to_rownames(var="barcode_shotgun"), t(species), by = 0) %>%
  mutate(Condition = paste0(genotype," ",timepoint))
colnames(s_met_and_taxa)[1] <- "sample_id"


interesting_taxa <- c("s__Adlercreutzia_caecicola","s__GGB24102_SGB40205","s__Turicibacter_sp_TS3",
                      "s__Parabacteroides_distasonis","s__Muribaculum_intestinale",
                      "s__Rikenellaceae_bacterium","s__GGB28849_SGB41516","s__GGB3171_SGB4185")


s_met_and_taxa %>% 
  select(sample_id,timepoint,genotype,interesting_taxa,Condition)

#for (i in 1:length(interesting_taxa)) {

interesting_taxa_long <- s_met_and_taxa %>% 
  select(sample_id,timepoint,genotype,interesting_taxa,Condition) %>%
  pivot_longer(.,cols=s__Adlercreutzia_caecicola:s__GGB3171_SGB4185, values_to = "RA", names_to = "taxa")

interesting_taxa_long$Condition = factor(interesting_taxa_long$Condition, levels = c("WT BL","KO BL",
                                                                                     "WT w3","KO w3",
                                                                                     "WT w6","KO w6"),
                                         ordered = T)

p <- ggplot(interesting_taxa_long, aes(x=Condition, y=RA, fill=Condition)) + geom_boxplot() + geom_point(pch=21, size=2)
p <- p + ylab("Relative Abundances")
p <- p + scale_color_manual(values = c("WT BL"="#5da1d4ff", "KO BL"="#c88cb2ff", 
                                       "WT w3"="#5da1d4ff","KO w3"="#c88cb2ff",
                                       "WT w6"="#5da1d4ff","KO w6"="#c88cb2ff"))
p <- p + scale_fill_manual(values = c("WT BL"="#5da1d4ff", "KO BL"="#c88cb2ff", 
                                      "WT w3"="#5da1d4ff","KO w3"="#c88cb2ff",
                                      "WT w6"="#5da1d4ff","KO w6"="#c88cb2ff"))
p <- p + stat_compare_means(label = "p.signif", comparisons = list(c("WT BL","KO BL"), 
                                                                   c("WT w3","KO w3"),
                                                                   c("WT w6","KO w6")
),  method = "wilcox.test", size = 2)
p <- p + theme_classic() + theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=7), legend.text = element_text(size=7), 
                                 legend.title = element_text(size=0), plot.title = element_text(size=7, face="bold", hjust = 0.5), axis.text.x = element_text(angle=45, hjust = 1))
p <- p + facet_wrap(~taxa, scales = "free")
p

ggsave(filename = "figures/Fig2D_longitudinal_interesting_species.pdf", width = 23, height = 20, units = "cm")

