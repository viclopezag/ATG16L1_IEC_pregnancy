
library(dplyr)
library(tidyverse)
library(variancePartition)
library(reshape2)
library(glue)
library(microbiome)
library(EnvStats)
library(ggpubr)
library(ggrepel)
library(data.table)


# Code for - Variance Partition Analyses: Shotgun metagenomics - HUMAnN 3.6 -------

# Setting Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# metadata raw #

metadata_raw <- read.table("metaInfo/metainformation_16S_and_shotgun_metagenomics.txt", 
                           sep = "\t", header=T) %>%
  filter(!is.na(shotgun_sample_id))

# Variance Partition Timepoint #

vp_time <- c("all","BL","w3","w6")

for (u in 1:length(vp_time)) {
# Loading metadata #

metadata <- metadata_raw %>%
  filter(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  mutate(mouse_id = as.character(mouse_id),
         timepoint = as.character(timepoint),
         genotype = as.character(genotype),
         delivery_time = as.character(delivery_time),
         female_cage_of_birth = as.character(female_cage_of_birth),
         male_cage_of_birth = as.character(male_cage_of_birth),
         male_genotype = as.character(male_genotype),
         weight_grams_scaled = scale(weight_grams),
         number_of_pups_scaled = scale(number_of_pups),
         lipocalin2_pg_ml_scaled = scale(lipocalin2_pg_ml))
  #column_to_rownames(var = "barcode_shotgun")
  #filter(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u])

# Gene families - EC #
EC_ids <- read_tsv("inputs/humann_merged_genefamilies-cpm-ec-named-no-taxa.tsv") %>% 
  filter(genefamily != "UNMAPPED" & genefamily != "UNGROUPED") %>%
  separate(.,col=genefamily,into=c("id","name"), sep = ": ") %>%
  filter(name != "NO_NAME")

sample_names <- substr(colnames(EC_ids),1,16)[-c(1:2)]

EC_ids <- EC_ids %>%
  select(id,name) %>%
  mutate(complete_name = paste0(id,"|",name))
  
EC_data <- read_tsv("inputs/humann_merged_genefamilies-cpm-ec-named-no-taxa.tsv") %>% 
  filter(genefamily != "UNMAPPED" & genefamily != "UNGROUPED") %>%
  separate(.,col=genefamily,into=c("id","name"), sep = ": ") %>%
  filter(name != "NO_NAME") %>% 
  column_to_rownames(var="id") %>%
  select(-name) %>%
  t(.)

rownames(EC_data) <- sample_names

EC_data <- EC_data[metadata$barcode_shotgun,]

colnames(EC_data) <- EC_ids$complete_name

# Gene families - GO #
GO_ids <- read_tsv("inputs/humann_merged_genefamilies-cpm-go-named-no-taxa.tsv") %>% 
  filter(genefamily != "UNMAPPED" & genefamily != "UNGROUPED") %>%
  separate(.,col=genefamily,into=c("id","name"), sep = ": ") %>%
  filter(name != "NO_NAME")

sample_names <- substr(colnames(GO_ids),1,16)[-c(1:2)]

GO_ids <- GO_ids %>%
  select(id,name) %>%
  mutate(complete_name = paste0(id,"|",name))

GO_data <- read_tsv("inputs/humann_merged_genefamilies-cpm-go-named-no-taxa.tsv") %>% 
  filter(genefamily != "UNMAPPED" & genefamily != "UNGROUPED") %>%
  separate(.,col=genefamily,into=c("id","name"), sep = ": ") %>%
  filter(name != "NO_NAME") %>% 
  column_to_rownames(var="id") %>%
  select(-name) %>%
  t(.)

rownames(GO_data) <- sample_names

GO_data <- GO_data[metadata$barcode_shotgun,]

colnames(GO_data) <- GO_ids$complete_name

# Gene families - KO #
KO_ids <- read_tsv("inputs/humann_merged_genefamilies-cpm-ko-named-no-taxa.tsv") %>% 
  filter(genefamily != "UNMAPPED" & genefamily != "UNGROUPED") %>%
  separate(.,col=genefamily,into=c("id","name"), sep = ": ") %>%
  filter(name != "NO_NAME")

sample_names <- substr(colnames(KO_ids),1,16)[-c(1:2)]

KO_ids <- KO_ids %>%
  select(id,name) %>%
  mutate(complete_name = paste0(id,"|",name))

KO_data <- read_tsv("inputs/humann_merged_genefamilies-cpm-ko-named-no-taxa.tsv") %>% 
  filter(genefamily != "UNMAPPED" & genefamily != "UNGROUPED") %>%
  separate(.,col=genefamily,into=c("id","name"), sep = ": ") %>%
  filter(name != "NO_NAME") %>% 
  column_to_rownames(var="id") %>%
  select(-name) %>%
  t(.)

rownames(KO_data) <- sample_names

KO_data <- KO_data[metadata$barcode_shotgun,]

colnames(KO_data) <- KO_ids$complete_name

# Gene families - RXN #
RXN_ids <- read_tsv("inputs/humann_merged_genefamilies-cpm-metacyc-rxn-named-no-taxa.tsv") %>% 
  filter(genefamily != "UNMAPPED" & genefamily != "UNGROUPED") %>%
  separate(.,col=genefamily,into=c("id","name"), sep = ": ") %>%
  filter(name != "NO_NAME")

sample_names <- substr(colnames(RXN_ids),1,16)[-c(1:2)]

RXN_ids <- RXN_ids %>%
  select(id,name) %>%
  mutate(complete_name = paste0(id,"|",name))

RXN_data <- read_tsv("inputs/humann_merged_genefamilies-cpm-metacyc-rxn-named-no-taxa.tsv") %>% 
  filter(genefamily != "UNMAPPED" & genefamily != "UNGROUPED") %>%
  separate(.,col=genefamily,into=c("id","name"), sep = ": ") %>%
  filter(name != "NO_NAME") %>% 
  column_to_rownames(var="id") %>%
  select(-name) %>%
  t(.)

rownames(RXN_data) <- sample_names

RXN_data <- RXN_data[metadata$barcode_shotgun,]

colnames(RXN_data) <- RXN_ids$complete_name

# Pathway Abundances - metacyc #
MCyc_ids <- read_tsv("inputs/humann_merged_pathabundance-cpm-no-taxa-metacyc-id.tsv") %>% 
  rename(id=metacyc_id,
         name=pathway) %>%
  filter(!is.na(name))

sample_names <- substr(colnames(MCyc_ids),1,16)[-c(1:2)]

MCyc_ids <- MCyc_ids %>%
  select(id,name) %>%
  mutate(complete_name = paste0(id,"|",name))

MCyc_data <- read_tsv("inputs/humann_merged_pathabundance-cpm-no-taxa-metacyc-id.tsv") %>% 
  rename(id=metacyc_id,
         name=pathway) %>%
  filter(!is.na(name)) %>%
  column_to_rownames(var="id") %>%
  select(-name) %>%
  t(.)

rownames(MCyc_data) <- sample_names

MCyc_data <- MCyc_data[metadata$barcode_shotgun,]

colnames(MCyc_data) <- MCyc_ids$complete_name

# Creating Directories for saving Variance Partition Outputs #
dir.create("outputs/variancePartition_hmn3.6")
#dir.create("plots/variancePartition_hmn3.6")

dir.create(glue("outputs/variancePartition_hmn3.6/{vp_time[u]}"))
#dir.create(glue("plots/variancePartition_no_dt/{vp_time[u]}"))

if (vp_time[u] == "all") {
formula_no_colinear <- list(formula_no_individual = ~ (1|genotype) + (1|timepoint) + number_of_pups_scaled + 
                                                       lipocalin2_pg_ml_scaled,
                            formula_individual = ~ (1|mouse_id) + (1|timepoint))

} else {
  
formula_no_colinear <- list(formula_no_individual = ~ (1|genotype)  + weight_grams_scaled + number_of_pups_scaled +
                                                      lipocalin2_pg_ml_scaled)
  
}  
                            
table_metadata_important <- metadata


level_abundances = list(EC_data,
                        GO_data,
                        KO_data,
                        RXN_data,
                        MCyc_data)


level_annotation = c("EC","GO","KO","RXN","MCyc")

#::::::::::::::::::::::::Variance Partition in all Taxa Levels :::::::::::::::::::::::::::::#
# Two Formulas will be evaluated 1) those variables non-colinear with group_g and 2) those variables non-colinear with R_NR #

non_colinear_formulas <- c("nc_no_individual","nc_individual")

set.seed(1989)


for (v in 1:length(formula_no_colinear)) {
  
  dir.create(glue("outputs/variancePartition_hmn3.6/{vp_time[u]}/{non_colinear_formulas[v]}")) 
  #dir.create(glue("plots/variancePartition_no_dt/{vp_time[u]}/{non_colinear_formulas[v]}")) 
  
  
  for (q in 1:length(level_annotation)){
    
    matrix_data <- level_abundances[[q]]

    # Delete Features not seen in at least 12% of the samples (80% of the samples per time/genotype group) 11 mouse/86 mouse = ~ 12% total samples , 
    # there are in average 15 mouse per time*genotype groups #
    
  if (vp_time[u]=="all"){
    # 80% rule : All timepoints : 9 samples/72 samples = ~ 12%
    matrix_data <- matrix_data[,colSums(matrix_data > 0) > 0.12*nrow(matrix_data)]
  } else if (vp_time[u]=="BL") {
    # 80% rule : BL : 12*0.8/24 mouse = ~ 0.37
    matrix_data <- matrix_data[,colSums(matrix_data > 0) > 0.37*nrow(matrix_data)]
    
  } else if (vp_time[u]=="w3") {
    # 80% rule : w3 : 12*0.8/24 mouse = 0.37
    matrix_data <- matrix_data[,colSums(matrix_data > 0) > 0.37*nrow(matrix_data)]
    
  } else {
    # 80% rule : w6 : 12*0.8/24 mouse = 0.37
    matrix_data <- matrix_data[,colSums(matrix_data > 0) > 0.37*nrow(matrix_data)] 
  }  
    # We are going to run Variance Partition on squared rooted transformed abundances #
    
    varPart <- fitExtractVarPartModel(t(matrix_data), formula_no_colinear[[v]], table_metadata_important, REML = F)
    
    # sort variables by mean fraction of varience explained
    
    vp <- sortCols(varPart)
    
    # Export results # 
    
    varMean <- colMeans(varPart)
    varPart_sorted <- varPart[, order(varMean, decreasing = TRUE)]
    
    write.table(varPart_sorted, glue("outputs/variancePartition_hmn3.6/{vp_time[u]}/{non_colinear_formulas[v]}/variance_partition_results_{level_annotation[q]}.txt"), sep = '\t', quote = FALSE)
    
    varMean <- as.data.frame(varMean)
    varMean$Parameter <- rownames(varMean)
    varMean <- varMean[order(varMean$varMean, decreasing = TRUE), ]
    
    write.table(varMean, glue("outputs/variancePartition_hmn3.6/{vp_time[u]}/{non_colinear_formulas[v]}/variance_partition_mean_{level_annotation[q]}.txt"), sep = '\t', quote = FALSE, row.names = FALSE)
    
    #### Violin plot of contribution of each variable to total variance ####
    
    varPart <- read.table(glue("outputs/variancePartition_hmn3.6/{vp_time[u]}/{non_colinear_formulas[v]}/variance_partition_results_{level_annotation[q]}.txt"), sep = "\t") %>%
      rownames_to_column(var = "feature") %>%
      as_tibble()
    
    covariates = structure(colnames(varPart))
    covariates = covariates[!covariates == "feature"] # Deleting feature from important variables
    covariates = covariates[!covariates == "Residuals"] # Deleting Residuals from important variables
    
    varPart_to_plot <- varPart %>%
      # select(-Residuals) %>%
      pivot_longer(2:ncol(.),names_to = "covariate")
    
    library(RColorBrewer)
    
    colourCount = ncol(varPart)-1
    nPalette = colorRampPalette(brewer.pal(8, "Set2"))(colourCount)
    #nPalette = c("#4A7BB7FF", "#DD3D2DFF", "#98CAE1FF", "#FDB366FF", "#000000")
    
    #pdf(glue("figures/variancePartition_no_dt/{vp_time[u]}/{non_colinear_formulas[v]}/violin_variance_partition_{level_annotation[q]}.pdf"), width = 7, height = 7)
    
    p <- ggplot(varPart_to_plot, aes(x = reorder(covariate,-value),  y = value, fill = covariate)) +
      geom_violin(trim = FALSE, alpha = 0.5, scale = "width") +
      geom_boxplot(width = 0.07) +
      labs(y = "Variance explained") +
      ylim(c(0, 1)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14), 
            axis.text.x = element_text(size = 14, angle = 45, hjust = 1.0),
            axis.title.x = element_blank(), 
            axis.title.y = element_text(size = 14),
            legend.position = "none",
            axis.text = element_text(size = 10, color = "black"), 
            axis.title = element_text(size = 10)) +
      theme(strip.placement = "outside", 
            strip.background  = element_blank(), 
            panel.border = element_blank()) + 
      scale_fill_manual(values = nPalette)
    
    plot(p)
    
    #dev.off()
    
   }
}

}

#------------------------------------------------------------------------------#
# Loading variance partition nc_no_individual - non-colinear variables #
#:::::::: All Timepoints ::::::::#
# EC #
varPart.mean.all.EC <- read.table("outputs/variancePartition_hmn3.6//all/nc_no_individual/variance_partition_mean_EC.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "EC") %>%
  mutate(time = "All Timepoints")

# GO 
varPart.mean.all.GO <- read.table("outputs/variancePartition_hmn3.6/all/nc_no_individual/variance_partition_mean_GO.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "GO") %>%
  mutate(time = "All Timepoints")
# KO
varPart.mean.all.KO <- read.table("outputs/variancePartition_hmn3.6/all/nc_no_individual/variance_partition_mean_KO.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "KO") %>%
  mutate(time = "All Timepoints")
# RXN
varPart.mean.all.RXN <- read.table("outputs/variancePartition_hmn3.6/all/nc_no_individual/variance_partition_mean_RXN.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "RXN") %>%
  mutate(time = "All Timepoints")
# MCyc
varPart.mean.all.MCyc <- read.table("outputs/variancePartition_hmn3.6/all/nc_no_individual/variance_partition_mean_MCyc.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "MCyc") %>%
  mutate(time = "All Timepoints")

#--------------------------------------------------------------------------------------------------------------#
#:::::::: BL: Baseline  ::::::::#
# EC 
varPart.mean.BL.EC <- read.table("outputs/variancePartition_hmn3.6/BL/nc_no_individual/variance_partition_mean_EC.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "EC") %>%
  mutate(time = "BL")

# GO 
varPart.mean.BL.GO <- read.table("outputs/variancePartition_hmn3.6/BL/nc_no_individual/variance_partition_mean_GO.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "GO") %>%
  mutate(time = "BL")
# KO
varPart.mean.BL.KO <- read.table("outputs/variancePartition_hmn3.6/BL/nc_no_individual/variance_partition_mean_KO.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "KO") %>%
  mutate(time = "BL")
# RXN
varPart.mean.BL.RXN <- read.table("outputs/variancePartition_hmn3.6/BL/nc_no_individual/variance_partition_mean_RXN.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "RXN") %>%
  mutate(time = "BL")
# MCyc
varPart.mean.BL.MCyc <- read.table("outputs/variancePartition_hmn3.6/BL/nc_no_individual/variance_partition_mean_MCyc.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "MCyc") %>%
  mutate(time = "BL")

#--------------------------------------------------------------------------------------------------------------#
#::::::::: w3 :::::::::::::#
# EC
varPart.mean.w3.EC <- read.table("outputs/variancePartition_hmn3.6/w3/nc_no_individual/variance_partition_mean_EC.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "EC") %>%
  mutate(time = "w3")
# GO 
varPart.mean.w3.GO <- read.table("outputs/variancePartition_hmn3.6/w3/nc_no_individual/variance_partition_mean_GO.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "GO") %>%
  mutate(time = "w3")
# KO
varPart.mean.w3.KO <- read.table("outputs/variancePartition_hmn3.6/w3/nc_no_individual/variance_partition_mean_KO.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "KO") %>%
  mutate(time = "w3")
# RXN
varPart.mean.w3.RXN <- read.table("outputs/variancePartition_hmn3.6/w3/nc_no_individual/variance_partition_mean_RXN.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "RXN") %>%
  mutate(time = "w3")
# MCyc
varPart.mean.w3.MCyc <- read.table("outputs/variancePartition_hmn3.6/w3/nc_no_individual/variance_partition_mean_MCyc.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "MCyc") %>%
  mutate(time = "w3")

#--------------------------------------------------------------------------------------------------------------#
#:::::::::::::: w6 ::::::::::::::::#
# EC 
varPart.mean.w6.EC <- read.table("outputs/variancePartition_hmn3.6/w6/nc_no_individual/variance_partition_mean_EC.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "EC") %>%
  mutate(time = "w6")
# GO 
varPart.mean.w6.GO <- read.table("outputs/variancePartition_hmn3.6/w6/nc_no_individual/variance_partition_mean_GO.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "GO") %>%
  mutate(time = "w6")
# KO
varPart.mean.w6.KO <- read.table("outputs/variancePartition_hmn3.6/w6/nc_no_individual/variance_partition_mean_KO.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "KO") %>%
  mutate(time = "w6")
# RXN
varPart.mean.w6.RXN <- read.table("outputs/variancePartition_hmn3.6/w6/nc_no_individual/variance_partition_mean_RXN.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "RXN") %>%
  mutate(time = "w6")
# MCyc
varPart.mean.w6.MCyc <- read.table("outputs/variancePartition_hmn3.6/w6/nc_no_individual/variance_partition_mean_MCyc.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "MCyc") %>%
  mutate(time = "w6")

#--------------------------------------------------------------------------------------------------------------#

varPart.mean <- rbind(varPart.mean.all.EC,
                      varPart.mean.all.GO,
                      varPart.mean.all.KO,
                      varPart.mean.all.RXN,
                      varPart.mean.all.MCyc,
                      varPart.mean.BL.EC,
                      varPart.mean.BL.GO,
                      varPart.mean.BL.KO,
                      varPart.mean.BL.RXN,
                      varPart.mean.BL.MCyc,
                      varPart.mean.w3.EC,
                      varPart.mean.w3.GO,
                      varPart.mean.w3.KO,
                      varPart.mean.w3.RXN,
                      varPart.mean.w3.MCyc,
                      varPart.mean.w6.EC,
                      varPart.mean.w6.GO,
                      varPart.mean.w6.KO,
                      varPart.mean.w6.RXN,
                      varPart.mean.w6.MCyc) %>%
  mutate(parameter_clean_name = ifelse(Parameter == "timepoint","Time",
                                       ifelse(Parameter == "mouse_id","Individual mice",
                                              ifelse(Parameter == "female_cage_of_birth","Maternal cage of birth",
                                                     ifelse(Parameter == "male_cage_of_birth","Paternal cage of birth",
                                                            ifelse(Parameter == "male_genotype","Paternal genotype",
                                                                   ifelse(Parameter == "number_of_pups_scaled", "Number of pups",
                                                                          ifelse(Parameter == "lipocalin2_pg_ml_scaled","Lipocalin (pg/ml)",
                                                                                 ifelse(Parameter == "genotype","Maternal genotype",
                                                                                        ifelse(Parameter == "delivery_time","Delivery time",
                                                                                               ifelse(Parameter == "weight_grams_scaled","Maternal weight (g)","Residuals")))))))))))


#parameters_color <- pal_jco("default", alpha = 1)(nrow(varPart.mean.all.MCyc))
palette_color = c("Individual mice"="#8c564b",
                  "Maternal genotype"="#9e0142",
                  "Maternal cage of birth"="#d53e4f",
                  "Maternal weight (g)"="#f46d43",
                  "Delivery time"="#fdae61",
                  "Number of pups"="#fee08b",
                  "Lipocalin (pg/ml)"="#abdda4",
                  "Paternal cage of birth"="#66c2a5",
                  "Paternal genotype" = "#3288bd",
                  "Time" = "#5e4fa2",
                  "Residuals"="darkgrey")


# Stacked Area Plot #

varPart.mean$level = factor(varPart.mean$level, levels = c("EC","GO","KO",
                                                           "RXN","MCyc"),
                            ordered = TRUE)
varPart.mean$Parameter = factor(varPart.mean$Parameter, levels = c("Residuals","female_cage_of_birth","genotype",
                                                                   "weight_scaled","number_of_pups_scaled","lipocalin2_pg_ml_scaled",
                                                                   "delivery_time","male_cage_of_birth","male_genotype","timepoint","mouse_id"
), ordered = TRUE)

varPart.mean %>%
  arrange(varMean) %>%
  ggbarplot(x = "level", y = "varMean", fill = "parameter_clean_name", add.params = list(size=0,fill="parameter_clean_name",alpha=1),
            ylab = "Mean Variance Explained (%)", xlab = "Functional Resolution", label = F) +
  theme_bw() + 
  theme(axis.text=element_text(size=7, color = NULL), 
        axis.title=element_text(size=7), 
        legend.text = element_text(size=7),
        legend.title = element_text(size=7), 
        plot.title = element_text(size=7, face="bold", hjust = 0.5),
        strip.background=element_rect(colour="black",
                                      fill="white")) +
  labs(fill = "Variables") + rotate_x_text(angle = 45) + 
  facet_grid(cols = vars(time)) + 
  scale_fill_manual(values = palette_color)

ggsave("figures/Fig3A_mean_variance_explained_humann3.6.pdf", width= 13, height=8, units = "cm")

write.table(varPart.mean,file = "outputs/variancePartition_hmn3.6/mean_variancePartition_nc_no_individual.txt", quote = F, sep = "\t", row.names = F)

# Computing mean +sd from variance_partition_means #

mean_all_functional_resolutions <- varPart.mean %>% 
  group_by(parameter_clean_name,time) %>%
  summarize(mean_vMean = mean(varMean),
            sd_vMean = sd(varMean))

write.table(mean_all_functional_resolutions,file = "outputs/variancePartition_hmn3.6/variancePartition_mean_all_functionalresolution.txt", quote = F, sep = "\t", row.names = F)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Making the same plot to check the effect of individual mice and time         #

# Loading variance partition nc_individual - non-colinear variables #
#:::::::: All Timepoints ::::::::#
# EC #
varPart.mean.all.EC <- read.table("outputs/variancePartition_hmn3.6/all/nc_individual/variance_partition_mean_EC.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "EC") %>%
  mutate(time = "All Timepoints")

# GO 
varPart.mean.all.GO <- read.table("outputs/variancePartition_hmn3.6/all/nc_individual/variance_partition_mean_GO.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "GO") %>%
  mutate(time = "All Timepoints")
# KO
varPart.mean.all.KO <- read.table("outputs/variancePartition_hmn3.6/all/nc_individual/variance_partition_mean_KO.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "KO") %>%
  mutate(time = "All Timepoints")
# RXN
varPart.mean.all.RXN <- read.table("outputs/variancePartition_hmn3.6/all/nc_individual/variance_partition_mean_RXN.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "RXN") %>%
  mutate(time = "All Timepoints")
# MCyc
varPart.mean.all.MCyc <- read.table("outputs/variancePartition_hmn3.6/all/nc_individual/variance_partition_mean_MCyc.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "MCyc") %>%
  mutate(time = "All Timepoints")

#--------------------------------------------------------------------------------------------------------------#

varPart.mean.all_individual <- rbind(varPart.mean.all.EC,
                                     varPart.mean.all.GO,
                                     varPart.mean.all.KO,
                                     varPart.mean.all.RXN,
                                     varPart.mean.all.MCyc) %>%
  mutate(parameter_clean_name = ifelse(Parameter == "timepoint","Time",
                                       ifelse(Parameter == "mouse_id","Individual mice",
                                              ifelse(Parameter == "female_cage_of_birth","Maternal cage of birth",
                                                     ifelse(Parameter == "male_cage_of_birth","Paternal cage of birth",
                                                            ifelse(Parameter == "male_genotype","Paternal genotype",
                                                                   ifelse(Parameter == "number_of_pups_scaled", "Number of pups",
                                                                          ifelse(Parameter == "lipocalin2_pg_ml_scaled","Lipocalin (pg/ml)",
                                                                                 ifelse(Parameter == "genotype","Maternal genotype",
                                                                                        ifelse(Parameter == "delivery_time","Delivery time",
                                                                                               ifelse(Parameter == "weight_grams_scaled","Maternal weight (g)","Residuals")))))))))))


#parameters_color <- pal_jco("default", alpha = 1)(nrow(varPart.mean.all.MCyc))
palette_color = c("Individual mice"="#8c564b",
                  "Maternal genotype"="#9e0142",
                  "Maternal cage of birth"="#d53e4f",
                  "Maternal weight (g)"="#f46d43",
                  "Delivery time"="#fdae61",
                  "Number of pups"="#fee08b",
                  "Lipocalin (pg/ml)"="#abdda4",
                  "Paternal cage of birth"="#66c2a5",
                  "Paternal genotype" = "#3288bd",
                  "Time" = "#5e4fa2",
                  "Residuals"="darkgrey")


# Stacked Area Plot #

varPart.mean.all_individual$level = factor(varPart.mean.all_individual$level, levels = c("EC","GO","KO",
                                                                                         "RXN","MCyc"),
                                           ordered = TRUE)
varPart.mean.all_individual$Parameter = factor(varPart.mean.all_individual$Parameter, levels = c("Residuals","female_cage_of_birth","genotype",
                                                                                                 "weight_scaled","number_of_pups_scaled","lipocalin2_pg_ml_scaled",
                                                                                                 "delivery_time","male_cage_of_birth","male_genotype","timepoint","mouse_id"
), ordered = TRUE)

varPart.mean.all_individual %>%
  arrange(varMean) %>%
  ggbarplot(x = "level", y = "varMean", fill = "parameter_clean_name", add.params = list(size=0,fill="parameter_clean_name",alpha=1),
            ylab = "Mean Variance Explained (%)", xlab = "Functional Resolution", label = F) +
  theme_bw() + 
  theme(axis.text=element_text(size=7, color = NULL), 
        axis.title=element_text(size=7), 
        legend.text = element_text(size=7),
        legend.title = element_text(size=7), 
        plot.title = element_text(size=7, face="bold", hjust = 0.5),
        strip.background=element_rect(colour="black",
                                      fill="white")) +
  labs(fill = "Variables") + rotate_x_text(angle = 45) + 
  facet_grid(cols = vars(time)) + 
  scale_fill_manual(values = palette_color)

#ggsave("figures/mean_variance_explained_nc_individual.pdf", width= 13, height=8, units = "cm")

#write.table(varPart.mean.all_individual,file = "outputs/variancePartition_hmn3.6/variancePartition_no_dt_mean_nc_individual.txt", quote = F, sep = "\t", row.names = F)

# Computing mean +sd from variance_partition_means #

mean_all_functional_resolutions_individual <- varPart.mean.all_individual %>% 
  group_by(parameter_clean_name,time) %>%
  summarize(mean_vMean = mean(varMean),
            sd_vMean = sd(varMean))

#write.table(mean_all_functional_resolutions_individual,file = "outputs/variancePartition_hmn3.6/variancePartition_no_dt_mean_all_functionalresolution_individual.txt", quote = F, sep = "\t", row.names = F)


#-------------------------------------------------------------------------------
#::::::::::::: Plot Variance Partition Features at week 3 ::::::::::::::: #
library(data.table)
library(magrittr)
library(phyloseq)
library(rstatix)
#------------------------------------------------------------------------------
# MCyc
    varPart.MCyc <- read.table(glue("outputs/variancePartition_hmn3.6/w3/nc_no_individual/variance_partition_results_MCyc.txt"), sep = "\t", header = T) %>%
      cbind(data.frame(feature = rownames(.)),.) %>%
      mutate(level = "MCyc") %>%
      data.table()
  
    top20_MCyc <- varPart.MCyc %>% 
        filter(feature != "PWY-241|C4 photosynthetic carbon assimilation cycle, NADP-ME type" & 
                 feature != "PWY-7117|C4 photosynthetic carbon assimilation cycle, PEPCK type" &
                 feature != "PWY-5723|Rubisco shunt") %>%
        select(feature,genotype) %>%
        top_n(20) %>%
        select(feature)
      
      # Variance partition long table at MCyc level
      varPart.MCyc.long <- varPart.MCyc %>%
        left_join(top20_MCyc,.,by = "feature") %>%
        pivot_longer(.,cols = colnames(select_if(., is.numeric)),names_to = "covariates", values_to = "variance") %>%
        mutate(covariates_clean = ifelse(covariates == "timepoint","Time",
                                         ifelse(covariates == "female_cage_of_birth","Maternal cage of birth",
                                                ifelse(covariates == "male_cage_of_birth","Paternal cage of birth",
                                                       ifelse(covariates == "male_genotype","Paternal genotype",
                                                              ifelse(covariates == "number_of_pups_scaled", "Number of pups",
                                                                     ifelse(covariates == "lipocalin2_pg_ml_scaled","Lipocalin (pg/ml)",
                                                                            ifelse(covariates == "genotype","Maternal genotype",
                                                                                   ifelse(covariates == "delivery_time","Delivery time",
                                                                                          ifelse(covariates == "weight_grams_scaled","Maternal weight (g)","Residuals"))))))))))
      
      
      # Color Palette
      parameters_color <- c("Maternal genotype"="#9e0142",
                            "Maternal cage of birth"="#d53e4f",
                            "Maternal weight (g)"="#f46d43",
                            "Delivery time"="#fdae61",
                            "Number of pups"="#fee08b",
                            "Lipocalin (pg/ml)"="#abdda4",
                            "Paternal cage of birth"="#66c2a5",
                            "Paternal genotype" = "#3288bd",
                            "Time" = "#5e4fa2",
                            "Residuals"="darkgrey")
      
      
      # Stacked Area Plot #
      # MCyc Features Plot #
      
      o <- varPart.MCyc.long %>% filter(covariates == "Residuals") %>% arrange(variance) %>% extract2("feature")
      
      varPart.MCyc.long %>%
        mutate(ValueName = factor(feature, o)) %>%
        ggplot(aes(fill=covariates_clean, y=variance, x=reorder(feature,-variance))) +
        geom_bar(position="stack", stat="identity") + theme_bw() +
        scale_fill_manual(values = parameters_color,
                          name = "Covariates") +
        labs( x = "MCyc", y = "Mean Variance Explained") + 
        theme_bw() + theme(axis.text=element_text(size=7, color = "black"), 
                           axis.title=element_text(size=9), 
                           legend.text = element_text(size=7), 
                           legend.title = element_text(size=0), 
                           plot.title = element_text(size=9, face="bold", hjust = 0.5)) +
        rotate_x_text(angle = 45) + coord_flip()
      
      ggsave("figures/Fig3B_top20_MCyc_features_week3_genotype.pdf", 
             width= 18, height=8, units = "cm")
      
 
  #----------------------------------- 
  # Plots top20 important features #

      # Loading Metadata #
      preg_metadata <- read.table("metaInfo/metainformation_16S_and_shotgun_metagenomics.txt", check.names = F, header = T, sep = "\t") %>% 
        filter(!is.na(barcode_shotgun)) %>%
        mutate(genotype = factor(genotype,levels = c("WT","KO"), ordered = T)) %>%
        arrange(timepoint,genotype)
      
      
      # Loading Abundance table #
      MCyc_ids <- read_tsv("inputs/humann_merged_pathabundance-cpm-no-taxa-metacyc-id.tsv") %>% 
        rename(id=metacyc_id,
               name=pathway) %>%
        filter(!is.na(name))
      
      sample_names <- substr(colnames(MCyc_ids),1,16)[-c(1:2)]
      
      MCyc_ids <- MCyc_ids %>%
        select(id,name) %>%
        mutate(complete_name = paste0(id,"|",name))
      
      MCyc_data <- read_tsv("inputs/humann_merged_pathabundance-cpm-no-taxa-metacyc-id.tsv") %>% 
        rename(id=metacyc_id,
               name=pathway) %>%
        filter(!is.na(name)) %>%
        column_to_rownames(var="id") %>%
        select(-name) %>%
        t(.)
      
      rownames(MCyc_data) <- sample_names
      
      MCyc_data <- MCyc_data[preg_metadata$barcode_shotgun,]
      
      colnames(MCyc_data) <- MCyc_ids$complete_name
      
      # Loading only significant Top20 metacyc features for genotype-variance partition at T3 #
      varPart.MCyc <- read.table("outputs/variancePartition_hmn3.6/w3/nc_no_individual/variance_partition_results_MCyc.txt", sep = "\t", header = T) %>%
        cbind(data.frame(feature = rownames(.)),.) %>%
        mutate(level = "MCyc") %>%
        data.table()
      
      top20_MCyc_genotype <- varPart.MCyc %>%
        select(feature,genotype) %>%
        filter(feature != "PWY-241|C4 photosynthetic carbon assimilation cycle, NADP-ME type" & 
                 feature != "PWY-7117|C4 photosynthetic carbon assimilation cycle, PEPCK type" &
                 feature != "PWY-5723|Rubisco shunt") %>%
        top_n(20)
      
      
      #selected_features <- top20_MCyc_genotype$feature
      selected_features <- c("POLYAMSYN-PWY|superpathway of polyamine biosynthesis I",
                             "PWY-6902|chitin degradation II (Vibrio)")
      
            # Pathway Annotation Table with Matched Order of samples metadata #
      top20_selected_pathways_abundances <- 
        cbind(preg_metadata,MCyc_data) %>%
        select(barcode_shotgun,timepoint,genotype,selected_features) %>%
        pivot_longer(.,cols=4:ncol(.), names_to = "pathways",values_to = "CPM")
      
      
      dir.create("figures/Fig3C_3D")
      
      # Computing mean of CPM abundances #
      
      for (i in 1:length(selected_features)){
        
        # Boxplots with statistics #
        
        p <-
          ggboxplot(top20_selected_pathways_abundances %>% filter(pathways == glue("{selected_features[i]}")),
                    x = "genotype", y = "CPM", color = "genotype",fill = "genotype", 
                    palette = c("WT"="#5da1d4ff", "KO"="#c88cb2ff"), ylab = "CPM", xlab = "Genotype/Timepoint",
                    facet.by = "timepoint", scales = "fixed", alpha = 0.5,
                    add = c("jitter")) +
          theme_bw() + 
          theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=7), legend.text = element_text(size=7),
                legend.title = element_text(size=7), plot.title = element_text(size=7, face="bold", hjust = 0.5),
                strip.background =element_rect(fill="white")) +
          ggtitle(glue("{selected_features[i]}")) + 
          geom_pwc(
            aes(group = genotype), tip.length = 0,
            method = "wilcox_test", label = "{p.adj.format}{p.adj.signif}",
            p.adjust.method = "BH", p.adjust.by = "group",
            hide.ns = FALSE, size = 0.3, pch = 0.2,   label.size = 2.5
          )
        
        plot(p)
        
        ggsave(file = glue("figures/Fig3C_3D/{selected_features[i]}_metacyc_pathway.pdf"), 
               width= 10, height=7, units = "cm")
        
      }
      
      
  #----------------------------------------------------------------------------     
    # KeggO    
        
      varPart.KeggO <- read.table(glue("outputs/variancePartition_hmn3.6/w3/nc_no_individual/variance_partition_results_KO.txt"), sep = "\t", header = T) %>%
            cbind(data.frame(feature = rownames(.)),.) %>%
            mutate(level = "KeggO") %>%
            data.table()
          
            
        top20_KeggO <- varPart.KeggO %>%
              filter(feature != "K07792|anaerobic C4-dicarboxylate transporter DcuB" & 
                       feature != "K11690|C4-dicarboxylate transporter, DctM subunit") %>%
              select(feature,genotype) %>%
              top_n(20) %>%
              select(feature)
            
      # Variance partition long table at KeggO level
        varPart.KeggO.long <- varPart.KeggO %>%
              left_join(top20_KeggO,.,by = "feature") %>%
              #arrange(.,covariates_to_plot[n]) %>%
              pivot_longer(.,cols = colnames(select_if(., is.numeric)),names_to = "covariates", values_to = "variance") %>%
              mutate(covariates_clean = ifelse(covariates == "timepoint","Time",
                                               ifelse(covariates == "female_cage_of_birth","Maternal cage of birth",
                                                      ifelse(covariates == "male_cage_of_birth","Paternal cage of birth",
                                                             ifelse(covariates == "male_genotype","Paternal genotype",
                                                                    ifelse(covariates == "number_of_pups_scaled", "Number of pups",
                                                                           ifelse(covariates == "lipocalin2_pg_ml_scaled","Lipocalin (pg/ml)",
                                                                                  ifelse(covariates == "genotype","Maternal genotype",
                                                                                         ifelse(covariates == "delivery_time","Delivery time",
                                                                                                ifelse(covariates == "weight_grams_scaled","Maternal weight (g)","Residuals"))))))))))
            
            
            # Color Palette
            parameters_color <- c("Maternal genotype"="#9e0142",
                                  "Maternal cage of birth"="#d53e4f",
                                  "Maternal weight (g)"="#f46d43",
                                  "Delivery time"="#fdae61",
                                  "Number of pups"="#fee08b",
                                  "Lipocalin (pg/ml)"="#abdda4",
                                  "Paternal cage of birth"="#66c2a5",
                                  "Paternal genotype" = "#3288bd",
                                  "Time" = "#5e4fa2",
                                  "Residuals"="darkgrey")
            
            
            # Stacked Area Plot #
            # KeggO Features Plot #
            
            o <- varPart.KeggO.long %>% filter(covariates == "Residuals") %>% arrange(variance) %>% extract2("feature")
            
            varPart.KeggO.long %>%
              mutate(ValueName = factor(feature, o)) %>%
              ggplot(aes(fill=covariates_clean, y=variance, x=reorder(feature,-variance))) +
              geom_bar(position="stack", stat="identity") + theme_bw() +
              scale_fill_manual(values = parameters_color,
                                name = "Covariates") +
              labs( x = "KO", y = "Mean Variance Explained") + 
              theme_bw() + theme(axis.text=element_text(size=7, color = "black"), 
                                 axis.title=element_text(size=9), 
                                 legend.text = element_text(size=7), 
                                 legend.title = element_text(size=0), 
                                 plot.title = element_text(size=9, face="bold", hjust = 0.5)) +
              rotate_x_text(angle = 45) + coord_flip()
            
            ggsave("figures/Fig3E_top20_KeggO_features_week3_genotype.pdf", 
                   width= 18, height=8, units = "cm")     
 
    #----------------------------------- 
    # Plots top20 important features #
            
            # Loading Abundance table #
            KO_ids <- read_tsv("inputs/humann_merged_genefamilies-cpm-ko-named-no-taxa.tsv") %>% 
              filter(genefamily != "UNMAPPED" & genefamily != "UNGROUPED") %>%
              separate(.,col=genefamily,into=c("id","name"), sep = ": ") %>%
              filter(name != "NO_NAME")
            
            sample_names <- substr(colnames(KO_ids),1,16)[-c(1:2)]
            
            KO_ids <- KO_ids %>%
              select(id,name) %>%
              mutate(complete_name = paste0(id,"|",name))
            
            KO_data <- read_tsv("inputs/humann_merged_genefamilies-cpm-ko-named-no-taxa.tsv") %>% 
              filter(genefamily != "UNMAPPED" & genefamily != "UNGROUPED") %>%
              separate(.,col=genefamily,into=c("id","name"), sep = ": ") %>%
              filter(name != "NO_NAME") %>% 
              column_to_rownames(var="id") %>%
              select(-name) %>%
              t(.)
            
            rownames(KO_data) <- sample_names
            
            KO_data <- KO_data[preg_metadata$barcode_shotgun,]
            
            colnames(KO_data) <- KO_ids$complete_name
            
            # Loading only significant Top20 metacyc features for genotype-variance partition at T3 #
            varPart.KO <- read.table("outputs/variancePartition_hmn3.6/w3/nc_no_individual/variance_partition_results_KO.txt", sep = "\t", header = T) %>%
              cbind(data.frame(feature = rownames(.)),.) %>%
              mutate(level = "KO") %>%
              data.table()
            
            top20_KeggO_genotype <- varPart.KO %>%
              select(feature,genotype) %>%
              filter(feature != "K07792|anaerobic C4-dicarboxylate transporter DcuB" & 
                       feature != "K11690|C4-dicarboxylate transporter, DctM subunit") %>%
              top_n(20)
            
            #selected_features <- top20_KeggO_genotype$feature
            selected_features <- c("K02377|GDP-L-fucose synthase [EC.1.1.1.271]")
            
            # Pathway Annotation Table with Matched Order of samples metadata #
            top20_selected_KeggO_abundances <- 
              cbind(preg_metadata,KO_data) %>%
              select(barcode_shotgun,timepoint,genotype,selected_features) %>%
              pivot_longer(.,cols=4:ncol(.), names_to = "KeggO",values_to = "CPM")
            
            
            dir.create("figures/Fig3F")
            
            
            # Computing mean of CPM abundances #
            
            for (i in 1:length(selected_features)){
              
              # Boxplots with statistics #
              
              p <-
                ggboxplot(top20_selected_KeggO_abundances %>% filter(KeggO == glue("{selected_features[i]}")),
                          x = "genotype", y = "CPM", color = "genotype",fill = "genotype", 
                          palette = c("WT"="#5da1d4ff", "KO"="#c88cb2ff"), ylab = "CPM", xlab = "Genotype/Timepoint",
                          facet.by = "timepoint", scales = "fixed", alpha = 0.5,
                          add = c("jitter")) +
                theme_bw() + 
                theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=7), legend.text = element_text(size=7),
                      legend.title = element_text(size=7), plot.title = element_text(size=7, face="bold", hjust = 0.5),
                      strip.background =element_rect(fill="white")) +
                ggtitle(glue("{selected_features[i]}")) + 
                geom_pwc(
                  aes(group = genotype), tip.length = 0,
                  method = "wilcox_test", label = "{p.adj.format}{p.adj.signif}",
                  p.adjust.method = "BH", p.adjust.by = "group",
                  hide.ns = FALSE, size = 0.3, pch = 0.2,   label.size = 2.5
                )
              
              plot(p)
              
              ggsave(file = glue("figures/Fig3F/{selected_features[i]}_KeggO_genefamily.pdf"), 
                     width= 10, height=7, units = "cm")
              
            }
            