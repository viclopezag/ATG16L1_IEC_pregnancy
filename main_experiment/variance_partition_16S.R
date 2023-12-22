
library(dplyr)
library(tidyverse)
library(variancePartition)
library(reshape2)
library(glue)
library(microbiome)
library(EnvStats)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(phyloseq)
library(rstatix)

# Code for - Variance Partition Analyses: 16S data -------

# Setting Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading Phyloseq Object #

ps <- readRDS("inputs/phyloseq_object.rds")

# Cutoff #

(ps_fr <- subset_samples(ps, sample_sums(ps) > 5000))

# Deleting taxa with relative abundance lowest than 1e-5 ()
minTotRelAbun <- 1e-5          
x <- taxa_sums(ps_fr)
keepTaxa <- (x / sum(x)) > minTotRelAbun
ps_f <- prune_taxa(keepTaxa, ps_fr)

# Variance Partition Timepoint #

vp_time <- c("all","BL","w3","w6")

for (u in 1:length(vp_time)) {
# Loading metadata #

metadata <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  sample_data() %>%
  data.frame() %>%
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

# Centered-log ration transformations #

# We will delete Unknown genus #

# ASVs #
asv_id <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_table(.) %>%
  data.frame(.) %>%
  cbind(data.frame(ASVs = rownames(.)),.) %>%
  mutate(complete_name = paste0(ASVs,"|",Phylum,"|",Family,"|",Species))

asv_data <- ps_f %>% 
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  microbiome::transform(., transform = "clr") %>%
  otu_table(.) %>%
  data.frame(.)

colnames(asv_data) <- asv_id$complete_name

# Species #
species_id <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Species") %>%
  tax_table(.) %>%
  data.frame(.) %>%
  cbind(data.frame(ASVs = rownames(.)),.) %>%
  mutate(complete_name = paste0(ASVs,"|",Phylum,"|",Family,"|",Species))

species_data <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Species") %>%
  microbiome::transform(., transform = "clr") %>%
  otu_table(.) %>%
  data.frame(.)

colnames(species_data) <- species_id$complete_name

# Genus #
genus_id <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Genus") %>%
  tax_table(.) %>%
  data.frame(.) %>%
  cbind(data.frame(ASVs = rownames(.)),.) %>%
  mutate(complete_name = paste0(ASVs,"|",Phylum,"|",Family,"|",Genus))

genus_data <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Genus") %>%
  microbiome::transform(., transform = "clr") %>%
  otu_table(.) %>%
  data.frame(.)

colnames(genus_data) <- genus_id$complete_name                       

# Family # 
family_id <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Family") %>%
  tax_table(.) %>%
  data.frame(.) %>%
  cbind(data.frame(ASVs = rownames(.)),.) %>%
  mutate(complete_name = paste0(ASVs,"|",Phylum,"|",Family))

family_data <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Family") %>%
  microbiome::transform(., transform = "clr") %>%
  otu_table(.) %>%
  data.frame(.)
  
colnames(family_data) <- family_id$complete_name 

# Order #
order_id <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Order") %>%
  tax_table(.) %>%
  data.frame(.) %>%
  cbind(data.frame(ASVs = rownames(.)),.) %>%
  mutate(complete_name = paste0(ASVs,"|",Phylum,"|",Order))

order_data <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Order") %>%
  microbiome::transform(., transform = "clr") %>%
  otu_table(.) %>%
  data.frame(.)

colnames(order_data) <- order_id$complete_name           

# Class #
class_id <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Class") %>%
  tax_table(.) %>%
  data.frame(.) %>%
  cbind(data.frame(ASVs = rownames(.)),.) %>%
  mutate(complete_name = paste0(ASVs,"|",Phylum,"|",Class))

class_data <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Class") %>%
  microbiome::transform(., transform = "clr") %>%
  otu_table(.) %>%
  data.frame(.)

colnames(class_data) <- class_id$complete_name                  

# Phylum #
phylum_id <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Phylum") %>%
  tax_table(.) %>%
  data.frame(.) %>%
  cbind(data.frame(ASVs = rownames(.)),.) %>%
  mutate(complete_name = paste0(Phylum))

phylum_data <- ps_f %>%
  subset_samples(if (vp_time[u] == "all") timepoint == "BL" | timepoint == "w3" | timepoint == "w6" else timepoint == vp_time[u]) %>%
  tax_glom(.,taxrank = "Phylum") %>%
  microbiome::transform(., transform = "clr") %>%
  otu_table(.) %>%
  data.frame(.)

colnames(phylum_data) <- phylum_id$complete_name                           

# Creating Directories for saving Variance Partition Outputs #
dir.create("outputs/variancePartition_16S")

dir.create(glue("outputs/variancePartition_16S/{vp_time[u]}"))

#::::::::::::::::::::::::::::::::::::::::::::::::#
# Identifying Co-linear Variables in the Metadata #
if (vp_time[u] == "all") {
formula <- ~ mouse_id + timepoint + genotype + weight_grams_scaled + 
  number_of_pups_scaled + lipocalin2_pg_ml_scaled

# Canonical Correlation Analysis #
C = canCorPairs(formula, metadata, showWarnings = TRUE)

#Correcting names 
rownames(C) <- c("Individual mice","Weeks","Maternal genotype","Maternal weight (g)", "Number of pups",
                 "Lipocalin (pg/ml)")
colnames(C) <- c("Individual mice","Weeks","Maternal genotype", "Maternal weight (g)", "Number of pups",
                 "Lipocalin (pg/ml)")
} else {

formula <- ~ mouse_id  + genotype + 
    female_cage_of_birth + male_cage_of_birth  + 
    male_genotype + weight_grams_scaled + number_of_pups_scaled + lipocalin2_pg_ml_scaled
  
  # Canonical Correlation Analysis #
  C = canCorPairs(formula, metadata, showWarnings = TRUE)
  
  #Correcting names 
  rownames(C) <- c("Individual mice","Maternal genotype", "Maternal cage of birth",
                   "Paternal cage of birth", "Paternal genotype", "Maternal weight (g)", "Number of pups",
                   "Lipocalin (pg/ml)")
  colnames(C) <- c("Individual mice","Maternal genotype", "Maternal cage of birth",
                   "Paternal cage of birth", "Paternal genotype", "Maternal weight (g)", "Number of pups",
                   "Lipocalin (pg/ml)")  
  
}  
# Make upper triangular matrix by setting NA to lower triangular part:

#C[lower.tri(C)] <- NA

# A Full Correlation Plot Using ggplot2
melted_C <- melt(C)

## Version One: Correlation Plot using ggplot2:
# Generate Color Gradient #
red_grey_gradient <- colorRampPalette(colors = c("lightgrey","darkred"))
n <- 4
red_grey_colors <- red_grey_gradient(n)

ggplot(data = melted_C, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = red_grey_colors, 
                       limits = c(0.0, 1.0)) +
  labs(title = "Canonical Correlation", 
       x = "", y = "", fill = "Correlation \n Measure") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, colour = "black"), 
        axis.title.x = element_text(face="bold", colour="black", size = 7),
        axis.title.y = element_text(face="bold", colour="black", size = 7),
        legend.title = element_text(face="bold", colour="black", size = 7),
        axis.text.x = element_text(face="bold",colour="black",size = 7),
        axis.text.y = element_text(face="bold",colour="black",size = 7)) +
  geom_text(aes(x = Var1, y = Var2, label = round(value, 2)), color = "black", 
            fontface = "bold", size = 4) + scale_x_discrete(guide = guide_axis(angle = 45))

if (vp_time[u] == "all"){
ggsave(file = glue("supp_figures/Supp_Fig3A_canonical_correlation_mice_covariates.pdf"), 
       width= 15, height=10, units = "cm")
}

# Colinear variables were selected if correlation coefficient is at least 0.8 #


if (vp_time[u] == "all") {
formula_no_colinear <- list(formula_no_individual = ~ (1|genotype) + (1|timepoint) + number_of_pups_scaled + 
                                                      lipocalin2_pg_ml_scaled,
                            formula_individual = ~ (1|mouse_id) + (1|timepoint))

} else {
  
formula_no_colinear <- list(formula_no_individual = ~ (1|genotype) + weight_grams_scaled + number_of_pups_scaled + 
                                                      lipocalin2_pg_ml_scaled)
  
}  
                            
table_metadata_important <- metadata


level_abundances = list(asv_data,
                        species_data,
                        genus_data,
                        family_data,
                        order_data,
                        class_data,
                        phylum_data)


level_annotation = c("ASVs","Species","Genus","Family","Order","Class","Phylum")

#::::::::::::::::::::::::Variance Partition in all Taxa Levels :::::::::::::::::::::::::::::#
# Two Formulas will be evaluated 1) those variables non-colinear with group_g and 2) those variables non-colinear with R_NR #

non_colinear_formulas <- c("nc_no_individual","nc_individual")

set.seed(1989)

for (v in 1:length(formula_no_colinear)) {
  
  dir.create(glue("outputs/variancePartition_16S/{vp_time[u]}/{non_colinear_formulas[v]}")) 
  
  for (q in 1:length(level_annotation)){
    
    matrix_data <- level_abundances[[q]]

  # Delete Features not seen in at least 12% of the samples (80% of the samples per time/genotype group) 11 mouse/86 mouse = ~ 12% total samples , 
  # there are in average 15 mouse per time*genotype groups #
    
  if (vp_time[u]=="all"){
    # 80% rule : All timepoints : 11 samples/86 samples = ~ 12%
    matrix_data <- matrix_data[,colSums(matrix_data > 0) > 0.12*nrow(matrix_data)]
  } else if (vp_time[u]=="BL") {
    # 80% rule : BL : 14*0.8/29 mouse = ~ 0.38
    matrix_data <- matrix_data[,colSums(matrix_data > 0) > 0.38*nrow(matrix_data)]
    
  } else if (vp_time[u]=="w3") {
    # 80% rule : w3 : 13*0.8/28 mouse = 0.37
    matrix_data <- matrix_data[,colSums(matrix_data > 0) > 0.37*nrow(matrix_data)]
    
  } else {
    # 80% rule : w6 : 14*0.8/29 mouse ~ 0.38
    matrix_data <- matrix_data[,colSums(matrix_data > 0) > 0.38*nrow(matrix_data)] 
  }  
    # We are going to run Variance Partition on squared rooted transformed abundances #
    
    varPart <- fitExtractVarPartModel(t(matrix_data), formula_no_colinear[[v]], table_metadata_important, REML = F)
    
    # sort variables by mean fraction of varience explained
    
    vp <- sortCols(varPart)
    
    # Export results # 
    
    varMean <- colMeans(varPart)
    varPart_sorted <- varPart[, order(varMean, decreasing = TRUE)]
    
    write.table(varPart_sorted, glue("outputs/variancePartition_16S/{vp_time[u]}/{non_colinear_formulas[v]}/variance_partition_results_{level_annotation[q]}.txt"), sep = '\t', quote = FALSE)
    
    varMean <- as.data.frame(varMean)
    varMean$Parameter <- rownames(varMean)
    varMean <- varMean[order(varMean$varMean, decreasing = TRUE), ]
    
    write.table(varMean, glue("outputs/variancePartition_16S/{vp_time[u]}/{non_colinear_formulas[v]}/variance_partition_mean_{level_annotation[q]}.txt"), sep = '\t', quote = FALSE, row.names = FALSE)
    
    #### Violin plot of contribution of each variable to total variance ####
    
    varPart <- read.table(glue("outputs/variancePartition_16S/{vp_time[u]}/{non_colinear_formulas[v]}/variance_partition_results_{level_annotation[q]}.txt"), sep = "\t") %>%
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
    
    if (level_annotation[q] == "Genus"){
    
      if(vp_time[u] == "all" & non_colinear_formulas[v] == "nc_no_individual"){
        pdf(glue("supp_figures/Supp_Fig3B_all_no_individual_violin_variance_partition_{level_annotation[q]}.pdf"), width = 7, height = 7)
    
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
    
          dev.off()
      }
      
      if(vp_time[u] == "all" & non_colinear_formulas[v] == "nc_individual"){
        pdf(glue("supp_figures/Supp_Fig3B_all_individual_violin_variance_partition_{level_annotation[q]}.pdf"), width = 7, height = 7)
        
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
        
        dev.off()
      }
      
      if(vp_time[u] == "BL"){
        pdf(glue("supp_figures/Supp_Fig3C_BL_violin_variance_partition_{level_annotation[q]}.pdf"), width = 7, height = 7)
        
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
        
        dev.off()
      }
      
      if(vp_time[u] == "w3"){
        pdf(glue("supp_figures/Supp_Fig3D_w3_violin_variance_partition_{level_annotation[q]}.pdf"), width = 7, height = 7)
        
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
        
        dev.off()
      }
      
      if(vp_time[u] == "w6"){
        pdf(glue("supp_figures/Supp_Fig3E_all_violin_variance_partition_{level_annotation[q]}.pdf"), width = 7, height = 7)
        
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
        
        dev.off()
      }
      
    } 
    
   }
}

}


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Stacked Area Plot - Variance partition - 16S data #

# Loading variance partition nc_no_individual - non-colinear variables #
#:::::::: All Timepoints ::::::::#
# ASVs #
varPart.mean.all.ASVs <- read.table("outputs/variancePartition_16S/all/nc_no_individual/variance_partition_mean_ASVs.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "ASVs") %>%
  mutate(time = "All Timepoints")

# Species 
varPart.mean.all.Species <- read.table("outputs/variancePartition_16S/all/nc_no_individual/variance_partition_mean_Species.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "Species") %>%
  mutate(time = "All Timepoints")
# Genus
varPart.mean.all.Genus <- read.table("outputs/variancePartition_16S/all/nc_no_individual/variance_partition_mean_Genus.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Genus") %>%
  mutate(time = "All Timepoints")
# Family
varPart.mean.all.Family <- read.table("outputs/variancePartition_16S/all/nc_no_individual/variance_partition_mean_Family.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Family") %>%
  mutate(time = "All Timepoints")
# Order
varPart.mean.all.Order <- read.table("outputs/variancePartition_16S/all/nc_no_individual/variance_partition_mean_Order.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Order") %>%
  mutate(time = "All Timepoints")
# Class
varPart.mean.all.Class <- read.table("outputs/variancePartition_16S/all/nc_no_individual/variance_partition_mean_Class.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Class") %>%
  mutate(time = "All Timepoints")

# Phylum
varPart.mean.all.Phylum <- read.table("outputs/variancePartition_16S/all/nc_no_individual/variance_partition_mean_Phylum.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Phylum") %>%
  mutate(time = "All Timepoints")
#--------------------------------------------------------------------------------------------------------------#

#:::::::: BL: Baseline  ::::::::#
# ASVs 
varPart.mean.BL.ASVs <- read.table("outputs/variancePartition_16S/BL/nc_no_individual/variance_partition_mean_ASVs.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "ASVs") %>%
  mutate(time = "BL")

# Species 
varPart.mean.BL.Species <- read.table("outputs/variancePartition_16S/BL/nc_no_individual/variance_partition_mean_Species.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "Species") %>%
  mutate(time = "BL")
# Genus
varPart.mean.BL.Genus <- read.table("outputs/variancePartition_16S/BL/nc_no_individual/variance_partition_mean_Genus.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Genus") %>%
  mutate(time = "BL")
# Family
varPart.mean.BL.Family <- read.table("outputs/variancePartition_16S/BL/nc_no_individual/variance_partition_mean_Family.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Family") %>%
  mutate(time = "BL")
# Order
varPart.mean.BL.Order <- read.table("outputs/variancePartition_16S/BL/nc_no_individual/variance_partition_mean_Order.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Order") %>%
  mutate(time = "BL")
# Class
varPart.mean.BL.Class <- read.table("outputs/variancePartition_16S/BL/nc_no_individual/variance_partition_mean_Class.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Class") %>%
  mutate(time = "BL")

# Phylum
varPart.mean.BL.Phylum <- read.table("outputs/variancePartition_16S/BL/nc_no_individual/variance_partition_mean_Phylum.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Phylum") %>%
  mutate(time = "BL")
#--------------------------------------------------------------------------------------------------------------#
#::::::::: w3 :::::::::::::#
# ASVs
varPart.mean.w3.ASVs <- read.table("outputs/variancePartition_16S/w3/nc_no_individual/variance_partition_mean_ASVs.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "ASVs") %>%
  mutate(time = "w3")
# Species 
varPart.mean.w3.Species <- read.table("outputs/variancePartition_16S/w3/nc_no_individual/variance_partition_mean_Species.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "Species") %>%
  mutate(time = "w3")
# Genus
varPart.mean.w3.Genus <- read.table("outputs/variancePartition_16S/w3/nc_no_individual/variance_partition_mean_Genus.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Genus") %>%
  mutate(time = "w3")
# Family
varPart.mean.w3.Family <- read.table("outputs/variancePartition_16S/w3/nc_no_individual/variance_partition_mean_Family.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Family") %>%
  mutate(time = "w3")
# Order
varPart.mean.w3.Order <- read.table("outputs/variancePartition_16S/w3/nc_no_individual/variance_partition_mean_Order.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Order") %>%
  mutate(time = "w3")
# Class
varPart.mean.w3.Class <- read.table("outputs/variancePartition_16S/w3/nc_no_individual/variance_partition_mean_Class.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Class") %>%
  mutate(time = "w3")

# Phylum
varPart.mean.w3.Phylum <- read.table("outputs/variancePartition_16S/w3/nc_no_individual/variance_partition_mean_Phylum.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Phylum") %>%
  mutate(time = "w3")

#--------------------------------------------------------------------------------------------------------------#
#:::::::::::::: w6 ::::::::::::::::#
# ASVs 
varPart.mean.w6.ASVs <- read.table("outputs/variancePartition_16S/w6/nc_no_individual/variance_partition_mean_ASVs.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "ASVs") %>%
  mutate(time = "w6")
# Species 
varPart.mean.w6.Species <- read.table("outputs/variancePartition_16S/w6/nc_no_individual/variance_partition_mean_Species.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "Species") %>%
  mutate(time = "w6")
# Genus
varPart.mean.w6.Genus <- read.table("outputs/variancePartition_16S/w6/nc_no_individual/variance_partition_mean_Genus.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Genus") %>%
  mutate(time = "w6")
# Family
varPart.mean.w6.Family <- read.table("outputs/variancePartition_16S/w6/nc_no_individual/variance_partition_mean_Family.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Family") %>%
  mutate(time = "w6")
# Order
varPart.mean.w6.Order <- read.table("outputs/variancePartition_16S/w6/nc_no_individual/variance_partition_mean_Order.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Order") %>%
  mutate(time = "w6")
# Class
varPart.mean.w6.Class <- read.table("outputs/variancePartition_16S/w6/nc_no_individual/variance_partition_mean_Class.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Class") %>%
  mutate(time = "w6")

# Phylum
varPart.mean.w6.Phylum <- read.table("outputs/variancePartition_16S/w6/nc_no_individual/variance_partition_mean_Phylum.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Phylum") %>%
  mutate(time = "w6")


#--------------------------------------------------------------------------------------------------------------#

varPart.mean <- rbind(varPart.mean.all.ASVs,
                      varPart.mean.all.Species,
                      varPart.mean.all.Genus,
                      varPart.mean.all.Family,
                      varPart.mean.all.Order,
                      varPart.mean.all.Class,
                      varPart.mean.all.Phylum,
                      varPart.mean.BL.ASVs,
                      varPart.mean.BL.Species,
                      varPart.mean.BL.Genus,
                      varPart.mean.BL.Family,
                      varPart.mean.BL.Order,
                      varPart.mean.BL.Class,
                      varPart.mean.BL.Phylum,
                      varPart.mean.w3.ASVs,
                      varPart.mean.w3.Species,
                      varPart.mean.w3.Genus,
                      varPart.mean.w3.Family,
                      varPart.mean.w3.Order,
                      varPart.mean.w3.Class,
                      varPart.mean.w3.Phylum,
                      varPart.mean.w6.ASVs,
                      varPart.mean.w6.Species,
                      varPart.mean.w6.Genus,
                      varPart.mean.w6.Family,
                      varPart.mean.w6.Order,
                      varPart.mean.w6.Class,
                      varPart.mean.w6.Phylum) %>%
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


#parameters_color <- pal_jco("default", alpha = 1)(nrow(varPart.mean.all.Order))
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

varPart.mean$level = factor(varPart.mean$level, levels = c("ASVs","Species","Genus",
                                                           "Family","Order","Class","Phylum"),
                            ordered = TRUE)
varPart.mean$Parameter = factor(varPart.mean$Parameter, levels = c("Residuals","female_cage_of_birth","genotype",
                                                                   "weight_grams_scaled","number_of_pups_scaled","lipocalin2_pg_ml_scaled",
                                                                   "male_cage_of_birth","male_genotype","timepoint","mouse_id"
), ordered = TRUE)

varPart.mean %>%
  arrange(varMean) %>%
  ggbarplot(x = "level", y = "varMean", fill = "parameter_clean_name", add.params = list(size=0,fill="parameter_clean_name",alpha=1),
            ylab = "Mean Variance Explained (%)", xlab = "Phylogenetic Resolution", label = F) +
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

ggsave("figures/Fig2A_mean_variance_explained.pdf", width= 13, height=8, units = "cm")

write.table(varPart.mean,file = "outputs/variancePartition_16S/mean_variancePartition_nc_no_individual.txt", quote = F, sep = "\t", row.names = F)

# Computing mean + sd from variance_partition_means #

mean_all_phylogenetic_resolutions <- varPart.mean %>% 
  group_by(time,parameter_clean_name) %>%
  summarize(mean_vMean = mean(varMean),
            sd_vMean = sd(varMean))

write.table(mean_all_phylogenetic_resolutions,file = "outputs/variancePartition_16S/variancePartition_mean_all_phyloresolution.txt", quote = F, sep = "\t", row.names = F)


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Making the same plot to check the effect of individual mice and time #

# Loading variance partition nc_no_individual - non-colinear variables #
#:::::::: All Timepoints ::::::::#
# ASVs #
varPart.mean.all.ASVs <- read.table("outputs/variancePartition_16S/all/nc_individual/variance_partition_mean_ASVs.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "ASVs") %>%
  mutate(time = "All Timepoints")

# Species 
varPart.mean.all.Species <- read.table("outputs/variancePartition_16S/all/nc_individual/variance_partition_mean_Species.txt", sep = "\t", header = T) %>%
  as_tibble() %>%
  mutate(level = "Species") %>%
  mutate(time = "All Timepoints")
# Genus
varPart.mean.all.Genus <- read.table("outputs/variancePartition_16S/all/nc_individual/variance_partition_mean_Genus.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Genus") %>%
  mutate(time = "All Timepoints")
# Family
varPart.mean.all.Family <- read.table("outputs/variancePartition_16S/all/nc_individual/variance_partition_mean_Family.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Family") %>%
  mutate(time = "All Timepoints")
# Order
varPart.mean.all.Order <- read.table("outputs/variancePartition_16S/all/nc_individual/variance_partition_mean_Order.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Order") %>%
  mutate(time = "All Timepoints")
# Class
varPart.mean.all.Class <- read.table("outputs/variancePartition_16S/all/nc_individual/variance_partition_mean_Class.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Class") %>%
  mutate(time = "All Timepoints")

# Phylum
varPart.mean.all.Phylum <- read.table("outputs/variancePartition_16S/all/nc_individual/variance_partition_mean_Phylum.txt", sep = "\t", header = T) %>%
  as_tibble()%>%
  mutate(level = "Phylum") %>%
  mutate(time = "All Timepoints")


varPart.mean.all.indv <- rbind(varPart.mean.all.ASVs,
                               varPart.mean.all.Species,
                               varPart.mean.all.Genus,
                               varPart.mean.all.Family,
                               varPart.mean.all.Order,
                               varPart.mean.all.Class,
                               varPart.mean.all.Phylum) %>%
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


#parameters_color <- pal_jco("default", alpha = 1)(nrow(varPart.mean.all.Order))
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

varPart.mean.all.indv$level = factor(varPart.mean.all.indv$level, levels = c("ASVs","Species","Genus",
                                                                             "Family","Order","Class","Phylum"),
                                     ordered = TRUE)
varPart.mean.all.indv$Parameter = factor(varPart.mean.all.indv$Parameter, levels = c("Residuals","female_cage_of_birth","genotype",
                                                                                     "weight_grams_scaled","number_of_pups_scaled","lipocalin2_pg_ml_scaled",
                                                                                     "male_cage_of_birth","male_genotype","timepoint","mouse_id"
), ordered = TRUE)

varPart.mean.all.indv %>%
  arrange(varMean) %>%
  ggbarplot(x = "level", y = "varMean", fill = "parameter_clean_name", add.params = list(size=0,fill="parameter_clean_name",alpha=1),
            ylab = "Mean Variance Explained (%)", xlab = "Phylogenetic Resolution", label = F) +
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

#ggsave("figures/Fig2A_mean_variance_explained_nc_individual.pdf", width= 5, height=8, units = "cm")

#write.table(varPart.mean.all.indv,file = "outputs/variancePartition_16S/variancePartition_mean_nc_individual.txt", quote = F, sep = "\t", row.names = F)


# Computing mean + sd from variance_partition_means #

mean_all_phylogenetic_resolutions_individual <- varPart.mean.all.indv %>% 
  group_by(time,parameter_clean_name) %>%
  summarize(mean_vMean = mean(varMean),
            sd_vMean = sd(varMean))

#write.table(mean_all_phylogenetic_resolutions_individual,file = "outputs/variancePartition_16S/variancePartition_no_dt_mean_all_phyloresolution_individual.txt", quote = F, sep = "\t", row.names = F)

#::::::::::::: Plot Variance Partition Features at week 3 ::::::::::::::: #


varPart.Genus <- read.table(glue("outputs/variancePartition_16S/w3/nc_no_individual/variance_partition_results_Genus.txt"), sep = "\t", header = T) %>%
      cbind(data.frame(feature = rownames(.)),.) %>%
      mutate(level = "Genus") %>%
      data.table()
    
  
  top20_genus <- varPart.Genus %>%
        select(feature,genotype) %>%
        top_n(20) %>%
        select(feature)
      
  # Variance partition long table at Genus level
      varPart.Genus.long <- varPart.Genus %>%
        left_join(top20_genus,.,by = "feature") %>%
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
      # Genus Features Plot #
      
      o <- varPart.Genus.long %>% filter(covariates == "Residuals") %>% arrange(variance) %>% extract2("feature")
      
      varPart.Genus.long %>%
        mutate(ValueName = factor(feature, o)) %>%
        ggplot(aes(fill=covariates_clean, y=variance, x=reorder(feature,-variance))) +
        geom_bar(position="stack", stat="identity") + theme_bw() +
        scale_fill_manual(values = parameters_color,
                          name = "Covariates") +
        labs( x = "Genus", y = "Mean Variance Explained") + 
        theme_bw() + theme(axis.text=element_text(size=7, color = "black"), 
                           axis.title=element_text(size=9), 
                           legend.text = element_text(size=7), 
                           legend.title = element_text(size=0), 
                           plot.title = element_text(size=9, face="bold", hjust = 0.5)) +
        rotate_x_text(angle = 45) + coord_flip() 
      
      ggsave("figures/Fig2B_top20_genus_features_variancePartition_week3.pdf", width= 15, height=8, units = "cm")
      
      
    