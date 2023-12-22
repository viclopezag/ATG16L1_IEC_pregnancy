
library(tidyverse)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(ggrepel)
library(vegan)
library(microbiome)
library(ggpubr)
library(dplyr)

# Code for - Stacked Area Plot: 16S data-------


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


################ genus Area Plot #################
# Extracting genus count data after rarefaction
genus_ps <- ps_f %>%
  rarefy_even_depth(., rngseed=1989) %>%
  tax_glom(., taxrank = "Genus")

# Extracting genus names from phyloseq object
genus_col_position <- 6
genus_names <- data.frame(genus_ps@tax_table[,genus_col_position])

# Extracting genus count table from phyloseq object
metadata <- data.frame(genus_ps@sam_data) %>%
  cbind(data.frame(samples=row.names(.)),.)

genus_table <- cbind(genus_names,t(data.frame(genus_ps@otu_table))) %>%
  pivot_longer(.,cols=2:ncol(.),names_to ="samples",values_to="abundances") %>%
  group_by(samples,Genus) %>%
  summarize(abundances = sum(abundances)) %>% # summing up all the "unknown" Genus
  filter(Genus != "unknown") # Deleting unknown Genus taxa.

# Identifying TopN (20) taxa #

topN_genera <- genus_table %>%
  group_by(Genus) %>%
  summarise(mean_ab = mean(abundances)) %>%
  top_n(20, mean_ab)

# Genus matrix with all "unknown" Genus summed up 

genus_matrix <- genus_table %>%
  pivot_wider(.,names_from = Genus, values_from = abundances)

others_matrix <- genus_matrix %>%
  as_tibble() %>%
  select(-topN_genera$Genus) %>%
  column_to_rownames(var="samples")

others_genera <- colnames(others_matrix)

others_sum <- others_matrix %>%
  rowSums() %>%
  as_tibble()

colnames(others_sum) <- "Others"

dir.create("outputs/stacked_area_genus_plot")
write.table(topN_genera, file = "outputs/stacked_area_genus_plot/top20_genera_main_experiment_16S.txt",sep="\t", row.names =F)

# Creating long table with metadata information #

genus_and_metadata <- genus_matrix %>%
  column_to_rownames(var="samples") %>%
  cbind(metadata,.,others_sum) %>%
  select(-others_genera) %>%
  pivot_longer(cols = A2:ncol(.), names_to = "taxa", values_to="abundances") %>%
  group_by(genotype,timepoint,taxa) %>%
  summarise(mean_ab = sum(abundances)) %>%
  group_by(genotype,timepoint) %>%
  mutate(percentage = mean_ab/ sum(mean_ab))


# Top 20 - Phylum to Genus assignments #
library(ggsci)

taxonomy_annotation <- data.frame(tax_table(genus_ps))

phylum_genus <- genus_and_metadata %>%
  left_join(.,taxonomy_annotation %>% select(Genus,Phylum), by = c("taxa"="Genus")) %>%
  mutate(Phylum = ifelse(is.na(Phylum), "Other",Phylum)) %>%
  select(genotype,timepoint,percentage,taxa,Phylum) %>%
  arrange(match(Phylum,c("Firmicutes","Bacteroidota",
                         "Proteobacteria","Actinobacteriota",
                         "Campylobacterota","Cyanobacteria",
                         "Verrucomicrobiota","Deferribacterota",
                         "Other")),desc(genotype),desc(timepoint),desc(percentage))
# arrange_at(2:6, desc) %>%
# arrange(match(Phylum, c("Firmicutes","Bacteroidota",
#                         "Actinobacteriota","Proteobacteria",
#                         "Verrucomicrobiota","Desulfobacterota",
#                         "Other"))) # https://stackoverflow.com/questions/46129322/arranging-rows-in-custom-order-using-dplyr


phylum_genus <- unique(phylum_genus[c("taxa","Phylum")])


# Assigning Colors # 

phylum_names <- unique(phylum_genus$Phylum)
phylum_color_list <- list(Firmicutes = c("darkred", "#BC3C2960"),
                          Bacteroidota = c("#0072B5ff", "#0072B560"),
                          Proteobacteria = c("#20854Eff","#20854E60"),
                         # Actinobacteriota = c("#FFDC91ff","#FFDC9160"),
                          #Deferribacterota= c("#7876B1ff","#7876B160"),
                          #Campylobacterota=c("#EE4C97ff","#EE4C9760"),
                          Other = c("lightgrey","darkgrey")
                          
)

library(glue)
taxonomy_color_assigned <- NULL

for (z in 1:length(phylum_names)){
  
  taxonomy_color <- phylum_genus %>%
    filter(Phylum == glue("{phylum_names[z]}"))
  
  mycolors <- cbind(data.frame(color_code = colorRampPalette(phylum_color_list[[z]],alpha = T)(nrow(taxonomy_color))), "color_order" = nrow(taxonomy_color):1)
  
  taxonomy_color_assigned <- rbind(taxonomy_color_assigned,cbind(taxonomy_color,mycolors))
  
}


# Join genus with color #
genus_mean_RA_top_colors <- genus_and_metadata %>%
  left_join(taxonomy_color_assigned,., by = "taxa")

#   
color_codes <- unique(genus_mean_RA_top_colors[c("taxa","Phylum","color_code")])

# Ordering Groups in my dataframe #
genus_mean_RA_top_colors$timepoint <- factor( genus_mean_RA_top_colors$timepoint, levels = c("BL", "w3","w6"), ordered = T)
genus_mean_RA_top_colors$taxa <- factor( genus_mean_RA_top_colors$taxa, levels = unique(genus_mean_RA_top_colors$taxa), ordered = T)
genus_mean_RA_top_colors$genotype <- factor( genus_mean_RA_top_colors$genotype, levels = c("WT","KO"), ordered = T)

genus_mean_RA_top_colors %>%
  ggbarplot(x = "genotype", y = "percentage", fill = "taxa", add.params = list(size=0,fill="taxa",alpha=1), sort.val = "none",
            palette = color_codes$color_code, ylab = "Mean Relative Abundances", xlab = NULL) +
  labs(fill = "Genus") + facet_wrap(vars(timepoint), scales = "fixed") + labs(x = NULL) + theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), 
        axis.title=element_text(size=9), 
        legend.text = element_text(size=7),
        legend.title = element_text(size=0), 
        plot.title = element_text(size=9, face="bold", hjust = 0.5),
        strip.background=element_rect(colour="black",
                                      fill="white")) +
  rotate_x_text(45)

ggsave(file = "figures/Fig1D_barplot_RA_stacked_area_genus.pdf", width= 14, height=8, units = "cm")

##########################################
# Testing the differences of top genera #

genus_RAs_ps <- transform_sample_counts(genus_ps, function(x) x / sum(x) )
genus_RAs <- data.frame(genus_RAs_ps@otu_table)

genus_RAs_table <- cbind(genus_names,t(genus_RAs)) %>%
  pivot_longer(.,cols=2:ncol(.),names_to ="samples",values_to="RAs") %>%
  group_by(samples,Genus) %>%
  summarize(RAs = sum(RAs)) %>% # summing up all the "unknown" Genus
  filter(Genus != "unknown") # Deleting unknown Genus taxa.

genus_RAs_matrix <- genus_RAs_table %>%
  pivot_wider(.,names_from = Genus, values_from = RAs) %>%
  column_to_rownames(var="samples")

genus_RA_to_plot <- genus_RAs_matrix %>%
  cbind(data.frame(genus_RAs_ps@sam_data),.) %>%
  select(-others_genera) %>%
  pivot_longer(cols = A2:ncol(.), names_to = "taxa", values_to="RAs") %>%
  mutate(Condition = paste0(genotype," ",timepoint),
         time_num = ifelse(timepoint == "BL",0,
                           ifelse(timepoint == "w3",3,6.2)))



genus_RA_to_plot$Condition <- factor(genus_RA_to_plot$Condition, levels= c("WT BL","KO BL",
                                                                           "WT w3","KO w3",                                                                           "WT T2","KO T2",
                                                                           "WT w6","KO w6"), ordered = T)


p <- ggplot(genus_RA_to_plot, aes(x=Condition, y=RAs, fill=Condition)) + geom_boxplot() + geom_point(pch=21, size=2)
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
p <- p + theme_classic() + theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=7), 
                                 legend.title = element_text(size=0), plot.title = element_text(size=9, face="bold", hjust = 0.5), axis.text.x = element_text(angle=45, hjust = 1))
p <- p + facet_wrap(~taxa, scales = "free")
p

# Longitudinal Plot #

ggline(genus_RA_to_plot,x = "time_num", y = "RAs", add = c("mean_se","jitter"), color = "genotype",
       add.params = list(size=2,alpha=0.4), facet.by = c("taxa"), scales = "free") +
  stat_compare_means(aes(group = genotype), label = "p.signif", method = "wilcox.test", size = 2) +
  stat_summary(mapping = aes(fill=genotype), fun = "mean", geom = "point", color = "black", pch = 21, size = 4)+
  ylab("Relative Abundances") + xlab("Weeks") +
  scale_fill_manual(values = c("WT"="#5da1d4ff", "KO"="#c88cb2ff")) +
  scale_color_manual(values = c("WT"="#5da1d4ff", "KO"="#c88cb2ff")) +
  theme_classic() + theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=7), legend.text = element_text(size=7), 
                          legend.title = element_text(size=0), plot.title = element_text(size=5, face="bold", hjust = 0.5)) 

