
library(vegan)
library(ggplot2)
library(tidyverse)
library(glue)
library(ggpubr)
library(EnvStats)
library(phyloseq)
library(data.table)
library(microbiome)
library(microViz)
library(MicEco)



# Code for - Alpha Diversity and Beta Diversity Analyses: 16S data : Validation Experiment -------

# Setting working directory #
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


ps <- readRDS("inputs/phyloseq_object_validation.rds")


# Cutoff #

#(ps_fr <- subset_samples(ps, sample_sums(ps) > 5000 & pregnancy_status == "P"))
(ps_fr <- subset_samples(ps, sample_sums(ps) > 5000))

# Deleting taxa with relative abundance lowest than 1e-5
minTotRelAbun <- 1e-5          
x <- taxa_sums(ps_fr)
keepTaxa <- (x / sum(x)) > minTotRelAbun
ps_f <- prune_taxa(keepTaxa, ps_fr)


# # Estimating Richness #
# 
estimated_richness_long <- estimate_richness(rarefy_even_depth(ps_f,rngseed=1989)) %>%
  cbind(data.frame(ps_f@sam_data),.) %>%
  data.table(.) %>%
  .[,c("animal_nr","genotype","timepoint","Observed","Shannon","Chao1",
       "Fisher","InvSimpson","Simpson")] %>%
  pivot_longer(.,
               cols=Observed:Simpson,
               names_to="diversity",
               values_to = "values")

estimated_richness_long$genotype <- factor(estimated_richness_long$genotype, levels = c("WT","KO"),
                                           ordered = T)

ggplot(estimated_richness_long %>% filter(diversity == "Shannon"), aes(x = genotype, y=values,fill = genotype, color = genotype)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(size = 2, alpha = 0.8) +
  scale_fill_manual(name=NULL,
                    values = c("WT"="#5da1d4ff","KO"="#c88cb2ff")) +
  scale_color_manual(name=NULL,
                    values = c("WT"="#5da1d4ff","KO"="#c88cb2ff")) +
  labs(x = "Genotype", y = "Shannon Diversity") +
  stat_compare_means(method = "wilcox", paired = FALSE, label = "p.signif", comparisons = list(c("WT","KO")), size = 2) + 
  facet_wrap(vars(timepoint), scales = "fixed", ncol = 4) + theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=9), 
        legend.title = element_text(size=9), plot.title = element_text(size=9, face="bold", hjust = 0.5), axis.text.x = element_text(size = 7,angle=45, hjust = 1),
        strip.background=element_rect(colour="black",
                                      fill="white"))

ggsave(file = "figures/Fig4B_shannon.pdf", width= 9, height=8, units = "cm")


#------------------------------------------------------------------------------
#####:::::::::::::::::::: Beta Diversity Plot t ::::::::::::::::::::::::::#####
#####

# Principal Coordinate Analysis #

ps_clr <- microbiome::transform(ps_f, transform = "clr")


ps.coord.pca <-
  ps_clr %>%
  ordinate(physeq=.,
           method="PCoA",
           distance="euclidean")

pca_plot <- plot_ordination(ps_clr, ps.coord.pca, color="genotype", axes = c(1,2)) +
  geom_point(size=2) + 
  scale_color_manual(name=NULL,
                     values = c("WT"="#5da1d4ff","KO"="#c88cb2ff")) +  
  stat_ellipse(aes(group = genotype, color = genotype), type = "norm", linetype = 2, level=0.80) +
  stat_ellipse(aes(group = genotype, color = genotype), type = "norm", linetype = 1,  level=0.70) +
  labs(color = "Genotype") + theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=9), 
        legend.title = element_text(size=9), plot.title = element_text(size=10, face="bold", hjust = 0.5), axis.text.x = element_text(size = 7,angle=45, hjust = 1),
        strip.background=element_rect(colour="black",
                                      fill="white")) +
  facet_wrap(vars(timepoint), scales = "fixed", ncol = 4) 


plot(pca_plot)


# Plotting Ordispide #

# computing centroids
centroids <- 
  cbind(data.frame(ps_clr@sam_data) %>% select(genotype,timepoint),data.frame(ps.coord.pca$vectors)) %>%
  group_by(genotype,timepoint) %>%
  summarize(Axis.1 = mean(Axis.1),
            Axis.2 = mean(Axis.2))

 # selecting scores
scrs <- cbind(data.frame(ps_clr@sam_data) %>% select(genotype,timepoint),data.frame(ps.coord.pca$vectors)[,c("Axis.1","Axis.2")])

# To draw the spider, we need to us geom_segment() which requires coordiates to draw the segment from and to. 
# Our to coordinates, the xend and yend aesthetics will be the centroids. So we need to replicate the group 
# centroid for each observation in the group. This we facilitate by a left join via merge:
segs <- scrs %>%
  inner_join(centroids %>% rename(oAxis.1 = Axis.1,
                                  oAxis.2 = Axis.2),
             by = c("genotype","timepoint"))

#::Spider Plot ::#
ggplot(scrs, aes(x = Axis.1, y = Axis.2, colour = genotype)) +
  geom_segment(data = segs,
               mapping = aes(xend = oAxis.1, yend = oAxis.2), size = 0.3) + # spiders
  geom_point(data = centroids, size = 4) +                         # centroids
  scale_color_manual(name=NULL,
                     values = c("WT"="#5da1d4ff","KO"="#c88cb2ff")) +  
  stat_ellipse(aes(group = genotype, color = genotype), type = "norm", linetype = 3, level=0.90) +
  geom_point() +                                              # sample scores
  facet_wrap(vars(timepoint), scales = "fixed", ncol=4) + theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=9), 
        legend.title = element_text(size=9), plot.title = element_text(size=10, face="bold", hjust = 0.5), axis.text.x = element_text(size = 7,angle=45, hjust = 1),
        strip.background=element_rect(colour="black",
                                      fill="white")) +
  xlab("PC1 [13%]") + ylab("PC2 [7.3%]")

ggsave(file = "figures/Fig4C_PCoA_ordispider_plot.pdf", width= 9, height=8, units = "cm")

#------------------------------------------------------------------------------
#:::::::::::::::::::::::::::::: PERMANOVA ANALYSIS :::::::::::::::::::::::::::::#

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA BL #

# Computing distance aitchison matrix #
micro.dis_BL <-
  subset_samples(ps_clr,timepoint == "BL") %>% #Only pregnancy BL
  phyloseq::distance("euclidean") # Aitchison Distances

# PERMANOVA
(permanova_CAPSCALE_BL <- adonis2(micro.dis_BL ~ genotype,
                               data = data.frame(subset_samples(ps_clr,timepoint == "BL")@sam_data),
                               #strata = data.frame(sample_data(subset_samples(ps_arrange(ps_clr, time_mouse), Timepoint != "d12" | Timepoint != "BL0")))[,"Mouse"],
                               permutations = 10000))

(FDR_corrected_BL <- data.frame(p.adjust(permanova_CAPSCALE_BL$`Pr(>F)`, method = "BH", n = length(permanova_CAPSCALE_BL$`Pr(>F)`))))

# Computing Effect size in PERMANOVA # 
# https://microucph.github.io/amplicon_data_analysis/html/omegasq.html

permanova_CAPSCALE_BL_table <- data.frame(adonis_OmegaSq(permanova_CAPSCALE_BL, partial = TRUE)) %>%
  cbind(.,data.frame(permanova_CAPSCALE_BL$R2),FDR_corrected_BL)

colnames(permanova_CAPSCALE_BL_table) <- c("Df","SumOfSqs","F","parOmegaSq","Pval","R2","FDR")

permanova_CAPSCALE_BL_table <-permanova_CAPSCALE_BL_table %>%
  select(Df,SumOfSqs,R2,`F`,Pval,parOmegaSq,FDR)

dir.create("outputs")
dir.create("outputs/permanova_Fig4C")

write.csv(permanova_CAPSCALE_BL_table, file = "outputs/permanova_Fig4C/aitchison_permanova_and_omegaSq_BL.csv", row.names = T)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA w3 #

# Computing distance aitchison matrix #
micro.dis_w3 <-
  subset_samples(ps_clr,timepoint == "w3") %>% #Only pregnancy w3
  phyloseq::distance("euclidean") # Aitchison Distances

# PERMANOVA
(permanova_CAPSCALE_w3 <- adonis2(micro.dis_w3 ~ genotype,
                                  data = data.frame(subset_samples(ps_clr,timepoint == "w3")@sam_data),
                                  permutations = 10000))

(FDR_corrected_w3 <- data.frame(p.adjust(permanova_CAPSCALE_w3$`Pr(>F)`, method = "BH", n = length(permanova_CAPSCALE_w3$`Pr(>F)`))))

# Computing Effect size in PERMANOVA # 
# https://microucph.github.io/amplicon_data_analysis/html/omegasq.html

permanova_CAPSCALE_w3_table <- data.frame(adonis_OmegaSq(permanova_CAPSCALE_w3, partial = TRUE)) %>%
  cbind(.,data.frame(permanova_CAPSCALE_w3$R2),FDR_corrected_w3)

colnames(permanova_CAPSCALE_w3_table) <- c("Df","SumOfSqs","F","parOmegaSq","Pval","R2","FDR")

permanova_CAPSCALE_w3_table <-permanova_CAPSCALE_w3_table %>%
  select(Df,SumOfSqs,R2,`F`,Pval,parOmegaSq,FDR)

write.csv(permanova_CAPSCALE_w3_table, file = "outputs/permanova_Fig4C/aitchison_permanova_and_omegaSq_w3.csv", row.names = T)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA w1 #

# Computing distance aitchison matrix #
micro.dis_w1 <-
  subset_samples(ps_clr,timepoint == "w1") %>% #Only pregnancy w1
  phyloseq::distance("euclidean") # Aitchison Distances

# PERMANOVA
(permanova_CAPSCALE_w1 <- adonis2(micro.dis_w1 ~ genotype,
                                  data = data.frame(subset_samples(ps_clr,timepoint == "w1")@sam_data),
                                  permutations = 10000))

(FDR_corrected_w1 <- data.frame(p.adjust(permanova_CAPSCALE_w1$`Pr(>F)`, method = "BH", n = length(permanova_CAPSCALE_w1$`Pr(>F)`))))

# Computing Effect size in PERMANOVA # 
# https://microucph.github.io/amplicon_data_analysis/html/omegasq.html

permanova_CAPSCALE_w1_table <- data.frame(adonis_OmegaSq(permanova_CAPSCALE_w1, partial = TRUE)) %>%
  cbind(.,data.frame(permanova_CAPSCALE_w1$R2),FDR_corrected_w1)

colnames(permanova_CAPSCALE_w1_table) <- c("Df","SumOfSqs","F","parOmegaSq","Pval","R2","FDR")

permanova_CAPSCALE_w1_table <-permanova_CAPSCALE_w1_table %>%
  select(Df,SumOfSqs,R2,`F`,Pval,parOmegaSq,FDR)

write.csv(permanova_CAPSCALE_w1_table, file = "outputs/permanova_Fig4C/aitchison_permanova_and_omegaSq_w1.csv", row.names = T)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA w2 #

# Computing distance aitchison matrix #
micro.dis_w2 <-
  subset_samples(ps_clr,timepoint == "w2") %>% #Only pregnancy w2
  phyloseq::distance("euclidean") # Aitchison Distances

# PERMANOVA
(permanova_CAPSCALE_w2 <- adonis2(micro.dis_w2 ~ genotype,
                                  data = data.frame(subset_samples(ps_clr,timepoint == "w2")@sam_data),
                                  permutations = 10000))

(FDR_corrected_w2 <- data.frame(p.adjust(permanova_CAPSCALE_w2$`Pr(>F)`, method = "BH", n = length(permanova_CAPSCALE_w2$`Pr(>F)`))))

# Computing Effect size in PERMANOVA # 
# https://microucph.github.io/amplicon_data_analysis/html/omegasq.html

permanova_CAPSCALE_w2_table <- data.frame(adonis_OmegaSq(permanova_CAPSCALE_w2, partial = TRUE)) %>%
  cbind(.,data.frame(permanova_CAPSCALE_w2$R2),FDR_corrected_w2)

colnames(permanova_CAPSCALE_w2_table) <- c("Df","SumOfSqs","F","parOmegaSq","Pval","R2","FDR")

permanova_CAPSCALE_w2_table <-permanova_CAPSCALE_w2_table %>%
  select(Df,SumOfSqs,R2,`F`,Pval,parOmegaSq,FDR)

write.csv(permanova_CAPSCALE_w2_table, file = "outputs/permanova_Fig4C/aitchison_permanova_and_omegaSq_w2.csv", row.names = T)

