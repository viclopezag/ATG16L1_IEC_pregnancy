
library(vegan)
library(ggplot2)
library(tidyverse)
library(glue)
library(ggpubr)
library(EnvStats)
library(phyloseq)
library(data.table)
library(microbiome)
library(rstudioapi)
library(microViz)
library(MicEco)

# Code for - Alpha Diversity and Beta Diversity Analyses: 16S data: Main Experiment -------

# Setting working directory #
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading phyloseq object #
ps <- readRDS("inputs/phyloseq_object.rds")

# Only working with samples containing more than 5000 high quality reads #

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
  .[,c("mouse_id","genotype","timepoint","Observed","Shannon","Chao1",
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
  facet_wrap(vars(timepoint), scales = "fixed") + theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=9), 
        legend.title = element_text(size=9), plot.title = element_text(size=9, face="bold", hjust = 0.5), axis.text.x = element_text(size = 7,angle=45, hjust = 1),
        strip.background=element_rect(colour="black",
                                      fill="white")) + rotate_x_text(45)

dir.create("figures")

ggsave(file = "figures/Fig1B_shannon.pdf", width= 9, height=8, units = "cm")

# All Diversity Measurements #
ggplot(estimated_richness_long, aes(x = genotype, y=values,fill = genotype, color = genotype)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(size = 2, alpha = 0.8) +
  scale_fill_manual(name=NULL,
                    values = c("WT"="#5da1d4ff","KO"="#c88cb2ff")) +
  scale_color_manual(name=NULL,
                     values = c("WT"="#5da1d4ff","KO"="#c88cb2ff")) +
  labs(x = "Genotype", y = "Shannon Diversity") +
  stat_compare_means(method = "wilcox", paired = FALSE, label = "p.signif", comparisons = list(c("WT","KO")), size = 2) + 
  facet_wrap(vars(timepoint), scales = "fixed") + theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=7), legend.text = element_text(size=7), 
        legend.title = element_text(size=7), plot.title = element_text(size=7, face="bold", hjust = 0.5), axis.text.x = element_text(size = 7,angle=45, hjust = 1),
        strip.background=element_rect(colour="black",
                                      fill="white")) + rotate_x_text(45) +
  facet_grid(diversity ~ timepoint, scales='free')

dir.create("supp_figures")
ggsave(file = "supp_figures/Supp_Fig1A_all_diversities.pdf", width= 12, height=25, units = "cm")

##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
#####:::::::::::::::::::: Beta Diversity Plot t ::::::::::::::::::::::::::#####
#####


# Principal Coordinate Analysis on Aitchison Distances #

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
  facet_wrap(vars(timepoint), scales = "fixed") 


plot(pca_plot)

# Plotting Ordispider #

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
               mapping = aes(xend = oAxis.1, yend = oAxis.2), size = 0.2) + # spiders
  geom_point(data = centroids, size = 2.5) +                         # centroids
  scale_color_manual(name=NULL,
                     values = c("WT"="#5da1d4ff","KO"="#c88cb2ff")) +  
  stat_ellipse(aes(group = genotype, color = genotype), type = "norm", linetype = 3, level=0.90) +
  geom_point(size=1) +                                              # sample scores
  facet_wrap(vars(timepoint), scales = "fixed") + theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=9), 
        legend.title = element_text(size=9), plot.title = element_text(size=9, face="bold", hjust = 0.5), axis.text.x = element_text(size = 7,angle=45, hjust = 1),
        strip.background=element_rect(colour="black",
                                      fill="white")) +
  xlab("PC1 [15.6%]") + ylab("PC2 [8.2%]") 


ggsave(file = "figures/Fig1C_PCoA_ordispider_plot.pdf", width= 9, height=8, units = "cm")

# :::::::::   Bray-Curtis :::::::: #

ps.coord.bray <-
  ps_f %>%
  rarefy_even_depth(.,rngseed=1989) %>% # let's rarefy the data to make sure that sequencing depth is similar for all samples
  ordinate(physeq=.,
           method="PCoA",
           distance="bray")

pcoa_plot_bray <- plot_ordination(ps_f, ps.coord.bray, color="genotype", axes = c(1,2)) +
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
  facet_wrap(vars(timepoint), scales = "fixed") 


plot(pcoa_plot_bray)


# Plotting Ordispide #

# computing centroids
centroids <- 
  cbind(data.frame(ps_f@sam_data) %>% select(genotype,timepoint),data.frame(ps.coord.bray$vectors)) %>%
  group_by(genotype,timepoint) %>%
  summarize(Axis.1 = mean(Axis.1),
            Axis.2 = mean(Axis.2))

# selecting scores
scrs <- cbind(data.frame(ps_f@sam_data) %>% select(genotype,timepoint),data.frame(ps.coord.bray$vectors)[,c("Axis.1","Axis.2")])

segs <- scrs %>%
  inner_join(centroids %>% rename(oAxis.1 = Axis.1,
                                  oAxis.2 = Axis.2),
             by = c("genotype","timepoint"))

#::Spider Plot ::#
ggplot(scrs, aes(x = Axis.1, y = Axis.2, colour = genotype)) +
  geom_segment(data = segs,
               mapping = aes(xend = oAxis.1, yend = oAxis.2), size = 0.2) + # spiders
  geom_point(data = centroids, size = 2.5) +                         # centroids
  scale_color_manual(name=NULL,
                     values = c("WT"="#5da1d4ff","KO"="#c88cb2ff")) +  
  stat_ellipse(aes(group = genotype, color = genotype), type = "norm", linetype = 3, level=0.90) +
  geom_point(size=1) +                                              # sample scores
  facet_wrap(vars(timepoint), scales = "fixed") + theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=9), 
        legend.title = element_text(size=9), plot.title = element_text(size=9, face="bold", hjust = 0.5), axis.text.x = element_text(size = 7,angle=45, hjust = 1),
        strip.background=element_rect(colour="black",
                                      fill="white")) +
  xlab("PCoA1 [26.1%]") + ylab("PCoA2 [14.4%]") 


ggsave(file = "supp_figures/Supp_Fig1B_PCoA_ordispider_plot_bray_curtis.pdf", width= 9, height=8, units = "cm")

# :::: Jaccard Index :::: #

ps.coord.jaccard <-
  ps_f %>%
  rarefy_even_depth(.,rngseed=1989) %>% # let's rarefy the data to make sure that sequencing depth is similar for all samples
  ordinate(physeq=.,
           method="PCoA",
           distance="jaccard")

pcoa_plot_jaccard <- plot_ordination(ps_f, ps.coord.jaccard, color="genotype", axes = c(1,2)) +
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
  facet_wrap(vars(timepoint), scales = "fixed") 


plot(pcoa_plot_jaccard)


# Plotting Ordispide #

# computing centroids
centroids <- 
  cbind(data.frame(ps_f@sam_data) %>% select(genotype,timepoint),data.frame(ps.coord.jaccard$vectors)) %>%
  group_by(genotype,timepoint) %>%
  summarize(Axis.1 = mean(Axis.1),
            Axis.2 = mean(Axis.2))

# selecting scores
scrs <- cbind(data.frame(ps_f@sam_data) %>% select(genotype,timepoint),data.frame(ps.coord.jaccard$vectors)[,c("Axis.1","Axis.2")])

segs <- scrs %>%
  inner_join(centroids %>% rename(oAxis.1 = Axis.1,
                                  oAxis.2 = Axis.2),
             by = c("genotype","timepoint"))

#::Spider Plot ::#
ggplot(scrs, aes(x = Axis.1, y = Axis.2, colour = genotype)) +
  geom_segment(data = segs,
               mapping = aes(xend = oAxis.1, yend = oAxis.2), size = 0.2) + # spiders
  geom_point(data = centroids, size = 2.5) +                         # centroids
  scale_color_manual(name=NULL,
                     values = c("WT"="#5da1d4ff","KO"="#c88cb2ff")) +  
  stat_ellipse(aes(group = genotype, color = genotype), type = "norm", linetype = 3, level=0.90) +
  geom_point(size=1) +                                              
  facet_wrap(vars(timepoint), scales = "fixed") + theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=9), legend.text = element_text(size=9), 
        legend.title = element_text(size=9), plot.title = element_text(size=9, face="bold", hjust = 0.5), axis.text.x = element_text(size = 7,angle=45, hjust = 1),
        strip.background=element_rect(colour="black",
                                      fill="white")) +
  xlab("PCoA1 [17.4%]") + ylab("PCoA2 [10.2%]") 


ggsave(file = "supp_figures/Supp_Fig1C_PCoA_ordispider_plot_jaccard.pdf", width= 9, height=8, units = "cm")


#::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Constraint Analysis of Principal Coordinates #

library(microViz)

arrow_colors <- data.frame(subset_samples(ps_arrange(ps_clr, time_mouse))@sam_data) %>%
  mutate(arrow_colors = if_else(genotype == "KO","#c88cb2ff","#5da1d4ff")) %>%
  select(arrow_colors) %>%
  pull()

ps.cap.aitchison <-
  subset_samples(ps_arrange(ps_clr, time_mouse)) %>%
  ordinate(physeq=.,
           method="CAP",
           distance="euclidean",
           formula= ~ timepoint+genotype)


cap_plot <- plot_ordination(subset_samples(ps_arrange(ps_clr, time_mouse)), ps.cap.aitchison, color="genotype_with_baseline",shape = "timepoint", axes = c(1,2)) +
  geom_point(size=2) + 
  scale_color_manual(values=c("BL"="darkgrey","KO"="#c88cb2ff", "WT"="#5da1d4ff")) +
  scale_shape_manual(values=c("BL"=19,"w3"=17,"w6"=18))+
  stat_ellipse(aes(group = genotype_with_baseline, color = genotype_with_baseline), type = "norm", linetype = 1, level=0.70) +
  stat_ellipse(aes(group = genotype_with_baseline, color = genotype_with_baseline), type = "norm", linetype = 2, level=0.80) +
  labs(color = "Genotype", shape = "Weeks") +
  theme_bw() +
  theme(axis.text=element_text(size=7, color = "black"), axis.title=element_text(size=7), legend.text = element_text(size=7), 
        legend.title = element_text(size=7), plot.title = element_text(size=7, face="bold", hjust = 0.5), axis.text.x = element_text(size = 7,angle=45, hjust = 1),
        strip.background=element_rect(colour="black",
                                      fill="white"))


# Add the environmental variables as arrows
arrowmat <- vegan::scores(ps.cap.aitchison, display="bp")

rownames(arrowmat) <- c("  w3","  w6","     WT")

labels = rownames(arrowmat)

# function to extract the last two letters of "labels"
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# Add labels, make a data.frame
arrowdf <- data.frame(labels = labels, arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1,
                 yend = CAP2,
                 x = 0,
                 y = 0,
                 shape = NULL,
                 color = NULL,
                 label = labels)

label_map <- aes(x = 1.5 * CAP1,
                 y = 1.5 * CAP2,
                 shape = NULL,
                 color = NULL,
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))


# Make a new graphic
cap_plot <- cap_plot +
  geom_segment(
    mapping = arrow_map,
    size = .5,
    data = arrowdf,
    color = "black",
    arrow = arrowhead
  ) +
  geom_text(
    mapping = label_map,
    size = 6,
    data = arrowdf,
    show.legend = FALSE) + theme_bw()

plot(cap_plot)

cap_plot <- cap_plot + geom_path(aes(group = subset_samples(ps_arrange(ps_clr, time_mouse))@sam_data$mouse_id),
                                 colour=arrow_colors,
                                 arrow = arrow(type = "closed",
                                               
                                               length=unit(0.1, "inches")), size = 1, alpha = 0.2)
print(cap_plot)

ggsave(file = "supp_figures/Supp_Fig1D_CAP_Aitchison-Genotype-Time-Arrows.pdf", width= 9, height=8, units = "cm")

#:::::::::::::::::::::::::::: PERMANOVA ANALYSIS ::::::::::::::::::::::::::::::#

# PERMANOVA on Aitchison Distances #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA BL #

# Computing distance aitchison matrix #
micro.dis_BL <-
  subset_samples(ps_clr,timepoint == "BL") %>% #Only pregnancy BL
  phyloseq::distance("euclidean") # Aitchison Distances

# PERMANOVA
(permanova_CAPSCALE_BL <- adonis2(micro.dis_BL ~ genotype,
                                  data = data.frame(subset_samples(ps_clr,timepoint == "BL")@sam_data),
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
dir.create("outputs/permanova_Fig1C_and_supp_Fig1BC")
write.csv(permanova_CAPSCALE_BL_table, file = "outputs/permanova_Fig1C_and_supp_Fig1BC/aitchison_permanova_and_omegaSq_BL.csv", row.names = T)

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

write.csv(permanova_CAPSCALE_w3_table, file = "outputs/permanova_Fig1C_and_supp_Fig1BC/aitchison_permanova_and_omegaSq_w3.csv", row.names = T)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA w6#

# Computing distance aitchison matrix #
micro.dis_w6 <-
  subset_samples(ps_clr,timepoint == "w6") %>% #Only pregnancy w6
  phyloseq::distance("euclidean") # Aitchison Distances

# PERMANOVA
(permanova_CAPSCALE_w6 <- adonis2(micro.dis_w6 ~ genotype,
                                    data = data.frame(subset_samples(ps_clr,timepoint == "w6")@sam_data),
                                    permutations = 10000))

(FDR_corrected_w6 <- data.frame(p.adjust(permanova_CAPSCALE_w6$`Pr(>F)`, method = "BH", n = length(permanova_CAPSCALE_w6$`Pr(>F)`))))

# Computing Effect size in PERMANOVA # 
# https://microucph.github.io/amplicon_data_analysis/html/omegasq.html

permanova_CAPSCALE_w6_table <- data.frame(adonis_OmegaSq(permanova_CAPSCALE_w6, partial = TRUE)) %>%
  cbind(.,data.frame(permanova_CAPSCALE_w6$R2),FDR_corrected_w6)

colnames(permanova_CAPSCALE_w6_table) <- c("Df","SumOfSqs","F","parOmegaSq","Pval","R2","FDR")

permanova_CAPSCALE_w6_table <-permanova_CAPSCALE_w6_table %>%
  select(Df,SumOfSqs,R2,`F`,Pval,parOmegaSq,FDR)

write.csv(permanova_CAPSCALE_w6_table, file = "outputs/permanova_Fig1C_and_supp_Fig1BC/aitchison_permanova_and_omegaSq_w6.csv", row.names = T)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA on Bray Curtis Distances #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA BL #

# Computing distance Bray Curtis matrix #
micro.dis_BL <-
  rarefy_even_depth(ps_f,rngseed=1989) %>% # let's rarefy the data to make sure that sequencing depth is similar for all samples
  subset_samples(.,timepoint == "BL") %>% #Only pregnancy BL
  phyloseq::distance("bray") # Bray Curtis Distances

# PERMANOVA
(permanova_CAPSCALE_BL <- adonis2(micro.dis_BL ~ genotype,
                                  data = data.frame(subset_samples(ps_clr,timepoint == "BL")@sam_data),
                                  permutations = 10000))

(FDR_corrected_BL <- data.frame(p.adjust(permanova_CAPSCALE_BL$`Pr(>F)`, method = "BH", n = length(permanova_CAPSCALE_BL$`Pr(>F)`))))

# Computing Effect size in PERMANOVA # 
# https://microucph.github.io/amplicon_data_analysis/html/omegasq.html

permanova_CAPSCALE_BL_table <- data.frame(adonis_OmegaSq(permanova_CAPSCALE_BL, partial = TRUE)) %>%
  cbind(.,data.frame(permanova_CAPSCALE_BL$R2),FDR_corrected_BL)

colnames(permanova_CAPSCALE_BL_table) <- c("Df","SumOfSqs","F","parOmegaSq","Pval","R2","FDR")

permanova_CAPSCALE_BL_table <-permanova_CAPSCALE_BL_table %>%
  select(Df,SumOfSqs,R2,`F`,Pval,parOmegaSq,FDR)

write.csv(permanova_CAPSCALE_BL_table, file = "outputs/permanova_Fig1C_and_supp_Fig1BC/bray_permanova_and_omegaSq_BL.csv", row.names = T)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA w3 #

# Computing distance Bray Curtis matrix #
micro.dis_w3 <-
  rarefy_even_depth(ps_f,rngseed=1989) %>% # let's rarefy the data to make sure that sequencing depth is similar for all samples
  subset_samples(.,timepoint == "w3") %>% #Only pregnancy w3
  phyloseq::distance("bray") # Bray Curtis Distances

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

write.csv(permanova_CAPSCALE_w3_table, file = "outputs/permanova_Fig1C_and_supp_Fig1BC/bray_permanova_and_omegaSq_w3.csv", row.names = T)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA w6 #

# Computing distance Bray Curtis matrix #
micro.dis_w6 <-
  rarefy_even_depth(ps_f,rngseed=1989) %>% # let's rarefy the data to make sure that sequencing depth is similar for all samples
  subset_samples(.,timepoint == "w6") %>% #Only pregnancy w6
  phyloseq::distance("bray") # Bray Curtis Distances

# PERMANOVA
(permanova_CAPSCALE_w6 <- adonis2(micro.dis_w6 ~ genotype,
                                    data = data.frame(subset_samples(ps_clr,timepoint == "w6")@sam_data),
                                    permutations = 10000))

(FDR_corrected_w6 <- data.frame(p.adjust(permanova_CAPSCALE_w6$`Pr(>F)`, method = "BH", n = length(permanova_CAPSCALE_w6$`Pr(>F)`))))

# Computing Effect size in PERMANOVA # 
# https://microucph.github.io/amplicon_data_analysis/html/omegasq.html

permanova_CAPSCALE_w6_table <- data.frame(adonis_OmegaSq(permanova_CAPSCALE_w6, partial = TRUE)) %>%
  cbind(.,data.frame(permanova_CAPSCALE_w6$R2),FDR_corrected_w6)

colnames(permanova_CAPSCALE_w6_table) <- c("Df","SumOfSqs","F","parOmegaSq","Pval","R2","FDR")

permanova_CAPSCALE_w6_table <-permanova_CAPSCALE_w6_table %>%
  select(Df,SumOfSqs,R2,`F`,Pval,parOmegaSq,FDR)

write.csv(permanova_CAPSCALE_w6_table, file = "outputs/permanova_Fig1C_and_supp_Fig1BC/bray_permanova_and_omegaSq_w6.csv", row.names = T)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA on Jaccard Distances #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA BL #

# Computing distance Jaccard matrix #
micro.dis_BL <-
  rarefy_even_depth(ps_f,rngseed=1989) %>% # let's rarefy the data to make sure that sequencing depth is similar for all samples
  subset_samples(.,timepoint == "BL") %>% #Only pregnancy BL
  phyloseq::distance("jaccard") # Jaccard Distances

# PERMANOVA
(permanova_CAPSCALE_BL <- adonis2(micro.dis_BL ~ genotype,
                                  data = data.frame(subset_samples(ps_clr,timepoint == "BL")@sam_data),
                                  permutations = 10000))

(FDR_corrected_BL <- data.frame(p.adjust(permanova_CAPSCALE_BL$`Pr(>F)`, method = "BH", n = length(permanova_CAPSCALE_BL$`Pr(>F)`))))

# Computing Effect size in PERMANOVA # 
# https://microucph.github.io/amplicon_data_analysis/html/omegasq.html

permanova_CAPSCALE_BL_table <- data.frame(adonis_OmegaSq(permanova_CAPSCALE_BL, partial = TRUE)) %>%
  cbind(.,data.frame(permanova_CAPSCALE_BL$R2),FDR_corrected_BL)

colnames(permanova_CAPSCALE_BL_table) <- c("Df","SumOfSqs","F","parOmegaSq","Pval","R2","FDR")

permanova_CAPSCALE_BL_table <-permanova_CAPSCALE_BL_table %>%
  select(Df,SumOfSqs,R2,`F`,Pval,parOmegaSq,FDR)

write.csv(permanova_CAPSCALE_BL_table, file = "outputs/permanova_Fig1C_and_supp_Fig1BC/jaccard_permanova_and_omegaSq_BL.csv", row.names = T)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA w3 #

# Computing distance Jaccard matrix #
micro.dis_w3 <-
  rarefy_even_depth(ps_f,rngseed=1989) %>% # let's rarefy the data to make sure that sequencing depth is similar for all samples
  subset_samples(.,timepoint == "w3") %>% #Only pregnancy T3
  phyloseq::distance("jaccard") # Jaccard Distances

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

write.csv(permanova_CAPSCALE_w3_table, file = "outputs/permanova_Fig1C_and_supp_Fig1BC/jaccard_permanova_and_omegaSq_w3.csv", row.names = T)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# PERMANOVA w6 #

# Computing distance Jaccard matrix #
micro.dis_w6 <-
  rarefy_even_depth(ps_f,rngseed=1989) %>% # let's rarefy the data to make sure that sequencing depth is similar for all samples
  subset_samples(.,timepoint == "w6") %>% #Only pregnancy w6
  phyloseq::distance("jaccard") # Jaccard Distances

# PERMANOVA
(permanova_CAPSCALE_w6 <- adonis2(micro.dis_w6 ~ genotype,
                                    data = data.frame(subset_samples(ps_clr,timepoint == "w6")@sam_data),
                                    permutations = 10000))

(FDR_corrected_w6 <- data.frame(p.adjust(permanova_CAPSCALE_w6$`Pr(>F)`, method = "BH", n = length(permanova_CAPSCALE_w6$`Pr(>F)`))))

# Computing Effect size in PERMANOVA # 
# https://microucph.github.io/amplicon_data_analysis/html/omegasq.html

permanova_CAPSCALE_w6_table <- data.frame(adonis_OmegaSq(permanova_CAPSCALE_w6, partial = TRUE)) %>%
  cbind(.,data.frame(permanova_CAPSCALE_w6$R2),FDR_corrected_w6)

colnames(permanova_CAPSCALE_w6_table) <- c("Df","SumOfSqs","F","parOmegaSq","Pval","R2","FDR")

permanova_CAPSCALE_w6_table <-permanova_CAPSCALE_w6_table %>%
  select(Df,SumOfSqs,R2,`F`,Pval,parOmegaSq,FDR)

write.csv(permanova_CAPSCALE_w6_table, file = "outputs/permanova_Fig1C_and_supp_Fig1BC/jaccard_permanova_and_omegaSq_w6.csv", row.names = T)
