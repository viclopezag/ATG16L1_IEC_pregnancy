
library(dplyr)
library(ggpubr)
library(tidyr)


# Code for - Plotting Physiology Differences in: Validation experiment -------

# Setting working directory #
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading Physiological Data #
taqman_expression <- read.csv("inputs/taqman_normalized_expression.csv", check.names = F)
tissue_measurements <- read.csv("inputs/tissue_measurements.csv", check.names = F)

taqman_expression_ileum <- subset(taqman_expression, Location == "Ileum")
taqman_expression_ileum_long <- gather(taqman_expression_ileum, cytokines, expression_value, gapdh:cxcl1, factor_key = TRUE) %>%
unite("genotype_pregnant_status", Status:Genotype, remove = F) # Create a concatenated column genotype_pregnant_status
taqman_expression_colon <- subset(taqman_expression, Location == "Colon")
taqman_expression_colon_long <- gather(taqman_expression_colon, cytokines, expression_value, gapdh:cxcl1, factor_key = TRUE) %>%
  unite("genotype_pregnant_status", Status:Genotype, remove = F) # Create a concatenated column genotype_pregnant_status


tissue_measurements_long <- gather(tissue_measurements, genotype, measurement, WT:KO, factor_key = TRUE) %>%
                            unite("genotype_pregnant_status", pregnant_status:genotype, remove = F) # Create a concatenated column genotype_pregnant_status


my_comparisons <- list( c("P_WT", "NP_WT"),
                        c("P_WT", "P_KO"),
                        c("P_WT", "NP_KO"),
                        c("NP_WT", "P_KO"),
                        c("NP_WT", "NP_KO"),
                        c("P_KO", "NP_KO"))


tissue_measurements_long <- tissue_measurements_long %>%
  mutate(genotype_pregnant_status = factor(genotype_pregnant_status, 
                                           levels = c("NP_WT","NP_KO","P_WT","P_KO"), ordered = T))

# Plotting tissue measurements
# Pregnant and non-pregnant mice

pdf(file="supp_figures/Supp_Fig6A_body_weight.pdf", width = 3, height = 3)

# Body Weight #
bxp <- ggboxplot(
  tissue_measurements_long %>% filter(feature=="body_weight"), 
  x = "genotype_pregnant_status", y = "measurement",
  color = "genotype",
  fill = "genotype",
  add = "jitter",
  alpha=0.5,
  add.params = list(size=3, color="genotype",alpha=0.5),
  palette = c("WT"="#5da1d4ff","KO"="#c88cb2ff"), notch = FALSE, ylab = "Body weight (g)", xlab = "Pregnancy/Genotype"
  #facet.by = "feature", scales = "free"
) +
  stat_compare_means(method = "wilcox", paired = FALSE, comparisons = my_comparisons, label = "p.signif",size=2) +
  labs(fill="genotype") +
  theme_bw() + theme(axis.text=element_text(size=7, color="black"), 
                          axis.title=element_text(size=9), 
                          legend.text = element_text(size=9), 
                          legend.title = element_text(size=0), 
                          plot.title = element_text(size=9, face="bold", hjust = 0.5)) + rotate_x_text(45)

bxp

dev.off()


# Spleen Weight #
pdf(file="supp_figures/Supp_Fig6B_spleen_weight.pdf", width = 3, height = 3)

bxp <- ggboxplot(
  tissue_measurements_long %>% filter(feature=="normalized_spleen_weight"), 
  x = "genotype_pregnant_status", y = "measurement",
  color = "genotype",
  fill = "genotype",
  add = "jitter",
  alpha=0.5,
  add.params = list(size=3, color="genotype",alpha=0.5),
  palette = c("WT"="#5da1d4ff","KO"="#c88cb2ff"), notch = FALSE, ylab = "Spleen weight (g)", xlab = "Pregnancy/Genotype"
  #facet.by = "feature", scales = "free"
) +
  stat_compare_means(method = "wilcox", paired = FALSE, comparisons = my_comparisons, label = "p.signif",size=2) +
  labs(fill="genotype") +
  theme_bw() + theme(axis.text=element_text(size=7, color="black"), 
                          axis.title=element_text(size=9), 
                          legend.text = element_text(size=9), 
                          legend.title = element_text(size=0), 
                          plot.title = element_text(size=9, face="bold", hjust = 0.5)) + rotate_x_text(45)

bxp

dev.off()

# Liver Weight #
pdf(file="supp_figures/Supp_Fig6C_liver_weight.pdf", width = 3, height = 3)

bxp <- ggboxplot(
  tissue_measurements_long %>% filter(feature=="normalized_liver_weight"), 
  x = "genotype_pregnant_status", y = "measurement",
  color = "genotype",
  fill = "genotype",
  add = "jitter",
  alpha=0.5,
  add.params = list(size=3, color="genotype",alpha=0.5),
  palette = c("WT"="#5da1d4ff","KO"="#c88cb2ff"), notch = FALSE, ylab = "Liver weight (g)", xlab = "Pregnancy/Genotype"
  #facet.by = "feature", scales = "free"
) +
  stat_compare_means(method = "wilcox", paired = FALSE, comparisons = my_comparisons, label = "p.signif",size=2) +
  labs(fill="genotype") +
  theme_bw() + theme(axis.text=element_text(size=7, color="black"), 
                          axis.title=element_text(size=9), 
                          legend.text = element_text(size=9), 
                          legend.title = element_text(size=0), 
                          plot.title = element_text(size=9, face="bold", hjust = 0.5)) + rotate_x_text(45)

bxp

dev.off()

pdf(file="supp_figures/Supp_Fig6D_caecum_weight.pdf", width = 3, height = 3)

# Caecum weight #

bxp <- ggboxplot(
  tissue_measurements_long %>% filter(feature=="caecum_weight"), 
  x = "genotype_pregnant_status", y = "measurement",
  color = "genotype",
  fill = "genotype",
  add = "jitter",
  alpha=0.5,
  add.params = list(size=3, color="genotype",alpha=0.5),
  palette = c("WT"="#5da1d4ff","KO"="#c88cb2ff"), notch = FALSE, ylab = "Caecum weight (g)", xlab = "Pregnancy/Genotype"
  #facet.by = "feature", scales = "free"
) +
  stat_compare_means(method = "wilcox", paired = FALSE, comparisons = my_comparisons, label = "p.signif",size=2) +
  labs(fill="genotype") +
  theme_bw() + theme(axis.text=element_text(size=7, color="black"), 
                          axis.title=element_text(size=9), 
                          legend.text = element_text(size=9), 
                          legend.title = element_text(size=0), 
                          plot.title = element_text(size=9, face="bold", hjust = 0.5)) + rotate_x_text(45)

bxp

dev.off()

pdf(file="supp_figures/Supp_Fig6E_colon_length.pdf", width = 3, height = 3)

# Colon length #
bxp <- ggboxplot(
  tissue_measurements_long %>% filter(feature=="colon_length"), 
  x = "genotype_pregnant_status", y = "measurement",
  color = "genotype",
  fill = "genotype",
  add = "jitter",
  alpha=0.5,
  add.params = list(size=3, color="genotype",alpha=0.5),
  palette = c("WT"="#5da1d4ff","KO"="#c88cb2ff"), notch = FALSE, ylab = "Colon length", xlab = "Pregnancy/Genotype"
  #facet.by = "feature", scales = "free"
) +
  stat_compare_means(method = "wilcox", paired = FALSE, comparisons = my_comparisons, label = "p.signif",size=2) +
  labs(fill="genotype") +
  theme_bw() + theme(axis.text=element_text(size=7, color="black"), 
                          axis.title=element_text(size=9), 
                          legend.text = element_text(size=9), 
                          legend.title = element_text(size=0), 
                          plot.title = element_text(size=9, face="bold", hjust = 0.5)) + rotate_x_text(45)

bxp

dev.off()


pdf(file="supp_figures/Supp_Fig6F_puppy_weight.pdf", width = 3, height = 3)

# Puppy Weight #

bxp <- ggboxplot(
  tissue_measurements_long %>% filter(feature=="puppy_weight"), 
  x = "genotype_pregnant_status", y = "measurement",
  color = "genotype",
  fill = "genotype",
  add = "jitter",
  alpha=0.5,
  add.params = list(size=3, color="genotype",alpha=0.5),
  palette = c("WT"="#5da1d4ff","KO"="#c88cb2ff"), notch = FALSE, ylab = "Puppy weight (g)", xlab = "Pregnancy/Genotype"
  #facet.by = "feature", scales = "free"
) +
  stat_compare_means(method = "wilcox", paired = FALSE, comparisons = list(c("P_WT","P_KO")), label = "p.signif",size=2) +
  labs(fill="genotype") +
  theme_bw() + theme(axis.text=element_text(size=7, color="black"), 
                     axis.title=element_text(size=9), 
                     legend.text = element_text(size=9), 
                     legend.title = element_text(size=0), 
                     plot.title = element_text(size=9, face="bold", hjust = 0.5)) + rotate_x_text(45)

bxp

dev.off()

pdf(file="supp_figures/Supp_Fig6G_puppy_length.pdf", width = 3, height = 3)

# Puppy length #

bxp <- ggboxplot(
  tissue_measurements_long %>% filter(feature=="puppy_length"), 
  x = "genotype_pregnant_status", y = "measurement",
  color = "genotype",
  fill = "genotype",
  add = "jitter",
  alpha=0.5,
  add.params = list(size=3, color="genotype",alpha=0.5),
  palette = c("WT"="#5da1d4ff","KO"="#c88cb2ff"), notch = FALSE, ylab = "Puppy length (cm)", xlab = "Pregnancy/Genotype"
  #facet.by = "feature", scales = "free"
) +
  stat_compare_means(method = "wilcox", paired = FALSE, comparisons = list(c("P_WT","P_KO")), label = "p.signif",size=2) +
  labs(fill="genotype") +
  theme_bw() + theme(axis.text=element_text(size=7, color="black"), 
                     axis.title=element_text(size=9), 
                     legend.text = element_text(size=9), 
                     legend.title = element_text(size=0), 
                     plot.title = element_text(size=9, face="bold", hjust = 0.5)) + rotate_x_text(45)

bxp



dev.off()

#-------------------------------------------------------------------------------
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Taqman expression -  Colonic mucose # 
taqman_expression_colon_long <- taqman_expression_colon_long %>%
  mutate(genotype_pregnant_status = ifelse(Status == "Pregnant" & Genotype == "WT","P_WT",
                                           ifelse(Status == "Pregnant" & Genotype == "KO","P_KO",
                                                  ifelse(Status == "Non-Pregnant" & Genotype == "WT","NP_WT","NP_KO")))) %>%
  filter(cytokines != "gapdh") %>%
  mutate(genotype_pregnant_status = factor(genotype_pregnant_status, 
                                           levels = c("NP_WT","NP_KO","P_WT","P_KO"), ordered = T))



grDevices::cairo_pdf(file = "figures/Fig4E_TNFa_expression.pdf", width = 3, height = 4,onefile=T)

# TNFalpha #
ggbarplot(
  taqman_expression_colon_long %>% filter(cytokines == "TNF"), 
  x = "genotype_pregnant_status", 
  y = "expression_value",
  add = c("mean_se","jitter"), 
  alpha = 0.5, color = "Genotype",fill = "Genotype",
  #facet.by = "cytokines", scale = "free",
  palette = c("WT"="#5da1d4ff","KO"="#c88cb2ff"), ylab = "Relative expression", xlab = "Pregnancy/Genotype"
) + ggtitle("TNF\u03b1") +
  stat_compare_means( method = "wilcox",
                      comparisons = list(c("P_WT", "P_KO"), c("NP_WT", "NP_KO"),
                                         c("P_WT", "NP_WT"),c("P_WT", "NP_KO"),
                                         c("P_KO", "NP_WT"),c("P_KO", "NP_KO")),
                      label = "p.signif"
  ) +
  theme_bw() + theme(axis.text=element_text(size=7, color="black"), 
                          axis.title=element_text(size=9), 
                          legend.text = element_text(size=9), 
                          legend.title = element_text(size=0), 
                          plot.title = element_text(size=9, face="bold", hjust = 0.5)) + 
  rotate_x_text(45) 

dev.off()


grDevices::cairo_pdf(file = "figures/Fig4F_CXCL1_expression.pdf", width = 3, height = 4,onefile=T)

# CXCL1 #
ggbarplot(
  taqman_expression_colon_long %>% filter(cytokines == "cxcl1"), 
  x = "genotype_pregnant_status", 
  y = "expression_value",
  add = c("mean_se","jitter"), 
  alpha = 0.5, color = "Genotype",fill = "Genotype",
  #facet.by = "cytokines", scale = "free",
  palette = c("WT"="#5da1d4ff","KO"="#c88cb2ff"), ylab = "Relative expression", xlab = "Pregnancy/Genotype"
) + ggtitle("CXCL1") +
  stat_compare_means( method = "wilcox",
                      comparisons = list(c("P_WT", "P_KO"), c("NP_WT", "NP_KO"),
                                         c("P_WT", "NP_WT"),c("P_WT", "NP_KO"),
                                         c("P_KO", "NP_WT"),c("P_KO", "NP_KO")),
                      label = "p.signif"
  ) +
  theme_bw() + theme(axis.text=element_text(size=7, color="black"), 
                          axis.title=element_text(size=9), 
                          legend.text = element_text(size=9), 
                          legend.title = element_text(size=0), 
                          plot.title = element_text(size=9, face="bold", hjust = 0.5)) + 
  rotate_x_text(45) 




dev.off()


