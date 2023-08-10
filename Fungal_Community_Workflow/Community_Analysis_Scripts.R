# set working directory
setwd("C:/Users/elinrg26/OneDrive - UBC/Research/GitHub/Cover Crop Study")

# Load `tidyverse` and requisite packages from Bioconductor
library(tidyverse)
library(BiocStyle)
library(BiocGenerics)
library(phyloseq)

# Load in the data
load(file = "RDS/results_covercrop_noGreenhouse.rda")
results.df <- readRDS(file = "RDS/results.df_noGreenhouse.rds")

# test <- readRDS(file = "RDS/metadata_covercrop.rds")
# load(file = "RDS/preliminary_results_covercrop.rda")
# 
# for (i in 1:length(RESULTS)){
#   
#   RESULTS[[i]]$results.samples <- RESULTS[[i]]$results.samples[1,1:10]
#   RESULTS[[i]]$results.samples <- cbind(RESULTS[[i]]$results.samples, as.data.frame(test[[i]]))
# }
# 
# # Remove Greenhouse samples [[44:141]]
# RESULTS <- RESULTS[c(1:43, 142:221)]
# # Remove turnip cover crop samples if needed
# RESULTS2 <- RESULTS[c(1:27, 32, 34:72, 77: 123)]
# # Save both objects
# save(RESULTS, file = "RDS/results_covercrop_noGreenhouse.rda")
# save(RESULTS2, file = "RDS/results_covercrop_noTurnip.rda")
# load(file = "RDS/results_covercrop_noGreenhouse.rda")
# load(file = "RDS/results_covercrop_noTurnip.rda")
# 
# RES <- list()
# # fix results.samples
# for (i in 1:length(RESULTS)){
#   
#   RES[[i]] <- RESULTS[[i]]$results.samples[1,1:10]
#   RES[[i]] <- cbind(RES[[i]], as.data.frame(test[[i]]))
# }
# 
# # Build combined data frame of metrics:
# results.df <- do.call(rbind, RES) %>%
#   dplyr::select("sample.ID", everything()) # move sample.ID column to leftmost
# 
# # Remove Greenhouse samples (rows 44:141)
# results.df <- results.df[c(1:43, 142:221),]
# # SuRDC is missing samples for the turnip cover crop which may need to be removed from the analysis
# results.df2 <- results.df[!results.df$cover_crop %in% "Turnip",]
# Save both objects
# saveRDS(results.df, file = "RDS/results.df_noGreenhouse.rds")
# saveRDS(results.df2, file = "RDS/results.df_noTurnip.rds")
# results.df <- readRDS(file = "RDS/results.df_noGreenhouse.rds")
# results.df2 <- readRDS(file = file = "RDS/results.df_noTurnip.rds")

# Load plotting packages
library(ggplot2)
library(cowplot)
library(rphylopic)


# Restructuring data for easier subsequent handling ----


# Extract all phyloseq objects from RESULTS into a list
ps.list <- list()
for (i in 1:length(RESULTS)) {
  ps.list[[i]] <- RESULTS[[i]]$ps
}

## All data ----

# Merge phyloseq OTU tables to return single combined OTU table
ps_C1 <- ps.list[[1]]
merge.otu <- merge_phyloseq(otu_table(ps_C1)) # Merging other OTU tables requires a base to add onto

for (i in 1:length(ps.list)) {
  # Extract each phyloseq object from the list
  ps_otu <- ps.list[[i]]
  merge.otu <- merge_phyloseq(otu_table(ps_otu), merge.otu)
}

# Merge taxonomy 
merge.tax <- merge_phyloseq(tax_table(ps_C1))

for (i in 1:length(ps.list)) {
  # Extract each phyloseq object from the list
  ps_otu <- ps.list[[i]]
  merge.tax <- merge_phyloseq(tax_table(ps_otu), merge.tax)
}

rownames(results.df) <- results.df$sample.ID  # Needed for to create the phyloseq-class object

# Create a new phyloseq-class object with all combined elements
# Note that the `merge_phyloseq()` does not work if phylogenetic trees have different numbers of tips

merge.ps <- phyloseq::phyloseq(otu_table(merge.otu),
                                tax_table(merge.tax),
                                sample_data(results.df))

# Create site-specific lists for phyloseq objects
ps.Cov <- list()
ps.Kal <- list()
ps.Su <- list()

for (i in 1:length(ps.list)) {
  # Sort ps.list by site
  if ("Covert" %in% sample_data(ps.list[[i]])$site) {
    ps.Cov[[i]] <- ps.list[[i]]
  } 
  if ("Kalala" %in% sample_data(ps.list[[i]])$site) {
    ps.Kal[[i]] <- ps.list[[i]]
  } 
  if ("SuRDC" %in% sample_data(ps.list[[i]])$site) {
    ps.Su[[i]] <- ps.list[[i]]
  }
}

# Remove "NULL" elements from lists
ps.Kal <- ps.Kal[c(44:87)]
ps.Su <- ps.Su[c(88:123)]

## Covert only ----

# Merge phyloseq OTU tables to return combined OTU table
ps_C1 <- ps.Cov[[1]]
otu.Cov <- merge_phyloseq(otu_table(ps_C1)) # Merging other OTU tables requires a base to add onto

for (i in 1:length(ps.Cov)) {
  # Extract each phyloseq object from the list
  ps_Cov <- ps.Cov[[i]]
  otu.Cov <- merge_phyloseq(otu_table(ps_Cov), otu.Cov)
}

# Merge taxonomy 
tax.Cov <- merge_phyloseq(tax_table(ps_C1))

for (i in 1:length(ps.Cov)) {
  # Extract each phyloseq object from the list
  ps_Cov <- ps.Cov[[i]]
  tax.Cov <- merge_phyloseq(tax_table(ps_Cov), tax.Cov)
}

# Create a new phyloseq-class object for Covert
merge.Cov <- phyloseq::phyloseq(otu_table(otu.Cov),
                               tax_table(tax.Cov),
                               sample_data(results.df))

## Kalala only ----

# Merge phyloseq OTU tables to return combined OTU table
ps_K1 <- ps.Kal[[1]]
otu.Kal <- merge_phyloseq(otu_table(ps_K1)) # Merging other OTU tables requires a base to add onto

for (i in 1:length(ps.Kal)) {
  # Extract each phyloseq object from the list
  ps_Kal <- ps.Kal[[i]]
  otu.Kal <- merge_phyloseq(otu_table(ps_Kal), otu.Kal)
}

# Merge taxonomy 
tax.Kal <- merge_phyloseq(tax_table(ps_K1))

for (i in 1:length(ps.Kal)) {
  # Extract each phyloseq object from the list
  ps_Kal <- ps.Kal[[i]]
  tax.Kal <- merge_phyloseq(tax_table(ps_Kal), tax.Kal)
}

# Create a new phyloseq-class object for Covert
merge.Kal <- phyloseq::phyloseq(otu_table(otu.Kal),
                                tax_table(tax.Kal),
                                sample_data(results.df))

## SuRDC only ----

# Merge phyloseq OTU tables to return combined OTU table
ps_S1 <- ps.Su[[1]]
otu.Su <- merge_phyloseq(otu_table(ps_S1)) # Merging other OTU tables requires a base to add onto

for (i in 1:length(ps.Su)) {
  # Extract each phyloseq object from the list
  ps_Su <- ps.Su[[i]]
  otu.Su <- merge_phyloseq(otu_table(ps_Su), otu.Su)
}

# Merge taxonomy 
tax.Su <- merge_phyloseq(tax_table(ps_S1))

for (i in 1:length(ps.Su)) {
  # Extract each phyloseq object from the list
  ps_Su <- ps.Su[[i]]
  tax.Su <- merge_phyloseq(tax_table(ps_Su), tax.Su)
}

# Create a new phyloseq-class object for Covert
merge.Su <- phyloseq::phyloseq(otu_table(otu.Su),
                                tax_table(tax.Su),
                                sample_data(results.df))


# Big picture descriptive analyses ----

palette <- c("#c44601", "#FCC9B5", "#E1B239", "#FCF2C7", "#A3D8C6", "#329973", "#7D99E6", "#E0D2EB", "#98669F", "#353A70", "#814B08", "gray60", "black")
# colours = c("#7D99E6", "#92A9EA", "#fcf5c7", "#B6CBE2", "#ACC8E5", "#8babf1", "#b390bc", "#94D2BD", "#E9C46A", "gray60", "#814B08", "#CE8188", "#bf0603", "#6E4007")
palette2 <- c("#005F73", "#0A9396", "#94D2BD", "#7A9EC6", "#E9D8A6", "#fbb13c", "#CA6702", "#9B2226", "#9F72AA", "#6d597a", "#355070", "#EE9B00")

## Boxplots of diversity measures across cover crops and sites ----

### Shannon diversity ----
png(file = "figures/Shannon_boxplot.png", height = 1200, width = 1400)
ggplot(data = results.df) +
  geom_boxplot(aes(y = Shannon, fill = cover_crop)) +  # y = c(PD, rich, Shannon, Simpson)
  labs(fill = "Cover Crop") +
  scale_y_log10() +
  scale_fill_manual(values = palette) +
  facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 20, family = "sans", face = "bold", vjust = 2),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.text.x  = element_text(colour = "transparent"),
        axis.ticks.x = element_line(colour = "transparent"),
        legend.text = element_text(size = 16, family = "sans", face = "bold"),
        legend.title = element_text(size = 20, family = "sans", face = "bold"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.3, fill = NA),
        strip.text.x = element_text(size = 20,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left
dev.off()

### Simpson diversity ----
plot3 <-
  ggplot(data = results.df) +
  geom_boxplot(aes(y = Simpson, fill = cover_crop), size = 0.2, outlier.size = 0.2) +  # y = c(PD, rich, Shannon, Simpson)
  labs(fill = NULL, y = "Simpson Index", tag = "C") +
  scale_fill_manual(values = palette) +
  scale_y_continuous(limits = c(0.77,1.03), expand = c(0,0)) +
  facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 8, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 6, family = "sans"),
        axis.text.x  = element_text(colour = "transparent"),
        axis.ticks.x = element_line(colour = "transparent"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("figures/Simpson_boxplot.png", width = 6.86, height = 5)

### Phylogenetic diversity ----
plot1 <- 
  ggplot(data = results.df) +
  geom_boxplot(aes(y = PD, fill = cover_crop), size = 0.2, outlier.size = 0.2) +  # y = c(PD, rich, Shannon, Simpson)
  scale_y_log10() +
  scale_fill_manual(values = palette) +
  labs(fill = NULL, y = "Phylogenetic Diversity", tag = "B") +
  facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 8, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 6, family = "sans"),
        axis.text.x  = element_text(colour = "transparent"),
        axis.ticks.x = element_line(colour = "transparent"),
        # axis.title.x = element_text(size = 10, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("figures/PD_boxplot2.png", width = 6.86, height = 5)

### Species Richness ----
plot2 <-
  ggplot(data = results.df) +
  geom_boxplot(aes(y = rich, fill = cover_crop), size = 0.2, outlier.size = 0.2) +  # y = c(PD, rich, Shannon, Simpson)
  labs(fill = NULL, y = "Species Richness", tag = "A") +
  scale_fill_manual(values = palette) +
  facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        plot.tag.position = c(0.05, 0.95), # horizontal, vertical
        axis.title.y = element_text(size = 8, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 6, family = "sans"),
        axis.text.x  = element_text(colour = "transparent"),
        axis.ticks.x = element_line(colour = "transparent"),
        legend.text = element_text(size = 4.2, family = "sans", face = "bold"),
        legend.position = c(0.83, 0.86), # horizontal, vertical
        legend.key.height = unit(0.15, "cm"),
        legend.key.width = unit(0.14, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("figures/richness_boxplot.png", width = 6.86, height = 4.5)

# CC_boxplots <- gridExtra::grid.arrange(nrow = 3, plot2, plot1, plot3)
# ggsave("figures/CC_boxplots2.png", plot = CC_boxplots, width = 3.23, height = 7)

CC_boxplots <-
  plot_grid(plot2, plot1, plot3, ncol = 1)
ggsave("figures/CC_boxplots3.png", plot = CC_boxplots, width = 3.23, height = 7)


png(file = "figures/richness_boxplot.png", height = 900, width = 1000)
ggplot(data = results.df, aes(y = rich)) +
  geom_boxplot(aes(fill = cover_crop)) +
  scale_fill_manual(values = palette) +
  facet_wrap(vars(site))
dev.off()


png(file = "figures/richness_vs_PD.png", height = 900, width = 1000)
ggplot(data = results.df, aes(x = rich, y = PD)) +
  geom_point(aes(fill = cover_crop), size = 4, shape = 21) +
  geom_smooth(method = "gam", method.args= list(family = "Gamma"), se = F,
              colour = "black") +
  scale_fill_manual(values = palette) +
  scale_colour_manual("black") +
  facet_wrap(vars(site))
dev.off()

## Beta Diversity/Co-Occurrence Affinity ----

# Co-Occurrence affinity was calculated using the alpha index and CRAN package formulated by Mainali et al. (2022).
# https://doi.org/10.1126/sciadv.abj9204

# Load package for computing co-occurrence affinity (alpha index)
library(CooccurrenceAffinity)

# Agglomerate by species for overall OTU table
ps.species <- tax_glom(merge.ps, taxrank = 'Species', NArm = FALSE)
# Merge samples by site
sitegroup <- merge_samples(otu_table(ps.species), group = sample_data(ps.species)$site, fun = sum)
sitegroup <- t(as.data.frame(otu_table(sitegroup)))
# Replace sequences with species names
taxaps <- as.data.frame(tax_table(ps.species))
rownames(sitegroup) <- taxaps$Species

# Prepare an occurrence matrix based on the OTU table
occur.site <- dataprep(data = sitegroup, row.or.col = "col", datatype = "abundance", threshold = 1, class0.rule = "less")

# Calculate co-occurrence affinity between sites
affinity.site <- CooccurrenceAffinity::affinity(occur.site, row.or.col = "col")

# View the occurrence matrix:
affinity.site$occur_mat
# view all computed outputs:
affinity.site$all

# Merge samples by cover crop
CCgroup <- merge_samples(otu_table(ps.species), group = sample_data(ps.species)$cover_crop, fun = sum)
CCgroup <- t(as.data.frame(otu_table(CCgroup)))
# Replace sequences with species names
rownames(CCgroup) <- taxaps$Species

# Prepare an occurrence matrix based on the OTU table
occur.CC <- dataprep(data = CCgroup, row.or.col = "col", datatype = "abundance", threshold = 1, class0.rule = "less")

# Calculate co-occurrence affinity between sites
affinity.CC <- CooccurrenceAffinity::affinity(occur.CC, row.or.col = "col", squarematrix = c("alpha_mle", "jaccard"))

# View the occurrence matrix:
affinity.CC$occur_mat
# view all computed outputs:
affinity.CC$all

# Add duplicates back into the resulting data frame
DATA.fix <- affinity.CC$all
# Switch the entities and their counts accordingly
DATA.fix <- relocate(DATA.fix, "entity_2", .before = "entity_1")
DATA.fix <- relocate(DATA.fix, "entity_2_count_mB", .before = "entity_1_count_mA")
colnames(DATA.fix)[1] = "entity_1"
colnames(DATA.fix)[2] = "entity_2"
colnames(DATA.fix)[3] = "entity_1_count_mA"
colnames(DATA.fix)[4] = "entity_2_count_mB"

# Remove periods added to cover crop names
DATA.fix$entity_1 <- gsub(".", " ", DATA.fix$entity_1, fixed = TRUE)
DATA.fix$entity_2 <- gsub(".", " ", DATA.fix$entity_2, fixed = TRUE)
affinity.CC$all$entity_1 <- gsub(".", " ", affinity.CC$all$entity_1, fixed = TRUE)
affinity.CC$all$entity_2 <- gsub(".", " ", affinity.CC$all$entity_2, fixed = TRUE)

# Combine the data frames
DATA <- rbind(affinity.CC$all, DATA.fix)

# Visualize alpha statistic by site
# plot(affinity.site$all$alpha_mle ~ affinity.site$all$obs_cooccur_X)

SUMCOUNT <- DATA$entity_1_count_mA + DATA$entity_2_count_mB
ALPHA <- DATA$alpha_mle
CC <- DATA$entity_1

# Random effects model for alpha
fit <- nlme::gls(ALPHA ~ SUMCOUNT, method = "ML", data = as.data.frame(DATA))
fit_random <- nlme::lme(ALPHA ~ SUMCOUNT, random = ~ 1|CC, data = as.data.frame(DATA))
fit_random2 <- nlme::lme(ALPHA ~ SUMCOUNT, random = ~ SUMCOUNT|CC, data = as.data.frame(DATA))

AIC(fit); AIC(fit_random); AIC(fit_random2)

# Visualize alpha by cover crop and overall trend
FIG3 <-
  ggplot(data = DATA, aes(x = SUMCOUNT, y = ALPHA)) +
  geom_point(aes(fill = CC), size = 1, shape = 21, stroke = 0.2) +
  geom_smooth(aes(col = CC), method = "glm", method.args = list(family = Gamma()), se = F, size = 0.5) +
  geom_smooth(aes(col = "All"), method = "glm", method.args = list(family = Gamma()), se = F, size = 1, colour = "black") +
  labs(x = "Species Count per Cover Crop Pair", y = "Co-Occurence Affinity (\u03b1)") +
  scale_x_continuous(limits = c(1600,2130), expand = c(0,0)) +
  scale_y_continuous(limits = c(1.73,2.45), expand = c(0,0)) +
  scale_fill_manual(values = palette) +
  scale_colour_manual(values = palette) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        # plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        axis.title.y = element_text(size = 8, family = "sans", face = "bold", vjust = 1.5),
        axis.text.y = element_text(size = 6, family = "sans"),
        axis.text.x  = element_text(size = 6, family = "sans"),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size = 8, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 4.5, family = "sans", face = "bold"),
        legend.position = "top", # horizontal, vertical
        legend.direction = "horizontal",
        legend.box.spacing = unit(c(0,0,0,0), "cm"),
        legend.key.height = unit(0.02, "cm"),
        legend.key.width = unit(0.02, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.2, fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

ggsave("figures/alpha_affinity_CC.png", plot = FIG3, width = 3.23, height = 4)


# Triangle heatmap to visualize alpha among pairs
FIG4 <-
  ggplot(data = affinity.CC$all, aes(x = entity_2, y = entity_1)) +
  theme_bw() +
  geom_tile(aes(fill = alpha_mle), colour = "gray60") +
  labs(fill = "\u03b1") +
  scale_x_discrete(limits = rev, labels = function(x) stringr::str_wrap(x, width = 10)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  scale_fill_gradient2(low = "#190087", mid = "#E72476", high = "#FAEC50", midpoint = 2.1, space = "Lab") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        # plot.tag.position = c(0.03, 0.93), # horizontal, vertical
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), #element_text(size = 4, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 3, family = "sans", face= "bold", angle = 40),
        axis.text.x  = element_text(size = 3, family = "sans", face = "bold", angle = 40),
        axis.ticks.y = element_line(colour = "transparent"),
        axis.ticks.x = element_line(colour = "transparent"),
        legend.title = element_text(size = 4, family = "sans", face = "bold"),
        legend.text = element_text(size = 3.5, family = "sans", face = "bold"),
        legend.position = c(0.38,0.95),#c(1, 0.45), # horizontal, vertical
        legend.direction = "horizontal",
        legend.margin = margin(t = 0, r = 0, b = 0, l =0, unit = "pt"),
        legend.box.spacing = unit(c(0,0,0,0), "cm"),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.18, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))
        # plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left)

ggsave("figures/alpha_matrix_CC.png", plot = FIG4, width = 3.5, height = 2)

# Inset FIG3 and FIG4
FIG5 <-
  ggdraw(FIG3) +
  draw_plot(FIG4,
    x = 0.098,
    y = 0.04,
    width = 0.6,
    height = 0.47
  )
ggsave("figures/combined_alpha_CC.png", plot = FIG5, width = 3.23, height = 3.85)


### IGNORE: Site-specific affinity between cover crops ----

# Agglomerate by species for Covert-specific OTU table
Cov.species <- tax_glom(merge.Cov, taxrank = 'Species', NArm = FALSE)
# Merge samples by cover crop
Cov.CCgroup <- merge_samples(otu_table(Cov.species), group = sample_data(Cov.species)$cover_crop, fun = sum)
Cov.CCgroup <- t(as.data.frame(otu_table(Cov.CCgroup)))
# Replace sequences with species names
taxaCov <- as.data.frame(tax_table(Cov.species))
rownames(Cov.CCgroup) <- taxaCov$Species

# Prepare an occurrence matrix based on the OTU table
occur.Cov <- dataprep(data = Cov.CCgroup, row.or.col = "col", datatype = "abundance", threshold = 1, class0.rule = "less")

# Calculate co-occurrence affinity (alpha index)
affinity.Cov <- CooccurrenceAffinity::affinity(occur.Cov, row.or.col = "col")

# View the occurrence matrix:
affinity.Cov$occur_mat
# view all computed outputs:
affinity.Cov$all


# Agglomerate by species for Kalala-specific OTU table
Kal.species <- tax_glom(merge.Kal, taxrank = 'Species', NArm = FALSE)
# Merge samples by cover crop
Kal.CCgroup <- merge_samples(otu_table(Kal.species), group = sample_data(Kal.species)$cover_crop, fun = sum)
Kal.CCgroup <- t(as.data.frame(otu_table(Kal.CCgroup)))
# Replace sequences with species names
taxaKal <- as.data.frame(tax_table(Kal.species))
rownames(Kal.CCgroup) <- taxaKal$Species

# Prepare an occurrence matrix based on the OTU table
occur.Kal <- dataprep(data = Kal.CCgroup, row.or.col = "col", datatype = "abundance", threshold = 1, class0.rule = "less")

# Calculate co-occurrence affinity (alpha index)
affinity.Kal <- CooccurrenceAffinity::affinity(occur.Kal, row.or.col = "col")

# View the occurrence matrix:
affinity.Kal$occur_mat
# view all computed outputs:
affinity.Kal$all


# Agglomerate by species for SuRDC-specific OTU table
Su.species <- tax_glom(merge.Su, taxrank = 'Species', NArm = FALSE)
# Merge samples by cover crop
Su.CCgroup <- merge_samples(otu_table(Su.species), group = sample_data(Su.species)$cover_crop, fun = sum)
Su.CCgroup <- t(as.data.frame(otu_table(Su.CCgroup)))
# Replace sequences with species names
taxaSu <- as.data.frame(tax_table(Su.species))
rownames(Su.CCgroup) <- taxaSu$Species

# Prepare an occurrence matrix based on the OTU table
occur.Su <- dataprep(data = Su.CCgroup, row.or.col = "col", datatype = "abundance", threshold = 1, class0.rule = "less")

# Calculate co-occurrence affinity (alpha index)
affinity.Su <- CooccurrenceAffinity::affinity(occur.Su, row.or.col = "col")

# View the occurrence matrix:
affinity.Su$occur_mat
# view all computed outputs:
affinity.Su$all


## Random effects model (gam) ----

# Load package for random effects glm modelling
library(lme4)

# Source the code from the betals.r script
source(file = "Fungal_Community_Workflow/betals.r")


results.df$site <- as.factor(results.df$site)
results.df$cover_crop <- as.factor(results.df$cover_crop)
# GLMM plot
model <- gam(list(PD ~ cover_crop + s(site, bs = "re"),  # random effects for sites, PD as a function of cover crop
               ~ cover_crop + s(site, bs = "re")),
              data = results.df,
              family = gammals())  # gamma location scale distribution for values from 0 to infinity (betals for 0-1)

summary(model)

library("multcomp")
g1 <- glht(model, linfct = mcp(cover_crop = "Tukey"))

# Error in linfct[[nm]] %*% C : requires numeric/complex matrix/vector arguments <-- follow through link on gam modeling



# Network plots ----

# load network plotting packages
library(phyloseq)
library(igraph)

png(file="figures/network_fungi_C28_taxa.png", height = 800, width = 720)
# Network plot (edge colour is not adjustable)
plot_network(RESULTS[[20]]$net, RESULTS[[20]]$ps, shape = "site",
             point_size = 10, line_weight = 0.6,
             line_alpha = 0.8, layout.method = layout.fruchterman.reingold,
             label = "sample.ID", title = NULL) +
  theme(legend.position = "right")
dev.off()

RESULTS[[1]]$net

C1_tree <- RESULTS[[1]]$phylo_tree[["fitGTR"]]

plot(C1_tree$tree, label.tips = NULL)


# Abundance bar plots ----

# Abundance plots were built following the workflow of Hui (2021):
# https://www.yanh.org/2021/01/01/microbiome-r/#build-phyloseq-project
ps.melt <- list()
for (i in 1:length(ps.list)) {
  # Calculate relative sample counts of taxa
  ps.rel <- list()
  ps.rel[[i]] <- transform_sample_counts(ps.list[[i]], function(x) x/sum(x)*100)
  # Agglomerate samples by taxon of choice
  agglomerated <- list()
  agglomerated[[i]] <- tax_glom(ps.rel[[i]], taxrank = 'Class', NArm = FALSE)
  # Melt into individual data frames
  ps.melt[[i]] <- psmelt(agglomerated[[i]])
}

ps.melt2 <- do.call(rbind, ps.melt)

# Organize the data for abundance plots
ps.melt2 <- ps.melt2 %>%
  group_by(Class) %>%
  mutate(median = median(Abundance))
rare <- unique(ps.melt2$Class[ps.melt2$median > 1])
ps.melt2$Class[!(ps.melt2$Class %in% rare)] <- "< 1%"
ps.sum2 <- ps.melt2 %>%
  group_by(sample.ID, cover_crop, Class) %>%
  summarise(Abundance = sum(Abundance))


# Load icons for cover crops
buckwheat <- get_phylopic("325ebaf5-9056-4f7c-9e31-5c9a9a82c755")
buffalo <- get_phylopic("461f7280-3636-42c0-98fd-4fca668460c5")
clover <- get_phylopic("4a86bfd2-2a14-4439-b70e-123a5b90f2ff")
fieldpea <- get_phylopic("a3d3f760-d5c1-45a1-8e6e-5b4797700014")
mustard <- get_phylopic("2b910601-b012-4c9e-896a-aea3ac04fab3")
phacelia <- get_phylopic("29fa5168-0fcc-432f-a37b-fe18da8887f9")
lentil <- get_phylopic("6295da8f-98b2-429a-8588-0e3e48012656")
turnip <- get_phylopic("208cb01d-fe92-4432-b6a6-5fac87deb9bf")
brassica <- get_phylopic("f20144d1-d243-4cca-aba2-24bce6c81d42")


# Abundance bar plot
FIG <-
  ggplot(ps.sum2, aes(x = sample.ID, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", aes(fill = Class)) + 
  labs(x = "", y= "Abundance (%)") +
  facet_wrap(~cover_crop, scales = "free_x", nrow = 1,
             labeller = labeller(cover_crop = label_wrap_gen(width = 10))) +
  theme_classic() +
  scale_fill_manual(values = palette2) + 
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(# plot.tag = element_text(size = 14, family = "sans", face = "bold"),
        # plot.tag.position = c(0.83, 0.85),
        axis.title.y = element_text(size = 10, family = "sans", face = "bold"), # vjust = 2.5),
        axis.text.y = element_text(size = 7, family = "sans"),
        axis.text.x  = element_text(size = 2, angle = 45),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        # legend.position = "none",
        legend.text = element_text(size = 7, family = "sans", face = "bold"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        # legend.position = c(0.82, 0.86), #horizontal, vertical
        legend.key.height = unit(0.32, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        # panel.border = element_rect(linewidth = 0.3, fill = NA),
        strip.text.x = element_text(size = 6, family = "sans", face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("figures/abundance_barplot_Class_CC.png", plot = FIG, width = 6.86, height = 4.5)

icons <-
  ggplot() +
  # geom_blank() +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) +
  add_phylopic(buckwheat, alpha = 1, x = 9.2, y = 88.5, ysize = 8, color = "#c44601") +
  add_phylopic(buffalo, alpha = 1, x = 16.5, y = 88.5, ysize = 8, color = "#FCC9B5") +
  add_phylopic(clover, alpha = 1, x = 25.5, y = 88, ysize = 7, color = "#E1B239") +
  add_phylopic(fieldpea, alpha = 1, x = 33.6, y = 88.5, ysize = 8, color = "#FCF2C7") +
  add_phylopic(mustard, alpha = 1, x = 42, y = 88.5, ysize = 7.8, color = "#A3D8C6") +
  add_phylopic(phacelia, alpha = 1, x = 50.4, y = 88.5, ysize = 7.8, color = "#329973") +
  add_phylopic(lentil, alpha = 1, x = 58, y = 88.5, ysize = 8, color = "#7D99E6") +
  add_phylopic(turnip, alpha = 1, x = 65.9, y = 88.5, ysize = 8, color = "#E0D2EB") +
  add_phylopic(brassica, alpha = 1, x = 74.3, y = 88.5, ysize = 6.5, color = "#98669F") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave("figures/phylopics.png", plot = icons, width = 6.86, height = 4.5)

# Layer cover crop icons on top of abundance bar plot
FIG6 <-
  ggdraw(FIG) +
  draw_plot(icons,
            x = 0,
            y = 0,
            width = 1,
            height = 1
  )
ggsave("figures/abundance_barplot_CC.png", plot = FIG6, width = 6.86, height = 4.5)


# Things to take a closer look at:
# Phylum: Phacelia S32 + C31, Buckwheat S12
# Class: Phacelia C31 + S31 + S32, Buckwheat S12, Spring Lentil S36, Brassica S10
# Genus: Phacelia C31 + S35, Mustard White C24, Brassica SuRDC, Buckwheat S12
# Overall: Phacelia, Buckwheat, and Brassica seem to have more variation in abundance


# Heatmaps ----

# Differential abundance analyses were conducted following the workflow of Hui (2021):
# https://www.yanh.org/2021/01/01/microbiome-r/#build-phyloseq-project

## BASED ON DESeq2 PACKAGE:

# Load `DESeq2` package
library(GenomeInfoDb)
library(DESeq2)

# Agglomerate data by species 
sample_data(merge.ps)$cover_crop <- as.factor(sample_data(merge.ps)$cover_crop) # factorize for DESeq2
ps.species <- tax_glom(merge.ps, taxrank = 'Species', NArm = FALSE)
# ps.class <- tax_glom(merge.ps, taxrank = 'Class', NArm = FALSE)http://127.0.0.1:36915/graphics/plot_zoom_png?width=1200&height=900

## MAYBE NOT NEEDED -- Subset samples by site(?) and cover crop
# ps.sub <- subset_samples(ps.species, cover_crop %in% c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover",
#                                                                  "Field Pea", "Mustard White", "Phacelia", "Spring Lentil",
#                                                                  "Turnip", "Winfred Brassica"))

# Filter features with abundances of zero > 90% of the sequences in a sample
# ps.sub2 <- prune_taxa(rowSums(otu_table(ps.sub) == 0) < ncol(otu_table(ps.sub)) * 0.9, ps.sub) <-- might not be needed
ps.ds = phyloseq_to_deseq2(ps.species, ~ 1)  ## `ps.sub` instead if subsampled by site/cover crop

# Use "poscounts": an alternative estimator on the condition that a gene has some zero counts.
# It calculates a modified geometric mean using the nth root of the product of non-zero counts.
ds <- estimateSizeFactors(ps.ds, type = "poscounts") 
ds = DESeq(ds, test = "Wald", fitType = "parametric")
alpha = 0.05 
res = results(ds, alpha = alpha)
res = res[order(res$padj, na.last = NA), ]

# Keeping only bottom 20 w/ lowest p.adj values, considered significant
# taxa_sig = rownames(res[1:20, ]) # select bottom 20
# ps.species.rel <- transform_sample_counts(ps.species, function(x) x/sum(x)*100)
# ps.rel.sig <- prune_taxa(taxa_sig, ps.species.rel)

# IF ALL KEPT
ps.species.rel <- transform_sample_counts(ps.species, function(x) x/sum(x)*100) # if all kept
# ps.class.rel <- transform_sample_counts(ps.class, function(x) x/sum(x)*100) # if all kept -- class
# ps.phylum <- tax_glom(ps.species.rel, taxrank = 'Phylum', NArm = FALSE)

matrix <- t(as.data.frame(otu_table(ps.rel.sig)))
matrix <- as.matrix(matrix)
rownames(matrix) <- as.character(tax_table(ps.rel.sig)[, "Species"])
results.df_sub <- data.frame(sample_data(ps.rel.sig))

# Define columns and rows by site and cover crops
heatmap_col = data.frame(
  Site = as.factor(results.df_sub$site), 
  `Cover crop` = as.factor(results.df_sub$cover_crop), 
  check.names = FALSE
)

rownames(heatmap_col) = rownames(results.df_sub)

# heatmap_row = data.frame(
  # Phylum = as.factor(tax_table(ps.rel.sig)[, "Phylum"])
# )

# rownames(heatmap_row) = rownames(matrix)

# Define colours to fill heatmap boxes
# phylum_col = RColorBrewer::brewer.pal(length(levels(heatmap_row$Phylum)), "Paired")
# names(phylum_col) = levels(heatmap_row$Phylum)
box_colours = list(
  Site = c(`Covert` = "red", `Kalala` = "blue", `SuRDC` = "yellow"),
  `Cover crop` = c(`Buckwheat` = "#c44601", `Buffalo Grass` = "#FCC9B5", `Crescendo Ladino Clover` = "#E1B239", 
                   `Field Pea` = "#FCF2C7", `Mustard White` = "#A3D8C6", `Phacelia` = "#329973",
                   `Spring Lentil` = "#7D99E6", `Turnip` = "#E0D2EB", `Winfred Brassica` =  "#98669F")
  # Phylum = phylum_col
)

# Heatmap of most significant 20 species grouped by cover crop and phylum
png(file = "figures/heatmap_20species_covercrop.png", height = 670, width = 1600, bg = "transparent")
ComplexHeatmap::pheatmap(matrix, scale = "row", border_color = "gray60",
                         cellwidth = 10, cellheight = 25,
                         fontsize = 13,
                         annotation_col = heatmap_col, 
                         # annotation_row = heatmap_row, 
                         annotation_colors = box_colours)
dev.off()

# Heatmap colun annotations
map_col <- data.frame(`Cover crop` = as.factor(results.df_sub$cover_crop), check.names = FALSE)
map_col2 <- data.frame(Site = as.factor(results.df_sub$site), check.names = FALSE)
col_fill = list(`Cover crop` = c(`Buckwheat` = "#c44601", `Buffalo Grass` = "#FCC9B5", `Crescendo Ladino Clover` = "#E1B239", 
                                  `Field Pea` = "#FCF2C7", `Mustard White` = "#A3D8C6", `Phacelia` = "#329973",
                                  `Spring Lentil` = "#7D99E6", `Turnip` = "#E0D2EB", `Winfred Brassica` =  "#98669F"))
col_fill2 <- list(Site = c(`Covert` = "red", `Kalala` = "blue", `SuRDC` = "yellow"))

# Fixing heatmap
annot_col <- ComplexHeatmap::HeatmapAnnotation(df = map_col, name = "Cover crop", col = col_fill, which = "column",
                                               annotation_name_gp = grid::gpar(fontsize = 15.5))
annot_col2 <- ComplexHeatmap::HeatmapAnnotation(df = map_col2, name = "Site", col = col_fill2, which = "column",
                                                annotation_name_gp = grid::gpar(fontsize = 15.5))

png(file = "figures/heatmap_20species_cc.png", height = 670, width = 1600, bg = "transparent")
# FIG2 <-
ComplexHeatmap::Heatmap(matrix, name = "Abundance", col = circlize::colorRamp2(c(0, 3, 6), c("#190087", "#E72476", "#FAEC50")),
                        border = "gray60",
                        show_row_names = TRUE,
                        row_names_gp = grid::gpar(fontsize = 15.5),
                        show_row_dend = TRUE,
                        row_dend_side = "left",
                        row_dend_width = unit(2, "cm"),
                        show_column_names = TRUE,
                        column_names_gp = grid::gpar(fontsize = 10),
                        show_column_dend = TRUE,
                        column_dend_side = "top",
                        column_dend_height = unit(2, "cm"),
                        # rect_gp = grid::gpar(col = "gray60", lwd = 0.5),
                        top_annotation = c(annot_col, annot_col2),
                        heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 12)))
dev.off()



## Heatmaps by site ----

### Covert ----

# Agglomerate data by species 
sample_data(merge.Cov)$cover_crop <- as.factor(sample_data(merge.Cov)$cover_crop) # factorize for DESeq2
Cov.species <- tax_glom(merge.Cov, taxrank = 'Species', NArm = FALSE)
# Cov.class <- tax_glom(merge.Cov, taxrank = 'Class', NArm = FALSE)

# Filter features with abundances of zero > 90% of the sequences in a sample
# Cov.sub2 <- prune_taxa(rowSums(otu_table(Cov.sub) == 0) < ncol(otu_table(Cov.sub)) * 0.9, Cov.sub) <-- might not be needed
Cov.ds = phyloseq_to_deseq2(Cov.species, ~ 1)  ## `Cov.sub` instead if subsampled by site/cover crop

# Use "poscounts": an alternative estimator on the condition that a gene has some zero counts.
# It calculates a modified geometric mean using the nth root of the product of non-zero counts.
ds.C <- estimateSizeFactors(Cov.ds, type = "poscounts") 
ds.C = DESeq(ds, test = "Wald", fitType = "parametric")
alpha = 0.05 
res.C = results(ds.C, alpha = alpha)
res.C = res.C[order(res.C$padj, na.last = NA), ]

# Keeping only bottom 20 w/ lowest p.adj values, considered significant
taxa_sig.C = rownames(res.C[1:20, ]) # select bottom 20
Cov.species.rel <- transform_sample_counts(Cov.species, function(x) x/sum(x)*100)
Cov.rel.sig <- prune_taxa(taxa_sig.C, Cov.species.rel)

## IF ALL KEPT
# Cov.species.rel <- transform_sample_counts(Cov.species, function(x) x/sum(x)*100) # if all kept
# Cov.class.rel <- transform_sample_counts(Cov.class, function(x) x/sum(x)*100) # if all kept -- class
# Cov.phylum <- tax_glom(Cov.species.rel, taxrank = 'Phylum', NArm = FALSE)

matrix.Cov <- t(as.data.frame(otu_table(Cov.rel.sig)))
matrix.Cov <- as.matrix(matrix.Cov)
rownames(matrix.Cov) <- as.character(tax_table(Cov.rel.sig)[, "Species"])
results.df_Cov <- data.frame(sample_data(Cov.rel.sig))

# Define columns and rows by site and cover crops
heatmap_col.C = data.frame(
  'Cover crop' = as.factor(results.df_Cov$cover_crop), 
  check.names = FALSE
)

rownames(heatmap_col.C) = rownames(results.df_Cov)

# heatmap_row.C = data.frame(
#   Phylum = as.factor(tax_table(Cov.rel.sig)[, "Phylum"])
# )
# 
# rownames(heatmap_row.C) = rownames(matrix.Cov)

# Define colours to fill heatmap boxes
# phylum_col.C = RColorBrewer::brewer.pal(length(levels(heatmap_row.C$Phylum)), "Paired")
# names(phylum_col.C) = levels(heatmap_row.C$Phylum)
box_colours.C = list(
  `Cover crop` = c(`Buckwheat` = "#c44601", `Buffalo Grass` = "#FCC9B5", `Crescendo Ladino Clover` = "#E1B239", 
                   `Field Pea` = "#FCF2C7", `Mustard White` = "#A3D8C6", `Phacelia` = "#329973",
                   `Spring Lentil` = "#7D99E6", `Turnip` = "#E0D2EB", `Winfred Brassica` =  "#98669F")#,
  # Phylum = phylum_col.C
)

# Heatmap of most significant 20 species grouped by cover crop and phylum
png(file = "figures/heatmap_20species_Cov.png", height = 500, width = 1400)
ComplexHeatmap::pheatmap(matrix.Cov, scale= "row", border_color = "gray60",
                         cellwidth = 20, cellheight = 20,
                         fontsize = 15,
                         annotation_col = heatmap_col.C,
                         # annotation_row = heatmap_row.C,
                         annotation_colors = box_colours)
dev.off()

### Kalala ----

# Agglomerate data by species 
sample_data(merge.Kal)$cover_crop <- as.factor(sample_data(merge.Kal)$cover_crop) # factorize for DESeq2
Kal.species <- tax_glom(merge.Kal, taxrank = 'Species', NArm = FALSE)
# Kal.class <- tax_glom(merge.Kal, taxrank = 'Class', NArm = FALSE)

# Filter features with abundances of zero > 90% of the sequences in a sample
# Kal.sub2 <- prune_taxa(rowSums(otu_table(Kal.sub) == 0) < ncol(otu_table(Kal.sub)) * 0.9, Kal.sub) <-- might not be needed
Kal.ds = phyloseq_to_deseq2(Kal.species, ~ 1)  ## `Kal.sub` instead if subsampled by site/cover crop

# Use "poscounts": an alternative estimator on the condition that a gene has some zero counts.
# It calculates a modified geometric mean using the nth root of the product of non-zero counts.
ds.K <- estimateSizeFactors(Kal.ds, type = "poscounts") 
ds.K = DESeq(ds, test = "Wald", fitType = "parametric")
alpha = 0.05 
res.K = results(ds.K, alpha = alpha)
res.K = res.K[order(res.K$padj, na.last = NA), ]

# Keeping only bottom 20 w/ lowest p.adj values, considered significant
taxa_sig.K = rownames(res.K[1:20, ]) # select bottom 20
Kal.species.rel <- transform_sample_counts(Kal.species, function(x) x/sum(x)*100)
Kal.rel.sig <- prune_taxa(taxa_sig.K, Kal.species.rel)

## IF ALL KEPT
# Kal.species.rel <- transform_sample_counts(Kal.species, function(x) x/sum(x)*100) # if all kept
# Kal.class.rel <- transform_sample_counts(Kal.class, function(x) x/sum(x)*100) # if all kept -- class
# Kal.phylum <- tax_glom(Kal.species.rel, taxrank = 'Phylum', NArm = FALSE)

matrix.Kal <- t(as.data.frame(otu_table(Kal.rel.sig)))
matrix.Kal <- as.matrix(matrix.Kal)
rownames(matrix.Kal) <- as.character(tax_table(Kal.rel.sig)[, "Species"])
results.df_Kal <- data.frame(sample_data(Kal.rel.sig))

# Define columns and rows by site and cover crops
heatmap_col.K = data.frame(
  'Cover crop' = as.factor(results.df_Kal$cover_crop), 
  check.names = FALSE
)

rownames(heatmap_col.K) = rownames(results.df_Kal)

# heatmap_row.K = data.frame(
#   Phylum = as.factor(tax_table(Kal.rel.sig)[, "Phylum"])
# )

# rownames(heatmap_row.K) = rownames(matrix.Kal)

# Define colours to fill heatmap boxes
# phylum_col.K = RColorBrewer::brewer.pal(length(levels(heatmap_row.K$Phylum)), "Paired")
# names(phylum_col.K) = levels(heatmap_row.K$Phylum)
box_colours.K = list(
  `Cover crop` = c(`Buckwheat` = "#c44601", `Buffalo Grass` = "#FCC9B5", `Crescendo Ladino Clover` = "#E1B239", 
                   `Field Pea` = "#FCF2C7", `Mustard White` = "#A3D8C6", `Phacelia` = "#329973",
                   `Spring Lentil` = "#7D99E6", `Turnip` = "#E0D2EB", `Winfred Brassica` =  "#98669F")
  # Phylum = phylum_col.K
)

# Heatmap of most significant 20 species grouped by cover crop and phylum
png(file = "figures/heatmap_20species_Kal.png", height = 500, width = 1200)
ComplexHeatmap::pheatmap(matrix.Kal, scale= "row", border_color = "gray60",
                         cellwidth = 15, cellheight = 15,
                         fontsize = 16,
                         annotation_col = heatmap_col.K, 
                         # annotation_row = heatmap_row.K, 
                         annotation_colors = box_colours)
dev.off()

### SuRDC ----

# Agglomerate data by species 
sample_data(merge.Su)$cover_crop <- as.factor(sample_data(merge.Su)$cover_crop) # factorize for DESeq2
Su.species <- tax_glom(merge.Su, taxrank = 'Species', NArm = FALSE)
# Su.class <- tax_glom(merge.Su, taxrank = 'Class', NArm = FALSE)

# Filter features with abundances of zero > 90% of the sequences in a sample
# Su.sub2 <- prune_taxa(rowSums(otu_table(Su.sub) == 0) < ncol(otu_table(Su.sub)) * 0.9, Su.sub) <-- might not be needed
Su.ds = phyloseq_to_deseq2(Su.species, ~ 1)  ## `Su.sub` instead if subsampled by site/cover crop

# Use "poscounts": an alternative estimator on the condition that a gene has some zero counts.
# It calculates a modified geometric mean using the nth root of the product of non-zero counts.
ds.S <- estimateSizeFactors(Su.ds, type = "poscounts") 
ds.S = DESeq(ds, test = "Wald", fitType = "parametric")
alpha = 0.05 
res.S = results(ds.S, alpha = alpha)
res.S = res.S[order(res.S$padj, na.last = NA), ]

# Keeping only bottom 20 w/ lowest p.adj values, considered significant
taxa_sig.S = rownames(res.S[1:20, ]) # select bottom 20
Su.species.rel <- transform_sample_counts(Su.species, function(x) x/sum(x)*100)
Su.rel.sig <- prune_taxa(taxa_sig.S, Su.species.rel)

## IF ALL KEPT
# Su.species.rel <- transform_sample_counts(Su.species, function(x) x/sum(x)*100) # if all kept
# Su.class.rel <- transform_sample_counts(Su.class, function(x) x/sum(x)*100) # if all kept -- class
# Su.phylum <- tax_glom(Su.species.rel, taxrank = 'Phylum', NArm = FALSE)

matrix.Su <- t(as.data.frame(otu_table(Su.rel.sig)))
matrix.Su <- as.matrix(matrix.Su)
rownames(matrix.Su) <- as.character(tax_table(Su.rel.sig)[, "Species"])
results.df_Su <- data.frame(sample_data(Su.rel.sig))

# Define columns and rows by site and cover crops
heatmap_col.S = data.frame(
  'Cover crop' = as.factor(results.df_Su$cover_crop), 
  check.names = FALSE
)

rownames(heatmap_col.S) = rownames(results.df_Su)

# heatmap_row.S = data.frame(
#   Phylum = as.factor(tax_table(Su.rel.sig)[, "Phylum"])
# )
# 
# rownames(heatmap_row.S) = rownames(matrix.Su)

# Define colours to fill heatmap boxes
# phylum_col.S = RColorBrewer::brewer.pal(length(levels(heatmap_row.S$Phylum)), "Paired")
# names(phylum_col.S) = levels(heatmap_row.S$Phylum)
box_colours.S = list(
  `Cover crop` = c(`Buckwheat` = "#c44601", `Buffalo Grass` = "#FCC9B5", `Crescendo Ladino Clover` = "#E1B239", 
                   `Field Pea` = "#FCF2C7", `Mustard White` = "#A3D8C6", `Phacelia` = "#329973",
                   `Spring Lentil` = "#7D99E6", `Turnip` = "#E0D2EB", `Winfred Brassica` =  "#98669F")
  # Phylum = phylum_col.S
)

# Heatmap of most significant 20 species grouped by cover crop and phylum
png(file = "figures/heatmap_20species_Su.png", height = 500, width = 1200)
ComplexHeatmap::pheatmap(matrix.Su, scale= "row", border_color = "gray60",
                         cellwidth = 15, cellheight = 15,
                         fontsize = 16,
                         annotation_col = heatmap_col.S, 
                         # annotation_row = heatmap_row.S, 
                         annotation_colors = box_colours)
dev.off()



