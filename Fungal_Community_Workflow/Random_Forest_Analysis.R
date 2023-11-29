# Load in required packages
library(ggplot2)
library(ggridges)
library(cowplot)
library(phyloseq)
library(randomForest)
library(ellipse)
library(FactoMineR)
library(grid)
library(gridExtra)

# Load in the data
load(file = "RDS/results_covercrop_noGreenhouse.rda")
results.df <- readRDS(file = "RDS/results.df_noGreenhouse.rds")

palette <- c("#c44601", "#FCC9B5", "#E1B239", "#FCF2C7", "#A3D8C6", "#329973", "#7D99E6", "#E0D2EB", "#98669F", "#353A70", "#814B08", "gray60", "black")

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

merge.ps <- phyloseq(otu_table(merge.otu),
                               tax_table(merge.tax),
                               sample_data(results.df))


# Random Forest ASVs ----

# Build matrix of species abundances to feed into random forest model
data <- as.data.frame(otu_table(merge.ps))
results.df$cover_crop <- as.factor(results.df$cover_crop)
results.df$site <- as.factor(results.df$site)

# Random forest models were constructed and visualized following the workflow of Milazzotto et al. (2022):
# https://doi.org/10.1016/j.isci.2022.103904

# Random forest model classifying ASVs by cover crop
# model.cc <- randomForest(y = results.df$cover_crop,
#                          x = data,
#                          mtry = 5,
#                          ntree = 20000,
#                          importance = TRUE,
#                          proximity = TRUE,
#                          keep.forest = TRUE,
#                          replace = TRUE)
# model.cc
# saveRDS(model.cc, file = "RDS/random_forest_CC.rds")
model.cc <- readRDS(file = "RDS/random_forest_CC.rds")

# Call:
#   randomForest(x = data, y = results.df$cover_crop, ntree = 20000, mtry = 5, replace = TRUE, importance = TRUE, proximity = TRUE, keep.forest = TRUE) 
# Type of random forest: classification
# Number of trees: 20000
# No. of variables tried at each split: 5
# 
# OOB estimate of  error rate: 98.37%
# Confusion matrix:
#                         Buckwheat Buffalo Grass Crescendo Ladino Clover Field Pea Mustard White Phacelia Spring Lentil Turnip Winfred Brassica class.error
# Buckwheat                   0           5                   2               0           4          1            2         0          1          1.0000000
# Buffalo Grass               0           1                   0               0           9          2            2         0          0          0.9285714
# Crescendo Ladino Clover     1           5                   0               1           4          1            1         0          0          1.0000000
# Field Pea                   3           3                   0               0           4          2            3         0          0          1.0000000
# Mustard White               2           8                   1               0           0          1            3         0          0          1.0000000
# Phacelia                    2           6                   1               2           3          0            0         0          1          1.0000000
# Spring Lentil               1           6                   1               1           2          3            1         0          0          0.9333333
# Turnip                      0           4                   0               0           4          1            0         0          0          1.0000000
# Winfred Brassica            1           2                   1               1           4          0            3         0          0          1.0000000

# Variable importance plot
varImpPlot(model.cc, type = 1, scale = FALSE, bg = "transparent")

# Random forest model classifying ASVs by site
# model.site <- randomForest(y = results.df$site,
#                          x = data,
#                          mtry = 5,
#                          ntree = 20000,
#                          importance = TRUE,
#                          proximity = TRUE,
#                          keep.forest = TRUE,
#                          replace = TRUE)
# model.site
# saveRDS(model.site, file = "RDS/random_forest_site.rds")
model.site <- readRDS(file = "RDS/random_forest_site.rds")

# Call:
#   randomForest(x = data, y = results.df$site, ntree = 20000, mtry = 5, replace = TRUE, importance = TRUE, proximity = TRUE, keep.forest = TRUE) 
# Type of random forest: classification
# Number of trees: 20000
# No. of variables tried at each split: 5
# 
# OOB estimate of  error rate: 13.01%
# Confusion matrix:
#         Covert Kalala SuRDC class.error
# Covert     36      7     0  0.16279070
# Kalala      1     36     7  0.18181818
# SuRDC       1      0    35  0.02777778

# Variable importance plot
varImpPlot(model.site, type = 1, scale = FALSE, bg = "transparent")


# Random Forest Species ----

# Agglomerate data by species
ps.species <- tax_glom(merge.ps, taxrank = 'Species', NArm = FALSE)
sample_data(ps.species)$cover_crop <- as.factor(sample_data(ps.species)$cover_crop)
sample_data(ps.species)$site <- as.factor(sample_data(ps.species)$site)

# Build matrix of species abundances to feed into random forest model
data2 <- as.data.frame(otu_table(ps.species))
colnames(data2) <- as.character(paste(tax_table(ps.species)[, "Genus"], tax_table(ps.species)[, "Species"]))

# Random forest models were constructed and visualized following the workflow of Milazzotto et al. (2022):
# https://doi.org/10.1016/j.isci.2022.103904

# Random forest model classifying species by cover crop
# model.cc2 <- randomForest(y = sample_data(ps.species)$cover_crop,
#                           x = data2,
#                           mtry = 5,
#                           ntree = 20000,
#                           importance = TRUE,
#                           proximity = TRUE,
#                           keep.forest = TRUE,
#                           replace = TRUE)
# model.cc2
# saveRDS(model.cc2, file = "RDS/random_forest_CC2.rds")
model.cc2 <- readRDS(file = "RDS/random_forest_CC2.rds")

# Call:
#   randomForest(x = data2, y = sample_data(ps.species)$cover_crop, ntree = 20000, mtry = 5, replace = TRUE, importance = TRUE, proximity = TRUE, keep.forest = TRUE)
# Type of random forest: classification
# Number of trees: 20000
# No. of variables tried at each split: 5
# 
# OOB estimate of  error rate: 98.37%
# Confusion matrix:
#                         Buckwheat Buffalo Grass Crescendo Ladino Clover Field Pea Mustard White Phacelia Spring Lentil Turnip Winfred Brassica class.error
# Buckwheat                   0           2                   2                2           3          2          3         0            1         1.0000000
# Buffalo Grass               4           1                   2                1           1          4          1         0            0         0.9285714
# Crescendo Ladino Clover     4           1                   0                1           6          0          1         0            0         1.0000000
# Field Pea                   3           2                   1                0           4          1          2         0            2         1.0000000
# Mustard White               6           0                   3                3           0          2          1         0            0         1.0000000
# Phacelia                    2           2                   1                2           3          0          3         0            2         1.0000000
# Spring Lentil               4           4                   0                1           1          0          1         3            1         0.9333333
# Turnip                      1           1                   0                1           0          2          4         0            0         1.0000000
# Winfred Brassica            3           0                   2                2           2          1          1         1            0         1.0000000

# Variable importance plot
varImpPlot(model.cc2, type = 1, scale = FALSE, bg = "transparent")

# Random forest model classifying species by site
# model.site2 <- randomForest(y = sample_data(ps.species)$site,
#                             x = data2,
#                             mtry = 5,
#                             ntree = 20000,
#                             importance = TRUE,
#                             proximity = TRUE,
#                             keep.forest = TRUE,
#                             replace = TRUE)
# model.site2
# saveRDS(model.site2, file = "RDS/random_forest_site2.rds")
model.site2 <- readRDS(file = "RDS/random_forest_site2.rds")

# Call:
#   randomForest(x = data2, y = sample_data(ps.species)$site, ntree = 20000, mtry = 5, replace = TRUE, importance = TRUE, proximity = TRUE, keep.forest = TRUE) 
# Type of random forest: classification
# Number of trees: 20000
# No. of variables tried at each split: 5
# 
# OOB estimate of  error rate: 0.81%
# Confusion matrix:
#         Covert Kalala SuRDC class.error
# Covert     43      0     0  0.00000000
# Kalala      0     43     1  0.02272727
# SuRDC       0      0    36  0.00000000

# Variable importance plot
varImpPlot(model.site2, type = 1, scale = FALSE, bg = "transparent")


# Principal Component Analysis ----

## ASVs ----
# Select top 20 ASVs for cover crop PCA
top_ASVs <- row.names(model.cc$importance)[order(model.cc$importance[,"MeanDecreaseAccuracy"],
                                                      decreasing = TRUE)][1:20]
top_ASVs

# Conduct a PCA on the proximity matrix from random forest
pca.cc <- PCA(model.cc$proximity, graph = FALSE) 
PC1.cc <- pca.cc$ind$coord[,1] # Store individual coordinates of PC1 as a vector
PC2.cc <- pca.cc$ind$coord[,2] # Store individual coordinates of PC2 as a vector
PCs.ID <- data.frame(cbind(PC1.cc, PC2.cc)) # Bind the coordinates together as a data frame
PCs.ID$cover_crop <- results.df$cover_crop # Add in cover crops to the data frame

# Define axis labels based on % data explained across each dimension of PCA
DIM_1.cc <- paste("Dim 1 (", round(pca.cc$eig[1,2], 1), "%)")
DIM_2.cc <- paste("Dim 2 (", round(pca.cc$eig[2,2], 1), "%)")

# PCA figure
PCA_FIG.cc <- 
  ggplot(data = PCs.ID, aes(x = PC1.cc, y = PC2.cc, colour = cover_crop, fill = cover_crop)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", lwd = 0.2, col = "grey60", alpha = 0.8) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", lwd = 0.2, col = "grey60", alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.2, segments = 200, linewidth = 0.2, show.legend = FALSE) +
  theme_bw() +
  ylab(DIM_2.cc) +
  xlab(DIM_1.cc) + 
  geom_point(size = 0.2, aes(colour = cover_crop)) +
  scale_color_manual(labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                              "Spring Lentil", "Turnip", "Winfred Brassica"),
                     values = palette) +
  scale_fill_manual(labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                             "Spring Lentil", "Turnip", "Winfred Brassica"),
                    values = palette,
                    guide = "none") +
  ylab(DIM_2.cc) +
  xlab(DIM_1.cc) + 
  labs(tag = "A") +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        axis.title.x  = element_text(size = 8, family = "sans", face = "bold"),
        axis.title.y  = element_text(size = 8, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 6, family = "sans"),
        axis.text.x  = element_text(size = 6, family = "sans"),
        legend.position = "top",  # horizontal, vertical
        legend.direction = "horizontal",
        legend.box.spacing = unit(c(0,0,0,0), "cm"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        legend.key.size = unit(0.2, "cm"),
        legend.key = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(file = "figures/PCA_CC.png", plot = PCA_FIG.cc, width = 3.23, height = 3.5)


# Select top 20 ASVs for site PCA
top_ASVs2 <- row.names(model.site$importance)[order(model.site$importance[,"MeanDecreaseAccuracy"],
                                                 decreasing = TRUE)][1:20]
top_ASVs2

# Conduct a PCA on the proximity matrix from random forest
pca.site <- PCA(model.site$proximity, graph = FALSE) 
PC1.site <- pca.site$ind$coord[,1] # Store individual coordinates of PC1 as a vector
PC2.site <- pca.site$ind$coord[,2] # Store individual coordinates of PC2 as a vector
PCs.ID2 <- data.frame(cbind(PC1.site, PC2.site)) # Bind the coordinates together as a data frame
PCs.ID2$site <- results.df$site # Add in sites to the data frame

# Define axis labels based on % data explained across each dimension of PCA
DIM_1.site <- paste("Dim 1 (", round(pca.site$eig[1,2], 1), "%)")
DIM_2.site <- paste("Dim 2 (", round(pca.site$eig[2,2], 1), "%)")

# PCA figure
PCA_FIG.site <- 
  ggplot(data = PCs.ID2, aes(x = PC1.site, y = PC2.site, colour = site, fill = site)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", lwd = 0.2, col = "grey60", alpha = 0.8) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", lwd = 0.2, col = "grey60", alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.2, segments = 200, linewidth = 0.2, show.legend = FALSE) +
  theme_bw() +
  ylab(DIM_2.site) +
  xlab(DIM_1.site) + 
  geom_point(size = 0.2, aes(colour = site)) +
  scale_color_manual(labels = c("Covert", "Kalala", "SuRDC"),
                     values = c("red", "blue","gold")) +
  scale_fill_manual(labels = c("Covert", "Kalala", "SuRDC"),
                    values = c("red", "blue","gold"),
                    guide = "none") +
  ylab(DIM_2.site) +
  xlab(DIM_1.site) + 
  labs(tag = "B") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.98), # horizontal, vertical
        axis.title.x  = element_text(size = 8, family = "sans", face = "bold"),
        axis.title.y  = element_text(size = 8, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 6, family = "sans"),
        axis.text.x  = element_text(size = 6, family = "sans"),
        legend.position = "top",  # horizontal, vertical
        legend.direction = "horizontal",
        legend.box.spacing = unit(c(0,0,0,0), "cm"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        legend.key.size = unit(0.2, "cm"),
        legend.key = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(file = "figures/PCA_site.png", plot = PCA_FIG.site, width = 3.23, height = 3.3)

PCA_ASVs <- plot_grid(PCA_FIG.cc, PCA_FIG.site, nrow = 2, rel_heights = c(1,0.9))
ggsave(file = "figures/PCA_ASVs.png", plot = PCA_ASVs, width = 3.23, height = 4.7)


## Species ----
# Select top 20 species for cover crop PCA
top_specs <- row.names(model.cc2$importance)[order(model.cc2$importance[,"MeanDecreaseAccuracy"],
                                                   decreasing = TRUE)][1:20]
top_specs

# Conduct a PCA on the proximity matrix from random forest
pca.cc2 <- PCA(model.cc2$proximity, graph = FALSE) 
PC1.cc2 <- pca.cc2$ind$coord[,1] # Store individual coordinates of PC1 as a vector
PC2.cc2 <- pca.cc2$ind$coord[,2] # Store individual coordinates of PC2 as a vector
PCs.ID3 <- data.frame(cbind(PC1.cc2, PC2.cc2)) # Bind the coordinates together as a data frame
PCs.ID3$cover_crop <- sample_data(ps.species)$cover_crop # Add in cover crops to the data frame

# Define axis labels based on % data explained across each dimension of PCA
DIM_1.cc2 <- paste("Dim 1 (", round(pca.cc2$eig[1,2], 1), "%)")
DIM_2.cc2 <- paste("Dim 2 (", round(pca.cc2$eig[2,2], 1), "%)")

# PCA figure
PCA_FIG.cc2 <- 
  ggplot(data = PCs.ID3, aes(x = PC1.cc2, y = PC2.cc2, colour = cover_crop, fill = cover_crop)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", lwd = 0.2, col = "grey60", alpha = 0.8) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", lwd = 0.2, col = "grey60", alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.2, segments = 200, linewidth = 0.2, show.legend = FALSE) +
  theme_bw() +
  ylab(DIM_2.cc2) +
  xlab(DIM_1.cc2) + 
  geom_point(size = 0.2, aes(colour = cover_crop)) +
  scale_color_manual(labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                                "Spring Lentil", "Turnip", "Winfred Brassica"),
                     values = palette) +
  scale_fill_manual(labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                               "Spring Lentil", "Turnip", "Winfred Brassica"),
                    values = palette,
                    guide = "none") +
  ylab(DIM_2.cc2) +
  xlab(DIM_1.cc2) + 
  labs(tag = "A") +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        plot.tag.position = c(0.025, 0.93), # horizontal, vertical
        axis.title.x  = element_text(size = 8, family = "sans", face = "bold"),
        axis.title.y  = element_text(size = 8, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 6, family = "sans"),
        axis.text.x  = element_text(size = 6, family = "sans"),
        legend.position = "top",  # horizontal, vertical
        legend.direction = "horizontal",
        legend.box.spacing = unit(c(0,0,0,0), "cm"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        legend.key.size = unit(0.2, "cm"),
        legend.key = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(file = "figures/PCA_CC2.png", plot = PCA_FIG.cc2, width = 3.23, height = 3.5)


# Select top 20 species for site PCA
top_specs2 <- row.names(model.site2$importance)[order(model.site2$importance[,"MeanDecreaseAccuracy"],
                                                     decreasing = TRUE)][1:20]
top_specs2

# Conduct a PCA on the proximity matrix from random forest
pca.site2 <- PCA(model.site2$proximity, graph = FALSE) 
PC1.site2 <- pca.site2$ind$coord[,1] # Store individual coordinates of PC1 as a vector
PC2.site2 <- pca.site2$ind$coord[,2] # Store individual coordinates of PC2 as a vector
PCs.ID4 <- data.frame(cbind(PC1.site2, PC2.site2)) # Bind the coordinates together as a data frame
PCs.ID4$site <- sample_data(ps.species)$site # Add in sites to the data frame

# Define axis labels based on % data explained across each dimension of PCA
DIM_1.site2 <- paste("Dim 1 (", round(pca.site2$eig[1,2], 1), "%)")
DIM_2.site2 <- paste("Dim 2 (", round(pca.site2$eig[2,2], 1), "%)")

# PCA figure
PCA_FIG.site2 <- 
  ggplot(data = PCs.ID4, aes(x = PC1.site2, y = PC2.site2, colour = site, fill = site)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", lwd = 0.2, col = "grey60", alpha = 0.8) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", lwd = 0.2, col = "grey60", alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.2, segments = 200, linewidth = 0.2, show.legend = FALSE) +
  theme_bw() +
  ylab(DIM_2.site2) +
  xlab(DIM_1.site2) + 
  geom_point(size = 0.2, aes(colour = site)) +
  scale_color_manual(labels = c("Covert", "Kalala", "SuRDC"),
                     values = c("red", "blue","gold")) +
  scale_fill_manual(labels = c("Covert", "Kalala", "SuRDC"),
                    values = c("red", "blue","gold"),
                    guide = "none") +
  ylab(DIM_2.site2) +
  xlab(DIM_1.site2) + 
  labs(tag = "C") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        plot.tag.position = c(0.025, 0.95), # horizontal, vertical
        axis.title.x  = element_text(size = 8, family = "sans", face = "bold"),
        axis.title.y  = element_text(size = 8, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 6, family = "sans"),
        axis.text.x  = element_text(size = 6, family = "sans"),
        legend.position = "top",  # horizontal, vertical
        legend.direction = "horizontal",
        legend.box.spacing = unit(c(0,0,0,0), "cm"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        legend.key.size = unit(0.2, "cm"),
        legend.key = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5, fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(file = "figures/PCA_site2.png", plot = PCA_FIG.site2, width = 3.23, height = 3.3)

PCA_species <- plot_grid(PCA_FIG.cc2, PCA_FIG.site2, nrow = 2, rel_heights = c(1,0.9))
ggsave(file = "figures/PCA_species.png", plot = PCA_species, width = 3.23, height = 4.8)


# Density Plots ----

## ASVs ----
# Density plots take a closer look at the top 3 ASVs for cover crop PCA
ASV_FIG1a <- 
  ggplot(data, aes(x = data[,top_ASVs[1]], y = results.df$cover_crop, fill = results.df$cover_crop)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                               "Spring Lentil", "Turnip", "Winfred Brassica"),
                    values = palette,
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                                                  "Spring Lentil", "Turnip", "Winfred Brassica")) +
  labs(x = "ASV 1 Expression", y = "Cover crop", tag = "i)") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 4, family = "sans", face = "bold"),
        axis.text.x  = element_text(size = 4, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

ASV_FIG1b <- 
  ggplot(data, aes(x = data[,top_ASVs[2]], y = results.df$cover_crop, fill = results.df$cover_crop)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                               "Spring Lentil", "Turnip", "Winfred Brassica"),
                    values = palette,
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                                                  "Spring Lentil", "Turnip", "Winfred Brassica")) +
  labs(x = "ASV 2 Expression", y = "Cover crop", tag = "ii)") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 4, family = "sans", face = "bold"),
        axis.text.x  = element_text(size = 4, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

ASV_FIG1c <- 
  ggplot(data, aes(x = data[,top_ASVs[3]], y = results.df$cover_crop, fill = results.df$cover_crop)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                               "Spring Lentil", "Turnip", "Winfred Brassica"),
                    values = palette,
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                                                  "Spring Lentil", "Turnip", "Winfred Brassica")) +
  labs(x = "ASV 3 Expression", y = "Cover crop", tag = "iii)") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 4, family = "sans", face = "bold"),
        axis.text.x  = element_text(size = 4, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

ASV_FIG1 <- plot_grid(ASV_FIG1a, ASV_FIG1b, ASV_FIG1c, ncol = 1)
ggsave(file = "figures/density_ASV1.png", plot = ASV_FIG1, width = 3.23, height = 4.8)


# Density plots take a closer look at the top 3 ASVs for site PCA
ASV_FIG2a <- 
  ggplot(data, aes(x = data[,top_ASVs2[1]], y = results.df$site, fill = results.df$site)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Covert", "Kalala", "SuRDC"),
                    values = c("red", "blue","gold"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Covert", "Kalala", "SuRDC")) +
  labs(x = "ASV 1 Expression", y = "Site", tag = "i)") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size = 5, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

ASV_FIG2b <- 
  ggplot(data, aes(x = data[,top_ASVs2[2]], y = results.df$site, fill = results.df$site)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Covert", "Kalala", "SuRDC"),
                    values = c("red", "blue","gold"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Covert", "Kalala", "SuRDC")) +
  labs(x = "ASV 2 Expression", y = "Site", tag = "ii)") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size = 5, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

ASV_FIG2c <- 
  ggplot(data, aes(x = data[,top_ASVs2[3]], y = results.df$site, fill = results.df$site)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Covert", "Kalala", "SuRDC"),
                    values = c("red", "blue","gold"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Covert", "Kalala", "SuRDC")) +
  labs(x = "ASV 3 Expression", y = "Site", tag = "iii)") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size = 5, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

ASV_FIG2 <- plot_grid(ASV_FIG2a, ASV_FIG2b, ASV_FIG2c, ncol = 1)
ggsave(file = "figures/density_ASV2.png", plot = ASV_FIG2, width = 3.23, height = 4.8)

ASV_FIG <- plot_grid(ASV_FIG1a, ASV_FIG1b, ASV_FIG2a, ASV_FIG2b, nrow = 4, rel_heights = c(1,1,0.9,0.9))
ggsave(file = "figures/density_ASV.png", plot = ASV_FIG, width = 3.23, height = 4.8)

combined_ASV <- plot_grid(PCA_ASVs, ASV_FIG, nrow = 1)
ggsave(file = "figures/combined_ASV.png", plot = combined_ASV, width = 6.86, height = 4.8)


## Species ----
# Density plots take a closer look at the top 3 species for cover crop PCA
spec_FIG1a <- 
  ggplot(data2, aes(x = data2[,top_specs[1]], y = sample_data(ps.species)$cover_crop, fill = sample_data(ps.species)$cover_crop)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                               "Spring Lentil", "Turnip", "Winfred Brassica"),
                    values = palette,
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                                                  "Spring Lentil", "Turnip", "Winfred Brassica")) +
  labs(x = paste(top_specs[1], "Abundance"), y = "Cover crop", tag = "B") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 4, family = "sans"),
        axis.text.x  = element_text(size = 4, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

spec_FIG1b <- 
  ggplot(data2, aes(x = data2[,top_specs[2]], y = sample_data(ps.species)$cover_crop, fill = sample_data(ps.species)$cover_crop)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                               "Spring Lentil", "Turnip", "Winfred Brassica"),
                    values = palette,
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                                                  "Spring Lentil", "Turnip", "Winfred Brassica")) +
  labs(x = paste("Unclassified_Phaeosphaeriaceae.1", "Abundance"), y = "Cover crop", tag = "C") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 4, family = "sans"),
        axis.text.x  = element_text(size = 4, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

spec_FIG1c <- 
  ggplot(data2, aes(x = data2[,top_specs[3]], y = sample_data(ps.species)$cover_crop, fill = sample_data(ps.species)$cover_crop)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                               "Spring Lentil", "Turnip", "Winfred Brassica"),
                    values = palette,
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Buckwheat", "Buffalo Grass", "Crescendo Ladino Clover", "Field Pea", "Mustard White", "Phacelia",
                                                  "Spring Lentil", "Turnip", "Winfred Brassica")) +
  labs(x = paste("Unclassified_Fungi.1", "Expression"), y = "Cover crop", tag = "iii)") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 4, family = "sans", face = "bold"),
        axis.text.x  = element_text(size = 4, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

spec_FIG1 <- plot_grid(spec_FIG1a, spec_FIG1b, spec_FIG1c, ncol = 1)
ggsave(file = "figures/density_spec1.png", plot = spec_FIG1, width = 3.23, height = 4.8)


# Density plots take a closer look at the top 3 species for site PCA
spec_FIG2a <- 
  ggplot(data2, aes(x = data2[,top_specs2[1]], y = sample_data(ps.species)$site, fill = sample_data(ps.species)$site)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Covert", "Kalala", "SuRDC"),
                    values = c("red", "blue","gold"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Covert", "Kalala", "SuRDC")) +
  labs(x = paste(top_specs2[1], "Expression"), y = "Site", tag = "E") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size = 5, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

spec_FIG2b <- 
  ggplot(data2, aes(x = data2[,top_specs2[2]], y = sample_data(ps.species)$site, fill = sample_data(ps.species)$site)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Covert", "Kalala", "SuRDC"),
                    values = c("red", "blue","gold"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Covert", "Kalala", "SuRDC")) +
  labs(x = paste("Unclassified_Microascales_gen_Incertae_sedis", "Expression"), y = "Site", tag = "F") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size = 5, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

spec_FIG2c <- 
  ggplot(data2, aes(x = data2[,top_specs2[3]], y = sample_data(ps.species)$site, fill = sample_data(ps.species)$site)) +
  geom_density_ridges(scale = 5, alpha = 0.6, size = 0.2) +
  theme_ridges() +
  scale_fill_manual(labels = c("Covert", "Kalala", "SuRDC"),
                    values = c("red", "blue","gold"),
                    guide = "none") +
  scale_y_discrete(expand = c(0.1, 0), labels = c("Covert", "Kalala", "SuRDC")) +
  labs(x = paste(top_specs2[3], "Expression"), y = "Site", tag = "iii)") +
  theme(plot.tag = element_text(size = 8, family = "sans", face = "bold"),
        plot.tag.position = c(0.09, 0.93), # horizontal, vertical
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.title.y  = element_text(hjust = 0.5, size = 6, family = "sans", face = "bold"),
        axis.text.y  = element_text(size = 5, family = "sans", face = "bold"),
        axis.text.x  = element_text(size = 5, family = "sans"),
        panel.border = element_rect(colour = "black", linewidth = 1, fill = "transparent"),
        legend.position = "none")

spec_FIG2 <- plot_grid(spec_FIG2a, spec_FIG2b, spec_FIG2c, ncol = 1)
ggsave(file = "figures/density_spec2.png", plot = spec_FIG2, width = 3.23, height = 4.8)

spec_FIG <- plot_grid(spec_FIG1a, spec_FIG1b, spec_FIG2a, spec_FIG2b, nrow = 4, rel_heights = c(1,1,0.9,0.9))
ggsave(file = "figures/density_spec.png", plot = spec_FIG, width = 3.23, height = 4.8)

combined_specs <- plot_grid(PCA_species, spec_FIG, nrow = 1)
ggsave(file = "figures/combined_specs.png", plot = combined_specs, width = 6.86, height = 4.8)


# Heatmaps ----

## ASVs ----
# Select top 50 ASVs
TOP50_ASVs <- row.names(model.cc$importance)[order(model.cc$importance[,"MeanDecreaseAccuracy"],
                                                     decreasing = TRUE)][1:50]
t50_ASVs <- t(data[,TOP50_ASVs])

# Define column annotations and colours
map_col <- data.frame(`Cover crop` = as.factor(results.df$cover_crop), check.names = FALSE)
map_col2 <- data.frame(Site = as.factor(results.df$site), check.names = FALSE)
col_fill = list(`Cover crop` = c(`Buckwheat` = "#c44601", `Buffalo Grass` = "#FCC9B5", `Crescendo Ladino Clover` = "#E1B239", 
                                 `Field Pea` = "#FCF2C7", `Mustard White` = "#A3D8C6", `Phacelia` = "#329973",
                                 `Spring Lentil` = "#7D99E6", `Turnip` = "#E0D2EB", `Winfred Brassica` =  "#98669F"))
col_fill2 <- list(Site = c(`Covert` = "red", `Kalala` = "blue", `SuRDC` = "yellow"))

# Heatmap annotations
annot_col <- ComplexHeatmap::HeatmapAnnotation(df = map_col, name = "Cover crop", col = col_fill, which = "column",
                                               annotation_name_gp = grid::gpar(fontsize = 5), simple_anno_size = unit(2, "mm"), show_legend = FALSE)
annot_col2 <- ComplexHeatmap::HeatmapAnnotation(df = map_col2, name = "Site", col = col_fill2, which = "column",
                                                annotation_name_gp = grid::gpar(fontsize = 5), simple_anno_size = unit(2, "mm"), show_legend = FALSE)

png(file = "figures/heatmap_50ASVs.png", width = 6.86, height = 4, units = "in", res = 600, bg = "transparent")
ComplexHeatmap::Heatmap(t50_ASVs, name = "Abundance", col = circlize::colorRamp2(c(0, 150, 300), c("#190087", "#E72476", "#FAEC50")),
                        show_row_names = FALSE,
                        row_names_gp = grid::gpar(fontsize = 4),
                        show_row_dend = TRUE,
                        row_dend_side = "left",
                        row_dend_width = unit(1, "cm"),
                        show_column_names = TRUE,
                        column_names_gp = grid::gpar(fontsize = 3),
                        show_column_dend = TRUE,
                        column_dend_side = "top",
                        column_dend_height = unit(1, "cm"),
                        # rect_gp = grid::gpar(col = "gray60", lwd = 0.5),
                        top_annotation = c(annot_col, annot_col2),
                        show_heatmap_legend = FALSE)
                        # heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 8)))
dev.off()


## Species ----
# Select top 50 species
TOP50_specs <- row.names(model.cc2$importance)[order(model.cc2$importance[,"MeanDecreaseAccuracy"],
                                                   decreasing = TRUE)][1:50]
t50_specs <- t(data2[,TOP50_specs])

# Fix duplicate words in row names
rows_spec <- as.data.frame(rownames(t50_specs))
rows_spec$Genus <- sapply(strsplit(rows_spec$`rownames(t50_specs)`, " "), `[`, 1)  # Separate genus and species
rows_spec$Species <- sapply(strsplit(rows_spec$`rownames(t50_specs)`, " "), `[`, 2)
# Only use species names if genus and species names contain duplicate words
for (i in 1:nrow(rows_spec)) {
  if (rows_spec$Genus[i] == rows_spec$Species[i]) {
    rows_spec$fixed[i] <- rows_spec$Species[i]
  } else if ("Fungi_gen_Incertae_sedis" %in% rows_spec$Genus[i]) {
    rows_spec$fixed[i] <- rows_spec$Species[i]
  } else if ("Chaetothyriales_gen_Incertae_sedis" %in% rows_spec$Genus[i]) {
    rows_spec$fixed[i] <- rows_spec$Species[i]
  } else if ("Mortierellales_gen_Incertae_sedis" %in% rows_spec$Genus[i]) {
    rows_spec$fixed[i] <- rows_spec$Species[i]
  } else if (rows_spec$Species[i] == paste("Unclassified_", rows_spec$Genus[i], sep = "")) {
    rows_spec$fixed[i] <- rows_spec$Species[i]
  } else if (rows_spec$Species[i] == paste("Unclassified_", rows_spec$Genus[i], ".1", sep = "")) {
    rows_spec$fixed[i] <- rows_spec$Species[i]
  } else {rows_spec$fixed[i] <- rows_spec$`rownames(t50_specs)`[i]
  }
}

# Reassign row names for top 50 species
rownames(t50_specs) <- rows_spec$fixed

# Define column annotations and colours
hmap_col <- data.frame(`Cover crop` = as.factor(sample_data(ps.species)$cover_crop), check.names = FALSE)
hmap_col2 <- data.frame(Site = as.factor(sample_data(ps.species)$site), check.names = FALSE)
hcol_fill = list(`Cover crop` = c(`Buckwheat` = "#c44601", `Buffalo Grass` = "#FCC9B5", `Crescendo Ladino Clover` = "#E1B239", 
                                 `Field Pea` = "#FCF2C7", `Mustard White` = "#A3D8C6", `Phacelia` = "#329973",
                                 `Spring Lentil` = "#7D99E6", `Turnip` = "#E0D2EB", `Winfred Brassica` =  "#98669F"))
hcol_fill2 <- list(Site = c(`Covert` = "red", `Kalala` = "blue", `SuRDC` = "yellow"))

# Heatmap annotations
hannot_col <- ComplexHeatmap::HeatmapAnnotation(df = hmap_col, name = "Cover crop", col = hcol_fill, which = "column",
                                               annotation_name_gp = grid::gpar(fontsize = 5), simple_anno_size = unit(2, "mm"), show_legend = FALSE)
hannot_col2 <- ComplexHeatmap::HeatmapAnnotation(df = hmap_col2, name = "Site", col = hcol_fill2, which = "column",
                                                annotation_name_gp = grid::gpar(fontsize = 5), simple_anno_size = unit(2, "mm"), show_legend = FALSE)

# png(file = "figures/heatmap_50specs.png", width = 6.86, height = 4, units = "in", res = 600, bg = "transparent")
hm_50specs <- grid.grabExpr(ComplexHeatmap::draw(
  ComplexHeatmap::Heatmap(t50_specs, name = "Abundance", col = circlize::colorRamp2(c(0, 150, 300), c("#190087", "#E72476", "#FAEC50")),
                        show_row_names = TRUE,
                        row_names_gp = grid::gpar(fontsize = 4),
                        show_row_dend = TRUE,
                        row_dend_side = "left",
                        row_dend_width = unit(0.75, "cm"),
                        show_column_names = TRUE,
                        column_names_gp = grid::gpar(fontsize = 3),
                        show_column_dend = TRUE,
                        column_dend_side = "top",
                        column_dend_height = unit(0.75, "cm"),
                        # rect_gp = grid::gpar(col = "gray60", lwd = 0.5),
                        top_annotation = c(hannot_col, hannot_col2),
                        show_heatmap_legend = FALSE,),
                        background = "transparent"))
                        # heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 5)))
# dev.off()


# Heatmap RF legend ----

a.scale <- data.frame("fill" = t50_specs[,"S38"], "species" = row.names(t50_specs))
hm.legend1 <-
  ggplot(data = a.scale) +
  theme_bw() +
  geom_tile(aes(x = a.scale$species, y = a.scale$fill, fill = a.scale$fill)) +
  labs(fill = "Abundance") +
  scale_fill_gradient2(low = "#190087", mid = "#E72476", high = "#FAEC50", midpoint = 150, space = "Lab") +
  theme(legend.title = element_text(size = 4, family = "sans", face = "bold"),
        legend.text = element_text(size = 3.5, family = "sans", face = "bold"),
        legend.position = c(0.9025,0.822), # horizontal, vertical
        legend.direction = "vertical",
        # legend.margin = margin(t = 0, r = 0, b = 0, l =0, unit = "pt"),
        legend.box.spacing = unit(c(0,0,0,0), "cm"),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.18, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"))

lgd1 <- get_legend(hm.legend1)
grid.newpage()
grid.draw(lgd1)

hm.legend2 <-
  ggplot(data = results.df) +
  geom_tile(aes(x = sample.ID, y = cover_crop, fill = cover_crop)) +
  scale_fill_manual(values = palette) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 3.5, family = "sans", face = "bold"),
        legend.position = c(0.94,0.71),#c(1, 0.45), # horizontal, vertical
        legend.direction = "vertical",
        # legend.margin = margin(t = 0, r = 0, b = 0, l =0, unit = "pt"),
        legend.box.spacing = unit(c(0,0,0,0), "cm"),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.18, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"))

lgd2 <- get_legend(hm.legend2)
grid.newpage()
grid.draw(lgd2)

hm.legend3 <-
  ggplot(data = results.df) +
  geom_tile(aes(x = sample.ID, y = site, fill = site)) +
  scale_fill_manual(values = c("red", "blue", "gold")) +
  labs(fill = "Site") +
  theme(legend.title = element_text(size = 4, family = "sans", face = "bold"),
        legend.text = element_text(size = 3.5, family = "sans", face = "bold"),
        legend.position = c(0.96,0.83),#c(1, 0.45), # horizontal, vertical
        legend.direction = "vertical",
        # legend.margin = margin(t = 0, r = 0, b = 0, l =0, unit = "pt"),
        legend.box.spacing = unit(c(0,0,0,0), "cm"),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.18, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"))

lgd3 <- get_legend(hm.legend3)
grid.newpage()
grid.draw(lgd3)

# Inset legends into heatmap
LGDS1 <-
  ggdraw(lgd1) +
  draw_plot(lgd2,
            x = 0,
            y = 0,
            width = 1,
            height = 1) +
  draw_plot(lgd3,
            x = 0,
            y = 0,
            width = 1,
            height = 1)
# ggsave("figures/hm_legend_test.png", plot = LGDS1, width = 6.86, height = 4.5)

hm_legend <-
  ggdraw(LGDS1) +
  draw_plot(hm_50specs,
            x = 0,
            y = 0,
            width = 0.94,
            height = 1) +
  draw_label("D", x = 0.025, y = 0.94,
             fontfamily = "sans",
             fontface = "bold",
             size = 12)


ggsave("figures/heatmap_fixlgd_cc.png", plot = hm_legend, width = 6.86, height = 4.25)

# Combine all random forest analysis figures
RF_species <- plot_grid(combined_specs, hm_legend, nrow = 2)
ggsave(file = "figures/RF_species.png", plot = RF_species, width = 6.86, height = 8.65)




# Figure 1 ----

# Species richness boxplot
plot2 <-
  ggplot(data = results.df) +
  geom_boxplot(aes(y = rich, fill = cover_crop), size = 0.2, outlier.size = 0.2) +  # y = c(PD, rich, Shannon, Simpson)
  labs(fill = NULL, y = "Species Richness", tag = "B") +
  scale_fill_manual(values = palette) +
  facet_wrap(vars(site)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        plot.tag.position = c(0.025, 0.97), # horizontal, vertical
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
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_text(size = 8,family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left

# Group boxplot and PCA plot
Fig1_BC <- plot_grid(plot2, PCA_FIG.site2, ncol = 1)

# Load package for venn diagram
library(ggvenn)

# Agglomerate by species for overall OTU table
ps.species <- tax_glom(merge.ps, taxrank = 'Species', NArm = FALSE)
# Merge samples by site
sitegroup <- merge_samples(otu_table(ps.species), group = sample_data(ps.species)$site, fun = sum)
sitegroup <- t(as.data.frame(otu_table(sitegroup)))
# Replace sequences with species names
venndata <- as.data.frame(sitegroup)
taxaps <- as.data.frame(tax_table(ps.species))
rownames(venndata) <- paste(taxaps$Genus, taxaps$Species, sep = ".")

for (i in 1:nrow(venndata)) {
  if (venndata[i,1] > 0) {  # TRUE if present
    venndata[i,1] <- row.names(venndata[i,])
  } else if (venndata[i,1] == 0) {  # FALSE if absent
    venndata[i,1] <- NA} 
  if (venndata[i,2] > 0) {  # TRUE if present
    venndata[i,2] <- row.names(venndata[i,])
  } else if (venndata[i,2] == 0) {  # FALSE if absent
    venndata[i,2] <- NA}
  if (venndata[i,3] > 0) {  # TRUE if present
    venndata[i,3] <- row.names(venndata[i,])
  } else if (venndata[i,3] == 0) {  # FALSE if absent
    venndata[i,3] <- NA}
}
venndata <- as.list(venndata)

# Venn diagram of species grouped by sites
venn <-
  ggvenn(venndata, columns = c("Covert", "Kalala", "SuRDC"),
         fill_color = c("red", "blue", "gold"),
         stroke_size = 0.5, set_name_size = 4, text_size = 2.4) +
  labs(tag = "A") +
  theme(plot.background = element_rect(colour = "transparent"),
        plot.tag = element_text(size = 12, family = "sans", face = "bold"),
        plot.tag.position = c(0.08, 0.97)) # horizontal, vertical

# Add Venn diagram to boxplot and PCA
Fig1 <- plot_grid(venn, Fig1_BC, ncol = 2)

ggsave(file = "figures/Figure_1_CC.png", plot = Fig1, width = 6.86, height = 4.7)


# Figure 2 ----

# Combine density plots with PCA
Fig4_BC <- plot_grid(spec_FIG1a, spec_FIG1b, ncol = 1)
Fig4_ABC <- plot_grid(PCA_FIG.cc2, Fig4_BC, ncol = 2)

# Add heatmap to density and PCA
Fig4 <- plot_grid(Fig4_ABC, hm_legend, nrow = 2, rel_heights = c(2,3))

ggsave(file = "figures/Figure_4_CC.png", plot = Fig4, width = 6.86, height = 7)
