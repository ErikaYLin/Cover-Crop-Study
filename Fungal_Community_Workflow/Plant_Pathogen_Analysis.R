# Pathogen Presence, Abundance, and Diversity

# Load in the data
load(file = "RDS/results_covercrop_noGreenhouse.rda")
results.df <- readRDS(file = "RDS/results.df_noGreenhouse.rds")
plant_pathogens <- read.csv(file = "data/Plant_Pathogens.csv")

# Load packages
library(phyloseq)
library(rphylopic)
library(dplyr)
library(ggplot2)
library(cowplot)

palette2 <- c("#9B2226", "#CA6702", "#fbb13c", "#E9D8A6", "#7A9EC6", "#94D2BD", "#0A9396")

# Extra colours: , "#005F73", "#EE9B00", "#355070", "#9F72AA","#6d597a", 

# Restructuring data for easier subsequent handling ----

# Extract all phyloseq objects from RESULTS into a list
ps.list <- list()
for (i in 1:length(RESULTS)) {
  ps.list[[i]] <- RESULTS[[i]]$ps
}

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

# Agglomerate by species for overall OTU table
ps.species <- tax_glom(merge.ps, taxrank = 'Species', NArm = FALSE)


# Identify vine pathogens ----
# View taxonomy table to identify vine pathogens
pathogen <- tax_table(ps.species)
pathogen <- as.data.frame(pathogen)
# Add genus to species name in Species column
pathogen$Species <- paste(pathogen$Genus, pathogen$Species, sep = " ")

# Add column to indicate possible pathogenicity (yes = 1, no = 0)
pathogen$pathogenic = 0
pathogen$pathogenic_taxon = NA

# Search for pathogens in taxonomy table
for (j in 1:nrow(plant_pathogens)) {
  # Search in Family
  if ("Family" %in% plant_pathogens$Taxon[j]){
    for (i in 1:nrow(pathogen)) {
      if (plant_pathogens$Pathogen[j] %in% pathogen$Family[i]){
        pathogen$pathogenic[i] = 1
        pathogen$pathogenic_taxon[i] = plant_pathogens$Pathogen[j]
      }}
  # Search in Genus
  } else if ("Genus" %in% plant_pathogens$Taxon[j]){
    for (i in 1:nrow(pathogen)) {
      if (plant_pathogens$Pathogen[j] %in% pathogen$Genus[i]){
        pathogen$pathogenic[i] = 1
        if (is.na(pathogen$pathogenic_taxon[i])){
        pathogen$pathogenic_taxon[i] = plant_pathogens$Pathogen[j]
        }}}
  # search in Species
  } else if ("Species" %in% plant_pathogens$Taxon[j]){
    for (i in 1:nrow(pathogen)) {
      if (plant_pathogens$Pathogen[j] %in% pathogen$Species[i]){
        pathogen$pathogenic[i] = 1
        if (is.na(pathogen$pathogenic_taxon[i])){
          pathogen$pathogenic_taxon[i] = plant_pathogens$Pathogen[j]
        }}}
  } else {
    pathogen$pathogenic[i] = 0
    pathogen$pathogenic_taxon[i] = NA}
}

# Export to .csv for others to easily update the table
# write.csv(pathogen, file = "Fungal_Community_Workflow/Plant_Pathogen_Classification.csv", row.names = TRUE)

# # Import updated table
# pathogen <- read.csv(file = "Fungal_Community_Workflow/Plant_Pathogen_Classification.csv")
# row.names(pathogen) <- pathogen$X
# pathogen <- pathogen[,-1]

tax = phyloseq::tax_table(as.matrix(pathogen))
tax_table(ps.species) <- tax

# Extract pathogens only
pathogen.only <- pathogen[pathogen$pathogenic == 1,]
path.names <- row.names(pathogen.only)

# Extract pathogens from ASV table
pathogen.ps <- prune_taxa(path.names, ps.species)


# Relative abundance bar plot ----
# Transform sample counts to relative abundances
pathogen.rel <- transform_sample_counts(pathogen.ps, function(x) x/sum(x)*100)
pathogen.melt <- psmelt(pathogen.rel)
pathogen.melt$pathogenic_taxon <- pathogen.only$pathogenic_taxon

# Organize the data for abundance plots
pathogen.melt1 <- pathogen.melt %>%
  group_by(sample.ID, cover_crop, pathogenic_taxon) %>%
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

# Relative abundance bar plot
FIG <-
  ggplot(pathogen.melt1, aes(x = sample.ID, y = Abundance, fill = pathogenic_taxon)) + 
  geom_bar(stat = "identity", aes(fill = pathogenic_taxon)) + 
  labs(x = "", y= "Abundance (%)", fill = "Taxon", tag = "A") +
  facet_wrap(~cover_crop, scales = "free_x", nrow = 1,
             labeller = labeller(cover_crop = label_wrap_gen(width = 10))) +
  theme_classic() +
  scale_fill_manual(values = palette2) + 
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(plot.tag = element_text(size = 14, family = "sans", face = "bold"),
        plot.tag.position = c(0.02, 0.95),
        axis.title.y = element_text(size = 10, family = "sans", face = "bold"), # vjust = 2.5),
        axis.text.y = element_text(size = 7, family = "sans"),
        axis.text.x  = element_text(size = 2, angle = 45),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        legend.text = element_text(size = 7, family = "sans", face = "bold"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.key.height = unit(0.32, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        strip.text.x = element_text(size = 6, family = "sans", face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("figures/pathogens_rel_barplot_CC.png", plot = FIG, width = 6.86, height = 4.5)

icons <-
  ggplot() +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) +
  add_phylopic(buckwheat, alpha = 1, x = 9.2, y = 88.5, ysize = 8, color = "#c44601") +
  add_phylopic(buffalo, alpha = 1, x = 16.5, y = 88.5, ysize = 8, color = "#FCC9B5") +
  add_phylopic(clover, alpha = 1, x = 25.5, y = 88, ysize = 7, color = "#E1B239") +
  add_phylopic(fieldpea, alpha = 1, x = 33.6, y = 88.5, ysize = 8, color = "#FCF2C7") +
  add_phylopic(mustard, alpha = 1, x = 41.7, y = 88.5, ysize = 7.8, color = "#A3D8C6") +
  add_phylopic(phacelia, alpha = 1, x = 50.2, y = 88.5, ysize = 7.8, color = "#329973") +
  add_phylopic(lentil, alpha = 1, x = 57.8, y = 88.5, ysize = 8, color = "#7D99E6") +
  add_phylopic(turnip, alpha = 1, x = 65.5, y = 88.5, ysize = 8, color = "#E0D2EB") +
  add_phylopic(brassica, alpha = 1, x = 74, y = 88.5, ysize = 6.5, color = "#98669F") +
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
        plot.background = element_blank())

# Layer cover crop icons on top of abundance bar plot
FIG2 <-
  ggdraw(FIG) +
  draw_plot(icons,
            x = 0,
            y = 0,
            width = 1,
            height = 1)
ggsave("figures/pathogens_relative_CC.png", plot = FIG2, width = 6.86, height = 4.5)

# Relative abundance of pathogen to non-pathogenic taxa
# Transform sample counts to relative abundances
pathogen.rel2 <- transform_sample_counts(ps.species, function(x) x/sum(x)*100)
pathogen.melt2 <- psmelt(pathogen.rel2)

# Organize the data for abundance plots
nonpathogenic <- unique(pathogen.melt2$Species[pathogen.melt2$pathogenic == 0])
pathogen.melt2$pathogenic[!(pathogen.melt2$Species %in% nonpathogenic)] <- "Pathogenic"
pathogen.melt2$pathogenic[pathogen.melt2$Species %in% nonpathogenic] <- "Non-pathogenic"
pathogen.melt2 <- pathogen.melt2 %>%
  group_by(sample.ID, cover_crop, pathogenic) %>%
  summarise(Abundance = sum(Abundance))

# Relative abundance bar plot
FIG6 <-
  ggplot(pathogen.melt2, aes(x = sample.ID, y = Abundance, fill = pathogenic)) + 
  geom_bar(stat = "identity", aes(fill = pathogenic)) + 
  labs(x = "", y= "Abundance (%)", tag = "B") +
  facet_wrap(~cover_crop, scales = "free_x", nrow = 1,
             labeller = labeller(cover_crop = label_wrap_gen(width = 10))) +
  theme_classic() +
  scale_fill_manual(values = c("#9EDAFA", "#FD8A8B")) + 
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(plot.tag = element_text(size = 14, family = "sans", face = "bold"),
    plot.tag.position = c(0.02, 0.95),
    axis.title.y = element_text(size = 10, family = "sans", face = "bold"), # vjust = 2.5),
    axis.text.y = element_text(size = 7, family = "sans"),
    axis.text.x  = element_text(size = 2, angle = 45),
    axis.line.y = element_line(linewidth = 0.2),
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks.x = element_line(linewidth = 0.2),
    legend.text = element_text(size = 7, family = "sans", face = "bold"),
    legend.title = element_blank(),
    legend.key.height = unit(0.32, "cm"),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    strip.text.x = element_text(size = 6, family = "sans", face = "bold"),
    strip.background = element_blank(),
    panel.spacing = unit(0.15, "lines"),
    plot.background = element_blank(),
    plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

icons2 <-
  ggplot() +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) +
  add_phylopic(buckwheat, alpha = 1, x = 9.2, y = 87.5, ysize = 8, color = "#c44601") +
  add_phylopic(buffalo, alpha = 1, x = 17, y = 87.5, ysize = 8, color = "#FCC9B5") +
  add_phylopic(clover, alpha = 1, x = 26, y = 87, ysize = 7, color = "#E1B239") +
  add_phylopic(fieldpea, alpha = 1, x = 34.6, y = 87.5, ysize = 8, color = "#FCF2C7") +
  add_phylopic(mustard, alpha = 1, x = 42.7, y = 87.5, ysize = 7.8, color = "#A3D8C6") +
  add_phylopic(phacelia, alpha = 1, x = 51.3, y = 87.5, ysize = 7.8, color = "#329973") +
  add_phylopic(lentil, alpha = 1, x = 59, y = 87.5, ysize = 8, color = "#7D99E6") +
  add_phylopic(turnip, alpha = 1, x = 67.4, y = 87.5, ysize = 8, color = "#E0D2EB") +
  add_phylopic(brassica, alpha = 1, x = 76, y = 87.5, ysize = 6.5, color = "#98669F") +
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

# Layer cover crop icons on top of abundance bar plot
FIG7 <-
  ggdraw(FIG6) +
  draw_plot(icons2,
            x = 0,
            y = 0,
            width = 1,
            height = 1)

ggsave("figures/pathogen-nonpathogen_CC.png", plot = FIG7, width = 6.86, height = 4.5)

# Arrange figures
FIG8 <- plot_grid(FIG2, FIG7, ncol = 1)
ggsave("figures/pathogens_barplots.png", plot = FIG8, width = 6.86, height = 7.8)

# Relative abundance grouped by site
# Represent site with a number (characters lost in merge_samples())
ps.pathogen <- pathogen.ps
sample_data(ps.pathogen)$treat <- paste(sample_data(ps.pathogen)$site, sample_data(ps.pathogen)$cover_crop, sep = ".")
# sample_data(ps.pathogen)$treatment <- gsub("Covert", 0, sample_data(ps.pathogen)$site)
# sample_data(ps.pathogen)$treatment <- gsub("Kalala", 1, sample_data(ps.pathogen)$site) 
# sample_data(ps.pathogen)$treatment <- gsub("SuRDC", 100, sample_data(ps.pathogen)$site)

# Organize the data for abundance plots
pathogen.CC <- merge_samples(ps.pathogen, "treat", fun = sum)
pathogen.rel3 <- transform_sample_counts(pathogen.CC, function(x) x/sum(x)*100)
pathogen.melt3 <- psmelt(pathogen.rel3)

# Add cover crop and site names back to data frame
for (n in 1:nrow(pathogen.melt3)) {
  if(grepl("Covert", pathogen.melt3$Sample[n])){
    pathogen.melt3$site[n] = "Covert"
  } else if (grepl("Kalala", pathogen.melt3$Sample[n])){
    pathogen.melt3$site[n] = "Kalala"
  } else if (grepl("SuRDC", pathogen.melt3$Sample[n])){
    pathogen.melt3$site[n] = "SuRDC"
  } 
  if (grepl("Buck", pathogen.melt3$Sample[n])){
    pathogen.melt3$cover_crop[n] = "Buckwheat"
  } else if (grepl("Buff", pathogen.melt3$Sample[n])){
    pathogen.melt3$cover_crop[n] = "Buffalo Grass"
  } else if (grepl("Cres", pathogen.melt3$Sample[n])){
    pathogen.melt3$cover_crop[n] = "Crescendo Ladino Clover"
  } else if (grepl("Pea", pathogen.melt3$Sample[n])){
    pathogen.melt3$cover_crop[n] = "Field Pea"
  } else if (grepl("Mustard", pathogen.melt3$Sample[n])){
    pathogen.melt3$cover_crop[n] = "Mustard White"
  } else if (grepl("Phac", pathogen.melt3$Sample[n])){
    pathogen.melt3$cover_crop[n] = "Phacelia"
  } else if (grepl("Spring", pathogen.melt3$Sample[n])){
    pathogen.melt3$cover_crop[n] = "Spring Lentil"
  } else if (grepl("Turnip", pathogen.melt3$Sample[n])){
    pathogen.melt3$cover_crop[n] = "Turnip"
  } else if (grepl("Brass", pathogen.melt3$Sample[n])){
    pathogen.melt3$cover_crop[n] = "Winfred Brassica"
}}

pathogen.melt3 <- pathogen.melt3 %>%
  group_by(cover_crop, site, pathogenic_taxon) %>%
  summarise(Abundance = sum(Abundance))

# Relative abundance bar plot
FIG9 <-
  ggplot(pathogen.melt3, aes(x = cover_crop, y = Abundance, fill = pathogenic_taxon)) + 
  geom_bar(stat = "identity", aes(fill = pathogenic_taxon)) + 
  labs(x = "", y= "Abundance (%)", fill = "Taxon") +
  facet_wrap(~site) + #, scales = "free_x", nrow = 1,
             #labeller = labeller(cover_crop = label_wrap_gen(width = 10))) +
  theme_classic() +
  scale_fill_manual(values = palette2) + 
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme(#plot.tag = element_text(size = 14, family = "sans", face = "bold"),
        #plot.tag.position = c(0.02, 0.95),
        axis.title = element_text(size = 10, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 7, family = "sans"),
        axis.text.x  = element_text(size = 4, angle = 90, vjust = 0.5),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        legend.text = element_text(size = 7, family = "sans", face = "bold"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.key.height = unit(0.32, "cm"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.3, fill = NA),
        strip.text.x = element_text(size = 10, family = "sans", face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("figures/pathogens_relative_site.png", plot = FIG9, width = 6.86, height = 4.5)


# Absolute abundance bar plot ----
# Combine abundance data and taxa
pathogen.melt3 <- psmelt(pathogen.ps)
pathogen.melt3$pathogenic_taxon <- pathogen.only$pathogenic_taxon

# Organize the data for abundance plots
pathogen.melt3 <- pathogen.melt3 %>%
  group_by(sample.ID, site, cover_crop, pathogenic_taxon) %>%
  summarise(Abundance)
  # summarise(Abundance)
pathogen.melt3 <- aggregate(pathogen.melt3, pathogen.melt3$Abundance, FUN = sum)

# Absolute abundance bar plot
FIG3 <-
  ggplot(pathogen.melt3, aes(x = sample.ID, y = Abundance, fill = pathogenic_taxon)) + 
  geom_bar(stat = "identity", aes(fill = pathogenic_taxon)) + 
  labs(x = "", y= "Abundance", fill = "Taxon") +
  facet_wrap(~cover_crop, scales = "free_x", nrow = 1,
             labeller = labeller(cover_crop = label_wrap_gen(width = 10))) +
  theme_classic() +
  scale_fill_manual(values = palette2) + 
  scale_y_continuous(limits = c(0,15000), expand = c(0,0)) +
  theme(# plot.tag = element_text(size = 14, family = "sans", face = "bold"),
    # plot.tag.position = c(0.83, 0.85),
    axis.title.y = element_text(size = 10, family = "sans", face = "bold"),
    axis.text.y = element_text(size = 7, family = "sans"),
    axis.text.x  = element_text(size = 2, angle = 45),
    axis.line.y = element_line(linewidth = 0.2),
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks.x = element_line(linewidth = 0.2),
    legend.text = element_text(size = 7, family = "sans", face = "bold"),
    legend.title = element_text(size = 8, family = "sans", face = "bold"),
    legend.key.height = unit(0.32, "cm"),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    strip.text.x = element_text(size = 6, family = "sans", face = "bold"),
    strip.background = element_blank(),
    panel.spacing = unit(0.15, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("figures/pathogens_abs_barplot_CC.png", plot = FIG3, width = 6.86, height = 4.5)

icons3 <-
  ggplot() +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) +
  add_phylopic(buckwheat, alpha = 1, x = 10.7, y = 86.5, ysize = 8, color = "#c44601") +
  add_phylopic(buffalo, alpha = 1, x = 18.2, y = 86.5, ysize = 8, color = "#FCC9B5") +
  add_phylopic(clover, alpha = 1, x = 26.6, y = 86, ysize = 7, color = "#E1B239") +
  add_phylopic(fieldpea, alpha = 1, x = 34.9, y = 86.5, ysize = 8, color = "#FCF2C7") +
  add_phylopic(mustard, alpha = 1, x = 42.4, y = 86.5, ysize = 7.8, color = "#A3D8C6") +
  add_phylopic(phacelia, alpha = 1, x = 50.4, y = 86.5, ysize = 7.8, color = "#329973") +
  add_phylopic(lentil, alpha = 1, x = 58, y = 86.5, ysize = 8, color = "#7D99E6") +
  add_phylopic(turnip, alpha = 1, x = 66, y = 86.5, ysize = 8, color = "#E0D2EB") +
  add_phylopic(brassica, alpha = 1, x = 74.1, y = 86.5, ysize = 6.5, color = "#98669F") +
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


# Layer cover crop icons on top of abundance bar plot
FIG4 <-
  ggdraw(FIG3) +
  draw_plot(icons3,
            x = 0,
            y = 0,
            width = 1,
            height = 1)
ggsave("figures/pathogens_absolute_CC.png", plot = FIG4, width = 6.86, height = 4.5)

# Abundance bar plot grouped by site and cover crop

pathogen.melt4 <- psmelt(pathogen.ps)
pathogen.melt4$pathogenic_taxon <- pathogen.only$pathogenic_taxon
pathogen.melt4 <- pathogen.melt4 %>%
  group_by(cover_crop, site, pathogenic_taxon, Abundance) %>%
  summarise(Abundance = mean(Abundance))

# Absolute abundance bar plot
FIG5 <-
  ggplot(pathogen.melt4, aes(x = cover_crop, y = Abundance, fill = pathogenic_taxon)) + 
  geom_bar(stat = "identity", aes(fill = pathogenic_taxon)) + 
  labs(x = "", y= "Abundance", fill = "Taxon") +
  facet_wrap(~site) +
  theme_classic() +
  scale_fill_manual(values = palette2) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  theme(#plot.tag = element_text(size = 14, family = "sans", face = "bold"),
    #plot.tag.position = c(0.02, 0.95),
    axis.title = element_text(size = 10, family = "sans", face = "bold"),
    axis.text.y = element_text(size = 7, family = "sans"),
    axis.text.x  = element_text(size = 4, angle = 90, vjust = 0.5),
    axis.line.y = element_line(linewidth = 0.2),
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks.x = element_line(linewidth = 0.2),
    legend.text = element_text(size = 7, family = "sans", face = "bold"),
    legend.title = element_text(size = 8, family = "sans", face = "bold"),
    legend.key.height = unit(0.32, "cm"),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(linewidth = 0.3, fill = NA),
    strip.text.x = element_text(size = 10, family = "sans", face = "bold"),
    strip.background = element_blank(),
    panel.spacing = unit(0.15, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) #top, right, bottom, left

ggsave("figures/pathogens_absolute_site.png", plot = FIG5, width = 6.86, height = 4.5)


