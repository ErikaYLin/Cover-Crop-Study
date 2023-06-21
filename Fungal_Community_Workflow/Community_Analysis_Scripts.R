# set working directory
setwd("C:/Users/elinrg26/OneDrive - UBC/Research/GitHub/Cover Crop Study")

# Load in the data
test <- readRDS(file = "RDS/metadata_covercrop.rds")
load(file = "RDS/preliminary_results_covercrop.rda")

for (i in 1:length(RESULTS)){
  
  RESULTS[[i]]$results.samples <- RESULTS[[i]]$results.samples[1,1:10]
  RESULTS[[i]]$results.samples <- cbind(RESULTS[[i]]$results.samples, as.data.frame(test[[i]]))
}


RES <- list()
# fix results.samples
for (i in 1:length(RESULTS)){
  
  RES[[i]] <- RESULTS[[i]]$results.samples[1,1:10]
  RES[[i]] <- cbind(RES[[i]], as.data.frame(test[[i]]))
}

# View data frame of metrics:
results.df <- do.call(rbind, RES) %>%
  dplyr::select("sample.ID", everything())

# Load ggplot2
library(ggplot2)

# Big picture descriptive analyses ----

boxplot(PD ~ cover_crop, data = results.df)

# boxplot of diversity measures across cover crops and sites
png(file = "figures/PD_boxplot.png", height = 1000, width = 1000)
ggplot(data = results.df) +
  geom_boxplot(aes(y = PD, fill = cover_crop)) +  # y = c(PD, rich, Shannon, Simpson)
  scale_y_log10() +
  facet_wrap(vars(site))
dev.off()

png(file = "figures/richness_boxplot.png", height = 1000, width = 1000)
ggplot(data = results.df, aes(y = rich)) +
  geom_boxplot(aes(fill = cover_crop)) +
  facet_wrap(vars(site))
dev.off()


png(file = "figures/connectivity_vs_connectance.png", height = 1000, width = 1000)
ggplot(data = results.df, aes(x = connectivity, y = connectance)) +
  geom_point(aes(colour = cover_crop), size = 4) +
  geom_smooth(method = "gam", method.args= list(family = "Gamma"), se = F) +
  facet_wrap(vars(site))
dev.off()

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

# Heatmaps ----
