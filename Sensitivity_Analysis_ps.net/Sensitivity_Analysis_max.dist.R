# set working directory
setwd("C:/Users/elinrg26/OneDrive - UBC/Research/GitHub/Cover Crop Study")

# Sensitivity Analysis on the `inv` argument in `ps.net` ----

# Load sample metadata
samdf <- read.csv("C:/Users/elinrg26/OneDrive - UBC/Research/GitHub/Cover Crop Study/Fungal_Community_Workflow/Mapping file for ITS sequencing.csv")
# Load reference FASTA file for taxonomic classification
fungi <- "./sh_general_release_dynamic_29.11.2022.fasta" # file name for fungal reference FASTA

# Define filter path for the forward and reverse reads of each sample fastq file
seq_path <- "./CC_Seq"  # Replace with the directory for fastq file after unzipping the folder

# Sort forward and reverse reads to be in the same order
fnFs <- sort(list.files(seq_path, pattern="_R1_001.fastq.gz"))
fnRs <- sort(list.files(seq_path, pattern="_R2_001.fastq.gz"))
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(seq_path, fnFs)
fnRs <- file.path(seq_path, fnRs)

# Extract sample names, assuming file paths have format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq
sNames <- sapply(strsplit(fnFs, "_"), `[`, 2)
sNames <- sapply(strsplit(sNames, "/"), `[`, 2)

# Define file names for filtered fastq.gz files and assign file paths
filt_path <- file.path(seq_path, "filtered") # Place filtered files in new subdirectory if one does not already exist
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sNames, "_R_filt.fastq.gz"))


filter.args = list(truncLen = c(240,160),
                   maxN = 0,
                   maxEE = c(2,2),
                   truncQ = 2,
                   rm.phix = TRUE, 
                   compress = TRUE,
                   multithread = FALSE)

out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                            truncLen = filter.args$truncLen,
                            maxN = filter.args$maxN,
                            maxEE = filter.args$maxEE,
                            truncQ = filter.args$truncQ,
                            rm.phix = filter.args$rm.phix,
                            compress = filter.args$compress,
                            multithread = filter.args$multithread)

source(file = "Microbiome_Analysis_Pipeline.R")
# `ps.only` function used to separate fastq processing from network construction to reduce redundancy which decreases computation time
  
# Empty list for storing the results
RESULTS <- list()

# Loop the analysis process for every file in the directory
for (i in 1:length(fnFs)){
  
  # Assign metadata to be input as a list for `ps.net` argument
  meta_data <- list()
  for (j in 1:length(fnFs)){
    
    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(samdf[samdf$sample.ID %in% ID,])
  }
  
  # Create a new row of results for each output from the `ps.net`
  RESULTS[[i]] <- ps.only(fastq.Fs = fnFs[i],
                         fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                         filtFs = filtFs[i],
                         filtRs = filtRs[i], # file paths for filtered fastq files
                         refFasta = fungi,
                         metadata = meta_data[[j]],
                         tree.args = list(k = 4, 
                                          inv = 0.2,
                                          model = "GTR",
                                          rearrangement = "stochastic"))
  
  # Save the results as they are produced, in case of crashes/errors
  save(RESULTS, file = "RDS/Sensitivity_Analysis/results_CC_ps.only.rda")
}

# Newly improved results.samples
RES <- list()
for (i in 1:length(RESULTS)){
  
  RES[[i]] <- RESULTS[[i]]$results.samples
}

# Build combined data frame of metrics
results.df <- do.call(rbind, RES)
results.df <- dplyr::relocate(results.df, "sample.ID", .before = "Shannon") # move sample.ID column to leftmost
# Save data frame of results
saveRDS(results.df, file = "RDS/Sensitivity_Analysis/results.df_ps.only.rds")


# List of max.dist values to be tested
max.dist_values <- list(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.2, 0.1, 0)

# Empty list for storing each iteration of the sensitivity analysis
SA_max.dist <- list()

# Loop for every value of max.dist
for (m in 1:length(max.dist_values)) {
  
  # Empty list for storing network results
  NETWORKS <- list()
  
  # Loop the analysis process for every file in the directory
  for (i in 1:length(RESULTS)){
  
  # Construct the network using `ps` from RESULTS
  net <- phyloseq::make_network(RESULTS[[i]]$ps,
                                type = "taxa",
                                distance = "jaccard",
                                max.dist = max.dist_values[[m]],
                                keep.isolates = TRUE)
  sampledata <- as.data.frame(phyloseq::sample_data(RESULTS[[i]]$ps)) # subset sample_info as a data frame
  
  # Count vertices and edges in site network
  v <- igraph::vcount(net)
  e <- igraph::ecount(net)
  
  connectivity <- e/v   # avg number associations between ASVs
  connectance <- e/v^2  # fraction of edges out of all possible edges
  
  # Data frame of sample network results
  results.net <- data.frame(vertices = v, 
                            edges = e, 
                            connectivity, 
                            connectance)
  
  results.net <- cbind(results.net, sampledata)
}  

  # Retrieve network results
  NET <- list()
  for (n in 1:length(NETWORKS)){
    
    NET[[n]] <- NETWORKS[[n]]$results.net
  }
  
  # Build combined data frame of metrics
  res.net[[m]] <- do.call(rbind, NET)
  results.df[[m]] <- cbind(results.df, res.net[[m]])
  
  # GLMM plot for Connectivity
  null_model[[m]] <- gam(list(connectivity ~ 1,  # Average connectivity
                          ~ 1),
                     data = results.df[[m]],
                     family = gammals())  
  
  model[[m]] <- gam(list(connectivity ~ cover_crop,  # Connectivity with the added effect of cover crop species
                     ~ cover_crop),
                data = results.df[[m]],
                family = gammals()) 
  
  # Calculate delta AIC for each iteration
  d_AIC[[m]] <- AIC(model)-AIC(null_model)
  
  
# Store function outputs for each iteration in a list
SA_max.dist[[m]] <- c(NETWORKS[[m]], results.df[[m]], model[[m]], d_AIC[[m]])
save(SA_max.dist, file = "RDS/Sensitivity_Analysis/Sensitivity_Results_max.dist.rda")
}

