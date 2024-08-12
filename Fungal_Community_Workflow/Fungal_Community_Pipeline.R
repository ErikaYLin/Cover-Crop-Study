# NEW CODE ----

# Install and load `mga` from GitHub
devtools::install_github("https://github.com/ErikaYLin/mga")
library(mga)

# Load sample metadata
samdf <- read.csv("data/Mapping file for ITS sequencing.csv")

# Load reference FASTA file for taxonomic classification
fungi <- "./sh_general_release_dynamic_29.11.2022.fasta" # UNITE General Release ITS reference FASTA

# Define filter path for forward and reverse reads of each sample fastq file
seq_path <- "./CC_Seq"  # directory containing extracted fastq files

# Sort forward and reverse reads to be in the same order
fnFs <- sort(list.files(seq_path, pattern="_R1_001.fastq.gz"))
fnRs <- sort(list.files(seq_path, pattern="_R2_001.fastq.gz"))
# Specify full file path to fnFs and fnRs
fnFs <- file.path(seq_path, fnFs)
fnRs <- file.path(seq_path, fnRs)

# Extract sample names
sNames <- sapply(strsplit(fnFs, "_"), `[`, 2)
sNames <- sapply(strsplit(sNames, "/"), `[`, 2)

# Define file names for filtered fastq.gz files
filt_path <- file.path(seq_path, "filtered")
# New subdirectory for filtered files if one does not already exist
if(!file_test("-d", filt_path)) dir.create(filt_path)
# Assign file paths
filtFs <- file.path(filt_path, paste0(sNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sNames, "_R_filt.fastq.gz"))

# Inspect the first 2 forward reads:
dada2::plotQualityProfile(fnFs[1:2])  ## forward reads maintain high throughput quality, trimmed at position 230
# Inspect the first 2 reverse reads:
dada2::plotQualityProfile(fnRs[1:2])  ## reverse read quality drops around position 180-200, trimmed at position 180
## first 10 bps also removed due to potential pathological issues

out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                            truncLen = c(230,180),
                            trimLeft = 10,
                            maxN = 0,
                            maxEE = c(2,2),
                            truncQ = 2,
                            rm.phix = TRUE, 
                            compress = TRUE,
                            multithread = FALSE)
# Inspect filter results
head(out)


# Empty list for storing the mga results
RESULTS <- list()

# Loop the analysis process for every file in the directory
for (i in 1:length(fnFs)){
  
  # Assign metadata to be input as a list for mga argument
  meta_data <- list()
  for (j in 1:length(fnFs)){
    
    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(samdf[samdf$sample.ID %in% ID,])
  }
  
  # Store results for each output from the `mga` in list
  RESULTS[[i]] <- mga(fastq.Fs = fnFs[i],
                      fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                      filtFs = filtFs[i],
                      filtRs = filtRs[i], # file paths for filtered fastq files
                      refFasta = fungi,
                      metadata = meta_data[[j]],
                      tree.args = list(k = 4, 
                                       inv = 0.2,
                                       model = "GTR",
                                       rearrangement = "stochastic"), 
                      network.args = list(type = "taxa",
                                          distance = "jaccard",
                                          max.dist = 0.35))
  
  # Save the results as they are produced, in case of crashes/errors
  save(RESULTS, file = "RDS/preliminary_results_CC_mga.rda")
}

# # Remove Greenhouse samples [[44:141]]
# RESULTS <- RESULTS[c(1:43, 142:221)]
# # Save modified `ps.net` output
# save(RESULTS, file = "RDS/results_CC_mga.rda")

# Extract data frame of diversity metrics from each mga object 
RES <- list()
for (i in 1:length(RESULTS)){
  
  RES[[i]] <- RESULTS[[i]]$results.samples
}

# Build combined data frame of metrics
results.df <- do.call(rbind, RES)
results.df <- dplyr::relocate(results.df, "sample.ID", .before = "Shannon") # move sample.ID column to leftmost
# Save data frame of results
saveRDS(results.df, file = "RDS/results.df_mga.rds")




# OLD CODE ----

## TEST: SuRDC Sample S1 ----

# Reference FASTA file paths:
fungi <- "./sh_general_release_dynamic_29.11.2022.fasta" # file name for fungal reference FASTA

bacteria <- "./silva_nr99_v138.1_train_set.fa.gz" # file name for bacterial reference FASTA

eukaryotes <- "./sh_general_release_dynamic_all_29.11.2022.fasta" # file name for all eukaryotes reference FASTA

# inputs for ps.net2
fastq.Fs <- "CC_seq/S1_S354_L001_R1_001.fastq.gz"
fastq.Rs <- "CC_seq/S1_S354_L001_R2_001.fastq.gz"
filtFs <- "filtered/S1_F_filt.fastq.gz"
filtRs <- "filtered/S1_R_filt.fastq.gz"

# inputs for ps.net
fastq.fwd <- "CC_seq/S1_S354_L001_R1_001.fastq.gz"
fastq.rev <- "CC_seq/S1_S354_L001_R2_001.fastq.gz"
filt.fwd <- "filtered/S1_F_filt.fastq.gz"
filt.rev <- "filtered/S1_R_filt.fastq.gz"

S1_test <- ps.net(fastq.fwd, fastq.rev, filt.fwd, filt.rev, 
                  refFasta = fungi, metadata = NULL, 
                  network.args = list(type = "taxa",
                                      distance = "jaccard",
                                      max.dist = 0.55,
                                      keep.isolates = TRUE))

list(sample.ID = "S1", site = "SuRDC", block = 1, plot = 27, trt = 16, cover.crop = "Crescendo Ladino Clover")

# TEST FILTERING

fastq.Fs <- sort(list.files(seq_path, pattern="_R1_001.fastq.gz"))
fastq.Rs <- sort(list.files(seq_path, pattern="_R2_001.fastq.gz"))
fastq.Fs <- file.path(seq_path, fastq.Fs)
fastq.Rs <- file.path(seq_path, fastq.Rs)



# TEST: All Samples ----

# Load sample metadata
samdf <- read.csv("C:/Users/elinrg26/OneDrive - UBC/Research/GitHub/Cover Crop Study/data/Mapping file for ITS sequencing.csv")
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

source("./Microbiome_Analysis_Pipeline.R")

# Empty list for storing the results
RESULTS <- list()

# Loop the analysis process for every file in the directory
for (i in 1:length(fnFs)){
  
  # Assign metadata to be input as a list for ps.net argument
  meta_data <- list()
  for (j in 1:length(fnFs)){
    
    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(samdf[samdf$sample.ID %in% ID,])
  }
  
  # Create a new row of results for each output from the `ps.net`
  RESULTS[[i]] <- ps.net(fastq.Fs = fnFs[i],
                         fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                         filtFs = filtFs[i],
                         filtRs = filtRs[i], # file paths for filtered fastq files
                         refFasta = fungi,
                         metadata = meta_data[[j]],
                         tree.args = list(k = 4, 
                                          inv = 0.2,
                                          model = "GTR",
                                          rearrangement = "stochastic"), 
                         network.args = list(type = "taxa",
                                             distance = "jaccard",
                                             max.dist = 0.35,
                                             keep.isolates = TRUE))
  
  # Save the results as they are produced, in case of crashes/errors
  save(RESULTS, file = "RDS/preliminary_results_CC_fix.rda")
}

# Remove Greenhouse samples [[44:141]]
RESULTS <- RESULTS[c(1:43, 142:221)]
# Save modified `ps.net` output
save(RESULTS, file = "RDS/results_covercrop_noGreenhouse_fix.rda")

# Newly improved results.samples
RES <- list()
for (i in 1:length(RESULTS)){
  
  RES[[i]] <- RESULTS[[i]]$results.samples
}

# Build combined data frame of metrics
results.df <- do.call(rbind, RES)
results.df <- dplyr::relocate(results.df, "sample.ID", .before = "Shannon") # move sample.ID column to leftmost
# Save data frame of results
saveRDS(results.df, file = "RDS/results.df_noGreenhouse_fix.rds")


# test <- readRDS(file = "RDS/metadata_covercrop.rds")
# load(file = "RDS/preliminary_results_covercrop.rda")
# 
# for (i in 1:length(RESULTS)){
#   
#   RESULTS[[i]]$results.samples <- RESULTS[[i]]$results.samples[1,1:10]
#   RESULTS[[i]]$results.samples <- cbind(RESULTS[[i]]$results.samples, as.data.frame(test[[i]]))
# }
# 
# 
# RES <- list()
# # fix results.samples
# for (i in 1:length(RESULTS)){
#   
#   RES[[i]] <- RESULTS[[i]]$results.samples[1,1:10]
#   RES[[i]] <- cbind(RES[[i]], as.data.frame(test[[i]]))
# }
# 
# # View data frame of metrics:
# results.df <- do.call(rbind, RES) %>%
#   dplyr::select("sample.ID", everything())



# TEST: C1 ----

# Load sample metadata
samdf <- read.csv("C:/Users/elinrg26/OneDrive - UBC/Research/GitHub/Cover Crop Study/data/Mapping file for ITS sequencing.csv")
# Load reference FASTA file for taxonomic classification
fungi <- "./sh_general_release_dynamic_29.11.2022.fasta" # file name for fungal reference FASTA

# Define filter path for the forward and reverse reads of each sample fastq file
seq_path <- "./CC_Seq1"  # Replace with the directory for fastq file after unzipping the folder

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

source("./Microbiome_Analysis_Pipeline.R")

# Empty list for storing the results
RESULTS_2 <- list()

# Loop the analysis process for every file in the directory
for (i in 1:length(fnFs)){
  
  # Assign metadata to be input as a list for ps.net argument
  meta_data <- list()
  for (j in 1:length(fnFs)){
    
    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(samdf[samdf$sample.ID %in% ID,])
  }
  
  # Create a new row of results for each output from the `ps.net`
  RESULTS_2[[i]] <- ps.net(fastq.Fs = fnFs[i],
                         fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                         filtFs = filtFs[i],
                         filtRs = filtRs[i], # file paths for filtered fastq files
                         refFasta = fungi,
                         metadata = meta_data[[j]],
                         tree.args = list(k = 4, 
                                          inv = 0.2,
                                          model = "GTR",
                                          rearrangement = "stochastic"), 
                         network.args = list(type = "taxa",
                                             distance = "jaccard",
                                             max.dist = 0.35,
                                             keep.isolates = TRUE))
  
  # Save the results as they are produced, in case of crashes/errors
  save(RESULTS_2, file = "RDS/C1_3/results_CC_test2.rda")
}

# Newly improved results.samples
RES <- list()
for (i in 1:length(RESULTS_2)){
  
  RES[[i]] <- RESULTS_2[[i]]$results.samples
}

# Build combined data frame of metrics
results.df <- do.call(rbind, RES)
results.df <- dplyr::relocate(results.df, "sample.ID", .before = "Shannon") # move sample.ID column to leftmost
# Save data frame of results
saveRDS(results.df, file = "RDS/C1_3/results.df_test2.rds")

