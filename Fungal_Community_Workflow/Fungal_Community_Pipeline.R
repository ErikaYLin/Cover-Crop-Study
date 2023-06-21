# set working directory
setwd("C:/Users/elinrg26/OneDrive - UBC/Research/GitHub/Cover Crop Study")

## TEST: SuRDC Sample S1 ----

# Reference FASTA file paths:
fungi <- "./sh_general_release_dynamic_29.11.2022.fasta" # file name for fungal reference FASTA

bacteria <- "./silva_nr99_v138.1_train_set.fa.gz" # file name for bacterial reference FASTA

eukaryotes <- "./sh_general_release_dynamic_all_29.11.2022.fasta" # file name for all eukaryotes reference FASTA

# inputs for ps.net2
fastq.Fs <- "SuRDC_seq/S1_S354_L001_R1_001.fastq.gz"
fastq.Rs <- "SuRDC_seq/S1_S354_L001_R2_001.fastq.gz"
filtFs <- "filtered/S1_F_filt.fastq.gz"
filtRs <- "filtered/S1_R_filt.fastq.gz"

# inputs for ps.net
fastq.fwd <- "SuRDC_seq/S1_S354_L001_R1_001.fastq.gz"
fastq.rev <- "SuRDC_seq/S1_S354_L001_R2_001.fastq.gz"
filt.fwd <- "filtered/S1_F_filt.fastq.gz"
filt.rev <- "filtered/S1_R_filt.fastq.gz"

S1_test <- ps.net(fastq.fwd, fastq.rev, filt.fwd, filt.rev, 
                  refFasta = fungi, metadata = NULL, 
                  network.args = list(type = "taxa",
                                      distance = "jaccard",
                                      max.dist = 0.55,
                                      keep.isolates = TRUE))

list(sample.ID = "S1", site = "SuRDC", block = 1, plot = 27, trt = 16, cover.crop = "Crescendo Ladino Clover")


# TEST: All Samples ----

# Load sample metadata
samdf <- as.list(read.csv("C:/Users/elinrg26/OneDrive - UBC/Research/GitHub/Cover Crop Study/Fungal Community Workflow/Mapping file for ITS sequencing.csv"))
# Load reference FASTA file for taxonomic classification
fungi <- "./sh_general_release_dynamic_29.11.2022.fasta" # file name for fungal reference FASTA

# Define filter path for the forward and reverse reads of each sample fastq file
seq_path <- "./Rosa_Seq"  # Replace with the directory for fastq file after unzipping the folder

# Sort forward and reverse reads to be in the same order
fnFs <- sort(list.files(seq_path, pattern="_R1_001.fastq.gz"))
fnRs <- sort(list.files(seq_path, pattern="_R2_001.fastq.gz"))
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(seq_path, fnFs)
fnRs <- file.path(seq_path, fnRs)

# Define file names for filtered fastq.gz files and assign file paths
filt_path <- file.path(seq_path, "filtered") # Place filtered files in new subdirectory if one does not already exist
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

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
                         filtFs,
                         filtRs, # file paths for filtered fastq files
                         refFasta = fungi,
                         metadata = meta_data[[j]],
                         network.args = list(type = "taxa",
                                             distance = "jaccard",
                                             max.dist = 0.55,
                                             keep.isolates = TRUE))
  
  # Save the results as they are produced, in case of crashes/errors
  save(RESULTS, file = "RDS/preliminary_results_covercrop.rda")
}


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


# Beta diversity ----
