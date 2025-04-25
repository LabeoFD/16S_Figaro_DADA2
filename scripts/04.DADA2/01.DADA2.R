# -------------------------
# METAGENOMIC WORKFLOW SCRIPT
# -------------------------

# Author: Hédia Tnani
# Created : 03/20/2025
# Last modified: 03/31/2025

########### SCRIPT ##############

### 1. SETUP & DATA IMPORT ###
# Load libraries
library(dada2)
library(phyloseq)
library(ggplot2)
library(ShortRead) # For readFastq() 
library(dplyr)
library(tidyverse)
library(Biostrings)  # For DNAString() etc.
library(patchwork)  # For combining plots
library(purrr)
library(furrr)
library(tictoc)  # For timing
library(kableExtra)
library(here)

# Check dada2 package version
packageVersion("dada2")
# [1] ‘1.36.0’

# List fastq files
list.files(here("data/processed-data/fully_trimmed_reads/figaro_300bp_reads"))
# [1] "Sample1_R1.fastq.gz"  "Sample1_R2.fastq.gz"  "Sample10_R1.fastq.gz"
# [4] "Sample10_R2.fastq.gz" "Sample11_R1.fastq.gz" "Sample11_R2.fastq.gz"
# [7] "Sample12_R1.fastq.gz" "Sample12_R2.fastq.gz" "Sample13_R1.fastq.gz"
# [10] "Sample13_R2.fastq.gz" "Sample14_R1.fastq.gz" "Sample14_R2.fastq.gz"
# [13] "Sample15_R1.fastq.gz" "Sample15_R2.fastq.gz" "Sample16_R1.fastq.gz"
# [16] "Sample16_R2.fastq.gz" "Sample17_R1.fastq.gz" "Sample17_R2.fastq.gz"
# [19] "Sample18_R1.fastq.gz" "Sample18_R2.fastq.gz" "Sample19_R1.fastq.gz"
# [22] "Sample19_R2.fastq.gz" "Sample2_R1.fastq.gz"  "Sample2_R2.fastq.gz" 
# [25] "Sample20_R1.fastq.gz" "Sample20_R2.fastq.gz" "Sample21_R1.fastq.gz"
# [28] "Sample21_R2.fastq.gz" "Sample22_R1.fastq.gz" "Sample22_R2.fastq.gz"
# [31] "Sample23_R1.fastq.gz" "Sample23_R2.fastq.gz" "Sample24_R1.fastq.gz"
# [34] "Sample24_R2.fastq.gz" "Sample25_R1.fastq.gz" "Sample25_R2.fastq.gz"
# [37] "Sample26_R1.fastq.gz" "Sample26_R2.fastq.gz" "Sample27_R1.fastq.gz"
# [40] "Sample27_R2.fastq.gz" "Sample28_R1.fastq.gz" "Sample28_R2.fastq.gz"
# [43] "Sample29_R1.fastq.gz" "Sample29_R2.fastq.gz" "Sample3_R1.fastq.gz" 
# [46] "Sample3_R2.fastq.gz"  "Sample30_R1.fastq.gz" "Sample30_R2.fastq.gz"
# [49] "Sample31_R1.fastq.gz" "Sample31_R2.fastq.gz" "Sample32_R1.fastq.gz"
# [52] "Sample32_R2.fastq.gz" "Sample33_R1.fastq.gz" "Sample33_R2.fastq.gz"
# [55] "Sample34_R1.fastq.gz" "Sample34_R2.fastq.gz" "Sample35_R1.fastq.gz"
# [58] "Sample35_R2.fastq.gz" "Sample36_R1.fastq.gz" "Sample36_R2.fastq.gz"
# [61] "Sample4_R1.fastq.gz"  "Sample4_R2.fastq.gz"  "Sample5_R1.fastq.gz" 
# [64] "Sample5_R2.fastq.gz"  "Sample6_R1.fastq.gz"  "Sample6_R2.fastq.gz" 
# [67] "Sample7_R1.fastq.gz"  "Sample7_R2.fastq.gz"  "Sample8_R1.fastq.gz" 
# [70] "Sample8_R2.fastq.gz"  "Sample9_R1.fastq.gz"  "Sample9_R2.fastq.gz"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(here("data/processed-data/fully_trimmed_reads/figaro_300bp_reads"), pattern="_R1.", full.names= TRUE))
# [1] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample1_R1.fastq.gz" 
# [2] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample10_R1.fastq.gz"
# [3] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample11_R1.fastq.gz"
# [4] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample12_R1.fastq.gz"
# [5] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample13_R1.fastq.gz"
# [6] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample14_R1.fastq.gz"
# [7] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample15_R1.fastq.gz"
# [8] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample16_R1.fastq.gz"
# [9] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample17_R1.fastq.gz"
# [10] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample18_R1.fastq.gz"
# [11] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample19_R1.fastq.gz"
# [12] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample2_R1.fastq.gz" 
# [13] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample20_R1.fastq.gz"
# [14] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample21_R1.fastq.gz"
# [15] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample22_R1.fastq.gz"
# [16] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample23_R1.fastq.gz"
# [17] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample24_R1.fastq.gz"
# [18] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample25_R1.fastq.gz"
# [19] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample26_R1.fastq.gz"
# [20] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample27_R1.fastq.gz"
# [21] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample28_R1.fastq.gz"
# [22] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample29_R1.fastq.gz"
# [23] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample3_R1.fastq.gz" 
# [24] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample30_R1.fastq.gz"
# [25] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample31_R1.fastq.gz"
# [26] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample32_R1.fastq.gz"
# [27] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample33_R1.fastq.gz"
# [28] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample34_R1.fastq.gz"
# [29] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample35_R1.fastq.gz"
# [30] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample36_R1.fastq.gz"
# [31] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample4_R1.fastq.gz" 
# [32] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample5_R1.fastq.gz" 
# [33] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample6_R1.fastq.gz" 
# [34] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample7_R1.fastq.gz" 
# [35] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample8_R1.fastq.gz" 
# [36] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample9_R1.fastq.gz" 

fnRs <- sort(list.files(here("data/processed-data/fully_trimmed_reads/figaro_300bp_reads"), pattern="_R2.", full.names= TRUE))
# [1] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample1_R2.fastq.gz" 
# [2] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample10_R2.fastq.gz"
# [3] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample11_R2.fastq.gz"
# [4] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample12_R2.fastq.gz"
# [5] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample13_R2.fastq.gz"
# [6] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample14_R2.fastq.gz"
# [7] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample15_R2.fastq.gz"
# [8] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample16_R2.fastq.gz"
# [9] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample17_R2.fastq.gz"
# [10] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample18_R2.fastq.gz"
# [11] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample19_R2.fastq.gz"
# [12] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample2_R2.fastq.gz" 
# [13] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample20_R2.fastq.gz"
# [14] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample21_R2.fastq.gz"
# [15] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample22_R2.fastq.gz"
# [16] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample23_R2.fastq.gz"
# [17] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample24_R2.fastq.gz"
# [18] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample25_R2.fastq.gz"
# [19] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample26_R2.fastq.gz"
# [20] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample27_R2.fastq.gz"
# [21] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample28_R2.fastq.gz"
# [22] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample29_R2.fastq.gz"
# [23] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample3_R2.fastq.gz" 
# [24] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample30_R2.fastq.gz"
# [25] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample31_R2.fastq.gz"
# [26] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample32_R2.fastq.gz"
# [27] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample33_R2.fastq.gz"
# [28] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample34_R2.fastq.gz"
# [29] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample35_R2.fastq.gz"
# [30] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample36_R2.fastq.gz"
# [31] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample4_R2.fastq.gz" 
# [32] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample5_R2.fastq.gz" 
# [33] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample6_R2.fastq.gz" 
# [34] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample7_R2.fastq.gz" 
# [35] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample8_R2.fastq.gz" 
# [36] "/home/hediatnani/16S_Figaro_DADA2/data/processed-data/fully_trimmed_reads/figaro_300bp_reads/Sample9_R2.fastq.gz" 

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) %>% str_sort(numeric = TRUE)
# [1] "Sample1"  "Sample2"  "Sample3"  "Sample4"  "Sample5"  "Sample6"  "Sample7" 
# [8] "Sample8"  "Sample9"  "Sample10" "Sample11" "Sample12" "Sample13" "Sample14"
# [15] "Sample15" "Sample16" "Sample17" "Sample18" "Sample19" "Sample20" "Sample21"
# [22] "Sample22" "Sample23" "Sample24" "Sample25" "Sample26" "Sample27" "Sample28"
# [29] "Sample29" "Sample30" "Sample31" "Sample32" "Sample33" "Sample34" "Sample35"
# [36] "Sample36"

# Inspect read quality profiles 
# plotQualityProfile` function displays the quality scores across each position in the reads
# Inspect read quality profiles of forward reads
plotQualityProfile(fnFs[1:length(sample.names)])

# Inspect read quality profiles of reverse reads
plotQualityProfile(fnRs[1:length(sample.names)])


##############################################################################
# 2) Define the "trimmed-data" directory
##############################################################################

# Define the path to the "trimmed-data" directory
trimmed_dir <- file.path(here("data/processed-data/"), "figaro-filtered-data")
# Check if the directory exists; if not, create it (including parent directories if needed)
if (!dir.exists(trimmed_dir)) {
  dir.create(trimmed_dir, recursive = TRUE)
}

##############################################################################
# 3) Define file paths for the trimmed output files
##############################################################################
# Trimmed forward reads
filtFs <- file.path(here("data/processed-data"), "figaro-filtered-data", paste0(sample.names, "_F_filt.fastq.gz"))
# Trimmed reverse reads
filtRs <- file.path(here("data/processed-data"), "figaro-filtered-data", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
# [1] "Sample1"  "Sample2"  "Sample3"  "Sample4"  "Sample5"  "Sample6"  "Sample7"  "Sample8" 
# [9] "Sample9"  "Sample10" "Sample11" "Sample12" "Sample13" "Sample14" "Sample15" "Sample16"
# [17] "Sample17" "Sample18" "Sample19" "Sample20" "Sample21" "Sample22" "Sample23" "Sample24"
# [25] "Sample25" "Sample26" "Sample27" "Sample28" "Sample29" "Sample30" "Sample31" "Sample32"
# [33] "Sample33" "Sample34" "Sample35" "Sample36"


names(filtRs) <- sample.names
# [1] "Sample1"  "Sample2"  "Sample3"  "Sample4"  "Sample5"  "Sample6"  "Sample7"  "Sample8" 
# [9] "Sample9"  "Sample10" "Sample11" "Sample12" "Sample13" "Sample14" "Sample15" "Sample16"
# [17] "Sample17" "Sample18" "Sample19" "Sample20" "Sample21" "Sample22" "Sample23" "Sample24"
# [25] "Sample25" "Sample26" "Sample27" "Sample28" "Sample29" "Sample30" "Sample31" "Sample32"
# [33] "Sample33" "Sample34" "Sample35" "Sample36"

# DEFINE FORWARD/REVERSE PRIMERS AND ADAPTERS
##############################################################################
# Primers = define and amplify the target region of interest.
# Adapters = enable the prepared library (PCR product) to be sequenced on a specific platform and allow sample identification (indexing).

# 16S primer
FWD_PRIMER  <- "CCTACGGGNGGCWGCAG"        # forward primer
REV_PRIMER  <- "GACTACHVGGGTATCTAATCC"    # reverse primer

# Check length
FWD_PRIMER_LEN <- nchar(FWD_PRIMER)
REV_PRIMER_LEN <- nchar(REV_PRIMER)

# Nextera-style adapter
FWD_ADPT <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"  # forward adapter
REV_ADPT <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" # reverse adapter

# Use Figaro-optimized parameters (high retention & good score)
# Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(298, 225),    # From Figaro optimal parameters
                     maxEE=c(4, 3),           # From Figaro optimal parameters
                     truncQ=2,                # Trim reads at first Q2 base
                     maxN=0,                  # Remove any reads with Ns
                     compress=TRUE,           # Compress output files
                     multithread=TRUE)        # Use multiple threads

# Look at the filtering results
out

#                        reads.in reads.out
# Sample1_R1.fastq.gz    148408    109921
# Sample10_R1.fastq.gz    28391     20892
# Sample11_R1.fastq.gz   147001    116337
# Sample12_R1.fastq.gz   123954    106476
# Sample13_R1.fastq.gz   230651    187618
# Sample14_R1.fastq.gz    72014     61112
# Sample15_R1.fastq.gz   130227    103363
# Sample16_R1.fastq.gz    83896     68406
# Sample17_R1.fastq.gz   245588    187716
# Sample18_R1.fastq.gz   303070    249760
# Sample19_R1.fastq.gz   214417    170319
# Sample2_R1.fastq.gz    297037    239274
# Sample20_R1.fastq.gz   124841     99810
# Sample21_R1.fastq.gz   119605     92843
# Sample22_R1.fastq.gz   255080    204508
# Sample23_R1.fastq.gz   135690    103664
# Sample24_R1.fastq.gz   279954    241258
# Sample25_R1.fastq.gz   180080    153475
# Sample26_R1.fastq.gz    28836     21287
# Sample27_R1.fastq.gz    72258     62755
# Sample28_R1.fastq.gz    60391     52541
# Sample29_R1.fastq.gz    80723     72204
# Sample3_R1.fastq.gz     42952     26725
# Sample30_R1.fastq.gz    98278     87746
# Sample31_R1.fastq.gz   103645     90279
# Sample32_R1.fastq.gz   120309    106973
# Sample33_R1.fastq.gz   162441    144222
# Sample34_R1.fastq.gz    45203     37793
# Sample35_R1.fastq.gz    52648     46260
# Sample36_R1.fastq.gz    30555     26284
# Sample4_R1.fastq.gz     87730     67832
# Sample5_R1.fastq.gz    126966    100372
# Sample6_R1.fastq.gz    259125    207230
# Sample7_R1.fastq.gz    226999    175727
# Sample8_R1.fastq.gz    124154     93428
# Sample9_R1.fastq.gz     22007     16870

# Calculates retention rates as a percentage, adds formatted retention rates with % symbol,
out_retained <- out %>% 
  as.data.frame() %>% 
  mutate(
    retention_rate = reads.out / reads.in * 100,
    retention_rate_formatted = paste(round(retention_rate, 2), "%"),
    sample_number = as.numeric(gsub("_R1.fastq.gz", "", rownames(.)))
  )

# Calculate summary statistics
summary_stats <- list(
  mean_retention = mean(out_retained$retention_rate),
  median_retention = median(out_retained$retention_rate),
  sd_retention = sd(out_retained$retention_rate),
  min_retention = min(out_retained$retention_rate),
  max_retention = max(out_retained$retention_rate),
  samples_below_75 = sum(out_retained$retention_rate < 75),
  samples_above_85 = sum(out_retained$retention_rate > 85)
)

# Calculate overall retention rate across all samples
overall_retention <- sum(out_retained$reads.out) / sum(out_retained$reads.in) * 100

# Output summary statistics
cat("Summary Statistics for Retention Rates:\n")
cat(sprintf("Overall retention rate: %.2f%%\n", overall_retention))
cat(sprintf("Mean retention rate: %.2f%%\n", summary_stats$mean_retention))
cat(sprintf("Median retention rate: %.2f%%\n", summary_stats$median_retention))
cat(sprintf("Range: %.2f%% - %.2f%%\n", summary_stats$min_retention, summary_stats$max_retention))
cat(sprintf("Samples with retention < 75%%: %d\n", summary_stats$samples_below_75))
cat(sprintf("Samples with retention > 85%%: %d\n", summary_stats$samples_above_85))


# Identify samples with retention rates below 75%
# Extract samples with retention rates below 75% using tidyverse approach
low_retention_samples <- out_retained %>%
  as_tibble(rownames = "sample_name") %>%  # Convert to tibble with rownames as a column
  filter(retention_rate < 75) %>%
  pull(sample_name) %>% 
   gsub("_R1.fastq.gz", "", .)  # Remove prefix and suffix # Extract just the sample names


# Identify samples with retention rates above 85%
high_retention_samples <- out_retained %>%
  as_tibble(rownames = "sample_name") %>%  # Convert to tibble with rownames as a column
  filter(retention_rate > 85) %>%
  pull(sample_name) %>% 
  gsub("_R1.fastq.gz", "", .) # Extract just the sample names

# Create a data frame for the summary statistics with sample names
stats_table <- data.frame(
  Statistic = c(
    "Overall retention rate",
    "Mean retention rate", 
    "Median retention rate",
    "Minimum retention rate",
    "Maximum retention rate",
    paste("Samples with retention < 75% (", length(low_retention_samples), ")", sep=""),
    paste("Samples with retention > 85% (", length(high_retention_samples), ")", sep="")
  ),
  Value = c(
    sprintf("%.2f%%", overall_retention),
    sprintf("%.2f%%", summary_stats$mean_retention),
    sprintf("%.2f%%", summary_stats$median_retention),
    sprintf("%.2f%%", summary_stats$min_retention),
    sprintf("%.2f%%", summary_stats$max_retention),
    paste(low_retention_samples, collapse = ", "),
    paste(high_retention_samples, collapse = ", ")
  )
)

# Generate the table as a string using kable 
table_string <- capture.output(
  kable(stats_table, 
        format = "pipe", 
        caption = "Summary Statistics for Retention Rates")
)

# Print the table using cat
cat(paste(table_string, collapse = "\n"), "\n")


# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
# 105380548 total bases in 353626 reads from 4 samples will be used for learning the error rates.

# When learnErrors() runs, it doesn't use all of your samples to learn the error model. 
# It randomly selects a subset of samples and reads that provide sufficient coverage to build a reliable error model.
# Only used 4 samples (randomly selected) out of your 36 samples
# From those 4 samples, it used 353,626 reads containing 105,380,548 bases

errR <- learnErrors(filtRs, multithread=TRUE)
# 121779900 total bases in 541244 reads from 5 samples will be used for learning the error rates.

# Visualize the estimated error rates for the forward reads
plotErrors(errF, nominalQ = TRUE)

# Visualize the estimated error rates for the reverse reads
plotErrors(errR, nominalQ = TRUE)


# Sample Inference
# Dereplicate filtered reads
derepFs <- derepFastq(filtFs, verbose = TRUE)


derepRs <- derepFastq(filtRs, verbose = TRUE)


# Compute dereplication stats using purrr
derep_stats <- tibble(
  Sample = names(derepFs),
  
  # Number of unique sequences
  UniqueF = map_int(derepFs, ~length(.$uniques)),
  UniqueR = map_int(derepRs, ~length(.$uniques)),
  
  # Total reads per sample (sum of abundance counts)
  TotalF = map_int(derepFs, ~sum(.$uniques)),
  TotalR = map_int(derepRs, ~sum(.$uniques)),
  
  # Compression = reads per unique sequence
  CompressionF = map2_dbl(TotalF, UniqueF, ~round(.x / .y, 1)),
  CompressionR = map2_dbl(TotalR, UniqueR, ~round(.x / .y, 1))
)

# Print summary
cat("Dereplication summary:\n")
cat("- Total unique forward sequences:", sum(derep_stats$UniqueF), "\n")
cat("- Total unique reverse sequences:", sum(derep_stats$UniqueR), "\n")
cat("- Average compression ratio (forward):", round(mean(derep_stats$CompressionF), 1), "reads per unique sequence\n")
cat("- Average compression ratio (reverse):", round(mean(derep_stats$CompressionR), 1), "reads per unique sequence\n")

# Denoise with DADA2
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
# Sample 1 - 109921 reads in 63358 unique sequences.
# Sample 2 - 20892 reads in 11894 unique sequences.
# Sample 3 - 116337 reads in 44563 unique sequences.
# Sample 4 - 106476 reads in 31945 unique sequences.
# Sample 5 - 187618 reads in 70869 unique sequences.
# Sample 6 - 61112 reads in 21326 unique sequences.
# Sample 7 - 103363 reads in 43172 unique sequences.
# Sample 8 - 68406 reads in 26317 unique sequences.
# Sample 9 - 187716 reads in 75832 unique sequences.
# Sample 10 - 249760 reads in 88350 unique sequences.
# Sample 11 - 170319 reads in 72737 unique sequences.
# Sample 12 - 239274 reads in 95266 unique sequences.
# Sample 13 - 99810 reads in 40595 unique sequences.
# Sample 14 - 92843 reads in 44478 unique sequences.
# Sample 15 - 204508 reads in 82624 unique sequences.
# Sample 16 - 103664 reads in 52107 unique sequences.
# Sample 17 - 241258 reads in 64920 unique sequences.
# Sample 18 - 153475 reads in 51095 unique sequences.
# Sample 19 - 21287 reads in 12619 unique sequences.
# Sample 20 - 62755 reads in 22499 unique sequences.
# Sample 21 - 52541 reads in 20384 unique sequences.
# Sample 22 - 72204 reads in 23254 unique sequences.
# Sample 23 - 26725 reads in 16147 unique sequences.
# Sample 24 - 87746 reads in 23013 unique sequences.
# Sample 25 - 90279 reads in 38237 unique sequences.
# Sample 26 - 106973 reads in 34298 unique sequences.
# Sample 27 - 144222 reads in 53830 unique sequences.
# Sample 28 - 37793 reads in 15804 unique sequences.
# Sample 29 - 46260 reads in 19703 unique sequences.
# Sample 30 - 26284 reads in 13227 unique sequences.
# Sample 31 - 67832 reads in 32300 unique sequences.
# Sample 32 - 100372 reads in 41221 unique sequences.
# Sample 33 - 207230 reads in 84445 unique sequences.
# Sample 34 - 175727 reads in 73034 unique sequences.
# Sample 35 - 93428 reads in 47635 unique sequences.
# Sample 36 - 16870 reads in 9410 unique sequences.

dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
# Sample 1 - 109921 reads in 54295 unique sequences.
# Sample 2 - 20892 reads in 11766 unique sequences.
# Sample 3 - 116337 reads in 39998 unique sequences.
# Sample 4 - 106476 reads in 29913 unique sequences.
# Sample 5 - 187618 reads in 61226 unique sequences.
# Sample 6 - 61112 reads in 19454 unique sequences.
# Sample 7 - 103363 reads in 38349 unique sequences.
# Sample 8 - 68406 reads in 24987 unique sequences.
# Sample 9 - 187716 reads in 64199 unique sequences.
# Sample 10 - 249760 reads in 71444 unique sequences.
# Sample 11 - 170319 reads in 63687 unique sequences.
# Sample 12 - 239274 reads in 77774 unique sequences.
# Sample 13 - 99810 reads in 29932 unique sequences.
# Sample 14 - 92843 reads in 40226 unique sequences.
# Sample 15 - 204508 reads in 72843 unique sequences.
# Sample 16 - 103664 reads in 44926 unique sequences.
# Sample 17 - 241258 reads in 56568 unique sequences.
# Sample 18 - 153475 reads in 46477 unique sequences.
# Sample 19 - 21287 reads in 11461 unique sequences.
# Sample 20 - 62755 reads in 22287 unique sequences.
# Sample 21 - 52541 reads in 19281 unique sequences.
# Sample 22 - 72204 reads in 23644 unique sequences.
# Sample 23 - 26725 reads in 15761 unique sequences.
# Sample 24 - 87746 reads in 24697 unique sequences.
# Sample 25 - 90279 reads in 37023 unique sequences.
# Sample 26 - 106973 reads in 33134 unique sequences.
# Sample 27 - 144222 reads in 50476 unique sequences.
# Sample 28 - 37793 reads in 16790 unique sequences.
# Sample 29 - 46260 reads in 20996 unique sequences.
# Sample 30 - 26284 reads in 12824 unique sequences.
# Sample 31 - 67832 reads in 30657 unique sequences.
# Sample 32 - 100372 reads in 35615 unique sequences.
# Sample 33 - 207230 reads in 71906 unique sequences.
# Sample 34 - 175727 reads in 63682 unique sequences.
# Sample 35 - 93428 reads in 40744 unique sequences.
# Sample 36 - 16870 reads in 8398 unique sequences.

# Robust denoising stats
denoising_stats <- tibble(
  Sample = names(dadaFs),
  
  DenoisedF = map_int(dadaFs, ~ sum(.x$denoised)),
  DenoisedR = map_int(dadaRs, ~ sum(.x$denoised))
)

# Print summary
cat("Denoising summary:\n")
cat("- Average denoised reads per sample (F):", round(mean(denoising_stats$DenoisedF), 0), "\n")
cat("- Average denoised reads per sample (R):", round(mean(denoising_stats$DenoisedR), 0), "\n")


# Merge Paired-End Reads
mergers <- mergePairs(
  dadaFs, filtFs, 
  dadaRs, filtRs, 
  minOverlap = 20,   # Require ≥20bp overlap
  maxMismatch = 0    # Allow 0 mismatches in overlap
) 


# Calculate average merged sequence length per sample
avg_seq_length <- tibble(
  Sample = names(mergers),
  
  # Weighted average length = sum(length * abundance) / sum(abundance)
  AvgMergedLength = map_dbl(mergers, ~ {
    seqs <- .$sequence
    abunds <- .$abundance
    weighted.mean(nchar(seqs), abunds)
  })
)

# Print
print(avg_seq_length)

# Summary
cat("Average merged sequence length across samples:", round(mean(avg_seq_length$AvgMergedLength), 1), "bp\n")

# Create sequence table
seqtab <- makeSequenceTable(mergers)

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab, 
  method = "consensus", # "consensus" - compares across samples
  multithread = TRUE, 
  verbose = TRUE
)

# Identified 17701 bimeras out of 37424 input sequences.

# Track reads at each step (optional)
getN <- function(x) sum(getUniques(x))
track.nbr.reads <- cbind(
  out, 
  sapply(dadaFs, getN), 
  sapply(dadaRs, getN), 
  sapply(mergers, getN), 
  rowSums(seqtab.nochim)
)
colnames(track.nbr.reads) <- c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "NonChim")
rownames(track.nbr.reads) <- sample.names
head(track.nbr.reads)
#         Input Filtered DenoisedF DenoisedR Merged NonChim
# Sample1 148408   109921    102776    103978  87176   57533
# Sample2  28391    20892     19012     19244  16921   11149
# Sample3 147001   116337    113988    114202 109819   70225
# Sample4 123954   106476    104956    105492 102798   72320
# Sample5 230651   187618    185241    185458 179518  102861
# Sample6  72014    61112     59852     60011  57877   37970

# Convert the matrix to a dataframe 
track_df <- as.data.frame(track.nbr.reads)

# ==== Calculate Percent Retained and Chimera Rate ====
# For each processing step, calculate the percentage of input reads retained
# Also calculate the chimera removal rate: (Merged - NonChim) / Merged
track_df <- track_df %>%
  mutate(
    Filtered_pct  = round(Filtered / Input * 100, 1),
    DenoisedF_pct = round(DenoisedF / Input * 100, 1),
    DenoisedR_pct = round(DenoisedR / Input * 100, 1),
    Merged_pct    = round(Merged / Input * 100, 1),
    NonChim_pct   = round(NonChim / Input * 100, 1),
    ChimeraRate   = round((Merged - NonChim) / Merged * 100, 1)
  )

# ==== Calculate Summary Statistics ====
# Compute averages and central tendencies for reporting
avg_filtered          <- mean(track_df$Filtered_pct, na.rm = TRUE)
avg_denoisedF         <- mean(track_df$DenoisedF_pct, na.rm = TRUE)
avg_denoisedR         <- mean(track_df$DenoisedR_pct, na.rm = TRUE)
avg_merged            <- mean(track_df$Merged_pct, na.rm = TRUE)
avg_nonchim           <- mean(track_df$NonChim_pct, na.rm = TRUE)
avg_chimera           <- mean(track_df$ChimeraRate, na.rm = TRUE)
mean_nonchim_reads    <- mean(track_df$NonChim, na.rm = TRUE)
median_nonchim_reads  <- median(track_df$NonChim, na.rm = TRUE)

# ==== Console Summary Output ====
# Use cat() to log the QC summary in a readable form
cat("===== DADA2 Read Tracking Summary (All Samples) =====\n")
cat("- Avg % retained after filtering:     ", avg_filtered, "%\n")
cat("- Avg % denoised (forward reads):     ", avg_denoisedF, "%\n")
cat("- Avg % denoised (reverse reads):     ", avg_denoisedR, "%\n")
cat("- Avg % merged:                       ", avg_merged, "%\n")
cat("- Avg % non-chimeric:                 ", avg_nonchim, "%\n")
cat("- Avg chimera removal rate:           ", avg_chimera, "%\n")
cat("- Mean non-chimeric reads/sample:     ", round(mean_nonchim_reads), "\n")
cat("- Median non-chimeric reads/sample:   ", round(median_nonchim_reads), "\n")

# If track_df doesn't have a Sample column, recreate it
if (!"Sample" %in% colnames(track_df)) {
  track_df <- track_df %>%
    mutate(Sample = rownames(track.nbr.reads)) %>%
    relocate(Sample)  # Move Sample to the front
}

# ==== Identify Best/Worst Retained Samples (by % NonChim) ====
top5_best <- track_df %>%
  arrange(desc(NonChim_pct)) %>%
  select(Sample, NonChim, NonChim_pct) %>%
  head(5)

bottom5_worst <- track_df %>%
  arrange(NonChim_pct) %>%
  select(Sample, NonChim, NonChim_pct) %>%
  head(5)

# ==== Print Summaries ====
cat("\n===== Top 5 Best Retained Samples (by % NonChim) =====\n")
print(top5_best)

cat("\n===== Bottom 5 Worst Retained Samples (by % NonChim) =====\n")
print(bottom5_worst)


# Convert the object to a tibble 
track_tibble <- as_tibble(track.nbr.reads, rownames = "Sample")

# Calculate percentages using purrr and dplyr
track_percentages <- track_tibble %>%
  mutate(across(-Sample, 
                ~ round(. / Input * 100, 1),
                .names = "{.col}_pct"))

# Format percentages with % symbol
track_formatted <- track_percentages %>%
  mutate(across(ends_with("_pct"), 
                ~ paste0(., "%"),
                .names = "{.col}_formatted"))

# Create a clean display table
result_table <- track_formatted %>%
  select(Sample, -Input, ends_with("_formatted")) %>%
  rename_with(~ gsub("_pct_formatted", " (%)", .), ends_with("_formatted"))


# Test different values of truncLen
# Define parameters
params <- list(
  truncLen = list(c(298, 225), c(291, 232), c(277, 246)),
  maxEE = list(c(4, 3), c(3, 3), c(3, 4)),
  suffix = c("298_225", "291_232", "277_246")
)

# #Apply filterAndTrim to each parameter set
# results <- pmap(params, function(truncLen, maxEE, suffix) {
#   filterAndTrim(
#     fnFs,
#     paste0(filtFs, "_", suffix),
#     fnRs,
#     paste0(filtRs, "_", suffix),
#     truncLen = truncLen,
#     maxEE = maxEE,
#     truncQ = 2,
#     maxN = 0,
#     compress = TRUE,
#     multithread = TRUE
#   )
# })
# 
# #Name the list elements according to the suffix values
# names(results) <- params$suffix

# Create a named list that will be returned with suffixes as names
results <- pmap(params, function(truncLen, maxEE, suffix) {
  filterAndTrim(
    fnFs,
    paste0(filtFs, "_", suffix),
    fnRs,
    paste0(filtRs, "_", suffix),
    truncLen = truncLen,
    maxEE = maxEE,
    truncQ = 2,
    maxN = 0,
    compress = TRUE,
    multithread = TRUE
  )
}) %>% 
  set_names(params$suffix)

# Create a function to process each parameter set through the DADA2 pipeline
process_dada2 <- function(filtFs_suffix, filtRs_suffix, suffix) {
  # Learn error rates
  errF <- learnErrors(filtFs_suffix, multithread=TRUE)
  errR <- learnErrors(filtRs_suffix, multithread=TRUE)
  
  # Dereplicate filtered reads
  derepFs <- derepFastq(filtFs_suffix, verbose = TRUE)
  derepRs <- derepFastq(filtRs_suffix, verbose = TRUE)
  
  # Denoise with DADA2
  dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
  dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
  
  # Merge Paired-End Reads
  mergers <- mergePairs(
    dadaFs, filtFs_suffix, 
    dadaRs, filtRs_suffix, 
    minOverlap = 20,
    maxMismatch = 0
  )
  
  # Create sequence table
  seqtab <- makeSequenceTable(mergers)
  
  # Remove Chimeras
  seqtab.nochim <- removeBimeraDenovo(
    seqtab, 
    method = "consensus",
    multithread = TRUE, 
    verbose = TRUE
  )
  
  # Track reads through the pipeline
  getN <- function(x) sum(getUniques(x))
  track.reads <- cbind(
    out,  # Make sure 'out' contains your initial read counts
    sapply(dadaFs, getN),
    sapply(dadaRs, getN),
    sapply(mergers, getN),
    rowSums(seqtab.nochim)
  )
  colnames(track.reads) <- c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "NonChim")
  rownames(track.reads) <- sample.names
  
  # Return a list of the results
  return(list(
    errF = errF,
    errR = errR,
    dadaFs = dadaFs,
    dadaRs = dadaRs,
    mergers = mergers,
    seqtab = seqtab,
    seqtab.nochim = seqtab.nochim,
    tracking = track.reads
  ))
}

# Apply the function to each parameter set
dada2_results <- map2(params$suffix, names(results), function(suffix, result_name) {
  filtFs_suffix <- paste0(filtFs, "_", suffix)
  filtRs_suffix <- paste0(filtRs, "_", suffix)
  
  process_dada2(filtFs_suffix, filtRs_suffix, suffix)
}) %>% 
  set_names(params$suffix)

# Extract all tracking tables 
combined_tables <- dada2_results %>% 
  map(~.$tracking) %>% 
  imap(~{colnames(.x) <- paste0(colnames(.x), "_", .y); as.data.frame(.x)}) %>% 
  purrr::reduce(cbind)

# Extract tracking data, reshape to long format, and plot with log scale
dada2_results %>% 
  map(~.$tracking) %>% 
  # Convert to dataframes with identifiers
  imap_dfr(~{
    as.data.frame(.x) %>% 
      mutate(Sample = rownames(.)) %>% 
      mutate(Parameter_Set = .y)
  }) %>% 
  # Reshape to long format for ggplot
  pivot_longer(
    cols = c(Input, Filtered, DenoisedF, DenoisedR, Merged, NonChim),
    names_to = "Stage",
    values_to = "Reads"
  ) %>% 
  # Convert Stage to factor to maintain order
  mutate(Stage = factor(Stage, 
                        levels = c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "NonChim"))) %>%
  # Create the plot with log scale
  ggplot(aes(x = Stage, y = Reads, color = Parameter_Set, group = interaction(Sample, Parameter_Set))) +
  geom_line(alpha = 0.3) +  # Individual sample lines
  geom_point(alpha = 0.3) +  # Points for each sample
  stat_summary(aes(group = Parameter_Set), fun = mean, geom = "line", size = 1.5) +  # Mean lines
  stat_summary(fun = mean, geom = "point", size = 3) +  # Mean points
  scale_y_log10() +  # Log scale for y-axis
  facet_wrap(~Parameter_Set) +  # Separate panels for each parameter set
  theme_bw() +
  labs(
    title = "Read counts through DADA2 pipeline steps (log scale)",
    subtitle = "Comparing different filtering parameters",
    x = "Pipeline Stage",
    y = "Number of Reads (log10 scale)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


dada2_results %>% 
  map(~.$tracking) %>% 
  imap_dfr(~{
    as.data.frame(.x) %>% 
      mutate(Sample = rownames(.)) %>% 
      mutate(Parameter_Set = .y)
  }) %>% 
  # Calculate retention percentage
  group_by(Sample, Parameter_Set) %>%
  mutate(across(c(Filtered, DenoisedF, DenoisedR, Merged, NonChim), 
                ~ 100 * . / Input)) %>%
  ungroup() %>%
  # Add Input as 100%
  mutate(Input = 100) %>%
  # Reshape to long format
  pivot_longer(
    cols = c(Input, Filtered, DenoisedF, DenoisedR, Merged, NonChim),
    names_to = "Stage",
    values_to = "Percentage"
  ) %>%
  mutate(Stage = factor(Stage, 
                        levels = c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "NonChim"))) %>%
  # Calculate mean percentage per parameter set and stage
  group_by(Parameter_Set, Stage) %>%
  summarize(Mean_Percentage = mean(Percentage, na.rm = TRUE),
            SE = sd(Percentage, na.rm = TRUE)/sqrt(n()),
            .groups = "drop") %>%
  # Plot mean percentages
  ggplot(aes(x = Stage, y = Mean_Percentage, color = Parameter_Set, group = Parameter_Set)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean_Percentage - SE, ymax = Mean_Percentage + SE), width = 0.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_bw() +
  labs(
    title = "Percentage of reads retained through DADA2 pipeline",
    subtitle = "Comparing different filtering parameters",
    x = "Pipeline Stage",
    y = "Percentage of Input Reads Retained"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create a two-panel figure with both absolute counts and percentages
p1 <- dada2_results %>% 
  map(~.$tracking) %>% 
  imap_dfr(~{
    as.data.frame(.x) %>% 
      mutate(Sample = rownames(.)) %>% 
      mutate(Parameter_Set = .y)
  }) %>% 
  pivot_longer(
    cols = c(Input, Filtered, DenoisedF, DenoisedR, Merged, NonChim),
    names_to = "Stage",
    values_to = "Reads"
  ) %>% 
  mutate(Stage = factor(Stage, 
                        levels = c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "NonChim"))) %>%
  # Calculate mean reads per parameter set and stage
  group_by(Parameter_Set, Stage) %>%
  summarize(Mean_Reads = mean(Reads, na.rm = TRUE), 
            .groups = "drop") %>%
  # Plot just the means on a single plot
  ggplot(aes(x = Stage, y = Mean_Reads, color = Parameter_Set, group = Parameter_Set)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "A) Average read counts",
    x = NULL,
    y = "Mean Reads (log scale)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p2 <- dada2_results %>% 
  map(~.$tracking) %>% 
  imap_dfr(~{
    as.data.frame(.x) %>% 
      mutate(Sample = rownames(.)) %>% 
      mutate(Parameter_Set = .y)
  }) %>% 
  # Calculate retention percentage
  group_by(Sample, Parameter_Set) %>%
  mutate(across(c(Filtered, DenoisedF, DenoisedR, Merged, NonChim), 
                ~ 100 * . / Input)) %>%
  ungroup() %>%
  # Add Input as 100%
  mutate(Input = 100) %>%
  # Reshape to long format
  pivot_longer(
    cols = c(Input, Filtered, DenoisedF, DenoisedR, Merged, NonChim),
    names_to = "Stage",
    values_to = "Percentage"
  ) %>%
  mutate(Stage = factor(Stage, 
                        levels = c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "NonChim"))) %>%
  # Calculate mean percentage per parameter set and stage
  group_by(Parameter_Set, Stage) %>%
  summarize(Mean_Percentage = mean(Percentage, na.rm = TRUE),
            .groups = "drop") %>%
  # Plot mean percentages
  ggplot(aes(x = Stage, y = Mean_Percentage, color = Parameter_Set, group = Parameter_Set)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  # Add percentage labels to the NonChim point
  geom_text(data = . %>% filter(Stage == "NonChim"), 
            aes(label = sprintf("%.1f%%", Mean_Percentage)),
            vjust = -0.8, size = 3.5) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_bw() +
  labs(
    title = "B) Percentage retention",
    x = "Pipeline Stage",
    y = "% of Input Reads Retained"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine the plots
combined_plot <- p1 / p2 + plot_layout(guides = "collect")

# Add a shared title
combined_plot <- combined_plot + 
  plot_annotation(
    title = "DADA2 Pipeline Performance Across Different Filtering Parameters",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Print combined plot
combined_plot


# Make the data in the long format
comparison_data <- dada2_results %>% 
  map(~.$tracking) %>% 
  imap_dfr(~{
    as.data.frame(.x) %>% 
      mutate(Sample = rownames(.)) %>% 
      mutate(Parameter_Set = .y)
  }) %>% 
  # Calculate retention percentage
  group_by(Sample, Parameter_Set) %>%
  mutate(across(c(Filtered, DenoisedF, DenoisedR, Merged, NonChim), 
                ~ 100 * . / Input)) %>%
  ungroup() %>%
  # Reshape to long format
  pivot_longer(
    cols = c(Filtered, DenoisedF, DenoisedR, Merged, NonChim),
    names_to = "Stage",
    values_to = "Percentage"
  ) %>%
  mutate(Stage = factor(Stage, 
                        levels = c("Filtered", "DenoisedF", "DenoisedR", "Merged", "NonChim")))

# Create side-by-side boxplots with jittered points
ggplot(comparison_data, aes(x = Stage, y = Percentage, fill = Parameter_Set, color = Parameter_Set)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5, outlier.shape = NA, width = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("277_246" = "#E41A1C", "291_232" = "#4DAF4A", "298_225" = "#377EB8")) +
  scale_color_manual(values = c("277_246" = "#E41A1C", "291_232" = "#4DAF4A", "298_225" = "#377EB8")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_bw() +
  labs(
    title = "Read retention through DADA2 pipeline stages",
    subtitle = "Comparison of different filtering parameters",
    x = "Pipeline Stage",
    y = "Percentage of Input Reads Retained"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.position = "top"
  )

# Focus on Final Retention 
dada2_results %>% 
  map(~.$tracking) %>% 
  imap_dfr(~{
    as.data.frame(.x) %>% 
      mutate(Sample = rownames(.)) %>% 
      mutate(Parameter_Set = .y)
  }) %>% 
  # Calculate final retention percentage
  mutate(Final_Retention_Pct = 100 * NonChim / Input) %>%
  # Plot boxplot of final retention by parameter set
  ggplot(aes(x = Parameter_Set, y = Final_Retention_Pct, fill = Parameter_Set)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_bw() +
  labs(
    title = "Final read retention after DADA2 pipeline",
    subtitle = "Comparing different filtering parameters",
    x = "Parameter Set",
    y = "Percentage of Input Reads Retained After Chimera Removal"
  )

### 3. TAXONOMIC ASSIGNATION ###

# for taxonomic reference data have a look at https://benjjneb.github.io/dada2/training.html
# Silva version 138.1: https://zenodo.org/records/4587955
# Silva version 138.2: https://zenodo.org/records/14169026
# GTDB r220: https://zenodo.org/records/13984843
# Refseq v16: https://zenodo.org/records/10403693
# Greengenes2: https://zenodo.org/records/14169078
# PR2: https://github.com/pr2database/pr2database/releases


# Downloaded SILVA file
silva_ref <- here("DB/Silva/v138_2/silva_nr99_v138.2_toGenus_trainset.fa.gz")
# RefSeq_ref  <- here("DB/RefSeq/v16/RefSeq_16S_6-11-20_RDPv16_Genus.fa")
# Greengenes2_ref  <- here("DB/Greengenes2/09/gg2_2024_09_toGenus_trainset.fa.gz")
# 

# Assign Taxonomy with SILVA
taxa <- assignTaxonomy(
  seqtab.nochim,
  refFasta = silva_ref,
  multithread = TRUE,
  tryRC = TRUE,  # Check reverse complements
)

# Special case
# Greengenes2_species_ref <- here("DB/Greengenes2/09/gg2_2024_09_toSpecies_trainset.fa.gz")
# # 
# taxa_Species <- assignTaxonomy(
#   seqtab.nochim, 
#   refFasta = Greengenes2_species_ref, 
#   multithread = TRUE, 
#   tryRC = TRUE,  # Check reverse complements
# )

# Add species assignment
silva_species_ref <- here("DB/Silva/v138_2/silva_v138.2_assignSpecies.fa.gz")
# RefSeq_species_ref <- here("DB/RefSeq/v16/RefSeq_16S_6-11-20_RDPv16_Species.fa")

# Add species-level assignment
taxa <- addSpecies(
  taxa, 
  silva_species_ref,
  verbose = TRUE 
)
# 4 out of 451 were assigned to the species level.
# Of which 4 had genera consistent with the input table.>  

# Count taxonomy
unname(taxa) %>%
  as.data.frame() %>% 
  dplyr::rename(Genus = V1) %>% 
  group_by(Genus) %>% 
  summarize(count= n())
# A tibble: 4 × 2
# Genus     count
#  <chr>     <int>
# 1 Archaea      30
# 2 Bacteria  18719
# 3 Eukaryota   868
# 4 NA          106
