## Quick script to generate the config.yaml file for pacbio-preprocessing pipeline.

## Rscript generate_config_run.r  bam_dir  out_dir phasing > config.yaml

args <- commandArgs(trailingOnly=TRUE)

## Check if the number of arguments is correct.
if (length(args) != 3) {
  stop("Usage: Rscript generate_config_run.r bam_dir out_dir phasing")
}

bam_dir <- args[1]
out_dir <- args[2]
phasing <- args[3]

## Check bam_dir exists. 
if (!file.exists(bam_dir)) {
  stop("Input directory does not exist.")
}

## Check if the in and out directories exist.
if (!file.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
## Check that phasing is either Yes or No. 
if (!phasing %in% c("Yes", "No")) {
  stop("Phasing must be either Yes or No.")
}

## Get all the bam files in the directory
bam_files <- list.files(bam_dir, pattern = ".bam$", full.names = TRUE)

## Check if bam_files is empty.
if (length(bam_files) == 0) {
  stop("No bam files found in the input directory.")
}

## Get sample names.
sample_names <- gsub(".bam|.aligned|.sort|.reads", "", basename(bam_files))

## Printing to stdout the config.yaml
cat("## Config file for pacbio-preprocessing pipeline\n\n")
cat(paste0('out_dir: "', out_dir, '"',"\n\n"))
cat(paste0('phasing: "', phasing, '"',"\n\n"))
cat("samples:\n")

for(i in 1:length(sample_names)){
  cat(paste0("   ", sample_names[i], ': "', bam_files[i], '"', "\n"))
  }
