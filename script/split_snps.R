#!/usr/bin/env Rscript


mylib <- "./R_lib"
dir.create(mylib, showWarnings = F)

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages('dplyr', repos = "https://cloud.r-project.org", lib = mylib)
if (!requireNamespace("parallel", quietly = TRUE)) install.packages('parallel', repos = "https://cloud.r-project.org", lib = mylib)
library(dplyr)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

snp_file <- args[1] #as.character(snakemake@input)

output <- args[2] #as.character(snakemake@output)

snp_info <- read.csv(snp_file, header = F, stringsAsFactors = F)
colnames(snp_info) <- c("id", "pos", "chr")

# random_snps <- snakemake@params$random_snps
random_snps <- TRUE

# snp_den <- snakemake@params$snp_den
snp_den <- c(args[3], args[4]) # c(snakemake@params$nb_snps, snakemake@params$chrom_length)

# batch <- as.numeric(snakemake@params$batch)
# size <- as.numeric(snakemake@params$batch_size)
batch <- 1
size <-  args[4] #snakemake@params$chrom_length

cores <-  args[5] #as.numeric(snakemake@threads)

# print(batch)
for (b in batch) {
  dir.create(dirname(output[b]), showWarnings = F)
  if (isTRUE(random_snps)) {
    if (is.null(snp_den) || length(snp_den) != 2) {
      snp_batch <- snp_info[sort(sample(1:nrow(snp_info), size)),]
    } else {
      n_snp <- snp_den[1]
      size <- snp_den[2]
      chr_seq <- seq(1, max(snp_info$pos), size)
      seq_ind <- seq_along(chr_seq)
      sample_pos <- mclapply(seq_ind, function(i) {
        if (i < max(seq_ind)) {
          pos <- snp_info$pos[snp_info$pos %in% chr_seq[i]:(chr_seq[i+1]-1)]
        } else {
          pos <- snp_info$pos[snp_info$pos > chr_seq[i]]
        }
        n_snp <- ifelse(n_snp > length(pos), length(pos), n_snp)
        s <- sample(pos, size = n_snp)
        return(s)
      }, mc.cores = cores, mc.preschedule = F) %>% unlist
      snp_batch <- snp_info %>% filter(pos %in% sample_pos)
    }

  } else {
    start <- size*(b-1) + 1
    end <- size*b
    if (end > nrow(snp_info)) end <- nrow(snp_info)

    snp_batch <- snp_info[start:end,]
  }

  write.table(snp_batch, output[b], sep = ",", row.names = F, col.names = F, quote = F)
}
