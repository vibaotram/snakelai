#!/usr/bin/env Rscript


mylib <- "/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.0"
.libPaths(mylib)
# dir.create(mylib, showWarnings = F)

# if (!requireNamespace("dplyr", quietly = TRUE)) install.packages('dplyr', repos = "https://cloud.r-project.org", lib = mylib, INSTALL_opts = "--no-lock")
# if (!requireNamespace("parallel", quietly = TRUE)) install.packages('parallel', repos = "https://cloud.r-project.org", lib = mylib, INSTALL_opts = "--no-lock")
library(dplyr, lib.loc = mylib)
library(parallel, lib.loc = mylib)
library(optparse, lib.loc = mylib)

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              help = "input file"),
  make_option(c("-o", "--output"),
              type = "character",
              help = "output file"),
  make_option(c("-n", "--nb_snps"),
              type = "character",
              help = "number of snps"),
  make_option(c("-w", "--window_size"),
              type = "integer",
              help = "window length"),
  make_option(c("-t", "--threads"),
              type = "integer",
              help = "number of threads")
)

myArgs <- parse_args(
  OptionParser(usage = "%prog [-i input] [-o output] [-n nb_snps] [-w window_size] [-t threads]",
  option_list = option_list)
  )

cores <- myArgs$threads

snp_file <- myArgs$input #as.character(snakemake@input)

output <- myArgs$output #as.character(snakemake@output)
outdir <- dirname(output)
dir.create(outdir, showWarnings = F)

snp_info <- read.csv(snp_file, header = F, stringsAsFactors = F)
colnames(snp_info) <- c("id", "pos", "chr")

n_snp <- myArgs$nb_snps
size <- as.numeric(myArgs$window_size)
if (n_snp == "all" && size != "chrom") { # subset all snps
  batch <- ceiling(nrow(snp_info)/size)
} else { # random snps
    n_snp <- as.numeric(n_snp)
    batch <- 1
}

# print(batch)
###--->>> fixing this code block
for (b in 1:batch) {
  dir.create(dirname(output[b]), showWarnings = F)
  if (is.numeric(n_snp) && size != "chrom") { # even sample
    chr_seq <- seq(1, max(snp_info$pos), size)
    seq_ind <- seq_along(chr_seq)
    sample_pos <- mclapply(seq_ind, function(i) {
      if (i < max(seq_ind)) {
        pos <- snp_info$pos[snp_info$pos %in% c(chr_seq[i]:(chr_seq[i+1]-1))]
      } else {
        pos <- snp_info$pos[snp_info$pos > chr_seq[i]]
      }
      n_snp <- ifelse(n_snp > length(pos), length(pos), n_snp)
      s <- sample(pos, size = n_snp, replace = F)
      return(s)
    }, mc.cores = cores, mc.preschedule = F) %>% unlist
    snp_batch <- snp_info %>% filter(pos %in% sample_pos)
  } else if (is.numeric(n_snp) && size == "chrom") { # random sample
    sample_pos <- sample(n_snp, snp_info$pos, replace = F)
    snp_batch <- snp_info %>% filter(pos %in% sample_pos)
  } else { # subset of all snps
    start <- size*(b-1) + 1
    end <- size*b
    if (end > nrow(snp_info)) end <- nrow(snp_info)
    snp_batch <- snp_info[start:end,]
  }

  write.table(snp_batch, file.path(outdir, paste0("test_", b, ".pos")), sep = ",", row.names = F, col.names = F, quote = F)
}

file.create(output)
