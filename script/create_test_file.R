#!/usr/bin/env Rscript

mylib <- "/home/baotram/R/x86_64-pc-linux-gnu-library/4.0"
.libPaths(mylib)
# dir.create(mylib, showWarnings = F)

# if (!requireNamespace("dplyr", quietly = TRUE)) install.packages('dplyr', repos = "https://cloud.r-project.org", lib = mylib, INSTALL_opts = "--no-lock")
# if (!requireNamespace("parallel", quietly = TRUE)) install.packages('parallel', repos = "https://cloud.r-project.org", lib = mylib, INSTALL_opts = "--no-lock")
# if (!requireNamespace("vcfR", quietly = TRUE)) install.packages('vcfR', repos = "https://cloud.r-project.org", lib = mylib, INSTALL_opts = "--no-lock")
# if (!requireNamespace("adegenet", quietly = TRUE)) install.packages('adegenet', repos = "https://cloud.r-project.org", dependencies = TRUE, lib = mylib, INSTALL_opts = "--no-lock")
library(dplyr, lib.loc = mylib)
library(parallel, lib.loc = mylib)
library(vcfR, lib.loc = mylib)
library(adegenet, lib.loc = mylib)
library(data.table, lib.loc = mylib)
library(optparse, lib.loc = mylib)

option_list <- list(
  make_option(c("-t", "--test"),
              type = "character",
              help = "test file"),
  make_option(c("-s", "--snp"),
              type = "character",
              help = "snp file"),
  make_option(c("-g", "--genotype_id"),
              type = "character",
              help = "genotype id"),
  make_option(c("-v", "--vcf"),
              type = "integer",
              help = "vcf file"),
  make_option(c("-t", "--threads"),
              type = "integer",
              help = "number of threads")
)

myArgs <- parse_args(
  OptionParser(usage = "%prog [-i input] [-o output] [-n nb_snps] [-w window_size] [-t threads]",
               option_list = option_list)
)

source("/data3/projects/vietcaf/baotram/scripts/robusta_vn/R/geno_functions.R", local = knitr::knit_global())

# args <- commandArgs(trailingOnly = TRUE)

test_file <- myArgs$test #snakemake@output$test_file
snp_file <- myArgs$snp #snakemake@output$snp_file
genotype_id <- myArgs$genotype_id #snakemake@params$genotype_id
vcf_file <- myArgs$vcf #snakemake@input
cores <- myArgs$threads



vcf_R <- read.vcfR(vcf_file)
# chr_gl <- vcfR2genlight(vcf_R, n.cores = cores)

# chr_geno <- as.matrix(chr_gl)
# chr_geno[is.na(chr_geno)] <- 9
#
# # test_geno <- t(chr_geno[genotype_id,])
# test_geno <- t(chr_geno)
# test_gt <- geno2elai_gt(test_geno, chr_snps)
# rownames(test_gt) <- genotype_id
# colnames(test_gt) <- chr_snps$ID
test_gt <- extract.gt(vcf_R, return.alleles = T)
test_gt <- gsub("/", "", test_gt, fixed = T)
test_gt <- gsub(".", "??", test_gt, fixed = T)
write_elai_geno(test_gt, test_file)


chr_snps <- as.data.frame(getFIX(vcf_R), stringsAsFactors = F)
chr_snps$ID <- paste(chr_snps$CHROM, chr_snps$POS, sep = "_")
snp_pos <- chr_snps %>%
  dplyr::mutate(CHROM = as.numeric(gsub(".*Chr|_\\d+", "", CHROM))) %>%
  dplyr::select("ID", "POS", "CHROM")
fwrite(snp_pos, snp_file, sep = ",", row.names = F, col.names = F, quote = F)
