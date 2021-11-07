#!/usr/bin/env Rscript

# mylib <- "/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.0"
# .libPaths(mylib)
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

source("script/geno_functions.R")

args <- commandArgs(trailingOnly = TRUE)

test_file <- args[1] #snakemake@output$test_file
snp_file <- args[2] #snakemake@output$snp_file
# genotype_id <- args[3] #snakemake@params$genotype_id
vcf_file <- args[3] #snakemake@input
# cores <- args[5]



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
