#!/usr/bin/env Rscript


mylib <- "/home/baotram/R/x86_64-pc-linux-gnu-library/4.0"

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages('dplyr', repos = "https://cloud.r-project.org", lib = mylib)
if (!requireNamespace("parallel", quietly = TRUE)) install.packages('parallel', repos = "https://cloud.r-project.org", lib = mylib)
if (!requireNamespace("vcfR", quietly = TRUE)) install.packages('vcfR', repos = "https://cloud.r-project.org", lib = mylib)
if (!requireNamespace("adegenet", quietly = TRUE)) install.packages('adegenet', repos = "https://cloud.r-project.org", dependencies = TRUE, lib = mylib)
if (!requireNamespace("optparse", quietly = TRUE)) install.packages('optparse', repos = "https://cloud.r-project.org", lib = mylib)
.libPaths(mylib)
suppressPackageStartupMessages(library("optparse"))
library(dplyr)
library(parallel)
library(vcfR)
library(adegenet)

source("/data3/projects/vietcaf/baotram/scripts/robusta_vn/R/geno_functions.R", local = knitr::knit_global())

option_list <- list(
  make_option(c("-o", "--output"),
              type = "character",
              help = "output files"),
  make_option(c("-i", "--input"),
              type = "character",
              help = "input files"),
  make_option(c("-g", "--nb_groups"),
              type = "integer",
              help = "number of groups"),
  make_option(c("-n", "--nb_genotypes"),
              type = "integer",
              help = "number of genotypes per group"),
  make_option(c("-t", "--threads"),
              type = "integer",
              help = "number of threads")
)

myArgs <- parse_args(
  OptionParser(usage = "%prog [-o output] [-i input] [-h nb_groups] [-n nb_genotypes] [-t threads]",
  option_list = option_list)
  )


out_files <- myArgs$output
outdir <- unique(dirname(out_files))
vcf_file <- myArgs$input
cores <- myArgs$threads
# nb_batches = nb_batches
nb_groups <- myArgs$nb_groups
nb_genotypes <- myArgs$nb_groups





vcf_R <- read.vcfR(vcf_file)
# small_vcfR <- vcf_R[sample(1:nrow(vcf_R), nb_snps)]



snmf_dir <- file.path(outdir, "snmf")
dir.create(snmf_dir, showWarnings = F)
geno_file <- file.path(snmf_dir, "chr.geno")



# get genotypes
# vcf_R@fix[,1] <- gsub(".", "_", vcf_R@fix[,1], fixed = T)
chr_gl <- vcfR2genlight(vcf_R, n.cores = cores)
chr_geno <- as.matrix(chr_gl)
chr_geno[is.na(chr_geno)] <- 9
# run snmf
write.geno(chr_geno, geno_file)
chr_snmf <- snmf(geno_file, K = nb_groups,
                 CPU = cores, repetitions = 10,
                 project = "new", entropy = T, seed = 987654321)

# chr1_snmf <- load.snmfProject(file.path(snmf_chr1_dir, "chr1.snmfProject"))


## generate genotypes based on ancestral genotypic freq (100 genotypes/group)


# get ancestral geno freq
chr_freq <- G(chr_snmf, K = nb_groups, run = best_run)
# colnames(chr_freq) <- c("C", "AG", "OB", "D", "ER")
chr_snps <- as.data.frame(getFIX(vcf_R), stringsAsFactors = F)
chr_snps$ID <- paste(chr_snps$CHROM, chr_snps$POS, sep = "_")
# chr_snps <- chr_snps[chr_snps$ID %in% chr1_gl@loc.names,]

# simulate 100 genotype each  group
n_inds <- nb_genotypes
n_snps <- nrow(chr_freq)/3

mclapply(1:ncol(chr_freq), function(v) {
  grp_elai_input <- out_files[v]
  ind_names <- paste0("group_", colnames(chr_freq)[v], "-", 1:n_inds)
  writeLines(c(n_inds, n_snps, paste(c("ID", ind_names), collapse = ",")), grp_elai_input)
  freq <- chr_freq[,v] # get allele freq of the group
  geno <- matrix(nrow=n_inds, ncol=n_snps) # genotype matrix of the group
  for (i in 1:n_snps) {
    # for locus i, frequencies of 0, 1, 2 genotype are in row 3i-2, 3i-1, 3i of freq matrix
    sml_geno <-  sample(c(0,1,2), n_inds, prob = c(freq[3*i-2], freq[3*i-1], freq[3*i]), replace = T)
    snp_id <- chr_gl@loc.names[i]
    sml_gt <- geno2elai_gt(sml_geno, chr_snps, snp_id)
    geno_li <- paste(c(snp_id, sml_gt), collapse = ",")
    cat(paste0(geno_li, "\n"), file = grp_elai_input, append = T)
  }
  return(v)
  }, mc.cores = 5, mc.preschedule = F)
