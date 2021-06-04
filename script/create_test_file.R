
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages('dplyr', repos = "https://cloud.r-project.org")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages('parallel', repos = "https://cloud.r-project.org")
if (!requireNamespace("vcfR", quietly = TRUE)) install.packages('vcfR', repos = "https://cloud.r-project.org")
if (!requireNamespace("adegenet", quietly = TRUE)) install.packages('adegenet', repos = "https://cloud.r-project.org")

library(dplyr)
library(parallel)
library(vcfR)

source("/data3/projects/vietcaf/baotram/scripts/robusta_vn/R/geno_functions.R", local = knitr::knit_global())


test_file <- snakemake@output$test_file
snp_file <- snakemake@output$snp_file
genotype_id <- snakemake@params$genotype_id
vcf_file <- snakemake@input




vcf_R <- read.vcfR(vcf_file)
chr_gl <- vcfR2genlight(vcf_R, n.cores = cores)
chr_snps <- <- as.data.frame(getFIX(vcf_R), stringsAsFactors = F)
chr_snps$ID <- paste(chr_snps$CHROM, chr_snps$POS, sep = "_")
chr_geno <- as.matrix(chr_gl)
chr_geno[is.na(chr_geno)] <- 9

test_geno <- t(chr_geno[genotype_id,])
test_gt <- geno2elai_gt(test_geno, chr_snps)
rownames(test_gt) <- genotype_id
colnames(test_gt) <- chr_snps$ID
write_elai_geno(test_gt, out_file) 



snp_pos <- chr_snps %>% dplyr::mutate(CHROM = as.numeric(gsub(".*Chr|_\\d+", "", CHROM))) %>% dplyr::select("ID", "POS", "CHROM")
fwrite(snp_pos, snp_file, sep = ",", row.names = F, col.names = F, quote = F)

