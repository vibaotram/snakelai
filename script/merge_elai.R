#!/usr/bin/env Rscript


mylib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
.libPaths(mylib)

library(LEA)
suppressPackageStartupMessages(library(Rsamtools))
library(tidyverse)
library(viridis)
library(Rcpp)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

chr_nb <- args[1]

outdir <- "/shared/projects/elai_most/data/snakelai_out_all/elai3runs_45VN_merged"
dir.create(outdir, showWarnings = FALSE)

elai_outdir <- "/shared/projects/elai_most/data/snakelai_out_all/elai_average"
source_dir <- "/shared/projects/elai_most/data/snakelai_out_all/source"


# function
sourceCpp("/shared/projects/vietcaf/scripts/snakelai/script/elai2GR.cpp")

plot_elai <- function(snakedir, chrom, groups, plot_path = NULL) {
  snp <- c("all", "even")
  all_dosage <- do.call(rbind, lapply(snp, function(s) {
    result_data <- file.path(snakedir, chrom, s, "c5_mg20", "local_dosage.csv")
    print(result_data)
    all_dosage <- fread(result_data)
    all_dosage1 <- data.frame(all_dosage, snp = s)
    return(all_dosage1)
  }))
  print("done1")
  all_dosage$ancestry <- factor(all_dosage$ancestry)
  print("done2")
  levels(all_dosage$ancestry) <- groups
  print("done3")
  all_dosage$ancestry <- factor(all_dosage$ancestry, levels = c("ER", "OB", "C", "AG", "D"))
  print("done4")

  if (!is.null(plot_path)) {
    cat(plot_path)
    # tiff(file = plot_path,
    # width = 20, height = 12, units = "in", res = 800)
    p <- ggplot(all_dosage) +
      facet_grid(rows = vars(individual)) +
      # geom_point(aes(x = pos, y = dosage, color = ancestry), size = 1.5, shape = 16, alpha = .5) +
      geom_line(aes(x = pos, y = dosage, color = ancestry, linetype = snp, alpha = snp)) +
      scale_linetype_manual(values=c("dashed", "solid")) +
      # scale_color_viridis_d(end = 0.9) +
      scale_color_manual(values = group_col) +
      scale_alpha_manual(values = c(.8, 1)) +
      xlab(paste("Position on chromosome", unique(all_dosage$chr))) +
      ylab("Ancestral dosage") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            legend.key.size = unit(5, "mm"),
            axis.ticks.y = element_line(size = .1),
            axis.text.y = element_text(size = 8),
            panel.spacing = unit(2, "mm"),
            strip.text.x = element_text(angle = 90, size = 10, hjust = 1))
    ggsave(filename = plot_path, plot = p, width = 10, height = 6, dpi = 1000)
    # dev.off()
  }

  return(invisible(all_dosage))
}

rounding <- function(x) {
  n <- sapply(x, function(i) {
    if (i < 0.1) {
      i <- 0
    } else if (i > 0.4 && i < .6) {
      i <- 0.5
    } else if (i > 0.9) {
      i <- 1
    } else {
      i <- NA
    }
  }, USE.NAMES = F, simplify = T)

  return(n)
}

elai_GR <- function(all_dosage, gap = 1e6){
  # rounding admixture
  all_dosage <- all_dosage %>%
    mutate(dosage = rounding(dosage)) %>%
    arrange(snp, pos)
  com_dosage <- reshape2::dcast(all_dosage, ... ~ snp,
                                value.var = "dosage")
  # remove ambiguous positions
  com_dosage <- na.omit(com_dosage)
  # for each individual and each ancestry,
  # retain positions that have equal admixture
  fin_dosage <- com_dosage %>%
    group_by(individual, ancestry) %>%
    arrange(pos) %>%
    rowwise() %>%
    filter(all == even)
  elai_GR_list <- dosage2GR(fin_dosage = fin_dosage, gap_length = gap)
  return(elai_GR_list)
}



# list of individuals
vcf_file <- list.files(source_dir, "CC1.8.Chr07.recode.vcf", recursive = T, full.names = T)
vcf_header <- scanBcfHeader(vcf_file)
ind_names <- vcf_header[[1]]@listData$Sample
info_file <- "/shared/projects/elai_most/data/snakelai_out_all/config/sequence_info_final.tsv"
ind_info <- read.delim(info_file)
ind_names <- sapply(ind_names, function(i)
  ifelse(any(grepl(i, ind_info$Label)),
         grep(i, ind_info$Label, value = T),
         i),
  USE.NAMES = F)


# get elai results
chrom <- paste0("CC1.8.Chr", chr_nb)
chr_snmf <- load.snmfProject(file.path(source_dir, chrom, "snmf/chr.snmfProject"))
chr_bestrun <- which.min(cross.entropy(chr_snmf))
chr_prop <- Q(chr_snmf, run = chr_bestrun)
rownames(chr_prop) <- ind_names
chr_groups <- chr_prop %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  filter(!grepl("S.+", id)) %>%
  mutate(prior = gsub(".+_", "", id)) %>%
  mutate(prior = case_when(prior %in% c("O", "B") ~ "OB",
                           prior %in% c("E", "R") ~ "ER",
                           prior %in% c("A", "G") ~ "AG",
                           TRUE ~ prior)) %>%
  pivot_longer(cols = V1:V5, names_to = "ancestry", values_to = "proportion") %>%
  group_by(ancestry) %>%
  arrange(desc(proportion), .by_group = T) %>%
  slice_head(n = 3) %>% slice_tail(n = 1) %>% pull(prior) # n=2 for chr08, n=3 for chr11, n=1 for others
print(chrom)
print(chr_groups)

chr_dose <- plot_elai(snakedir = elai_outdir, chrom = chrom, groups = chr_groups)

chr_GR <- elai_GR(chr_dose)

saveRDS(chr_GR, file = file.path(outdir, paste0("chrom_", chr_nb, ".RDS")))
