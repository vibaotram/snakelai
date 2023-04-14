#!/usr/bin/env Rscript

mylib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
.libPaths(mylib)

library(data.table)
library(tidyverse)
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-e", "--elaidir"),
              type = "character",
              default = NULL,
              help = "Path of elai dir"))


average_elai <- function(files) {
  outfile <- gsub("elai_results_\\d", "elai_average", files[1])
  dir.create(dirname(outfile), recursive = T, showWarnings = F)
  
  dosage_list <- lapply(files, fread)
  dosage <- dosage_list[[1]]
  colnames(dosage)[colnames(dosage) == "dosage"] <- "dosage_1"
  dosage$dosage_2 <- dosage_list[[2]]$dosage
  dosage$dosage_3 <- dosage_list[[3]]$dosage
  dosage <- dosage %>%
    rowwise() %>%
    mutate(dosage = mean(c(dosage_1, dosage_2, dosage_3))) %>%
    select(!starts_with("dosage_"))

  fwrite(dosage, outfile,
         sep = ",", row.names = F, col.names = T,
         append = F)
}


myArgs <- parse_args(
  OptionParser(usage = "%prog [options]", option_list = option_list,
               description = "average elai results"))

ld1 <- myArgs$elaidir
fname <- gsub(".+/elai_results_1", "", ld1)
ldlist <- list.files(gsub("/elai_result.+", "", ld1), pattern = "elai_result", full.names = T)
ldlist <- paste0(ldlist, fname)
average_elai(ldlist)

# ev <- dirname(gsub("elai_output", "elai_average", er))
# average_elai(er, ev)

# snakedir <- "/shared/projects/elai_most/data/snakelai_out_all/elai_results_1"
#
# elai_raw <- list.files(snakedir, pattern = "local_dosage",
#                        recursive = T, full.names = T)
#
# elai_raw <- unique(dirname(dirname(dirname(elai_raw))))
#
#
#
# writeLines(elai_raw, "/shared/projects/elai_most/data/snakelai_out_all/local_dosage_dirs.txt")
