read_elai <- function(elai_outdir) {
  adm_file <- file.path(elai_outdir, "elai_r.admix.txt")
  dosage_file <- file.path(elai_outdir, "elai_r.ps21.txt")
  snpout_file <- file.path(elai_outdir, "elai_r.snpinfo.txt")
  log_file <- file.path(elai_outdir, "elai_r.log.txt")

  log_cmd <- readLines(log_file, n = 3)
  log_cmd <- log_cmd[grep("COMMAND", log_cmd)]
  # log_cmd <- stringr::str_trim(gsub("## COMMAND: ", "", log_cmd))
  # log_cmd <- unlist(strsplit(log_cmd, " "))
  test_file <- stringr::str_split(log_cmd, " -p 1 ", simplify = T)[1]
  test_file <- stringr::str_split(test_file, " -g ", simplify = T)
  test_file <- test_file[length(test_file)]
  stopifnot(file.exists(test_file))

  sources <- stringr::str_split(log_cmd, " -p 1\\d+ ", simplify = T)
  sources <- sources[!grepl("-p 1", sources)]
  sources <- gsub(".*-g ", "", sources) # **
  stopifnot(all(file.exists(sources)))

  ancestries <- basename(sources)

  admixture <- data.table::fread(adm_file, header = F, strip.white = T, sep = " ", stringsAsFactors = F, data.table = F)
  colnames(admixture) <- ancestries
  test_info <- read.csv(test_file, skip = 2, nrows = 5, header = T)
  ind_names <- colnames(test_info)[-1]
  rownames(admixture) <- ind_names


  ## read ancestry dosage
  ## snp_out_file <- file.path(elai_outdir, "output/vn_var.snpinfo.txt")
  snp_out <- data.table::fread(snpout_file, header = T, stringsAsFactors = F, data.table = F)
  n_snp <- nrow(snp_out) # number of snps in snpinfo.txt
  n_ind <- ncol(test_info)-1 # number of test individuals in the test genotype file
  n_src <- length(sources) # number of source populations

  dosage <- scan(dosage_file)
  dim(dosage) <- c(n_src, n_snp, n_ind)

  # get dosage info for each individual in a list
  dosage <- sapply(ind_names, function(ind) {
    # get dosage matrix of the individual
    i <- which(ind_names == ind)
    ind_dosage <- as.data.frame(t(dosage[, , i]))
    # ind_dosage has [n. sources] columns and [n. snps] rows
    colnames(ind_dosage) <- ancestries
    ind_dosage$rs <- snp_out$rs
    # merge dosage table with snp info table
    ind_dosage <- merge(snp_out, ind_dosage, by = "rs")
    ind_dosage <- reshape2::melt(ind_dosage,
                                 measure.vars = ancestries,
                                 variable.name = "ancestry",
                                 value.name = "dosage")
    return(ind_dosage)
  }, USE.NAMES = T, simplify = F)

  return(list("oveall_admixture" = admixture, "ancestral_dosage" = dosage))
}


elai_dir = unique(dirname(as.character(snakemake@input)))

#elai_dir_lists <- list.files(elai_dir, "output_.+", full.names = T)
print(elai_dir)

elai_results <- lapply(elai_dir, function(d) {
  elai_res <- read_elai(d)
  adm <- data.frame(elai_res$oveall_admixture, "batch" = basename(dirname(d)), stringsAsFactors = F)
  dosage <- elai_res$ancestral_dosage
  return(list("oveall_admixture" = adm, "ancestral_dosage" = dosage))
  })

all_adm <- do.call(rbind, lapply(elai_results, function(x) {
  a <- x$oveall_admixture
  a$id <- rownames(a)
  return(a)
  }))
write.table(all_adm, as.character(snakemake@params$adm_file), sep = "\t", row.names = F, col.names = T)

all_dosage <- do.call(rbind, lapply(elai_results, function(x) {
  a <- x$ancestral_dosage
  d <- do.call(rbind, lapply(names(a), function(n) data.frame(a[[n]], "individual" = n)))
  return(d)
  }))
write.table(all_dosage, as.character(snakemake@params$dosage_file), sep = "\t", row.names = F, col.names = T)
