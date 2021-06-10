vcfFile2geno <- function(vcf_file, geno_file = NULL, vcftools = "/usr/itrop/vcftools-0.1.16/bin/vcftools", threads = 20) {
  # convert vcf to 012 by vcftools
  temp_geno_file <- basename(tempfile())
  convert_cmd <- glue::glue("{vcftools} --vcf {test_file} --012 --out {temp_geno_file}") # vcf_file
  system(convert_cmd)
  # transpose and remove spaces 
  message("transposing geno table")
  geno <- data.table::fread(paste0(temp_geno_file, ".012"), nThread = threads, header = F, drop = 1)
  geno <- t(geno)
  message("Converting -1 to 9")
  geno[geno == -1] <- 9
  
  ind_names <- fread(paste0(temp_geno_file, ".012.indv"), nThread = 20, header = F)
  colnames(geno) <- ind_names$V1
  loci <- fread(paste0(temp_geno_file, ".012.pos"), nThread = 20, header = F)
  rownames(geno) <- paste(loci$V1, loci$V2, sep = "_")
  
  if (!is.null(geno_file)) {
    message("Removing spacing between columns")
    data.table::fwrite(geno, temp_geno_file, sep = " ", row.names = F, col.names = F, nThread = threads)
    sed_cmd <- glue::glue("sed -r 's/\\s+//g' {temp_geno_file} > {geno_file}")
    system(sed_cmd)
  }
  
  message("Removing temp files")
  temp_files <- list.files(getwd(), pattern = temp_geno_file, full.names = T)
  unlink(temp_files)
  
  return(geno)
}
# vcfFile2geno(test_file, "test_fun.geno")



write_elai_geno <- function(genotypes, file) {
  ## genotypes is a matrix/dataframe, with inds in columns and SNPs in rows, values are in ATGC format
  file.create(file)
  write(ncol(genotypes), file)
  write(nrow(genotypes), file, append = T)
  genotypes <- data.frame("ID" = rownames(genotypes), genotypes)
  data.table::fwrite(genotypes, file, append = T, sep = ",", quote = F, col.names = T, row.names = F)
}
# write_elai_geno(t(sml_genotypes), "test.gen")


geno2elai_gt <- function(geno, snps_vcf, snp_id = NULL) {
  # geno is a matrix/df in geno format, snps in rows, indv in columns
  stopifnot(all(geno %in% c(0:2, 9)))
  gt <- geno
  if (is.vector(geno)) {
    ref <- snps_vcf$REF[snps_vcf$ID == snp_id]
    alt <- snps_vcf$ALT[snps_vcf$ID == snp_id]
    for (i in 1:length(geno)) {
      if (geno[i] == 0) {
        gt[i] <- paste0(ref, ref)
      } else if (geno[i] == 1) {
        gt[i] <- paste0(ref, alt)
      } else if (geno[i] == 2) {
        gt[i] <- paste0(alt, alt)
      } else if (geno[i] == 9) {
        gt[i] <- "??"
      }
    }
  } else if (is.data.frame(geno) || is.matrix(geno) || is.data.table(geno)) {
    for (i in 1:nrow(geno)) {
      ref <- snps_vcf$REF[snps_vcf$ID == rownames(geno)[i]]
      alt <- snps_vcf$ALT[snps_vcf$ID == rownames(geno)[i]]
      for (j in 1:ncol(geno)) {
        if (geno[i,j] == 0) {
          gt[i,j] <- paste0(ref, ref)
        } else if (geno[i,j] == 1) {
          gt[i,j] <- paste0(ref, alt)
        } else if (geno[i,j] == 2) {
          gt[i,j] <- paste0(alt, alt)
        } else if (geno[i,j] == 9) {
          gt[i,j] <- "??"
        }
      }
    }
  }
  return(gt)
}


# /data3/projects/vietcaf/baotram/scripts/robusta_vn/elai/elai-lin -g E_R.geno -p 10 -g A_G.geno -p 11 -g VN.geno -p 1 -pos snp.pos -s 20 -C 2 -c 10 -o vn_var -mixgen 30 --exclude-nopos --exclude-miss1
elai <- function(source_files, test_file, snp_file, outdir, 
                 lower_layers, gens, options = "-s 20 --exclude-nopos --exclude-miss1", 
                 elai_src = "/data3/projects/vietcaf/baotram/scripts/robusta_vn/elai/elai-lin",
                 rerun = FALSE) {
  
  
  dir.create(outdir, recursive = T, showWarnings = F)
  # outdir <- normalizePath(outdir)
  adm_file <- file.path(outdir, "elai_r.admix.txt")
  dosage_file <- file.path(outdir, "elai_r.ps21.txt")
  snpout_file <- file.path(outdir, "elai_r.snpinfo.txt")
  log_file <- file.path(outdir, "elai_r.log.txt")
  
  
  src_opt <- do.call(paste, lapply(source_files, function(s) {
    paste0("-g ", normalizePath(s), " -p 1", which(source_files == s)-1)
  }))
  
  test_opt <- paste("-g", normalizePath(test_file), "-p 1")
  
  snp_opt <- paste("-pos", normalizePath(snp_file))
  
  other_opt <- paste("-C", length(source_files), "-c", lower_layers, "-mixgen", gens, "-o elai_r", options)
  
  elai_cmd <- paste(elai_src, src_opt, test_opt, snp_opt, other_opt)
  
  temp_outdir <- tempfile()
  dir.create(temp_outdir)
  
  
  
  if (all(file.exists(adm_file, dosage_file, snpout_file, log_file)) && !rerun) {
    # old_cmd <- readLines(log_file, n = 3)
    # old_cmd <- old_cmd[grep("COMMAND", old_cmd)]
    # old_cmd <- stringr::str_trim(gsub("## COMMAND: ", "", old_cmd))
    # 
    # if (old_cmd == elai_cmd) {
    #   message("same elai cmd has runned, elai outputs exists ... read them")
    # } else {
    #   run_cmd <- paste("cd", temp_outdir, "&&", elai_cmd)
    #   print(elai_cmd)
    #   
    #   # temp_outdir <- tempfile()
    #   # dir.create(temp_outdir)
    #   # setwd(temp_outdir)
    #   system(run_cmd)
    #   
    #   file.copy(from = list.files(file.path(temp_outdir, "output"), full.names = T), to = outdir, recursive = T)
    #   
    #   unlink(temp_outdir, recursive = T)
    # }
    message("skip re-running elai")
  } else {
    run_cmd <- paste("cd", temp_outdir, "&&", elai_cmd)
    print(elai_cmd)
    system(run_cmd)
    
    file.copy(from = list.files(file.path(temp_outdir, "output"), full.names = T), to = outdir, recursive = T)
    
    unlink(temp_outdir, recursive = T)
  }
  
  
  ## todo
  ## read overall admixture
  admixture <- data.table::fread(adm_file, header = F, strip.white = T, sep = " ", stringsAsFactors = F, data.table = F)
  test_info <- read.csv(test_file, skip = 2, nrows = 5, header = T)
  ind_names <- colnames(test_info)[-1]
  rownames(admixture) <- ind_names
  log_cmd <- readLines(log_file, n = 3)
  log_cmd <- log_cmd[grep("COMMAND", log_cmd)]
  # log_cmd <- stringr::str_trim(gsub("## COMMAND: ", "", log_cmd))
  log_cmd <- unlist(strsplit(log_cmd, " "))
  sources <- log_cmd[grep("-g", log_cmd)+1]
  sources <- sources[sapply(sources, function(s) grepl(basename(s), paste(source_files, collapse = " ")))] 
  ancestries <- basename(sources)
  colnames(admixture) <- ancestries
  
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

read_elai <- function(elai_outdir, elai_indir) {
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
  test_file <- basename(test_file[length(test_file)])
  test_file <- file.path(elai_indir, test_file)
  stopifnot(file.exists(test_file))
  
  sources <- stringr::str_split(log_cmd, " -p 1\\d+ ", simplify = T)
  sources <- sources[!grepl("-p 1", sources)]
  sources <- gsub(".*-g ", "", sources)
  sources <- file.path(elai_indir, basename(sources))
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
