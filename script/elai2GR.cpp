// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

/*** R
 dosage2MT <- function(dosage, p, matrix_smp = NULL) {
   df <- dosage %>%
     dplyr::filter(pos == p) %>%
     tibble::column_to_rownames("ancestry") %>%
     dplyr::select(even)
  df <- t(df)
  if (!is.null(matrix_smp) && setequal(colnames(matrix_smp), colnames(df))) {
    # coln <- colnames(matrix_smp)[colnames(matrix_smp) %in% colnames(df)]
    df <- subset(df, select = colnames(matrix_smp))
  }
  return(df)
}

filter_ind <- function(dosage, ind) {
  dosage_ind <- dosage %>% dplyr::filter(individual == ind)
  return(dosage_ind)
}

filter_chr <- function(dosage, chr) {
  dosage_chr <- dosage %>% dplyr::filter(chr == chr)
  return(dosage_chr)
}

get_chr_name <- function(rs) {
  chr <- gsub("_\\d+", "", rs[1])
  return(chr)
}

combine_GR <- function(...) {
  gr <- c(...)
  return(gr)
}

dosage2GR_R <- function(fin_dosage) {
  elai_GR_list <- GRanges()
  for (id in unique(fin_dosage$individual)) {
    fin_dosage_ind <- fin_dosage %>% filter(individual == id)
    elai_GR_ind <- GRanges()
    for (chr in sort(unique(fin_dosage_ind$chr))) {
      elai_GR_chr <- GRanges()
      fin_dosage_ind_chr <- fin_dosage_ind %>% filter(chr == chr)
      chr_name <- gsub("_\\d+", "", fin_dosage_ind_chr$rs[1])
      R_start <- sort(unique(fin_dosage_ind_chr$pos))[1]
      R_stop <- R_start
      adm_df <- fin_dosage_ind_chr %>%
        filter(pos == R_start) %>% 
        column_to_rownames("ancestry") %>% 
        dplyr::select(even)
      adm_df <- t(adm_df)
      for (p in sort(unique(fin_dosage_ind_chr$pos))) {
        if (p == R_start) {
          next
        } else if (p > R_stop + 10e3) {
          elai_GR_pos <- GRanges(seqnames = chr_name, strand = "+", ranges = IRanges(R_start, R_stop), adm_df, individual = id)
          elai_GR_chr <- c(elai_GR_chr, elai_GR_pos)
          R_start <- p
        } else if (p > R_stop && p < R_stop + 10e3) {
          adm_df_p <- fin_dosage_ind_chr %>%
            filter(pos == p) %>% 
            column_to_rownames("ancestry") %>% 
            dplyr::select(even)
          adm_df_p <- t(adm_df_p)
          adm_df_p <- adm_df_p[,colnames(adm_df)]
          if (all(adm_df_p == adm_df)) {
            R_stop <- p
          } else {
            elai_GR_pos <- GRanges(seqnames = chr_name, strand = "+", ranges = IRanges(R_start, R_stop), adm_df, individual = id)
            elai_GR_chr <- c(elai_GR_chr, elai_GR_pos)
            R_start <- p
            adm_df <- adm_df_p
          }
        }
        if (R_stop == max(unique(fin_dosage_ind_chr$pos))) {
          elai_GR_pos <- GRanges(seqnames = chr_name, strand = "+", ranges = IRanges(R_start, R_stop), adm_df, individual = id)
          elai_GR_chr <- c(elai_GR_chr, elai_GR_pos)
        }
      }
      elai_GR_ind <- c(elai_GR_ind, elai_GR_chr)
    }
    elai_GR_list <- c(elai_GR_list, elai_GR_ind)
  }
}
*/

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
RObject dosage2GR(DataFrame fin_dosage, int gap_length) {
  // should do: aplly dplyr functions directly to subset the dataframe;
  // dplyr, GenomicRanges, tibble
  Environment dplyr_pkg = Environment::namespace_env("dplyr");
  Function dplyr_filter = dplyr_pkg["filter"];
  Environment granges_pkg = Environment::namespace_env("GenomicRanges");
  Function GRanges = granges_pkg["GRanges"];
  Environment iranges_pkg = Environment::namespace_env("IRanges");
  Function IRanges = iranges_pkg["IRanges"];
  Function dosage2MT("dosage2MT");
  Function filter_ind("filter_ind");
  Function filter_chr("filter_chr");
  Function get_chr_name("get_chr_name");
  Function combine_GR("combine_GR");
  
  CharacterVector ind_raw_list = fin_dosage["individual"];
  CharacterVector ind_list = unique(ind_raw_list);
  RObject elai_GR_list = GRanges();
  

  for (auto& id : ind_list) {
  // parallelForEach(ind_list, [] (auto& id){
    String ind(id);
    Rcout << "individual: " << id << endl;
    DataFrame fin_dosage_ind = filter_ind(fin_dosage, ind);
    RObject elai_GR_ind = GRanges();
    CharacterVector chr_raw_list = fin_dosage["chr"];
    CharacterVector chr_list = sort_unique(chr_raw_list);
    for (auto& chr : chr_list) {
      Rcout << "chromosome: " << chr << endl;
      RObject elai_GR_chr = GRanges();
      DataFrame fin_dosage_ind_chr = filter_chr(fin_dosage_ind, chr);
      String chr_name = get_chr_name(fin_dosage_ind_chr["rs"]);
      IntegerVector pos_raw_list = fin_dosage_ind_chr["pos"];
      IntegerVector pos_list = sort_unique(pos_raw_list);
      int R_start = pos_list[0];
      int R_stop = pos_list[0];
      NumericMatrix adm_df = dosage2MT(fin_dosage_ind_chr, R_start);
      CharacterVector ancestry_raw_list = fin_dosage_ind_chr["ancestry"];
      CharacterVector ancestry_list = sort_unique(ancestry_raw_list);
      Progress progress(pos_list.size());
      // Rcout << pos_list << endl;
      for (auto& p : pos_list) { // first position
        // Rcout << p << endl;
        if (p == R_start) continue;
        else if (p == max(pos_list)) { // last position
          // Rcout << "end" << endl;
          R_stop = p;
          RObject ranges = IRanges(R_start, R_stop);
          DataFrame individual = DataFrame::create(Named("individual") = ind);
          RObject elai_GR_pos = GRanges(chr_name, ranges, "+", adm_df, individual);
          elai_GR_chr = combine_GR(elai_GR_chr, elai_GR_pos);
        } else if (p > (R_stop + gap_length)) { // gap with previous R_stop -> generate GR for previous ranges -> renew R_start, R_stop, adm_df
          // Rcout << "new range because of gap" << endl;
          RObject ranges = IRanges(R_start, R_stop);
          DataFrame individual = DataFrame::create(Named("individual") = ind);
          RObject elai_GR_pos = GRanges(chr_name, ranges, "+", adm_df, individual);
          elai_GR_chr = combine_GR(elai_GR_chr,elai_GR_pos);
          R_start = p;
          R_stop = p;
          adm_df = dosage2MT(fin_dosage_ind_chr, R_start);
        } else { // next position, no gap -> check if there is a switch in admixture
          // Rcout << "same range" << endl;
          NumericMatrix adm_df_p = dosage2MT(fin_dosage_ind_chr, p, adm_df);
          // Rcout << adm_df_p << endl;
          NumericVector adm = adm_df(0, _);
          // Rcout << "adm: " << adm << endl;
          NumericVector adm_p = adm_df_p(0, _);
          // Rcout << "adm_p: " << adm_p << endl;
          // bool check = is_true(all(adm == adm_p));
          // Rcout << check << endl;
          if (adm_p.size() < ancestry_list.size()) { // if no switch -> update R_stop
            continue;
          } else if ((adm.size() == adm_p.size()) & is_true(all(adm == adm_p))) {
            R_stop = p;
          } else { // if switch -> the same as gap
            // Rcout << "new range" << endl;
            RObject ranges = IRanges(R_start, R_stop);
            DataFrame individual = DataFrame::create(Named("individual") = ind);
            RObject elai_GR_pos = GRanges(chr_name, ranges, "+", adm_df, individual);
            elai_GR_chr = combine_GR(elai_GR_chr, elai_GR_pos);
            R_start = p;
            R_stop = p;
            adm_df = dosage2MT(fin_dosage_ind_chr, R_start);
            adm = adm_df(0, _);
            // Rcout << "new adm: " << adm << endl;
          }
        }
        progress.increment(); 
        // Rcout << p << " done" << endl;
      }
      elai_GR_ind = combine_GR(elai_GR_ind, elai_GR_chr);
      // Rcout << "done" << endl;
    }
    elai_GR_list = combine_GR(elai_GR_list, elai_GR_ind);
    // LogicalVector ind_row = ind_list == id;
    // IntegerVector row_id = seq(0, fin_dosage.nrow()-1);
    // IntegerVector ind_row_id = row_id[ind_row];
    // NumericMatrix::Sub fin_dosage_ind = fin_dosage(ind_row_id, _);
  };
  return elai_GR_list;
};


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


/*** R
# test_chr7_gr <- elai_GR(dose_chr7 %>% filter(pos > 4.0e7))
*/