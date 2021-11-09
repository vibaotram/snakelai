#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
void simul_geno(NumericMatrix chr_freq, int n_inds, CharacterVector snp_id, DataFrame chr_snp, CharacterVector out_files, int core){
  Environment pkg = Environment::namespace_env("data.table");
  Function fread = pkg["fread"];
  Function geno2elai_gt("geno2elai_gt");
  
  #if defined(_OPENMP)
  #pragma omp parallel for num_threads(core)
  #endif
  for (int i = 0; i < chr_freq.ncol(); ++i){ // for each source file
    String grp_elai_input(out_files[i]);
    ofstream source_input;
    source_input.open(grp_elai_input);
    // write number of inds
    source_input << n_inds << endl;
    // write number of snps
    int n_snps = chr_freq.nrow()/3;
    source_input << n_snps << endl;
    CharacterVector ind_names(n_inds);
    CharacterVector source_names = colnames(chr_freq);
    // write header containing individual ids
    String header("ID");
    for (int j = 0; j < n_inds; ++j){ // each simulated genotype
      String src(source_names[i]);
      String s("group_");
      s += src;
      s += "-";
      s += j+1;
      ind_names[j] = s;
      header += ",";
      header += s;
      }
    // Rcout << header.get_cstring() << endl;
    
    source_input << header.get_cstring() << endl;
    NumericVector src_freq = chr_freq(_, i);
    for (int l = 0; l < n_snps; ++l){ // for each locus, simulate genotypes
      String loc = snp_id[l];
      NumericVector freq = src_freq[Range(3*l, 3*l + 2)];
      NumericVector sml_geno = sample(NumericVector {0, 1, 2}, n_inds, true, freq);
      CharacterVector sml_gt = geno2elai_gt(sml_geno, chr_snp, loc);
      // Rcout << loc.get_cstring() << " ";
      // Rcout << sml_gt << endl;
      String geno_ploc(loc);
      for (CharacterVector::iterator g = sml_gt.begin(); g != sml_gt.end(); ++g){
        geno_ploc += ",";
        geno_ploc += *g;
        }
      
      source_input << geno_ploc.get_cstring() << endl; // write genotypes to file
    }
    
    
    source_input.close();
  }
}

