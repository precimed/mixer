/*
  bgmg - tool to calculate log likelihood of BGMG and UGMG mixture models
  Copyright (C) 2018 Oleksandr Frei 

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <stdint.h>

#include <unordered_map>
#include <memory>

#include <iostream>
#include <sstream>
#include <vector>
#include <valarray>

#include "boost/utility.hpp"

#include "bgmg_log.h"
#include "bgmg_parse.h"
#include "ld_matrix_csr.h"

enum CostCalculator {
  CostCalculator_Sampling = 0,
  CostCalculator_Gaussian = 1,
  CostCalculator_Convolve = 2,
}; 

// Singleton class to manage a collection of objects, identifiable with some integer ID.
template<class Type>
class TemplateManager : boost::noncopyable {
 public:
  static TemplateManager<Type>& singleton() {
    static TemplateManager<Type> manager;
    return manager;
  }

  std::shared_ptr<Type> Get(int id) {
    if (map_.find(id) == map_.end()) {
      LOG << " Create new context (id=" << id << ")";
      map_[id] = std::make_shared<Type>();
    }
    return map_[id];
  }

  void Erase(int id) {
    LOG << " Dispose context (id=" << id << ")";
    map_.erase(id);
  }

  void Clear() {
    LOG << " Dispose all context ids";
    map_.clear();
  }

 private:
  TemplateManager() { }  // Singleton (make constructor private)
  std::unordered_map<int, std::shared_ptr<Type>> map_;
};

template<typename T>
class DenseMatrix {
 public:
  explicit DenseMatrix(size_t no_rows = 0, size_t no_columns = 0, bool store_by_rows = true)
    : no_rows_(no_rows),
    no_columns_(no_columns),
    store_by_rows_(store_by_rows),
    data_(nullptr) {
    if (no_rows > 0 && no_columns > 0) {
      data_ = new T[no_rows_ * no_columns_];
    }
  }

  DenseMatrix(const DenseMatrix<T>& src_matrix) {
    no_rows_ = src_matrix.no_rows();
    no_columns_ = src_matrix.no_columns();
    store_by_rows_ = src_matrix.store_by_rows_;
    if (no_columns_ >0 && no_rows_ > 0) {
      data_ = new T[no_rows_ * no_columns_];

      for (size_t i = 0; i < no_rows_ * no_columns_; ++i) {
        data_[i] = src_matrix.get_data()[i];
      }
    } else {
      data_ = nullptr;
    }
  }

  virtual ~DenseMatrix() {
    delete[] data_;
  }

  void InitializeZeros() {
    memset(data_, 0, sizeof(T)* no_rows_ * no_columns_);
  }

  T& operator() (size_t index_row, size_t index_col) {
    if (store_by_rows_) {
      return data_[index_row * no_columns_ + index_col];
    }
    return data_[index_col * no_rows_ + index_row];
  }

  const T& operator() (size_t index_row, size_t index_col) const {
    assert(index_row < no_rows_);
    assert(index_col < no_columns_);
    if (store_by_rows_) {
      return data_[index_row * no_columns_ + index_col];
    }
    return data_[index_col * no_rows_ + index_row];
  }

  DenseMatrix<T>& operator= (const DenseMatrix<T>& src_matrix) {
    no_rows_ = src_matrix.no_rows();
    no_columns_ = src_matrix.no_columns();
    store_by_rows_ = src_matrix.store_by_rows_;
    if (data_ != nullptr) {
      delete[] data_;
    }
    if (no_columns_ >0 && no_rows_ > 0) {
      data_ = new  T[no_rows_ * no_columns_];

      for (size_t i = 0; i < no_rows_ * no_columns_; ++i) {
        data_[i] = src_matrix.get_data()[i];
      }
    } else {
      data_ = nullptr;
    }

    return *this;
  }

  size_t no_rows() const { return no_rows_; }
  size_t no_columns() const { return no_columns_; }
  size_t size() const { return no_rows_ * no_columns_; }
  bool is_equal_size(const DenseMatrix<T>& rhs) const {
    return no_rows_ == rhs.no_rows_ && no_columns_ == rhs.no_columns_;
  }

  std::string to_str() {
    int rows_to_str = std::min<int>(5, no_rows_ - 1);
    int cols_to_str = std::min<int>(5, no_columns_ - 1);
    std::stringstream ss;
    ss << "[";
    for (int i = 0; i < rows_to_str; i++) {
      bool last_row = (i == (rows_to_str - 1));
      for (int j = 0; j < cols_to_str; j++) {
        bool last_col = (j == (cols_to_str - 1));
        ss << (*this)(i, j);
        if (last_col) ss << ", ...";
        else ss << ", ";
      }
      if (last_row) ss << "; ...";
      else ss << "; ";
    }
    ss << "]";

    size_t nnz = 0;
    for (int i = 0; i < no_rows_; i++)
      for (int j = 0; j < no_columns_; j++)
        if ((*this)(i, j) != 0) nnz++;
    ss << ", nnz=" << nnz;

    return ss.str();
  }


  T* get_data() {
    return data_;
  }

  const T* get_data() const {
    return data_;
  }

 private:
  size_t no_rows_;
  size_t no_columns_;
  bool store_by_rows_;
  T* data_;
};

struct LoglikeCacheElem {
 public:
  LoglikeCacheElem() {}
  LoglikeCacheElem(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float  rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, double cost);
  void get(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float* rho_beta, int sig2_zero_len, float* sig2_zero, float* rho_zero, double* cost);

 private:
  float pi_vec1_;
  float pi_vec2_;
  float pi_vec3_;
  float sig2_beta1_;
  float sig2_beta2_;
  float sig2_zero1_;
  float sig2_zero2_;
  float rho_beta_;
  float rho_zero_;
  double cost_;
};

// History of log likelihood calculations
class LoglikeCache {
public:
 void    add_entry(float pi_vec, float sig2_zero, float sig2_beta, double cost);
 int64_t get_entry(int entry_index, float* pi_vec, float* sig2_zero, float* sig2_beta, double* cost);
 void    add_entry(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float  rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, double cost);
 int64_t get_entry(int entry_index, int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float* rho_beta, int sig2_zero_len, float* sig2_zero, float* rho_zero, double* cost);

 void clear() { cache_.clear(); }
 int num_entries() { return cache_.size(); }
private:
  std::vector<LoglikeCacheElem> cache_;
};

class BgmgCalculator : public TagToSnpMapping {
 public:
  BgmgCalculator();
  virtual ~BgmgCalculator() {}

  int64_t init(std::string bim_file, std::string frq_file, std::string chr_labels, std::string trait1_file, std::string trait2_file, std::string exclude, std::string extract);
  int64_t convert_plink_ld(std::string plink_ld_gz, std::string plink_ld_bin);  // require init() to be called first, e.i. doesn't work after set_tag_indices.

  // num_snp = total size of the reference (e.i. the total number of genotyped variants)
  // num_tag = number of tag variants to include in the inference (must be a subset of the reference)
  // indices = array of size num_tag, containing indices from 0 to num_snp-1
  // NB: all tag variants must have defined zvec, nvec, mafvec and weights.
  int64_t set_tag_indices(int num_snp, int num_tag, int* tag_indices);
  
  int64_t set_chrnumvec(int num_snp, const int* chrlabel);
  int64_t retrieve_chrnumvec(int length, int* buffer);

  int64_t set_snp_order(int component_id, int64_t length, const int* buffer);
  int64_t retrieve_snp_order(int component_id, int64_t length, int* buffer);
  int64_t retrieve_k_pdf(int length, double* buffer);

  // consume input in plink format, e.i. lower triangular LD r2 matrix
  // - snp_index is less than tag_index;
  // - does not contain r2=1.0 of the variant with itself
  // must be called after set_tag_indices
  // must be called one for each chromosome, sequentially, starting from lowest chromosome number
  // non-tag variants will be ignored
  int64_t set_ld_r2_coo(int chr_label, int64_t length, int* snp_index, int* tag_index, float* r);
  int64_t set_ld_r2_coo(int chr_label, const std::string& filename);
  int64_t set_ld_r2_csr(int chr_label = -1);  // finalize

  int64_t num_ld_r2_snp(int snp_index);
  int64_t retrieve_ld_r2_snp(int snp_index, int length, int* tag_index, float* r2);
  int64_t num_ld_r2_chr(int chr_label);
  int64_t retrieve_ld_r2_chr(int chr_label, int64_t length, int* snp_index, int* tag_index, float* r2);
  int64_t num_ld_r2_snp_range(int snp_index_from, int snp_index_to);
  int64_t retrieve_ld_r2_snp_range(int snp_index_from, int snp_index_to, int length, int* snp_index, int* tag_index, float* r2);

  // must be called after set_ld_r2, as it adjusts r2 matrix
  // one value for each snp (tag and non-tag)
  int64_t set_mafvec(int length, float* values);
  
  // zvec, nvec, weights for tag variants
  // all values must be defined
  int64_t set_zvec(int trait, int length, float* values);
  int64_t set_nvec(int trait, int length, float* values);
  int64_t set_causalbetavec(int trait, int length, float* values);
  int64_t set_weights(int length, float* values);
  int64_t set_weights_randprune(int n, float r2);   // alternative to set_weights; calculates weights based on random pruning from LD matrix

  int64_t retrieve_zvec(int trait, int length, float* buffer);
  int64_t retrieve_nvec(int trait, int length, float* buffer);
  int64_t retrieve_causalbetavec(int trait, int length, float* buffer);
  int64_t retrieve_mafvec(int length, float* buffer);
  int64_t retrieve_weights(int length, float* buffer);
  int64_t retrieve_tag_indices(int num_tag, int* tag_indices);

  int64_t find_snp_order();  // private - only for testing
  int64_t find_tag_r2sum(int component_id, float num_causals);  // private - only for testing
  void find_tag_r2sum_no_cache(int component_id, float num_causal, int k_index, std::vector<float>* buffer); // private - only for testing

  int64_t set_option(char* option, double value);
  
  void clear_state();
  void clear_tag_r2sum(int component_id);

  int64_t retrieve_tag_r2_sum(int component_id, float num_causal, int length, float* buffer);
  int64_t retrieve_ld_tag_r2_sum(int length, float* buffer);
  int64_t retrieve_ld_tag_r4_sum(int length, float* buffer);
  int64_t retrieve_fixed_effect_delta(int trait_index, int length, float* delta);
  
  // Calculate and retrieve weighted_causal_r2, wcr2[i], i runs from 0 to num_snp
  // wcr2[i] = sum_j w_j r^2_ij, where the sum runs across all tag variants, and 
  // w_j is random-pruning-based weight of j-th tag variant.
  // This function may help us to make a better selection of informative tag variants.
  int64_t retrieve_weighted_causal_r2(int length, float* buffer);

  double calc_univariate_cost(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);
  double calc_univariate_cost_cache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);
  double calc_univariate_cost_nocache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);        // default precision (see FLOAT_TYPE in bgmg_calculator.cc)
  double calc_univariate_cost_nocache_float(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);  // for testing single vs double precision
  double calc_univariate_cost_nocache_double(int trait_index, float pi_vec, float sig2_zero, float sig2_beta); // for testing single vs double precision
  int64_t calc_univariate_pdf(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* zvec, float* pdf);
  int64_t calc_univariate_power(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, float zthresh, int length, float* nvec, float* svec);
  int64_t calc_univariate_delta_posterior(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* c0, float* c1, float* c2);
  int64_t calc_bivariate_delta_posterior(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero,
                                         int length, float* c00, float* c10, float* c01, float* c20, float* c11, float* c02);

  double calc_bivariate_cost(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);
  double calc_bivariate_cost_nocache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);
  double calc_bivariate_cost_cache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);
  int64_t calc_bivariate_pdf(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, int length, float* zvec1, float* zvec2, float* pdf);
  void log_diagnostics();
 
  double calc_unified_univariate_cost(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* Ezvec2);
  double calc_unified_univariate_cost_gaussian(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* Ezvec2);
  double calc_unified_univariate_cost_convolve(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL);
  int64_t calc_unified_univariate_pdf(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, int length, float* zvec, float* pdf);
  int64_t calc_unified_univariate_power(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float zthresh, int length, float* nvec, float* svec);
  int64_t calc_unified_univariate_delta_posterior(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, int length, float* c0, float* c1, float* c2);

  int64_t seed() { return seed_; }
  void set_seed(int64_t seed) { seed_ = seed; }

  int max_causals() const { return max_causals_; }
  int k_max() const { return k_max_; }

  virtual int num_snp() { return num_snp_; }
  virtual int num_tag() { return num_tag_; }
  virtual const std::vector<int>& tag_to_snp() { return  tag_to_snp_; }
  virtual const std::vector<int>& snp_to_tag() { return snp_to_tag_; }
  virtual const std::vector<char>& is_tag() { return is_tag_; }
  virtual const std::vector<int>& chrnumvec() { return chrnumvec_; }
  virtual const std::vector<float>& mafvec() { return mafvec_; }

  int64_t clear_loglike_cache() { loglike_cache_.clear(); return 0; }
  int64_t get_loglike_cache_size() { return loglike_cache_.num_entries(); }
  int64_t get_loglike_cache_univariate_entry(int entry_index, float* pi_vec, float* sig2_zero, float* sig2_beta, double* cost) {
    return loglike_cache_.get_entry(entry_index, pi_vec, sig2_zero, sig2_beta, cost);
  }
  int64_t get_loglike_cache_bivariate_entry(int entry_index, int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float* rho_beta, int sig2_zero_len, float* sig2_zero, float* rho_zero, double* cost) {
    return loglike_cache_.get_entry(entry_index, pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, cost);
  }

 private:
  int num_snp_;
  int num_tag_;
  std::vector<int> tag_to_snp_; // 0..num_snp_-1, size=num_tag_
  std::vector<int> snp_to_tag_; // 0..num_tag-1,  size=num_snp_, -1 indicate non-tag snp
  std::vector<char> is_tag_;    // true or false, size=num_snp_, is_tag_[i] == (snp_to_tag_[i] != -1)
  std::vector<int> chrnumvec_;  // vector of chromosome labels, one per snp

  LdMatrixCsr ld_matrix_csr_;

  // all stored for for tag variants (only)
  std::vector<float> zvec1_;
  std::vector<float> nvec1_;
  std::vector<float> zvec2_;
  std::vector<float> nvec2_;
  std::vector<float> weights_;
  std::vector<float> mafvec_;

  std::vector<float> causalbetavec1_;  // assumed causal betas, added as fixed effects component in the model
  std::vector<float> causalbetavec2_;  // e.g. delta_j = \sqrt N_j \sum_i \sqrt H_i r_ij \beta_i   <- here "beta_i" is causalbetavec
  
  std::vector<float>* get_zvec(int trait_index);
  std::vector<float>* get_nvec(int trait_index);
  std::vector<float>* get_causalbetavec(int trait_index);

  LoglikeCache loglike_cache_;

  // vectors with one value for each component in the mixture
  // snp_order_ gives the order of how SNPs are considered to be causal 
  // tag_r2_sum_ gives cumulated r2 across causal SNPs, according to snp_order, where last_num_causals_ define the actual number of causal variants.
  std::vector<std::shared_ptr<DenseMatrix<int>>> snp_order_;  // permutation matrix; #rows = pimax*num_snp; #cols=k_max_
  std::vector<std::shared_ptr<DenseMatrix<float>>> tag_r2sum_;
  std::vector<float>                               last_num_causals_;

  // options, and what do they affect
  int k_max_;
  int max_causals_;
  int num_components_;
  int64_t seed_;
  bool use_complete_tag_indices_;  // an option that indicates that all SNPs are TAG (i.e. num_snp_ == num_tag_).
  float r2_min_;
  float z1max_;
  float z2max_;
  float max_chisq_trait2_;
  CostCalculator cost_calculator_;
  bool cache_tag_r2sum_;
  double cubature_abs_error_;
  double cubature_rel_error_;
  int cubature_max_evals_;
  int ld_format_version_;      // overwrite format version for LD matrix files. Default -1. Set this to 0 to read from MiXeR v1.0 LD files.
  std::vector<double> k_pdf_;  // the log-likelihood cost calculated independently for each of 0...k_max-1 selections of causal variants.            
  bool calc_k_pdf_;            // a calc_fixed_effect_delta_from_causalbetavecflag indicating whether we should calculate k_pdf_
  void check_num_snp(int length);
  void check_num_tag(int length);
  double calc_univariate_cost_fast(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);
  double calc_bivariate_cost_fast(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);
  double calc_univariate_cost_convolve(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);
  double calc_bivariate_cost_convolve(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);
  void calc_fixed_effect_delta_from_causalbetavec(int trait_index, std::valarray<float>* delta);
  void find_unified_tag_delta_sampling(int num_components, float* pi_vec, float* sig2_vec, float sig2_zeroC, int k_index, const float* nvec, const float* hvec, std::vector<float>* tag_delta2);

  BimFile bim_file_;

  template<typename T>
  friend double calc_univariate_cost_nocache_template(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, BgmgCalculator& rhs);
};

typedef TemplateManager<BgmgCalculator> BgmgCalculatorManager;
