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

#include "bgmg_calculator.h"

#include <assert.h>

#include <random>
#include <limits>
#include <algorithm>
#include <vector>
#include <valarray>
#include <cmath>

#include "boost/throw_exception.hpp"

#include "bgmg_log.h"

#define OMP_CHUNK 1000

void BgmgCalculator::check_num_snp(int length) {
  if (num_snp_ == -1) BOOST_THROW_EXCEPTION(::std::runtime_error("call set_tag_indices first"));
  if (num_snp_ != length) BOOST_THROW_EXCEPTION(::std::runtime_error("length != num_snps_"));
}

void BgmgCalculator::check_num_tag(int length) {
  if (num_tag_ == -1) BOOST_THROW_EXCEPTION(::std::runtime_error("call set_tag_indices first"));
  if (num_tag_ != length) BOOST_THROW_EXCEPTION(::std::runtime_error("length != num_snps_"));
}

int64_t BgmgCalculator::set_zvec(int trait, int length, float* values) {
  LOG << "set_zvec(trait=" << trait << "); ";
  if (trait != 1) BOOST_THROW_EXCEPTION(::std::runtime_error("trait != 1"));
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BOOST_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  check_num_tag(length);
  zvec1_.assign(values, values+length);
  return 0;
}

int64_t BgmgCalculator::set_nvec(int trait, int length, float* values) {
  if (trait != 1) BOOST_THROW_EXCEPTION(::std::runtime_error("trait != 1"));
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BOOST_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  LOG << "set_nvec(trait=" << trait << "); ";
  check_num_tag(length);
  nvec1_.assign(values, values+length);
  return 0;
}


int64_t BgmgCalculator::set_weights(int length, float* values) {
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BOOST_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  LOG << "set_weights; ";
  check_num_tag(length);
  weights_.assign(values, values+length);
  return 0;
}

int64_t BgmgCalculator::set_option(char* option, double value) {
  LOG << "set_option(" << option << "=" << value << "); ";

  if (!strcmp(option, "kmax")) {
    clear_state(); k_max_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "r2min")) {
    clear_state(); r2_min_ = value; return 0;
  } else if (!strcmp(option, "max_causals")) {
    if (!last_num_causals_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("can't change max_causals after find_snp_order"));
    clear_state(); max_causals_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "num_components")) {
    if (!last_num_causals_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("can't change num_components after find_snp_order"));
    clear_state(); num_components_ = static_cast<int>(value); return 0;
  }

  BOOST_THROW_EXCEPTION(::std::runtime_error("unknown option"));
  return 0;
}

#define CHECK_SNP_INDEX(i) if (i < 0 || i >= num_snp_) BOOST_THROW_EXCEPTION(::std::runtime_error("CHECK_SNP_INDEX failed"));
#define CHECK_TAG_INDEX(i) if (i < 0 || i >= num_tag_) BOOST_THROW_EXCEPTION(::std::runtime_error("CHECK_TAG_INDEX failed"));

int64_t BgmgCalculator::set_tag_indices(int num_snp, int num_tag, int* tag_indices) {
  if (num_snp_ != -1 || num_tag_ != -1) BOOST_THROW_EXCEPTION(::std::runtime_error("can not call set_tag_indices twice"));

  LOG << "set_tag_indices(num_snp=" << num_snp << ", num_tag=" << num_tag << "); ";
  num_snp_ = num_snp;
  num_tag_ = num_tag;

  is_tag_.resize(num_snp, 0);
  snp_to_tag_.resize(num_snp, -1);
  tag_to_snp_.assign(tag_indices, tag_indices + num_tag);
  for (int i = 0; i < tag_to_snp_.size(); i++) {
    CHECK_SNP_INDEX(tag_to_snp_[i]);
    is_tag_[tag_to_snp_[i]] = 1;
    snp_to_tag_[tag_to_snp_[i]] = i;
  }
  return 0;
}


int64_t BgmgCalculator::set_ld_r2_coo(int length, int* snp_index, int* tag_index, float* r2) {
  if (!csr_ld_r2_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("can't call set_ld_r2_coo after set_ld_r2_csr"));
  LOG << "set_ld_r2_coo(length=" << length << "); ";

  if (last_num_causals_.empty()) find_snp_order();

  for (int i = 0; i < length; i++) {
    if (!std::isfinite(r2[i])) BOOST_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  int was = coo_ld_.size();
  for (int i = 0; i < length; i++) {
    CHECK_SNP_INDEX(snp_index[i]); CHECK_SNP_INDEX(tag_index[i]);
    if (r2[i] < r2_min_) continue;
    // tricky part here is that we take into account snp_can_be_causal_
    // there is no reason to keep LD information about certain causal SNP if we never selecting it as causal
    // (see how snp_can_be_causal_ is created during find_snp_order() call)
    if (snp_can_be_causal_[snp_index[i]] && is_tag_[tag_index[i]]) coo_ld_.push_back(std::make_tuple(snp_index[i], snp_to_tag_[tag_index[i]], r2[i]));
    if (snp_can_be_causal_[tag_index[i]] && is_tag_[snp_index[i]]) coo_ld_.push_back(std::make_tuple(tag_index[i], snp_to_tag_[snp_index[i]], r2[i]));
  }
  LOG << "set_ld_r2_coo: coo_ld_.size()=" << coo_ld_.size() << " (new: " << coo_ld_.size() - was << ")";
  return 0;
}

int64_t BgmgCalculator::set_ld_r2_csr() {
  if (coo_ld_.empty()) 
    BOOST_THROW_EXCEPTION(::std::runtime_error("coo_ld_ is empty"));

  LOG << "set_ld_r2_csr (coo_ld_.size()==" << coo_ld_.size() << "); ";

  for (int i = 0; i < tag_to_snp_.size(); i++)
    coo_ld_.push_back(std::make_tuple(tag_to_snp_[i], i, 1.0f));
  
  std::sort(coo_ld_.begin(), coo_ld_.end());

  csr_ld_tag_index_.reserve(coo_ld_.size());
  csr_ld_r2_.reserve(coo_ld_.size());

  for (int i = 0; i < coo_ld_.size(); i++) {
    csr_ld_tag_index_.push_back(std::get<1>(coo_ld_[i]));
    csr_ld_r2_.push_back(std::get<2>(coo_ld_[i]));
  }

  // find starting position for each snp
  csr_ld_snp_index_.resize(snp_to_tag_.size() + 1, coo_ld_.size());
  for (int i = (coo_ld_.size() - 1); i >= 0; i--) {
    int snp_index = std::get<0>(coo_ld_[i]);
    csr_ld_snp_index_[snp_index] = i;
  }

  for (int i = (csr_ld_snp_index_.size() - 2); i >= 0; i--)
    if (csr_ld_snp_index_[i] > csr_ld_snp_index_[i + 1])
      csr_ld_snp_index_[i] = csr_ld_snp_index_[i + 1];

  coo_ld_.clear();
  return 0;
}

class xorshf96  //period 2^96-1
{
public:
  using result_type = unsigned long;
  static constexpr result_type min() { return 0; }
  static constexpr result_type max() { return std::numeric_limits<unsigned long>::max(); }
  result_type operator()() {
    static unsigned long x = 123456789, y = 362436069, z = 521288629;
    unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;

    return z;
  }
};

// TBD: Can be speed up by partial random shuffle 
int64_t BgmgCalculator::find_snp_order() {
  LOG << "find_snp_order(num_components_=" << num_components_ << ", k_max_=" << k_max_ << ", max_causals_=" << max_causals_ << ")";

  assert(max_causals_ <= num_snp_);
  assert(num_components_ > 0 && num_components_ <= 3); // not supported?
  
  if (last_num_causals_.size() > 0) BOOST_THROW_EXCEPTION(::std::runtime_error("find_snp_order called twice"));

  snp_can_be_causal_.resize(num_snp_, 0);

  xorshf96 random_engine;
  std::vector<int> perm(num_snp_, 0);

  for (int component_index = 0; component_index < num_components_; component_index++) {
    snp_order_.push_back(std::make_shared<DenseMatrix<int>>(max_causals_, k_max_));
    tag_r2sum_.push_back(std::make_shared<DenseMatrix<float>>(num_tag_, k_max_));
    
    tag_r2sum_[component_index]->InitializeZeros();
    last_num_causals_.push_back(0);
    
    for (int k = 0; k < k_max_; k++) {
      for (int i = 0; i < num_snp_; i++) perm[i] = i;
      std::shuffle(perm.begin(), perm.end(), random_engine);
      for (int i = 0; i < max_causals_; i++) {
        (*snp_order_[component_index])(i, k) = perm[i];
        snp_can_be_causal_[perm[i]] = 1;
      }
    }
  }

  int num_can_be_causal = 0;
  for (int i = 0; i < num_snp_; i++) num_can_be_causal += snp_can_be_causal_[i];
  LOG << "num_can_be_causal = " << num_can_be_causal;
  return 0;
}

int64_t BgmgCalculator::find_tag_r2sum(int component_id, int num_causals) {
  assert(component_id >= 0 && component_id < num_components_); 
  assert(num_causals >= 0 && num_causals < max_causals_);

  if (last_num_causals_.empty()) find_snp_order();

  int last_num_causals = last_num_causals_[component_id]; 
  
  LOG << "find_tag_r2sum(component_id=" << component_id << ", num_causals=" << num_causals << ", last_num_causals=" << last_num_causals << ")";

  if (num_causals == last_num_causals) return 0;
  float sign;              // +1.0f if we should increase r2sum, -1.0f to decrease
  int scan_from, scan_to;  // inclusive 0-based indices of causal variants
  if (num_causals > last_num_causals) {
    sign = 1.0f;
    scan_from = last_num_causals;
    scan_to = num_causals - 1;
  } else {
    sign = -1.0f;
    scan_from = num_causals;
    scan_to = last_num_causals - 1;
  }
  
  for (int scan_index = scan_from; scan_index <= scan_to; scan_index++) {
    for (int k_index = 0; k_index < k_max_; k_index++) {
      int snp_index = (*snp_order_[component_id])(scan_index, k_index);
      int r2_index_from = csr_ld_snp_index_[snp_index];
      int r2_index_to = csr_ld_snp_index_[snp_index + 1];
      for (int r2_index = r2_index_from; r2_index < r2_index_to; r2_index++) {
        int tag_index = csr_ld_tag_index_[r2_index];
        float r2 = csr_ld_r2_[r2_index];
        (*tag_r2sum_[component_id])(tag_index, k_index) += (sign * r2);
      }
    }
  }

  last_num_causals_[component_id] = num_causals;
  return 0;
}

int64_t BgmgCalculator::set_hvec(int length, float* values) {
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BOOST_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  if (!hvec_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("can not set hvec twice"));

  LOG << "set_hvec(" << length << "); ";
  check_num_snp(length);
  hvec_.assign(values, values + length);

  for (int i = 0; i < csr_ld_r2_.size(); i++) {
    int tag = csr_ld_tag_index_[i];
    int snp = tag_to_snp_[tag];
    csr_ld_r2_[i] *= values[snp];
  }
  return 0;
}


int64_t BgmgCalculator::retrieve_tag_r2_sum(int component_id, int num_causal, int length, float* buffer) {
  if (length != k_max_ * num_tag_) BOOST_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));

  LOG << "retrieve_tag_r2_sum(component_id=" << component_id << ", num_causal=" << num_causal << ")";

  find_tag_r2sum(component_id, num_causal);
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    for (int k_index = 0; k_index < k_max_; k_index++) {
      float tag_r2sum = (*tag_r2sum_[component_id])(tag_index, k_index);
      buffer[tag_index * k_max_ + k_index] = tag_r2sum;
    }
  }
}

int64_t BgmgCalculator::calc_univariate_pdf(float pi_vec, float sig2_zero, float sig2_beta, int length, float* zvec, float* pdf) {
  // input buffer "zvec" contains z scores (presumably an equally spaced grid)
  // output buffer contains pdf(z), aggregated across all SNPs with corresponding weights
  if (nvec1_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("nvec1 is not set"));
  if (weights_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  assert(num_components_ == 1);
  int num_causals = static_cast<int>(round(pi_vec * num_snp_));
  if (num_causals >= max_causals_) BOOST_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  const int component_id = 0;   // univariate is always component 0.

  LOG << "calc_univariate_pdf(pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";

  find_tag_r2sum(component_id, num_causals);

  const float pi_k = 1. / static_cast<float>(k_max_);
  static const float inv_sqrt_2pi = 0.3989422804014327f;

// omp reduction on std::vector ( https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector ) - did not work for microsoft compiler
// #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
//                               std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
//                     initializer(omp_priv = omp_orig)
// Final solution is to do a the reduction with omp critical (see here http://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-reduction.html )

  // we accumulate crazy many small values - each of them is OK as float; the sum is also OK as float;  
  // but accumulation must be done with double precision.
  // std::vector<double> pdf_double(length, 0.0);
  std::valarray<double> pdf_double(0.0, length);

#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, length);
#pragma omp for schedule(static)
    for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
      if (weights_[tag_index] == 0) continue;
      for (int k_index = 0; k_index < k_max_; k_index++) {
        float tag_r2sum = (*tag_r2sum_[component_id])(tag_index, k_index);
        float sig2eff = tag_r2sum * nvec1_[tag_index] * sig2_beta + sig2_zero;
        float s = sqrt(sig2eff);

        for (int z_index = 0; z_index < length; z_index++) {
          float a = zvec[z_index] / s;
          float pdf_tmp = pi_k * inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
          pdf_double_local[z_index] += static_cast<double>(pdf_tmp * weights_[tag_index]);
        }
      }
    }
#pragma omp critical
    pdf_double += pdf_double_local;
  }

  for (int i = 0; i < length; i++) pdf[i] = static_cast<float>(pdf_double[i]);

  return 0;
}

double BgmgCalculator::calc_univariate_cost(float pi_vec, float sig2_zero, float sig2_beta) {
  if (zvec1_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("zvec1 is not set"));
  if (nvec1_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("nvec1 is not set"));
  if (weights_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  assert(num_components_ == 1);
  int num_causals = static_cast<int>(round(pi_vec * num_snp_));
  if (num_causals >= max_causals_) return 1e100; // too large pi_vec
  const int component_id = 0;   // univariate is always component 0.

  LOG << "calc_univariate_cost(pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";

  find_tag_r2sum(component_id, num_causals);

  const float pi_k = 1. / static_cast<float>(k_max_);
  static const float inv_sqrt_2pi = 0.3989422804014327f;

  double log_pdf_total = 0.0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;

    float pdf_tag = 0.0f;
    for (int k_index = 0; k_index < k_max_; k_index++) {
      float tag_r2sum = (*tag_r2sum_[component_id])(tag_index, k_index);
      float sig2eff = tag_r2sum * nvec1_[tag_index] * sig2_beta + sig2_zero;

      float s = sqrt(sig2eff);
      float a = zvec1_[tag_index] / s;
      float pdf = pi_k * inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
      pdf_tag += pdf;
    }
    log_pdf_total += -log(pdf_tag) * weights_[tag_index];
  }

  return log_pdf_total;
}

double BgmgCalculator::calc_bivariate_cost(int num_components, double* pi_vec, double* sig2_beta, double* rho_beta, double* sig2_zero, double rho_zero) {
  LOG << "calc_bivariate_cost(!!!not implemented!!!)";
  return 0.0;
}

// TBD: validate if input contains undefined values (or other out of range values)