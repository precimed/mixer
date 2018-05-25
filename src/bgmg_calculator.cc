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
#include <cmath>

#include "boost/throw_exception.hpp"

void BgmgCalculator::set_num_snp(int length) {
  if ((num_snp_ != -1) && (num_snp_ != length)) BOOST_THROW_EXCEPTION(::std::runtime_error("length != num_snps_"));
  num_snp_ = length;
}

void BgmgCalculator::set_num_tag(int length) {
  if ((num_tag_ != -1) && (num_tag_ != length)) BOOST_THROW_EXCEPTION(::std::runtime_error("length != num_tag_"));
  num_tag_ = length;
}

int64_t BgmgCalculator::set_zvec(int trait, int length, float* values) {
  status_.back() << "set_zvec(trait=" << trait << "); ";
  if (trait != 1) BOOST_THROW_EXCEPTION(::std::runtime_error("trait != 1"));
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BOOST_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  set_num_tag(length);
  zvec1_.assign(values, values+length);
  return 0;
}

int64_t BgmgCalculator::set_nvec(int trait, int length, float* values) {
  if (trait != 1) BOOST_THROW_EXCEPTION(::std::runtime_error("trait != 1"));
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BOOST_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  status_.back() << "set_nvec(trait=" << trait << "); ";
  set_num_tag(length);
  nvec1_.assign(values, values+length);
  return 0;
}


int64_t BgmgCalculator::set_weights(int length, float* values) {
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BOOST_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  status_.back() << "set_weights; ";
  set_num_tag(length);
  weights_.assign(values, values+length);
  return 0;
}

int64_t BgmgCalculator::set_option(char* option, double value) {
  status_.back() << "set_option(" << option << "=" << value << "); ";

  if (!strcmp(option, "kmax")) {
    if (!last_num_causals_.empty()) {
      last_num_causals_.clear();
      snp_order_.clear();
      tag_r2sum_.clear();
    }
    k_max_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "r2min")) {
    r2_min_ = value; return 0;
  } else if (!strcmp(option, "max_causals")) {
    if (!last_num_causals_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("can't change max_causals after find_snp_order"));
    max_causals_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "num_components")) {
    if (!last_num_causals_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("can't change num_components after find_snp_order"));
    num_components_ = static_cast<int>(value); return 0;
  }

  BOOST_THROW_EXCEPTION(::std::runtime_error("unknown option"));
  return 0;
}

#define CHECK_SNP_INDEX(i) if (i < 0 || i >= num_snp_) BOOST_THROW_EXCEPTION(::std::runtime_error("CHECK_SNP_INDEX failed"));
#define CHECK_TAG_INDEX(i) if (i < 0 || i >= num_tag_) BOOST_THROW_EXCEPTION(::std::runtime_error("CHECK_TAG_INDEX failed"));

int64_t BgmgCalculator::set_tag_indices(int num_snp, int num_tag, int* tag_indices) {
  status_.back() << "set_tag_indices(num_snp=" << num_snp << ", num_tag=" << num_tag << "); ";
  set_num_snp(num_snp);
  set_num_tag(num_tag);
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
  status_.back() << "set_ld_r2_coo(length=" << length << "); ";

  for (int i = 0; i < length; i++) {
    if (!std::isfinite(r2[i])) BOOST_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  int was = coo_ld_.size();
  for (int i = 0; i < length; i++) {
    CHECK_SNP_INDEX(snp_index[i]); CHECK_SNP_INDEX(tag_index[i]);
    if (r2[i] < r2_min_) continue;
    if (is_tag_[tag_index[i]]) coo_ld_.push_back(std::make_tuple(snp_index[i], snp_to_tag_[tag_index[i]], r2[i]));
    if (is_tag_[snp_index[i]]) coo_ld_.push_back(std::make_tuple(tag_index[i], snp_to_tag_[snp_index[i]], r2[i]));
  }
  status_.back() << "(add " << coo_ld_.size() - was << "); ";
  return 0;
}

int64_t BgmgCalculator::set_ld_r2_csr() {
  if (coo_ld_.empty()) 
    BOOST_THROW_EXCEPTION(::std::runtime_error("coo_ld_ is empty"));

  status_.back() << "set_ld_r2_csr (coo_ld_.size()==" << coo_ld_.size() << "); ";

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
  assert(max_causals_ <= num_snp_);
  assert(num_components_ > 0 && num_components_ <= 3); // not supported?
  
  if (last_num_causals_.size() > 0) BOOST_THROW_EXCEPTION(::std::runtime_error("find_snp_order called twice"));

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
      for (int i = 0; i < max_causals_; i++) (*snp_order_[component_index])(i, k) = perm[i];
    }
  }
  return 0;
}

int64_t BgmgCalculator::find_tag_r2sum(int component_id, int num_causals) {
  assert(component_id >= 0 && component_id < num_components_); 
  assert(num_causals >= 0 && num_causals < max_causals_);

  if (last_num_causals_.empty()) find_snp_order();

  int last_num_causals = last_num_causals_[component_id]; 
  
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

  status_.back() << "set_hvec(" << length << "); ";
  set_num_snp(length);
  for (int i = 0; i < csr_ld_r2_.size(); i++) {
    int tag = csr_ld_tag_index_[i];
    int snp = tag_to_snp_[tag];
    csr_ld_r2_[i] *= values[snp];
  }
  return 0;
}


int64_t BgmgCalculator::retrieve_tag_r2_sum(int component_id, int num_causal, int length, float* buffer) {
  if (length != k_max_ * num_tag_) BOOST_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));

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
  find_tag_r2sum(component_id, num_causals);

  const float pi_k = 1. / static_cast<float>(k_max_);
  static const float inv_sqrt_2pi = 0.3989422804014327f;

  // we accumulate crazy many small values - each of them is OK as float; the sum is also OK as float;  
  // but accumulation must be done with double precision.
  std::vector<double> pdf_double(length, 0.0);

  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    for (int k_index = 0; k_index < k_max_; k_index++) {
      float tag_r2sum = (*tag_r2sum_[component_id])(tag_index, k_index);
      float sig2eff = tag_r2sum * nvec1_[tag_index] * sig2_beta + sig2_zero;
      float s = sqrt(sig2eff);

      for (int z_index = 0; z_index < length; z_index++) {
        float a = zvec[z_index] / s;
        float pdf_tmp = pi_k * inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
        pdf_double[z_index] += static_cast<double>(pdf_tmp * weights_[tag_index]);
      }
    }
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
  find_tag_r2sum(component_id, num_causals);

  const float pi_k = 1. / static_cast<float>(k_max_);
  static const float inv_sqrt_2pi = 0.3989422804014327f;

  double log_pdf_total = 0.0;
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
  return 0.0;
}

const char* BgmgCalculator::status() {
  status_.push_back(std::stringstream());
  static std::string str;
  str = status_[status_.size() - 2].str();
  return str.c_str();
}

// TBD: validate if input contains undefined values (or other out of range values)