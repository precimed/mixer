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

#include <random>
#include <limits>
#include <algorithm>

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
  set_num_tag(length);
  zvec1_.assign(values, values+length);
  return 0;
}

int64_t BgmgCalculator::set_nvec(int trait, int length, float* values) {
  if (trait != 1) BOOST_THROW_EXCEPTION(::std::runtime_error("trait != 1"));
  status_.back() << "set_nvec(trait=" << trait << "); ";
  set_num_tag(length);
  nvec1_.assign(values, values+length);
  return 0;
}


int64_t BgmgCalculator::set_weights(int length, float* values) {
  status_.back() << "set_weights; ";
  set_num_tag(length);
  weights_.assign(values, values+length);
  return 0;
}

int64_t BgmgCalculator::set_option(char* option, double value) {
  status_.back() << "set_option(" << option << "=" << value << "); ";

  if (!strcmp(option, "kmax")) {
    k_max_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "r2min")) {
    r2_min_ = value; return 0;
  } else if (!strcmp(option, "max_causals")) {
    max_causals_ = static_cast<int>(value); return 0;
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

  int was = coo_ld_.size();
  for (int i = 0; i < length; i++) {
    CHECK_SNP_INDEX(snp_index[i]); CHECK_SNP_INDEX(tag_index[i]);
    if (r2[i] < r2_min_) continue;
    if (is_tag_[tag_index[i]]) coo_ld_.push_back(std::make_tuple(snp_index[i], snp_to_tag_[tag_index[i]], r2[i]));
    if (is_tag_[snp_index[i]]) coo_ld_.push_back(std::make_tuple(tag_index[i], snp_to_tag_[snp_index[i]], r2[i]));
  }
  status_.back() << "(add " << coo_ld_.size() - was << "); ";
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
int64_t BgmgCalculator::find_snp_sample() {
  snp_sample_ = std::make_shared<DenseMatrix<int>>(max_causals_, k_max_);
  
  std::vector<int> perm(num_snp_, 0);
  
  xorshf96 random_engine;

  for (int k = 0; k < k_max_; k++) {
    for (int i = 0; i < num_snp_; i++) perm[i] = i;
    std::shuffle(perm.begin(), perm.end(), random_engine);
    for (int i = 0; i < max_causals_; i++) (*snp_sample_)(i, k) = perm[i];
  }
  return 0;
}

int64_t BgmgCalculator::set_hvec(int length, float* values) {
  status_.back() << "set_hvec(" << length << "); ";
  set_num_snp(length);
  for (int i = 0; i < csr_ld_r2_.size(); i++) {
    int tag = csr_ld_tag_index_[i];
    int snp = tag_to_snp_[tag];
    csr_ld_r2_[i] *= values[snp];
  }
  return 0;
}



double BgmgCalculator::calc_univariate_cost(float pi_vec, float sig2_zero, float sig2_beta) {
  /*
  if (zvec1_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("zvec1 is not set"));
  if (nvec1_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("nvec1 is not set"));
  if (w_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
  if (r2_.empty() || r2_hist_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("ref_ld is not set"));
  
  //std::default_random_engine generator;
  xorshf96 generator;


  const float pi_k = 1. / static_cast<float>(k_max_);
  std::vector<float> sig2eff(k_max_, 0.0);
   
  float pdf_total = 0.0;
  for (int i = 0; i < num_snps_; i++) {
    sig2eff.assign(k_max_, sig2_zero);
    for (int r2bini = 0; r2bini < r2bins_; r2bini++) {
      std::binomial_distribution<int> distribution(r2_hist_[i + num_snps_ * r2bini], pi_vec);
      float r2eff = r2_[i + num_snps_ * r2bini] * sig2_beta * nvec1_[i] * hvec_[i];
      for (int k = 0; k < k_max_; k++) {
        int num_causals = distribution(generator);
        sig2eff[k] += static_cast<float>(num_causals) * r2eff;
      }
    }
    
    static const float inv_sqrt_2pi = 0.3989422804014327f;

    float pdf = 0.0f;
    // #pragma omp parallel for schedule(static) reduction(+ : pdf) 
    for (int k = 0; k < k_max_; k++) {
      float s = sqrt(sig2eff[k]);
      float a = zvec1_[i] / s;
      float tmp = pi_k * inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
      pdf += tmp;
    }

    pdf_total += -log(pdf) * w_[i];
  }

  return pdf_total;*/
  return 0.0;
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