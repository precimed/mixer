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

#include "boost/throw_exception.hpp"

void BgmgCalculator::set_num_snps(int length) {
  if ((num_snps_ != -1) && (num_snps_ != length)) BOOST_THROW_EXCEPTION(::std::runtime_error("length != num_snps_"));
  num_snps_ = length;
}

int64_t BgmgCalculator::set_zvec(int trait, int length, double* values) {
  if (trait != 1) BOOST_THROW_EXCEPTION(::std::runtime_error("trait != 1"));
  set_num_snps(length);

  zvec1_.assign(length, 0.0f);
  for (int i = 0; i < length; i++)
    zvec1_[i] = static_cast<float>(values[i]);
  return 0;
}

int64_t BgmgCalculator::set_nvec(int trait, int length, double* values) {
  if (trait != 1) BOOST_THROW_EXCEPTION(::std::runtime_error("trait != 1"));
  set_num_snps(length);
  nvec1_.assign(length, 0.0f);
  for (int i = 0; i < length; i++)
    nvec1_[i] = static_cast<float>(values[i]);
  return 0;
}

int64_t BgmgCalculator::set_hvec(int length, double* values) {
  set_num_snps(length);
  hvec_.assign(length, 0.0f);
  for (int i = 0; i < length; i++)
    hvec_[i] = static_cast<float>(values[i]);
  return 0;
}

int64_t BgmgCalculator::set_weights(int length, double* values) {
  set_num_snps(length);
  w_.assign(length, 0.0f);
  for (int i = 0; i < length; i++)
    w_[i] = static_cast<float>(values[i]);
  return 0;
}

int64_t BgmgCalculator::set_ref_ld(int length, int r2bins, double* sum_r2, double* sum_r4) {
  set_num_snps(length);
  r2bins_ = r2bins;
  r2_.assign(length * r2bins, 0);
  r2_hist_.assign(length * r2bins, 0);
  for (int i = 0; i < length * r2bins; i++) {
    if (sum_r4[i] == 0) continue;
    double n = round(sum_r2[i] * sum_r2[i] / sum_r4[i]);
    r2_hist_[i] = static_cast<int>(n);
    r2_[i] = static_cast<float>(sum_r2[i] / n);
  }
  return 0;
}

int64_t BgmgCalculator::set_option(char* option, int value) {
  if (!strcmp(option, "kmax")) {
    k_max_ = value; return 0;
  }

  BOOST_THROW_EXCEPTION(::std::runtime_error("unknown option"));
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

double BgmgCalculator::calc_univariate_cost(float pi_vec, float sig2_zero, float sig2_beta) {
  if (zvec1_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("zvec1 is not set"));
  if (nvec1_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("nvec1 is not set"));
  if (hvec_.empty()) BOOST_THROW_EXCEPTION(::std::runtime_error("hvec is not set"));
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
        sig2eff[k] += static_cast<float>(distribution(generator)) * r2eff;
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

  return pdf_total;
}

double BgmgCalculator::calc_bivariate_cost(int num_components, double* pi_vec, double* sig2_beta, double* rho_beta, double* sig2_zero, double rho_zero) {
  return 0.0;
}

// TBD: validate if input contains undefined values (or other out of range values)