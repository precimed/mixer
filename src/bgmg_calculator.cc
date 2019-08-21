/*ma
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

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_set_num_threads(i)
#endif

#include <chrono>
#include <random>
#include <limits>
#include <algorithm>
#include <vector>
#include <valarray>
#include <cmath>
#include <sstream>
#include <fstream>
#include <set>
#include <numeric>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include "cubature/cubature.h"

// Include namespace SEMT & global operators.
#define SEMT_DISABLE_PRINT 0
#include "semt/Semt.h"
// Include macros: INT, DINT, VAR, DVAR, PARAM, DPARAM
#include "semt/Shortcuts.h"

#include "bgmg_log.h"
#include "bgmg_parse.h"
#include "bgmg_math.h"
#include "fmath.hpp"

#include <immintrin.h>  // _mm_setcsr, _mm_getcsr

#define FLOAT_TYPE float

std::vector<float>* BgmgCalculator::get_zvec(int trait_index) {
  if ((trait_index != 1) && (trait_index != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  return (trait_index == 1) ? &zvec1_ : &zvec2_;
}

std::vector<float>* BgmgCalculator::get_nvec(int trait_index) {
  if ((trait_index != 1) && (trait_index != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  return (trait_index == 1) ? &nvec1_ : &nvec2_;
}

std::vector<float>* BgmgCalculator::get_causalbetavec(int trait_index) {
  if ((trait_index != 1) && (trait_index != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  return (trait_index == 1) ? &causalbetavec1_ : &causalbetavec2_;
}

BgmgCalculator::BgmgCalculator() : num_snp_(-1), num_tag_(-1), k_max_(100), seed_(0), 
    use_complete_tag_indices_(false), r2_min_(0.0), z1max_(1e10), z2max_(1e10), ld_format_version_(-1), num_components_(1), 
    max_causals_(100000), cost_calculator_(CostCalculator_Sampling), cache_tag_r2sum_(false), ld_matrix_csr_(*this),
    cubature_abs_error_(0), cubature_rel_error_(1e-4), cubature_max_evals_(0), calc_k_pdf_(false) {
  boost::posix_time::ptime const time_epoch(boost::gregorian::date(1970, 1, 1));
  seed_ = (boost::posix_time::microsec_clock::local_time() - time_epoch).ticks();

  // flush denormals to zero --- implemented only for GCC. Need the same for clang and MS VS.
  // https://stackoverflow.com/questions/9314534/why-does-changing-0-1f-to-0-slow-down-performance-by-10x
  // https://carlh.net/plugins/denormals.php
  // #if defined(__clang__)
  // #include <fenv.h>
  // fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV);
  // #elif  defined(_MSC_VER)
  // #include <immintrin.h>
  // _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  // _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

  _mm_setcsr( _mm_getcsr() | (1<<15) | (1<<6));
}

void BgmgCalculator::check_num_snp(int length) {
  if (num_snp_ == -1) BGMG_THROW_EXCEPTION(::std::runtime_error("call set_tag_indices first"));
  if (num_snp_ != length) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_snps_"));
}

void BgmgCalculator::check_num_tag(int length) {
  if (num_tag_ == -1) BGMG_THROW_EXCEPTION(::std::runtime_error("call set_tag_indices first"));
  if (num_tag_ != length) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_snps_"));
}

int64_t BgmgCalculator::set_zvec(int trait, int length, float* values) {
  if ((trait != 1) && (trait != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));

  int num_undef = 0;
  for (int i = 0; i < length; i++) if (!std::isfinite(values[i])) num_undef++;
  LOG << " set_zvec(trait=" << trait << "); num_undef=" << num_undef;
  check_num_tag(length);
  get_zvec(trait)->assign(values, values + length);
  return 0;
}

int64_t BgmgCalculator::set_nvec(int trait, int length, float* values) {
  if ((trait != 1) && (trait != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  
  int num_undef = 0;
  for (int i = 0; i < length; i++) if (!std::isfinite(values[i])) num_undef++;
  LOG << " set_nvec(trait=" << trait << "); num_undef=" << num_undef;
  check_num_tag(length);
  get_nvec(trait)->assign(values, values + length);
  return 0;
}

int64_t BgmgCalculator::set_causalbetavec(int trait, int length, float* values) {
  if ((trait != 1) && (trait != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  
  int num_undef = 0;
  for (int i = 0; i < length; i++) if (!std::isfinite(values[i])) num_undef++;
  if (num_undef > 0) BGMG_THROW_EXCEPTION(::std::runtime_error("undefined values not allowed in causalbetavec, use zero instead"));
  LOG << " set_causalbetavec(trait=" << trait << "); num_undef=" << num_undef;
  check_num_snp(length);
  get_causalbetavec(trait)->assign(values, values + length);
  return 0;
}

int64_t BgmgCalculator::set_weights(int length, float* values) {
  int nnz = 0;
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
    if (values[i] != 0) nnz++;
  }

  LOG << " set_weights(length=" << length << "), nnz=" << nnz << "; ";
  check_num_tag(length);
  weights_.assign(values, values+length);
  return 0;
}

int64_t BgmgCalculator::set_option(char* option, double value) {
  LOG << " set_option(" << option << "=" << value << "); ";

  if (!strcmp(option, "diag")) {
    log_diagnostics(); return 0;
  } else if (!strcmp(option, "kmax")) {
    clear_state(); k_max_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "r2min")) {
    clear_state(); r2_min_ = value; return 0;
  } else if (!strcmp(option, "max_causals")) {
    if (!snp_order_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't change max_causals after find_snp_order"));
    clear_state(); max_causals_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "num_components")) {
    if (!snp_order_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't change num_components after find_snp_order"));
    clear_state(); num_components_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "seed")) {
    seed_ = static_cast<int64_t>(value); return 0;
  } else if (!strcmp(option, "calc_k_pdf")) {
    calc_k_pdf_ = (value != 0); return 0;
  } else if (!strcmp(option, "cubature_max_evals")) {
    cubature_max_evals_ = static_cast<int64_t>(value); return 0;
  } else if (!strcmp(option, "cubature_abs_error")) {
    cubature_abs_error_ = value; return 0;
  } else if (!strcmp(option, "cubature_rel_error")) {
    cubature_rel_error_ = value; return 0;
  } else if (!strcmp(option, "fast_cost")) {
    cost_calculator_ = (value != 0) ? CostCalculator_Sampling : CostCalculator_Gaussian; return 0;
  } else if (!strcmp(option, "cost_calculator")) {
    int int_value = (int)value;
    if (int_value < 0 || int_value > 2) BGMG_THROW_EXCEPTION(::std::runtime_error("cost_calculator value must be 0 (Sampling), 1 (Gaussian) or 2 (Convolve)"));
    cost_calculator_ = (CostCalculator)int_value; return 0;
  } else if (!strcmp(option, "z1max")) {
    if (value <= 0) BGMG_THROW_EXCEPTION(::std::runtime_error("zmax must be positive"));
    z1max_ = value; return 0;
  } else if (!strcmp(option, "z2max")) {
    if (value <= 0) BGMG_THROW_EXCEPTION(::std::runtime_error("zmax must be positive"));
    z2max_ = value; return 0;
  } else if (!strcmp(option, "ld_format_version")) {
    ld_format_version_ = int(value); return 0;
  } else if (!strcmp(option, "use_complete_tag_indices")) {
    use_complete_tag_indices_ = (value != 0); return 0;
  } else if (!strcmp(option, "threads")) {
    if (value > 0) {
      LOG << " omp_set_num_threads(" << static_cast<int>(value) << ")";
      omp_set_num_threads(static_cast<int>(value));
    }
    return 0;
  } else if (!strcmp(option, "cache_tag_r2sum")) {
    cache_tag_r2sum_ = (value != 0);
    for (int component_id = 0; component_id < num_components_; component_id++) clear_tag_r2sum(component_id);
    return 0;
  }

  BGMG_THROW_EXCEPTION(::std::runtime_error("unknown option"));
  return 0;
}

int64_t BgmgCalculator::set_tag_indices(int num_snp, int num_tag, int* tag_indices) {
  if (num_snp_ != -1 || num_tag_ != -1) BGMG_THROW_EXCEPTION(::std::runtime_error("can not call set_tag_indices twice"));

  LOG << " set_tag_indices(num_snp=" << num_snp << ", num_tag=" << num_tag << "); ";
  num_snp_ = num_snp;
  num_tag_ = num_tag;

  is_tag_.resize(num_snp, 0);
  snp_to_tag_.resize(num_snp, -1);
  tag_to_snp_.assign(tag_indices, tag_indices + num_tag);
  for (int i = 0; i < tag_to_snp_.size(); i++) {
    CHECK_SNP_INDEX((*this), tag_to_snp_[i]);
    is_tag_[tag_to_snp_[i]] = 1;
    snp_to_tag_[tag_to_snp_[i]] = i;
  }

  return 0;
}

int64_t BgmgCalculator::set_ld_r2_coo(int chr_label, int64_t length, int* snp_index, int* tag_index, float* r) {
  return ld_matrix_csr_.set_ld_r2_coo(chr_label, length, snp_index, tag_index, r, r2_min_);
}

int64_t BgmgCalculator::set_ld_r2_coo(int chr_label, const std::string& filename) {
  if (ld_format_version_ == 0)
    return ld_matrix_csr_.set_ld_r2_coo_version0(chr_label, filename, r2_min_);

  return ld_matrix_csr_.set_ld_r2_coo(chr_label, filename, r2_min_);
}

int64_t BgmgCalculator::set_ld_r2_csr(int chr_label) {
  int64_t retval = ld_matrix_csr_.set_ld_r2_csr(r2_min_, chr_label);
  return retval;
}


/*
// Nice trick, but not so important for performance.
// We use std::mt19937_64.
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
*/

int64_t BgmgCalculator::find_snp_order() {
  if (max_causals_ <= 0 || max_causals_ > num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("find_snp_order: max_causals_ <= 0 || max_causals_ > num_snp_"));
  if (num_components_ <= 0 || num_components_ > 3) BGMG_THROW_EXCEPTION(::std::runtime_error("find_snp_order: num_components_ must be between 1 and 3"));
  if (snp_order_.size() > 0) BGMG_THROW_EXCEPTION(::std::runtime_error("find_snp_order: called twice"));

  LOG << ">find_snp_order(num_components_=" << num_components_ << ", k_max_=" << k_max_ << ", max_causals_=" << max_causals_ << ")";
  SimpleTimer timer(-1);

  // Right now all SNPs must be included in snp_can_be_causal_.
  // Remember that the only purpose of snp_can_be_causal_ is to limit the information that store about the LD matrix.
  // At some point we use LD structure only to calculate tag_r2, so we only store r2 for SNPs that are selected as causal by find_snp_order.
  // Later we've started to use LD structure to 
  //  - perform random pruning (=> LD must be stored for all tag variants)
  //    "for (int i = 0; i < num_tag_; i++) snp_can_be_causal_[tag_to_snp_[i]] = 1;"
  //  - calculate sum_r2, sum_r4 (=> LD must be stored for all variants, OR we need to change the logic and load hvec so that sum_r2 and sum_r4 are calculated on the fly in set_ld_r2_coo.
  // For now we simply set snp_can_be_causal_ to 1 and store LD structure for all variants.
  std::vector<char> snp_can_be_causal(num_snp_, 0);  // mask of SNPs that may be causal (e.i. included in snp_order array)

  SimpleTimer log_timer(10000); // log some message each 10 seconds
  for (int component_index = 0; component_index < num_components_; component_index++) {
    if (log_timer.fire())
      LOG << " find_snp_order still working, component_id=" << component_index;

    snp_order_.push_back(std::make_shared<DenseMatrix<int>>(max_causals_, k_max_));

    if (cache_tag_r2sum_) clear_tag_r2sum(component_index);
    
#pragma omp parallel
    {
      std::vector<int> perm(num_snp_, 0);

#pragma omp for schedule(static)
      for (int k = 0; k < k_max_; k++) {
        for (int i = 0; i < num_snp_; i++) perm[i] = i;

        std::mt19937_64 random_engine;
        random_engine.seed(seed_ + component_index * k_max_ + k);  // ensure each k in each component starts with its own seed.

        // perform partial Fisher Yates shuffle (must faster than full std::shuffle)
        // swap_offset is a random integer, with max of n-1, n-2, n-3, ..., n-max_causals
        for (int i = 0; i < max_causals_; i++) {
          const int swap_offset = std::uniform_int_distribution<int>(0, num_snp_ - i - 1)(random_engine);
          std::iter_swap(perm.begin() + i, perm.begin() + i + swap_offset);
        }

        for (int i = 0; i < max_causals_; i++) {
          (*snp_order_[component_index])(i, k) = perm[i];
        }
      }
    }

    // Fill in snp_can_be_causal
    for (int k = 0; k < k_max_; k++) {
      for (int i = 0; i < max_causals_; i++) {
        snp_can_be_causal[(*snp_order_[component_index])(i, k)] = 1;
      }
    }
  }

  int num_can_be_causal = 0;
  for (int i = 0; i < num_snp_; i++) num_can_be_causal += snp_can_be_causal[i];
  LOG << "<find_snp_order: num_can_be_causal = " << num_can_be_causal << ", elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::find_tag_r2sum(int component_id, float num_causals) {
  if (!cache_tag_r2sum_) BGMG_THROW_EXCEPTION(::std::runtime_error("find_tag_r2sum can be used only with cache_tag_r2sum==true"));
  if (num_causals < 0 || num_causals >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("find_tag_r2sum: num_causals < 0 || num_causals >= max_causals_"));
  if (component_id < 0 || component_id >= num_components_) BGMG_THROW_EXCEPTION(::std::runtime_error("find_tag_r2sum: component_id must be between 0 and num_components_"));

  const float num_causals_original = num_causals;
  if (snp_order_.empty()) find_snp_order();

  float last_num_causals = last_num_causals_[component_id]; 
  const float last_num_causals_original = last_num_causals;
  
  LOG << ">find_tag_r2sum(component_id=" << component_id << ", num_causals=" << num_causals << ", last_num_causals=" << last_num_causals << ")";
  SimpleTimer timer(-1);

  // if num_causal is more than twice lower than last_num_causals we should re-calculate tag_r2sum from scratch.
  if (num_causals < (last_num_causals / 2)) {
    clear_tag_r2sum(component_id);
    last_num_causals = 0.0f;
  }

  // changeset contains a list of indices with corresponding weight
  // indices apply to snp_order_[component_id] array.
  // weights are typicaly +1 (to increase by r2) or -1 (to decrease by r2).
  // First and last weights is float-point number between 1 and -1,
  // to handle cases when num_causals is float-point number (derived from pivec).
  // This is important for fminsearch which get's confused if cost is a stepwise of pivec.
  std::vector<std::pair<int, float>> changeset;
  
  // Decreasing number of causals from B to A has an opposite effect to increasing from A to B.
  // To handle decreasing case we just swap num_causals and last_num_causals, and set sign to -1.0f.
  float sign = 1.0f;
  if (num_causals < last_num_causals) {
    float tmp = num_causals; num_causals = last_num_causals; last_num_causals = tmp;
    sign = -1.0f;
  }

  // There are 3 cases
  // 1. floor(num_causals) == floor(last_num_causals)
  // 2. floor(num_causals) == floor(last_num_causals) + 1
  // 3. floor(num_causals) >= floor(last_num_causals) + 2

  float floor_num_causals = floor(num_causals);
  float floor_last_num_causals = floor(last_num_causals);
  if ((int)floor_num_causals == (int)floor_last_num_causals) {
    changeset.push_back(std::make_pair((int)floor_last_num_causals, sign * (num_causals - last_num_causals)));
  }
  else if ((int)floor_num_causals >= ((int)floor_last_num_causals + 1)) {
    // handle case 2 and case 3 - lower boundary
    changeset.push_back(std::make_pair((int)floor_last_num_causals, sign * (floor_last_num_causals + 1.0f - last_num_causals)));

    // happends for the case 3 - bulk change (empty loop in case 2)
    for (int i = ((int)floor_last_num_causals + 1); i < (int)floor_num_causals; i++) {
      changeset.push_back(std::make_pair(i, sign));
    }

    // handle case 2 and case 3 - upper boundary
    changeset.push_back(std::make_pair((int)floor_num_causals, sign * (num_causals - floor_num_causals)));
  }
  else {
    BGMG_THROW_EXCEPTION(::std::runtime_error("floor_num_causals < floor_last_num_causals"));
  }

  // apply infinitesimal model to adjust tag_r2sum for all r2 that are below r2min (and thus do not contribute via resampling)
  const std::vector<float>& tag_sum_r2_below_r2min = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2(LD_TAG_COMPONENT_BELOW_R2MIN);
  const float pival_delta = (num_causals_original - last_num_causals_original) / static_cast<float>(num_snp_);

  std::vector<float> hvec;
  find_hvec(*this, &hvec);

  // it is OK to parallelize the following loop on k_index, because:
  // - all structures here are readonly, except tag_r2sum_ that we are accumulating
  // - two threads will never touch the same memory location (that's why we choose k_index as an outer loop)
#pragma omp parallel
{
  LdMatrixRow ld_matrix_row;

#pragma omp for schedule(static)
  for (int k_index = 0; k_index < k_max_; k_index++) {
    for (auto change : changeset) {
      int scan_index = change.first;
      float scan_weight = change.second;
      int snp_index = (*snp_order_[component_id])(scan_index, k_index);  // index of a causal snp
      ld_matrix_csr_.extract_row(snp_index, &ld_matrix_row);
      auto iter_end = ld_matrix_row.end();
      for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
        int tag_index = iter.tag_index();
        float r2 = iter.r2();
        float hval = hvec[snp_index];
        (*tag_r2sum_[component_id])(tag_index, k_index) += (scan_weight * r2 * hval);
      }
    }
  }
}

{
  SimpleTimer timer2(-1);
  int num_tag_inf_adjusted = 0;  // number of snps adjusted according to infinitesimal model
#pragma omp parallel for schedule(static) reduction(+: num_tag_inf_adjusted)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    float inf_adj = pival_delta * tag_sum_r2_below_r2min[tag_index];
    if (inf_adj != 0) {
      num_tag_inf_adjusted++;
      for (int k_index = 0; k_index < k_max_; k_index++) {
        (*tag_r2sum_[component_id])(tag_index, k_index) += inf_adj;
      }
    }
  }
  if (num_tag_inf_adjusted > 0)
    LOG << " apply infinitesimal model to " << num_tag_inf_adjusted << " tag SNPs, to adjust tag_r2sum for all r2 that are below r2min, elapsed time " << timer2.elapsed_ms() << "ms";
}

  LOG << "<find_tag_r2sum(component_id=" << component_id << ", num_causals=" << num_causals_original << ", last_num_causals=" << last_num_causals << "), elapsed time " << timer.elapsed_ms() << "ms";

  last_num_causals_[component_id] = num_causals_original;
  return 0;
}

int64_t BgmgCalculator::set_mafvec(int length, float* values) {
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  if (!mafvec_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can not set mafvec twice"));

  LOG << ">set_mafvec(" << length << "); ";
  check_num_snp(length);
  mafvec_.assign(values, values + length);
  LOG << "<set_mafvec(" << length << "); ";
  return 0;
}

int64_t BgmgCalculator::retrieve_ld_tag_r2_sum(int length, float* buffer) {
  check_num_tag(length);
  LOG << " retrieve_ld_tag_r2_sum()";
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    buffer[tag_index] = ld_matrix_csr_.ld_tag_sum()->ld_tag_sum_r2()[tag_index];
  }
  return 0;
}

int64_t BgmgCalculator::retrieve_ld_tag_r4_sum(int length, float* buffer) {
  check_num_tag(length);
  LOG << " retrieve_ld_tag_r4_sum()";
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    buffer[tag_index] = ld_matrix_csr_.ld_tag_sum()->ld_tag_sum_r4()[tag_index];
  }
  return 0;
}

int64_t BgmgCalculator::retrieve_tag_r2_sum(int component_id, float num_causal, int length, float* buffer) {
  if (length != (k_max_ * num_tag_)) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (num_causal < 0 && !cache_tag_r2sum_) BGMG_THROW_EXCEPTION(::std::runtime_error("retrieve_tag_r2sum with num_causal<0 is meant for cache_tag_r2sum==true"));
  if (component_id < 0 || component_id >= num_components_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong component_id"));

  LOG << " retrieve_tag_r2_sum(component_id=" << component_id << ", num_causal=" << num_causal << ")";

  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) {
    // use negative to retrieve tag_r2_sum for last_num_causal (for debugging purpose)
    if (num_causal >= 0) {
      find_tag_r2sum(component_id, num_causal);
    }

    for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
      for (int k_index = 0; k_index < k_max_; k_index++) {
        float tag_r2sum = (*tag_r2sum_[component_id])(tag_index, k_index);
        buffer[tag_index * k_max_ + k_index] = tag_r2sum;
      }
    }
  } else {
#pragma omp parallel
    {
      std::vector<float> tag_r2sum(num_tag_, 0.0f);

#pragma omp for schedule(static)
      for (int k_index = 0; k_index < k_max_; k_index++) {
        find_tag_r2sum_no_cache(0, num_causal, k_index, &tag_r2sum);
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          buffer[tag_index * k_max_ + k_index] = tag_r2sum[tag_index];
        }
      }
    }
  }
  return 0;
}

// pdf of gaussian (normal) distribution with 0 mean and std error of s
// z is the point where to calculate pdf
template<typename T>
inline T gaussian_pdf(const T z, const T s) {
  static const T inv_sqrt_2pi = static_cast<T>(0.3989422804014327);
  const T a = z / s;
  const T pdf = inv_sqrt_2pi / s * std::exp(static_cast<T>(-0.5) * a * a);
  return pdf + std::numeric_limits<T>::min();
}

// censored_cdf gives cumulated probability of |z|>zmax, where z is normal random variable
// http://www.public.iastate.edu/~stat415/meeker/ml_estimation_chapter.pdf
// Principles of Maximum Likelihood Estimation and The Analysis of Censored Data
template<typename T>
inline T censored_cdf(const T z, const T s) {
  assert(z >= 0);
  static const T inv_sqrt_2 = static_cast<T>(0.7071067811865475);
  return std::erfc((z / s) * inv_sqrt_2) + std::numeric_limits<T>::min();
}


/*
// partial specification for float, to use fmath::exp instead of std::exp
template<>
inline float gaussian_pdf<float>(const float z, const float s) {
  static const float inv_sqrt_2pi = static_cast<float>(0.3989422804014327);
  const float a = z / s;
  const float pdf = inv_sqrt_2pi / s * fmath::exp(static_cast<float>(-0.5) * a * a);
  return pdf;
}
*/

template<typename T>
inline T gaussian2_pdf(const T z1, const T z2, const T a11, const T a12, const T a22) {
  static const T log_pi = static_cast<T>(-1.0 * log(2.0 * 3.14159265358979323846));

  // Calculation of log - likelihood and pdf, specific to bivariate normal
  // distribution with zero mean.It takes into account an explicit formula
  // for inverse 2x2 matrix, S = [a b; c d], => S^-1 = [d - b; -c a] . / det(S)
  const T dt = a11 * a22 - a12 * a12;  // det(S)

  const T log_exp = -0.5 * (a22*z1*z1 + a11*z2*z2 - 2.0*a12*z1*z2) / dt;
  const T log_dt = -0.5 * std::log(dt);

  const T pdf = std::exp(log_pi + log_dt + log_exp);
  return pdf + std::numeric_limits<T>::min();
}

/*
template<>
inline float gaussian2_pdf(const float z1, const float z2, const float a11, const float a12, const float a22) {
  static const float log_pi = static_cast<float>(-1.0 * log(2.0 * 3.14159265358979323846));

  // Calculation of log - likelihood and pdf, specific to bivariate normal
  // distribution with zero mean.It takes into account an explicit formula
  // for inverse 2x2 matrix, S = [a b; c d], => S^-1 = [d - b; -c a] . / det(S)
  const float dt = a11 * a22 - a12 * a12;  // det(S)

  const float log_exp = -0.5 * (a22*z1*z1 + a11*z2*z2 - 2.0*a12*z1*z2) / dt;
  const float log_dt = -0.5 * fmath::log(dt);

  const float pdf = fmath::exp(log_pi + log_dt + log_exp) + std::numeric_limits<float>::min();
  return pdf;
}
*/

int64_t BgmgCalculator::calc_univariate_pdf(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* zvec, float* pdf) {
  // input buffer "zvec" contains z scores (presumably an equally spaced grid)
  // output buffer contains pdf(z), aggregated across all SNPs with corresponding weights
  
  std::vector<float>& nvec(*get_nvec(trait_index));
  if (nvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  float num_causals = pi_vec * static_cast<float>(num_snp_);
  if ((int)num_causals >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  const int component_id = 0;   // univariate is always component 0.

  LOG << ">calc_univariate_pdf(trait_index="<< trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ", length(zvec)=" << length << ")";
  SimpleTimer timer(-1);

  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) {
    find_tag_r2sum(component_id, num_causals);
  }

  const double pi_k = 1.0 / static_cast<double>(k_max_);

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
    std::vector<float> tag_r2sum(num_tag_, 0.0f);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {

      if (cache_tag_r2sum_) {
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          tag_r2sum[tag_index] = (*tag_r2sum_[0])(tag_index, k_index);
        }
      } else {
        find_tag_r2sum_no_cache(0, num_causals, k_index, &tag_r2sum);
      }

      for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
        if (weights_[tag_index] == 0) continue;
        if (!std::isfinite(nvec[tag_index])) continue;

        const double tag_weight = static_cast<double>(weights_[tag_index]);

        float tag_r2sum_value = tag_r2sum[tag_index];
        float sig2eff = tag_r2sum_value * nvec[tag_index] * sig2_beta + sig2_zero;
        float s = sqrt(sig2eff);

        for (int z_index = 0; z_index < length; z_index++) {
          double pdf_tmp = static_cast<double>(gaussian_pdf<FLOAT_TYPE>(zvec[z_index], s));
          pdf_double_local[z_index] += pi_k * pdf_tmp * tag_weight;
        }
      }
    }
#pragma omp critical
    pdf_double += pdf_double_local;
  }

  for (int i = 0; i < length; i++) pdf[i] = static_cast<float>(pdf_double[i]);
  LOG << "<calc_univariate_pdf(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::calc_univariate_power(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, float zthresh, int length, float* nvec, float* svec) {
  // input buffer "nvec" contains a set of sample sizes (N) to calculate power
  // output buffer svec(n), e.i. a fraction of heritability explained by genome-wide significant SNPs, aggregated across all SNPs with corresponding weights.
  //
  // NB! This function uses analytical formula for the following double integral:
  // C(z) = E(\delta^2 | z) * P(z) = \int_{z : |z| > zt} \int_{delta} p(z|delta) p(delta) delta^2 d[delta] d[z]
  // As long as we fix the set of causal variants p(delta) is just a normal distribution with zero mean variance "delta2eff" (see code below).
  // In this case the integral can be taken analytically, and benefit from the fact that both numerator and denominator in S(N) formula are additive across
  // tag SNPs and across resampling iterations (1..kmax).

  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  float num_causals = pi_vec * static_cast<float>(num_snp_);
  if ((int)num_causals >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  const int component_id = 0;   // univariate is always component 0.

  LOG << ">calc_univariate_power(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ", zthresh=" << zthresh << ", length(nvec)=" << length << ")";
  SimpleTimer timer(-1);

  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) {
    find_tag_r2sum(component_id, num_causals);
  }

  const double pi_k = 1.0 / static_cast<double>(k_max_);

  std::valarray<double> s_numerator_global(0.0, length);
  std::valarray<double> s_denominator_global(0.0, length);

#pragma omp parallel
  {
    std::valarray<double> s_numerator_local(0.0, length);
    std::valarray<double> s_denominator_local(0.0, length);
    std::vector<float> tag_r2sum(num_tag_, 0.0f);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {

      if (cache_tag_r2sum_) {
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          tag_r2sum[tag_index] = (*tag_r2sum_[0])(tag_index, k_index);
        }
      }
      else {
        find_tag_r2sum_no_cache(0, num_causals, k_index, &tag_r2sum);
      }

      for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
        if (weights_[tag_index] == 0) continue;
        const double tag_weight = static_cast<double>(weights_[tag_index]);

        float tag_r2sum_value = tag_r2sum[tag_index];
        for (int n_index = 0; n_index < length; n_index++) {
          float delta2eff = tag_r2sum_value * nvec[n_index] * sig2_beta;
          float sig2eff = delta2eff + sig2_zero;
          float sqrt_sig2eff = sqrt(sig2eff);
          static const float sqrt_2 = sqrtf(2.0);
          float numerator1 = gaussian_pdf<FLOAT_TYPE>(zthresh, sqrt_sig2eff) * 2 * delta2eff * delta2eff * zthresh / sig2eff;
          float numerator2 = std::erfcf(zthresh / (sqrt_2 * sqrt_sig2eff)) * delta2eff;
          s_numerator_local[n_index] += numerator1 + numerator2;
          s_denominator_local[n_index] += delta2eff;
        }
      }
    }
#pragma omp critical
    {
      s_numerator_global += s_numerator_local;
      s_denominator_global += s_denominator_local;
    }
  }

  for (int i = 0; i < length; i++) svec[i] = static_cast<float>(s_numerator_global[i] / s_denominator_global[i]);
  LOG << "<calc_univariate_power(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ", zthresh=" << zthresh << ", length(nvec)=" << length << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::calc_univariate_delta_posterior(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* c0, float* c1, float* c2) {
  // c0 = c(0), c1=c(1), c2=c(2), where c(q) = \int_\delta \delta^q P(z|delta) P(delta)
  // c(q) is define so that:
  //  E(\delta^2|z_j) = c2[j]/c0[j];
  //  E(\delta  |z_j) = c1[j]/c0[j];

  if ((length == 0) || (length != num_tag_)) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_tag_"));

  float num_causals = pi_vec * static_cast<float>(num_snp_);
  if ((int)num_causals >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  const int component_id = 0;   // univariate is always component 0.

  std::vector<float>& zvec(*get_zvec(trait_index));
  std::vector<float>& nvec(*get_nvec(trait_index));
  if (zvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec is not set"));
  if (nvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec is not set"));

  LOG << ">calc_univariate_delta_posterior(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ", length(nvec)=" << length << ")";
  SimpleTimer timer(-1);

  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) {
    find_tag_r2sum(component_id, num_causals);
  }

  std::valarray<double> c0_global(0.0f, num_tag_);
  std::valarray<double> c1_global(0.0f, num_tag_);
  std::valarray<double> c2_global(0.0f, num_tag_);

#pragma omp parallel
  {
    std::vector<float> tag_r2sum(num_tag_, 0.0f);
    std::valarray<double> c0_local(0.0f, num_tag_);
    std::valarray<double> c1_local(0.0f, num_tag_);
    std::valarray<double> c2_local(0.0f, num_tag_);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {
      if (cache_tag_r2sum_) {
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          tag_r2sum[tag_index] = (*tag_r2sum_[0])(tag_index, k_index);
        }
      }
      else {
        find_tag_r2sum_no_cache(0, num_causals, k_index, &tag_r2sum);
      }

      for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
        if (!std::isfinite(zvec[tag_index]) || !std::isfinite(nvec[tag_index])) continue;

        const float tag_r2sum_value = tag_r2sum[tag_index];
        const float delta2eff = tag_r2sum_value * nvec[tag_index] * sig2_beta;  // S^2_kj
        const float sig2eff = delta2eff + sig2_zero;
        const float sig2eff_1_2 = sqrt(sig2eff);
        const float sig2eff_3_2 = sig2eff_1_2 * sig2eff;
        const float sig2eff_5_2 = sig2eff_3_2 * sig2eff;

        const float z = zvec[tag_index];
        const float exp_common = std::exp(-0.5f*z*z / sig2eff);

        c0_local[tag_index] += (exp_common / sig2eff_1_2);
        c1_local[tag_index] += (exp_common / sig2eff_3_2) * z * delta2eff;
        c2_local[tag_index] += (exp_common / sig2eff_5_2) *     delta2eff * (sig2_zero*sig2_zero + sig2_zero*delta2eff + z*z*delta2eff);
      }
    }

#pragma omp critical
    {
      c0_global += c0_local;
      c1_global += c1_local;
      c2_global += c2_local;
    }
  }

  // save results to output buffers
  const double pi_k = 1.0 / static_cast<double>(k_max_);
  static const double inv_sqrt_2pi = 0.3989422804014327;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (!std::isfinite(zvec[tag_index]) || !std::isfinite(nvec[tag_index])) continue;
    c0[tag_index] = pi_k * inv_sqrt_2pi * c0_global[tag_index];
    c1[tag_index] = pi_k * inv_sqrt_2pi * c1_global[tag_index];
    c2[tag_index] = pi_k * inv_sqrt_2pi * c2_global[tag_index];
  }

  LOG << "<calc_univariate_delta_posterior(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ", length(nvec)=" << length << "), elapsed time " << timer.elapsed_ms() << "ms";
}

double BgmgCalculator::calc_univariate_cost(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<float>& zvec(*get_zvec(trait_index));
  if (zvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec is not set"));
  if (nvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  double cost;
  if (cost_calculator_ == CostCalculator_Gaussian) cost = calc_univariate_cost_fast(trait_index, pi_vec, sig2_zero, sig2_beta);
  else if (cost_calculator_ == CostCalculator_Convolve) cost = calc_univariate_cost_convolve(trait_index, pi_vec, sig2_zero, sig2_beta);
  else if (!cache_tag_r2sum_) cost = calc_univariate_cost_nocache(trait_index, pi_vec, sig2_zero, sig2_beta);
  else cost = calc_univariate_cost_cache(trait_index, pi_vec, sig2_zero, sig2_beta);

  loglike_cache_.add_entry(pi_vec, sig2_zero, sig2_beta, cost);
  return cost;
}

double BgmgCalculator::calc_univariate_cost_cache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<float>& zvec(*get_zvec(trait_index));

  float num_causals = pi_vec * static_cast<float>(num_snp_);
  if ((int)num_causals >= max_causals_) return 1e100; // too large pi_vec
  const int component_id = 0;   // univariate is always component 0.
    
  LOG << ">calc_univariate_cost(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";
  SimpleTimer timer(-1);

  find_tag_r2sum(component_id, num_causals);

  std::valarray<float> fixed_effect_delta(0.0, num_tag_);
  calc_fixed_effect_delta_from_causalbetavec(trait_index, &fixed_effect_delta);

  const double pi_k = 1.0 / static_cast<double>(k_max_);
  
  double log_pdf_total = 0.0;
  int num_infinite = 0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_infinite)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec[tag_index]) || !std::isfinite(nvec[tag_index])) continue;

    double pdf_tag = 0.0f;
    for (int k_index = 0; k_index < k_max_; k_index++) {
      float tag_r2sum = (*tag_r2sum_[component_id])(tag_index, k_index);
      float sig2eff = tag_r2sum * nvec[tag_index] * sig2_beta + sig2_zero;

      const float tag_z = zvec[tag_index] - fixed_effect_delta[tag_index];  // apply causalbetavec;
      float s = sqrt(sig2eff);
      const bool censoring = std::abs(tag_z) > z1max_;

      double pdf = static_cast<double>(censoring ? censored_cdf<FLOAT_TYPE>(z1max_, s) : gaussian_pdf<FLOAT_TYPE>(tag_z, s));
      pdf_tag += pi_k * pdf;
    }
    double increment = -std::log(pdf_tag) * static_cast<double>(weights_[tag_index]);
    if (!std::isfinite(increment)) num_infinite++;
    log_pdf_total += increment;
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_univariate_cost(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

struct pi_struct
{
  constexpr static SEMT_PRECISION value = 3.14159265358979323846;
};
struct inv_sqrt_2pi_struct
{
  constexpr static SEMT_PRECISION value = 0.3989422804014327;
};

SEMT::Expr<SEMT::Literal<pi_struct>> pi_value;
SEMT::Expr<SEMT::Literal<inv_sqrt_2pi_struct>> inv_sqrt_2pi_value;


double BgmgCalculator::calc_univariate_cost_cache_deriv(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int deriv_length, double* deriv) {
  // NB! censoring is not implemented in calc_univariate_cost_cache_deriv
  // NB! causalbetavec is not supported in calc_univariate_cost_cache_deriv
  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<float>& zvec(*get_zvec(trait_index));

  if (deriv_length != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("deriv_length != 3"));
  if (zvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec is not set"));
  if (nvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
  if (!cache_tag_r2sum_) BGMG_THROW_EXCEPTION(::std::runtime_error("bgmg_calc_univariate_cost_with_deriv only works with cache_tag_r2sum")); 

  float num_causals = pi_vec * static_cast<float>(num_snp_);
  if ((int)num_causals >= max_causals_) return 1e100; // too large pi_vec
  const int component_id = 0;   // univariate is always component 0.

  LOG << ">calc_univariate_cost_deriv(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";
  SimpleTimer timer(-1);

  find_tag_r2sum(component_id, num_causals);

  const float pi_k = 1. / static_cast<float>(k_max_);

  std::valarray<double> pdf_double(0.0, num_tag_);
  std::valarray<double> pdf_deriv_sig2zero(0.0, num_tag_);
  std::valarray<double> pdf_deriv_sig2beta(0.0, num_tag_);
  std::valarray<double> pdf_deriv_pivec(0.0, num_tag_);

#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, num_tag_);
    std::valarray<double> pdf_deriv_sig2zero_local(0.0, num_tag_);
    std::valarray<double> pdf_deriv_sig2beta_local(0.0, num_tag_);
    std::valarray<double> pdf_deriv_pivec_local(0.0, num_tag_);

#pragma omp for schedule(static)
    for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
      if (weights_[tag_index] == 0) continue;
      if (!std::isfinite(zvec[tag_index]) || !std::isfinite(nvec[tag_index])) continue;

      for (int k_index = 0; k_index < k_max_; k_index++) {
        float tag_r2sum = (*tag_r2sum_[component_id])(tag_index, k_index);

        DVAR(semt_sig2zero, 0); DVAR(semt_sig2beta, 1); DVAR(semt_pivec, 2); DVAR(zvec2, 3); DVAR(r2eff_times_n, 4);
        auto s = semt_sig2zero + r2eff_times_n * semt_sig2beta * semt_pivec;
        auto f = inv_sqrt_2pi_value * pow(s, RAT(-1, 2)) * exp(RAT(-1, 2) * zvec2 / s);
        auto f_sig2zero = SEMT::deriv_t(f, semt_sig2zero);
        auto f_sig2beta = SEMT::deriv_t(f, semt_sig2beta);
        auto f_pivec = SEMT::deriv_t(f, semt_pivec);
        SEMT::CAR semt_params = { sig2_zero, sig2_beta, pi_vec, zvec[tag_index] * zvec[tag_index], tag_r2sum * nvec[tag_index] / pi_vec };
        pdf_double_local[tag_index] += pi_k * f.apply(semt_params);
        pdf_deriv_sig2zero_local[tag_index] += pi_k * f_sig2zero.apply(semt_params);
        pdf_deriv_sig2beta_local[tag_index] += pi_k * f_sig2beta.apply(semt_params);
        pdf_deriv_pivec_local[tag_index] += pi_k * f_pivec.apply(semt_params);
      }
    }

#pragma omp critical
    {
      pdf_double += pdf_double_local;
      pdf_deriv_sig2zero += pdf_deriv_sig2zero_local;
      pdf_deriv_sig2beta += pdf_deriv_sig2beta_local;
      pdf_deriv_pivec += pdf_deriv_pivec_local;
    }
  }

  double log_pdf_total = 0.0;
  double& pi_vec_io = deriv[0]; pi_vec_io = 0;
  double& sig2_zero_io = deriv[1]; sig2_zero_io = 0;
  double& sig2_beta_io = deriv[2]; sig2_beta_io = 0;
  int num_infinite = 0;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec[tag_index]) || !std::isfinite(nvec[tag_index])) continue;
    double increment = -std::log(pdf_double[tag_index]) * weights_[tag_index];
    if (!std::isfinite(increment)) num_infinite++;
    log_pdf_total += increment;
    pi_vec_io += (-pdf_deriv_pivec[tag_index] / pdf_double[tag_index]) * weights_[tag_index];
    sig2_zero_io += (-pdf_deriv_sig2zero[tag_index] / pdf_double[tag_index]) * weights_[tag_index];
    sig2_beta_io += (-pdf_deriv_sig2beta[tag_index] / pdf_double[tag_index]) * weights_[tag_index];
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_univariate_cost_deriv(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

template<typename T>
double calc_univariate_cost_nocache_template(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, BgmgCalculator& rhs) {
  std::vector<float>& nvec(*rhs.get_nvec(trait_index));
  std::vector<float>& zvec(*rhs.get_zvec(trait_index));

  float num_causals = pi_vec * static_cast<float>(rhs.num_snp_);
  if ((int)num_causals >= rhs.max_causals_) return 1e100; // too large pi_vec
  const int component_id = 0;   // univariate is always component 0.
  if (rhs.snp_order_.empty()) rhs.find_snp_order();

  LOG << ">calc_univariate_cost_nocache(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";
  SimpleTimer timer(-1);

  const double pi_k = 1.0 / static_cast<double>(rhs.k_max_);

  rhs.k_pdf_.assign(rhs.k_max_, 0.0f);

  std::valarray<float> fixed_effect_delta(0.0, rhs.num_tag_);;
  rhs.calc_fixed_effect_delta_from_causalbetavec(trait_index, &fixed_effect_delta);

  std::valarray<double> pdf_double(0.0, rhs.num_tag_);
#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, rhs.num_tag_);
    std::vector<float> tag_r2sum(rhs.num_tag_, 0.0f);

#pragma omp for schedule(static)
      for (int k_index = 0; k_index < rhs.k_max_; k_index++) {

        rhs.find_tag_r2sum_no_cache(component_id, num_causals, k_index, &tag_r2sum);
        for (int tag_index = 0; tag_index < rhs.num_tag_; tag_index++) {
          if (rhs.weights_[tag_index] == 0) continue;
          if (!std::isfinite(zvec[tag_index]) || !std::isfinite(nvec[tag_index])) continue;

          float tag_r2sum_value = tag_r2sum[tag_index];
          float sig2eff = tag_r2sum_value * nvec[tag_index] * sig2_beta + sig2_zero;

          const float tag_z = zvec[tag_index] - fixed_effect_delta[tag_index];  // apply causalbetavec
          float s = sqrt(sig2eff);
          const bool censoring = std::abs(tag_z) > rhs.z1max_;
          double pdf = static_cast<double>(censoring ? censored_cdf<T>(rhs.z1max_, s) : gaussian_pdf<T>(tag_z, s));
          pdf_double_local[tag_index] += pdf * pi_k;

          if (rhs.calc_k_pdf_) rhs.k_pdf_[k_index] += (-std::log(pdf) * rhs.weights_[tag_index]);
        }
      }
#pragma omp critical
      pdf_double += pdf_double_local;
  }

  double log_pdf_total = 0.0;
  int num_infinite = 0;
  for (int tag_index = 0; tag_index < rhs.num_tag_; tag_index++) {
    if (rhs.weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec[tag_index])) continue;
    double increment = -std::log(pdf_double[tag_index]) * rhs.weights_[tag_index];
    if (!std::isfinite(increment)) num_infinite++;
    log_pdf_total += increment;
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_univariate_cost_nocache(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_univariate_cost_nocache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  return calc_univariate_cost_nocache_template<FLOAT_TYPE>(trait_index, pi_vec, sig2_zero, sig2_beta, *this);
}
double BgmgCalculator::calc_univariate_cost_nocache_float(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  return calc_univariate_cost_nocache_template<float>(trait_index, pi_vec, sig2_zero, sig2_beta, *this);
}
double BgmgCalculator::calc_univariate_cost_nocache_double(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  return calc_univariate_cost_nocache_template<double>(trait_index, pi_vec, sig2_zero, sig2_beta, *this);
}

std::string calc_bivariate_params_to_str(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, int length) {
  std::stringstream ss;
  ss << "pi_vec=[" << pi_vec[0] << ", " << pi_vec[1] << ", " << pi_vec[2] << "], "
     << "sig2_beta=[" << sig2_beta[0] << ", " << sig2_beta[1] << "], "
     << "rho_beta=" << rho_beta << ", "
     << "sig2_zero=[" << sig2_zero[0] << ", " << sig2_zero[1] << "], "
     << "rho_zero=" << rho_zero;
  if (length >= 0) ss << ", length(zvec)=" << length;
  return ss.str();
}

double BgmgCalculator::calc_bivariate_cost(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  if (zvec1_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec1 is not set"));
  if (nvec1_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec1 is not set"));
  if (zvec2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec2 is not set"));
  if (nvec2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec2 is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
  if (num_components_ != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: require num_components == 3. Remember to call set_option('num_components', 3)."));
  if (sig2_beta_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_beta_len != 2"));
  if (sig2_zero_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_zero_len != 2"));
  if (pi_vec_len != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: pi_vec_len != 3"));

  double cost;
  if (cost_calculator_==CostCalculator_Gaussian) cost = calc_bivariate_cost_fast(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);
  else if (cost_calculator_==CostCalculator_Convolve) cost = calc_bivariate_cost_convolve(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);
  else if (!cache_tag_r2sum_) cost = calc_bivariate_cost_nocache(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);
  else cost = calc_bivariate_cost_cache(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);

  loglike_cache_.add_entry(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, cost);
  return cost;
}

double BgmgCalculator::calc_bivariate_cost_cache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {

  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost(" << ss << ")";
  SimpleTimer timer(-1);

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) return 1e100; // too large pi_vec
  }

  for (int component_id = 0; component_id < 3; component_id++) {
    find_tag_r2sum(component_id, num_causals[component_id]);
  }

  std::valarray<float> fixed_effect_delta1(0.0, num_tag_), fixed_effect_delta2(0.0, num_tag_);
  calc_fixed_effect_delta_from_causalbetavec(1, &fixed_effect_delta1);
  calc_fixed_effect_delta_from_causalbetavec(2, &fixed_effect_delta2);

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  // pi_k is mixture weight
  const double pi_k = 1.0 / static_cast<double>(k_max_);

  double log_pdf_total = 0.0;
  int num_infinite = 0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_infinite)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec1_[tag_index]) || !std::isfinite(nvec1_[tag_index])) continue;
    if (!std::isfinite(zvec2_[tag_index]) || !std::isfinite(nvec2_[tag_index])) continue;

    const float z1 = zvec1_[tag_index] - fixed_effect_delta1[tag_index];  // apply causalbetavec;
    const float z2 = zvec2_[tag_index] - fixed_effect_delta2[tag_index];  // apply causalbetavec;
    const float n1 = nvec1_[tag_index];
    const float n2 = nvec2_[tag_index];

    double pdf_tag = 0.0f;
    for (int k_index = 0; k_index < k_max_; k_index++) {
      const float tag_r2sum_c1 = (*tag_r2sum_[0])(tag_index, k_index);
      const float tag_r2sum_c2 = (*tag_r2sum_[1])(tag_index, k_index);
      const float tag_r2sum_c3 = (*tag_r2sum_[2])(tag_index, k_index);

      // Sigma  = [A1+A3  B3;  B3  C2+C3] + Sigma0 = ...
      //        = [a11    a12; a12   a22]
      const float A1 = tag_r2sum_c1 * n1 * sig2_beta[0];
      const float C2 = tag_r2sum_c2 * n2 * sig2_beta[1];
      const float A3 = tag_r2sum_c3 * n1 * sig2_beta[0];
      const float C3 = tag_r2sum_c3 * n2 * sig2_beta[1];
      const float B3 = sqrt(A3*C3) * rho_beta;

      const float a11 = A1 + A3 + a0;
      const float a22 = C2 + C3 + c0;
      const float a12 =      B3 + b0;

      const bool censoring = (std::abs(z1) > z1max_) || (std::abs(z2) > z2max_);
      const double pdf = static_cast<double>(censoring ? censored2_cdf<FLOAT_TYPE>(z1max_, z2max_, a11, a12, a22) : gaussian2_pdf<FLOAT_TYPE>(z1, z2, a11, a12, a22));
      pdf_tag += pi_k * pdf;
    }

    double increment = static_cast<double>(-std::log(pdf_tag) * weights_[tag_index]);
    if (!std::isfinite(increment)) num_infinite++;
    log_pdf_total += increment;
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_bivariate_cost(" << ss << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_bivariate_cost_nocache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost_nocache(" << ss << ")";
  SimpleTimer timer(-1);

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) return 1e100; // too large pi_vec
  }

  if (snp_order_.empty()) find_snp_order();

  std::valarray<float> fixed_effect_delta1(0.0, num_tag_), fixed_effect_delta2(0.0, num_tag_);
  calc_fixed_effect_delta_from_causalbetavec(1, &fixed_effect_delta1);
  calc_fixed_effect_delta_from_causalbetavec(2, &fixed_effect_delta2);

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  // pi_k is mixture weight
  const double pi_k = 1.0 / static_cast<double>(k_max_);

  std::valarray<double> pdf_double(0.0, num_tag_);
#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, num_tag_);
    std::vector<float> tag_r2sum0(num_tag_, 0.0f);
    std::vector<float> tag_r2sum1(num_tag_, 0.0f);
    std::vector<float> tag_r2sum2(num_tag_, 0.0f);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {

      find_tag_r2sum_no_cache(0, num_causals[0], k_index, &tag_r2sum0);
      find_tag_r2sum_no_cache(1, num_causals[1], k_index, &tag_r2sum1);
      find_tag_r2sum_no_cache(2, num_causals[2], k_index, &tag_r2sum2);

      for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
        if (weights_[tag_index] == 0) continue;
        if (!std::isfinite(zvec1_[tag_index]) || !std::isfinite(nvec1_[tag_index])) continue;
        if (!std::isfinite(zvec2_[tag_index]) || !std::isfinite(nvec2_[tag_index])) continue;

        const float z1 = zvec1_[tag_index] - fixed_effect_delta1[tag_index];  // apply causalbetavec;
        const float z2 = zvec2_[tag_index] - fixed_effect_delta2[tag_index];  // apply causalbetavec;
        const float n1 = nvec1_[tag_index];
        const float n2 = nvec2_[tag_index];

        float pdf_tag = 0.0f;

        const float tag_r2sum_c1 = tag_r2sum0[tag_index];
        const float tag_r2sum_c2 = tag_r2sum1[tag_index];
        const float tag_r2sum_c3 = tag_r2sum2[tag_index];

        // Sigma  = [A1+A3  B3;  B3  C2+C3] + Sigma0 = ...
        //        = [a11    a12; a12   a22]
        const float A1 = tag_r2sum_c1 * n1 * sig2_beta[0];
        const float C2 = tag_r2sum_c2 * n2 * sig2_beta[1];
        const float A3 = tag_r2sum_c3 * n1 * sig2_beta[0];
        const float C3 = tag_r2sum_c3 * n2 * sig2_beta[1];
        const float B3 = sqrt(A3*C3) * rho_beta;

        const float a11 = A1 + A3 + a0;
        const float a22 = C2 + C3 + c0;
        const float a12 = B3 + b0;

        const bool censoring = (std::abs(z1) > z1max_) || (std::abs(z2) > z2max_);
        const double pdf = static_cast<double>(censoring ? censored2_cdf<FLOAT_TYPE>(z1max_, z2max_, a11, a12, a22) : gaussian2_pdf<FLOAT_TYPE>(z1, z2, a11, a12, a22));
        pdf_double_local[tag_index] += pi_k * pdf;
      }
    }
#pragma omp critical
    pdf_double += pdf_double_local;
  }

  double log_pdf_total = 0.0;
  int num_infinite = 0;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec1_[tag_index]) || !std::isfinite(nvec1_[tag_index])) continue;
    if (!std::isfinite(zvec2_[tag_index]) || !std::isfinite(nvec2_[tag_index])) continue;
    double increment = -std::log(pdf_double[tag_index]) * weights_[tag_index];
    if (!std::isfinite(increment)) num_infinite++;
    log_pdf_total += increment;
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_bivariate_cost_nocache(" << ss << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

int64_t BgmgCalculator::calc_bivariate_pdf(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, int length, float* zvec1, float* zvec2, float* pdf) {
  // input buffer "zvec1" and "zvec2" contains z scores (presumably an equally spaced grid)
  // output buffer contains pdf(z), aggregated across all SNPs with corresponding weights
  if (nvec1_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec1 is not set"));
  if (nvec2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec2 is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
  if (num_components_ != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_pdf: require num_components == 3. Remember to call set_option('num_components', 3)."));
  if (sig2_beta_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_beta_len != 2"));
  if (sig2_zero_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_zero_len != 2"));
  if (pi_vec_len != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: pi_vec_len != 3"));

  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, length);
  LOG << ">calc_bivariate_pdf(" << ss << ")";
  SimpleTimer timer(-1);

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  }

  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) {
    for (int component_id = 0; component_id < 3; component_id++) {
      find_tag_r2sum(component_id, num_causals[component_id]);
    }
  }

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  // pi_k is mixture weight
  const double pi_k = 1.0 / static_cast<double > (k_max_);

  std::valarray<double> pdf_double(0.0, length);

#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, length);
    std::vector<float> tag_r2sum0(num_tag_, 0.0f);
    std::vector<float> tag_r2sum1(num_tag_, 0.0f);
    std::vector<float> tag_r2sum2(num_tag_, 0.0f);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {
      if (cache_tag_r2sum_) {
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          tag_r2sum0[tag_index] = (*tag_r2sum_[0])(tag_index, k_index);
          tag_r2sum1[tag_index] = (*tag_r2sum_[1])(tag_index, k_index);
          tag_r2sum2[tag_index] = (*tag_r2sum_[2])(tag_index, k_index);
        }
      } else {
        find_tag_r2sum_no_cache(0, num_causals[0], k_index, &tag_r2sum0);
        find_tag_r2sum_no_cache(1, num_causals[1], k_index, &tag_r2sum1);
        find_tag_r2sum_no_cache(2, num_causals[2], k_index, &tag_r2sum2);
      }

      for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
        if (weights_[tag_index] == 0) continue;
        if (!std::isfinite(nvec1_[tag_index])) continue;
        if (!std::isfinite(nvec2_[tag_index])) continue;

        double tag_weight = static_cast<double>(weights_[tag_index]);

        const float n1 = nvec1_[tag_index];
        const float n2 = nvec2_[tag_index];

        const float tag_r2sum_c1 = tag_r2sum0[tag_index];
        const float tag_r2sum_c2 = tag_r2sum1[tag_index];
        const float tag_r2sum_c3 = tag_r2sum2[tag_index];

        // Sigma  = [A1+A3  B3;  B3  C2+C3] + Sigma0 = ...
        //        = [a11    a12; a12   a22]
        const float A1 = tag_r2sum_c1 * n1 * sig2_beta[0];
        const float C2 = tag_r2sum_c2 * n2 * sig2_beta[1];
        const float A3 = tag_r2sum_c3 * n1 * sig2_beta[0];
        const float C3 = tag_r2sum_c3 * n2 * sig2_beta[1];
        const float B3 = sqrt(A3*C3) * rho_beta;

        const float a11 = A1 + A3 + a0;
        const float a22 = C2 + C3 + c0;
        const float a12 = B3 + b0;
        
        for (int z_index = 0; z_index < length; z_index++) {
          double pdf_tmp = static_cast<double>(gaussian2_pdf<FLOAT_TYPE>(zvec1[z_index], zvec2[z_index], a11, a12, a22));
          pdf_double_local[z_index] += pi_k * pdf_tmp * tag_weight;
        }
      }
    }
#pragma omp critical
    pdf_double += pdf_double_local;
  }

  for (int i = 0; i < length; i++) pdf[i] = static_cast<float>(pdf_double[i]);
  LOG << "<calc_bivariate_pdf(" << ss << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

void calc_bivariate_delta_posterior_integrals(float a, float b, float c, float i, float j, float k, float z1, float z2,
                                              float* c00, float* c10, float* c01, float* c20, float* c11, float* c02) {
  // [a b; b c] is variance-covariance matrix of eps in "z=delta+eps"
  // [i j; j k] is variance-covariance matrix of delta for a specific choice of causal variants in sampling
  // here we are calculating moments cPQ = E[delta1^P delta2^Q | z1, z2].
  // the formulas are pretty tricky to derive - loads of integrals involved.

  static const float inv_sqrt2_pi = 0.2250790790392765f;
  const float z12 = z1*z1;
  const float z22 = z2*z2;
  const float z1z2 = z1*z2;

  const float ci = c*i;
  const float bj = b*j;
  const float j2 = j*j;
  const float ik = i*k;
  const float bi = b*i;
  const float aj = a*j;
  const float cj = c*j;
  const float bk = b*k;
  const float ak = a*k;
  const float ij = i*j;
  const float i2 = i*i;
  const float k2 = k*k;
  const float jk = j*k;

  const float ci2 = ci*i;
  const float bij = bi*j;
  const float aj2 = aj*j;
  const float ij2 = ij*j;
  const float i2k = i2*k;
  const float cj2 = cj*j;
  const float bjk = bj*k;
  const float j2k = j2*k;
  const float ak2 = ak*k;
  const float ik2 = ik*k;
  const float cij = ci*j;
  const float bj2 = bj*j;
  const float bik = bi*k;
  const float ajk = aj*k;
  const float ijk = ij*k;

  const float j3 = j2*j;

  const float eN = (c + k) * z12 - 2 * (b + j) * z1z2 + (a + i) * z22;
  const float eD = -2.0f * (b + j)*(b + j) + 2.0f * (a + i) * (c + k);

  const float eD2 = eD*eD;
  const float inv_eD = 1.0f / eD;
  const float inv_eD2 = inv_eD*inv_eD;
  const float x2eD = 2.0f*eD;
  const float x2eD_x4eN = x2eD - 4.0f*eN;

  const float c00_tmp = inv_sqrt2_pi * exp(-eN * inv_eD) * sqrt(inv_eD);

  const float c10_pol = ci * z1 - bj * z1 - j2 * z1 + ik * z1 - bi * z2 + aj * z2;
  const float c01_pol = cj * z1 - bk * z1 - bj * z2 - j2 * z2 + ak * z2 + ik * z2;
  const float c20_pol = eD2*i + x2eD_x4eN*(bij + bij - ci2 - aj2 + ij2 - i2k) - x2eD*(j2*z12 - 2.0f*ij*z1z2           + i2*z22);
  const float c02_pol = eD2*k + x2eD_x4eN*(bjk + bjk - cj2 - ak2 + j2k - ik2) - x2eD*(k2*z12 - 2.0f*jk*z1z2           + j2*z22);
  const float c11_pol = eD2*j + x2eD_x4eN*(bj2 + j3  - cij - ajk + bik - ijk) - x2eD*(jk*z12 -      j2*z1z2 - ik*z1z2 + ij*z22);

  (*c00) = c00_tmp;
  (*c10) = c00_tmp * 2.0f * inv_eD  * c10_pol;
  (*c01) = c00_tmp * 2.0f * inv_eD  * c01_pol;
  (*c20) = c00_tmp *        inv_eD2 * c20_pol;
  (*c11) = c00_tmp *        inv_eD2 * c11_pol;
  (*c02) = c00_tmp *        inv_eD2 * c02_pol;  
}

int64_t BgmgCalculator::calc_bivariate_delta_posterior(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero,
                                                       int length, float* c00, float* c10, float* c01, float* c20, float* c11, float* c02) {
  // where c(i,j) = \int_{\delta1, \delta2} \delta1^i \delta2^j P(z1, z2 | delta1, delta2) P(delta1, delta2)
  if (zvec1_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec1 is not set"));
  if (nvec1_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec1 is not set"));
  if (zvec2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec2 is not set"));
  if (nvec2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec2 is not set"));
  if (num_components_ != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: require num_components == 3. Remember to call set_option('num_components', 3)."));
  if (sig2_beta_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_beta_len != 2"));
  if (sig2_zero_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_zero_len != 2"));
  if (pi_vec_len != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: pi_vec_len != 3"));
  if ((length == 0) || (length != num_tag_)) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_tag_"));

  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, length);
  LOG << ">calc_bivariate_delta_posterior(" << ss << ")";
  SimpleTimer timer(-1);

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  }

  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) {
    for (int component_id = 0; component_id < 3; component_id++) {
      find_tag_r2sum(component_id, num_causals[component_id]);
    }
  }

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  std::valarray<double> c00_global(0.0f, num_tag_);
  std::valarray<double> c10_global(0.0f, num_tag_);
  std::valarray<double> c01_global(0.0f, num_tag_);
  std::valarray<double> c20_global(0.0f, num_tag_);
  std::valarray<double> c11_global(0.0f, num_tag_);
  std::valarray<double> c02_global(0.0f, num_tag_);

#pragma omp parallel
  {
    std::vector<float> tag_r2sum0(num_tag_, 0.0f);
    std::vector<float> tag_r2sum1(num_tag_, 0.0f);
    std::vector<float> tag_r2sum2(num_tag_, 0.0f);
    std::valarray<double> c00_local(0.0f, num_tag_);
    std::valarray<double> c10_local(0.0f, num_tag_);
    std::valarray<double> c01_local(0.0f, num_tag_);
    std::valarray<double> c20_local(0.0f, num_tag_);
    std::valarray<double> c11_local(0.0f, num_tag_);
    std::valarray<double> c02_local(0.0f, num_tag_);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {
      if (cache_tag_r2sum_) {
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          tag_r2sum0[tag_index] = (*tag_r2sum_[0])(tag_index, k_index);
          tag_r2sum1[tag_index] = (*tag_r2sum_[1])(tag_index, k_index);
          tag_r2sum2[tag_index] = (*tag_r2sum_[2])(tag_index, k_index);
        }
      } else {
        find_tag_r2sum_no_cache(0, num_causals[0], k_index, &tag_r2sum0);
        find_tag_r2sum_no_cache(1, num_causals[1], k_index, &tag_r2sum1);
        find_tag_r2sum_no_cache(2, num_causals[2], k_index, &tag_r2sum2);
      }

      for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
        if (!std::isfinite(zvec1_[tag_index]) || !std::isfinite(nvec1_[tag_index])) continue;
        if (!std::isfinite(zvec2_[tag_index]) || !std::isfinite(nvec2_[tag_index])) continue;

        const float n1 = nvec1_[tag_index];
        const float n2 = nvec2_[tag_index];
        const float z1 = zvec1_[tag_index];
        const float z2 = zvec2_[tag_index];

        const float tag_r2sum_c1 = tag_r2sum0[tag_index];
        const float tag_r2sum_c2 = tag_r2sum1[tag_index];
        const float tag_r2sum_c3 = tag_r2sum2[tag_index];

        // Sigma  = [A1+A3  B3;  B3  C2+C3] + Sigma0 = ...
        //        = [a11    a12; a12   a22]
        const float A1 = tag_r2sum_c1 * n1 * sig2_beta[0];
        const float C2 = tag_r2sum_c2 * n2 * sig2_beta[1];
        const float A3 = tag_r2sum_c3 * n1 * sig2_beta[0];
        const float C3 = tag_r2sum_c3 * n2 * sig2_beta[1];
        const float B3 = sqrt(A3*C3) * rho_beta;

        const float A = A1 + A3;
        const float C = C2 + C3;
        const float B = B3;

        float c00buf, c10buf, c01buf, c20buf, c11buf, c02buf;
        calc_bivariate_delta_posterior_integrals(a0, b0, c0, A, B, C, z1, z2, &c00buf, &c10buf, &c01buf, &c20buf, &c11buf, &c02buf);
        c00_local[tag_index] += static_cast<double>(c00buf);
        c10_local[tag_index] += static_cast<double>(c10buf);
        c01_local[tag_index] += static_cast<double>(c01buf);
        c20_local[tag_index] += static_cast<double>(c20buf);
        c11_local[tag_index] += static_cast<double>(c11buf);
        c02_local[tag_index] += static_cast<double>(c02buf);
      }
    }
#pragma omp critical
    {
      c00_global += c00_local;
      c10_global += c10_local;
      c01_global += c01_local;
      c20_global += c20_local;
      c11_global += c11_local;
      c02_global += c02_local;
    }
  }

  // save results to output buffers
  const double pi_k = 1.0 / static_cast<double>(k_max_);
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (!std::isfinite(zvec1_[tag_index]) || !std::isfinite(nvec1_[tag_index])) continue;
    if (!std::isfinite(zvec2_[tag_index]) || !std::isfinite(nvec2_[tag_index])) continue;
    c00[tag_index] = pi_k * c00_global[tag_index];
    c10[tag_index] = pi_k * c10_global[tag_index];
    c01[tag_index] = pi_k * c01_global[tag_index];    
    c20[tag_index] = pi_k * c20_global[tag_index];
    c11[tag_index] = pi_k * c11_global[tag_index];
      c02[tag_index] = pi_k * c02_global[tag_index];    
  }

  LOG << "<calc_bivariate_delta_posterior(" << ss << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

template<typename T>
std::string std_vector_to_str(const std::vector<T>& vec) {
  std::stringstream ss;
  int max_el = std::min<int>(5, vec.size() - 1);
  ss << "[";
  for (int i = 0; i < max_el; i++) {
    bool last = (i == (max_el - 1));
    ss << vec[i];
    if (last) ss << ", ...";
    else ss << ", ";
  }
  ss << "]";

  size_t nnz = 0;
  for (size_t i = 0; i < vec.size(); i++) if (vec[i] != 0) nnz++;
  ss << ", nnz=" << nnz;
  return ss.str();
}

void BgmgCalculator::log_diagnostics() {
  size_t mem_bytes = 0, mem_bytes_total = 0;
  LOG << " diag: num_snp_=" << num_snp_;
  LOG << " diag: num_tag_=" << num_tag_;
  mem_bytes_total += ld_matrix_csr_.log_diagnostics();
  LOG << " diag: zvec1_.size()=" << zvec1_.size();
  LOG << " diag: zvec1_=" << std_vector_to_str(zvec1_);
  LOG << " diag: nvec1_.size()=" << nvec1_.size();
  LOG << " diag: nvec1_=" << std_vector_to_str(nvec1_);
  LOG << " diag: causalbetavec1_.size()=" << causalbetavec1_.size();
  LOG << " diag: causalbetavec1_=" << std_vector_to_str(causalbetavec1_);
  LOG << " diag: zvec2_.size()=" << zvec2_.size();
  LOG << " diag: zvec2_=" << std_vector_to_str(zvec2_);
  LOG << " diag: nvec2_.size()=" << nvec2_.size();
  LOG << " diag: nvec2_=" << std_vector_to_str(nvec2_);
  LOG << " diag: causalbetavec2_.size()=" << causalbetavec2_.size();
  LOG << " diag: causalbetavec2_=" << std_vector_to_str(causalbetavec2_);
  LOG << " diag: weights_.size()=" << weights_.size();
  LOG << " diag: weights_=" << std_vector_to_str(weights_);
  LOG << " diag: mafvec_.size()=" << mafvec_.size();
  LOG << " diag: mafvec_=" << std_vector_to_str(mafvec_);
  for (int i = 0; i < snp_order_.size(); i++) {
    mem_bytes = snp_order_[i]->size() * sizeof(int); mem_bytes_total += mem_bytes;
    LOG << " diag: snp_order_[" << i << "].shape=[" << snp_order_[i]->no_rows() << ", " << snp_order_[i]->no_columns() << "]" << " (mem usage = " << mem_bytes << " bytes)";
    LOG << " diag: snp_order_[" << i << "]=" << snp_order_[i]->to_str();
  }
  for (int i = 0; i < tag_r2sum_.size(); i++) {
    mem_bytes = tag_r2sum_[i]->size() * sizeof(float); mem_bytes_total += mem_bytes;
    LOG << " diag: tag_r2sum_[" << i << "].shape=[" << tag_r2sum_[i]->no_rows() << ", " << tag_r2sum_[i]->no_columns() << "]" << " (mem usage = " << mem_bytes << " bytes)";
    LOG << " diag: tag_r2sum_[" << i << "]=" << tag_r2sum_[i]->to_str();
  }
  for (int i = 0; i < last_num_causals_.size(); i++) 
    LOG << " diag: last_num_causals_[" << i << "]=" << last_num_causals_[i];
  LOG << " diag: options.k_max_=" << k_max_;
  LOG << " diag: options.use_complete_tag_indices_=" << use_complete_tag_indices_;
  LOG << " diag: options.max_causals_=" << max_causals_;
  LOG << " diag: options.num_components_=" << num_components_;
  LOG << " diag: options.r2_min_=" << r2_min_;
  LOG << " diag: options.z1max_=" << z1max_;
  LOG << " diag: options.z2max_=" << z2max_;
  LOG << " diag: options.cost_calculator_=" << ((int)cost_calculator_) <<
    ((cost_calculator_==CostCalculator_Sampling) ? " (Sampling)" :
     (cost_calculator_==CostCalculator_Gaussian) ? " (Gaussian)" :
     (cost_calculator_==CostCalculator_Convolve) ? " (Convolve)" : " (Unknown)");
  LOG << " diag: options.cache_tag_r2sum_=" << (cache_tag_r2sum_ ? "yes" : "no");
  LOG << " diag: options.seed_=" << (seed_);
  LOG << " diag: options.cubature_abs_error_=" << (cubature_abs_error_);
  LOG << " diag: options.cubature_rel_error_=" << (cubature_rel_error_);
  LOG << " diag: options.cubature_max_evals_=" << (cubature_max_evals_);
  LOG << " diag: options.calc_k_pdf_=" << (calc_k_pdf_);
  LOG << " diag: options.ld_format_version_=" << (ld_format_version_);
  LOG << " diag: Estimated memory usage (total): " << mem_bytes_total << " bytes";
}

double BgmgCalculator::calc_univariate_cost_fast(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<float>& zvec(*get_zvec(trait_index));

  // Use an approximation that preserves variance and kurtosis.
  // This gives a robust cost function that scales up to a very high pivec, including infinitesimal model pi==1.

  std::stringstream ss;
  ss << "calc_univariate_cost_fast(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";
  LOG << ">" << ss.str();

  double log_pdf_total = 0.0;
  SimpleTimer timer(-1);

  std::valarray<float> fixed_effect_delta(0.0, num_tag_);
  calc_fixed_effect_delta_from_causalbetavec(trait_index, &fixed_effect_delta);

  int num_zero_tag_r2 = 0;
  int num_infinite = 0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_zero_tag_r2, num_infinite)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec[tag_index]) || !std::isfinite(nvec[tag_index])) continue;
    double tag_weight = static_cast<double>(weights_[tag_index]);
    
    const float tag_r2 = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2()[tag_index];
    const float tag_r4 = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r4()[tag_index];

    if (tag_r2 == 0 || tag_r4 == 0) {
      num_zero_tag_r2++; continue;
    }

    const float tag_chi = tag_r4 / tag_r2;

    const float tag_eta_factor = pi_vec * tag_r2 + (1.0f - pi_vec) * tag_chi;
    const double tag_pi1 = static_cast<double>(pi_vec * tag_r2 / tag_eta_factor);
    const double tag_pi0 = 1.0 - tag_pi1;
    const float tag_sig2beta = sig2_beta * tag_eta_factor;

    const float tag_z = zvec[tag_index] - fixed_effect_delta[tag_index];  // apply causalbetavec;
    const float tag_n = nvec[tag_index];

    const bool censoring = std::abs(tag_z) > z1max_;
    const float s1 = sqrt(sig2_zero);
    const float s2 = sqrt(sig2_zero + tag_n *tag_sig2beta);

    const double tag_pdf0 = static_cast<double>(censoring ? censored_cdf<FLOAT_TYPE>(z1max_, s1) : gaussian_pdf<FLOAT_TYPE>(tag_z, s1));
    const double tag_pdf1 = static_cast<double>(censoring ? censored_cdf<FLOAT_TYPE>(z1max_, s2) : gaussian_pdf<FLOAT_TYPE>(tag_z, s2));
    const double tag_pdf = tag_pi0 * tag_pdf0 + tag_pi1 * tag_pdf1;
    const double increment = (-std::log(tag_pdf) * tag_weight);
    if (!std::isfinite(increment)) num_infinite++;
    log_pdf_total += increment;
  }

  if (num_zero_tag_r2 > 0)
    LOG << " warning: zero tag_r2 encountered " << num_zero_tag_r2 << " times";
  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<" << ss.str() << ", cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_bivariate_cost_fast(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost_fast(" << ss << ")";

  double log_pdf_total = 0.0;
  SimpleTimer timer(-1);

  std::valarray<float> fixed_effect_delta1(0.0, num_tag_), fixed_effect_delta2(0.0, num_tag_);
  calc_fixed_effect_delta_from_causalbetavec(1, &fixed_effect_delta1);
  calc_fixed_effect_delta_from_causalbetavec(2, &fixed_effect_delta2);

  int num_zero_tag_r2 = 0;
  int num_infinite = 0;

  const float s0_a11 = sig2_zero[0];
  const float s0_a22 = sig2_zero[1];
  const float s0_a12 = sqrt(sig2_zero[0] * sig2_zero[1]) * rho_zero;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_zero_tag_r2, num_infinite)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec1_[tag_index]) || !std::isfinite(nvec1_[tag_index])) continue;
    if (!std::isfinite(zvec2_[tag_index]) || !std::isfinite(nvec2_[tag_index])) continue;

    const float z1 = zvec1_[tag_index] - fixed_effect_delta1[tag_index];
    const float n1 = nvec1_[tag_index];
    const float z2 = zvec2_[tag_index] - fixed_effect_delta2[tag_index];
    const float n2 = nvec2_[tag_index];

    const float tag_r2 = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2()[tag_index];
    const float tag_r4 = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r4()[tag_index];

    if (tag_r2 == 0 || tag_r4 == 0) {
      num_zero_tag_r2++; continue;
    }

    const float tag_chi = tag_r4 / tag_r2;

    const float tag_eta_factor[3] = {
      pi_vec[0] * tag_r2 + (1.0f - pi_vec[0]) * tag_chi,
      pi_vec[1] * tag_r2 + (1.0f - pi_vec[1]) * tag_chi,
      pi_vec[2] * tag_r2 + (1.0f - pi_vec[2]) * tag_chi
    };

    const float tag_pi1[3] = {
      pi_vec[0] * tag_r2 / tag_eta_factor[0],
      pi_vec[1] * tag_r2 / tag_eta_factor[1],
      pi_vec[2] * tag_r2 / tag_eta_factor[2]
    };

    const float tag_pi0[3] = {
      1.0f - tag_pi1[0],
      1.0f - tag_pi1[1],
      1.0f - tag_pi1[2]
    };

    const float a11[3] = { tag_eta_factor[0] * n1 * sig2_beta[0], 0,                                     tag_eta_factor[2] * n1 * sig2_beta[0] };
    const float a22[3] = { 0,                                     tag_eta_factor[1] * n2 * sig2_beta[1], tag_eta_factor[2] * n2 * sig2_beta[1] };
    const float a12[3] = { 0,                                     0,                                     rho_beta * sqrt(a11[2] * a22[2]) };

    const float f0[8] = { 0,0,0,0,1,1,1,1 };
    const float f1[8] = { 0,0,1,1,0,0,1,1 };
    const float f2[8] = { 0,1,0,1,0,1,0,1 };

    const bool censoring = (std::abs(z1) > z1max_) || (std::abs(z2) > z2max_);

    FLOAT_TYPE tag_pdf = 0.0f;
    for (int i = 0; i < 8; i++) {
      const float pi1 = (f0[i] ? tag_pi1[0] : tag_pi0[0]);
      const float pi2 = (f1[i] ? tag_pi1[1] : tag_pi0[1]);
      const float pi3 = (f2[i] ? tag_pi1[2] : tag_pi0[2]);
      const float a11i = s0_a11 + f0[i] * a11[0] + f1[i] * a11[1] + f2[i] * a11[2];
      const float a22i = s0_a22 + f0[i] * a22[0] + f1[i] * a22[1] + f2[i] * a22[2];
      const float a12i = s0_a12 + f0[i] * a12[0] + f1[i] * a12[1] + f2[i] * a12[2];
      tag_pdf += static_cast<double>(pi1*pi2*pi3) * static_cast<double>(censoring ? censored2_cdf<FLOAT_TYPE>(z1max_, z2max_, a11i, a12i, a22i) : gaussian2_pdf<FLOAT_TYPE>(z1, z2, a11i, a12i, a22i));
    }

    if (tag_pdf <= 0)
      tag_pdf = 1e-100;

    double increment = static_cast<double>(-std::log(tag_pdf) * weights_[tag_index]);
    if (!std::isfinite(increment)) num_infinite++;
    log_pdf_total += increment;
  }

  if (num_zero_tag_r2 > 0)
    LOG << " warning: zero tag_r2 encountered " << num_zero_tag_r2 << " times";
  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_bivariate_cost_fast(" << ss << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

class UnivariateCharacteristicFunctionData {
 public:
  float pi_vec;
  float sig2_zero;
  float sig2_beta;
  int tag_index;  // for which SNP to calculate the characteristic function
                  // note that use_complete_tag_indices_ must be enabled for "convolve" calculator, 
                  // so "tag" and "snp" are in the same indexing
  LdMatrixRow* ld_matrix_row;
  const std::vector<float>* hvec;
  const std::vector<float>* zvec;
  const std::vector<float>* nvec;
  const std::valarray<float>* fixed_effect_delta;
  const std::vector<float>* ld_tag_sum_r2_below_r2min_adjust_for_hvec;
  int func_evals;
};

int calc_univariate_characteristic_function_times_cosinus(unsigned ndim, const double *x, void *raw_data, unsigned fdim, double* fval) {
  assert(ndim == 1);
  assert(fdim == 1);
  const float t = (float)x[0];
  const float minus_tsqr_half = -t*t/2.0;
  UnivariateCharacteristicFunctionData* data = (UnivariateCharacteristicFunctionData *)raw_data;
  const float pi1 = data->pi_vec;
  const float pi0 = 1.0 - pi1;
  const float sig2beta_times_nval = (*data->nvec)[data->tag_index] * data->sig2_beta;
  const float zval = (*data->zvec)[data->tag_index] - (*data->fixed_effect_delta)[data->tag_index];  // apply causalbetavec
  
    // apply infinitesimal model to adjust tag_r2sum for all r2 that are below r2min (and thus do not contribute via resampling)
  const float inf_adj = pi1 * (*data->ld_tag_sum_r2_below_r2min_adjust_for_hvec)[data->tag_index] * sig2beta_times_nval;
  const float m_1_pi = M_1_PI;

  double result = m_1_pi * cos(t * zval) * std::exp(minus_tsqr_half * (data->sig2_zero + inf_adj));
  
  auto iter_end = data->ld_matrix_row->end();
  for (auto iter = data->ld_matrix_row->begin(); iter < iter_end; iter++) {
    int snp_index = iter.tag_index();  // yes, this is correct - snp_index on the LHS, tag_index on the RHS.
                                       // ld_matrix was designed to work with "sampling" calculator which require 
                                       // a mapping from a causal SNP "i" to all tag SNPs "j" that are in LD with "i".
                                       // however, for for "convolve" we are interested in mapping from a tag SNP "j"
                                       // to all causal SNPs "i" in LD with "j". To make this work we store a complete
                                       // LD matrix (e.g. _num_snp == _num_tag), and explore symmetry of the matrix. 
    float r2 = iter.r2();
    float hval = (*data->hvec)[snp_index];
    result *= (double)(pi0 + pi1 * std::exp(minus_tsqr_half * sig2beta_times_nval * r2 * hval));
  }

  data->func_evals++;
  *fval = result;
  return 0; 
}

int calc_univariate_characteristic_function_for_integration(unsigned ndim, const double *x, void *raw_data, unsigned fdim, double* fval) {
  const double t = x[0];
  const double inv_1_minus_t = 1.0 / (1.0 - t);
  const double x_transformed = t * inv_1_minus_t;
  const double jacob = inv_1_minus_t * inv_1_minus_t;
  int retval = calc_univariate_characteristic_function_times_cosinus(ndim, &x_transformed, raw_data, fdim, fval);
  (*fval) *= jacob;
  return retval;
}

double BgmgCalculator::calc_univariate_cost_convolve(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  if (!use_complete_tag_indices_) BGMG_THROW_EXCEPTION(::std::runtime_error("Convolve calculator require 'use_complete_tag_indices' option"));

  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<float>& zvec(*get_zvec(trait_index));

  std::stringstream ss;
  ss << "trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta;
  LOG << ">calc_univariate_cost_convolve(" << ss.str() << ")";

  double log_pdf_total = 0.0;
  int num_snp_failed = 0;
  int num_infinite = 0;
  double func_evals = 0.0;
  SimpleTimer timer(-1);

  std::valarray<float> fixed_effect_delta(0.0, num_tag_);
  calc_fixed_effect_delta_from_causalbetavec(trait_index, &fixed_effect_delta);

  std::vector<float> hvec;
  find_hvec(*this, &hvec);

  // in use_complete_tag_indices_ many tag indices are undefined
  // we compute deftag to avoid unbalanced load across OpenMP threads (static scheduler is still the way to go here.)
  std::vector<int> deftag_indices;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec[tag_index]) || !std::isfinite(nvec[tag_index])) continue;
    deftag_indices.push_back(tag_index);
  }
  const int num_deftag = deftag_indices.size();

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    UnivariateCharacteristicFunctionData data;
    data.pi_vec = pi_vec;
    data.sig2_zero = sig2_zero;
    data.sig2_beta = sig2_beta;
    data.ld_matrix_row = &ld_matrix_row;
    data.hvec = &hvec;
    data.zvec = &zvec;
    data.nvec = &nvec;
    data.fixed_effect_delta = &fixed_effect_delta;
    data.ld_tag_sum_r2_below_r2min_adjust_for_hvec = &ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2(LD_TAG_COMPONENT_BELOW_R2MIN);

#pragma omp for schedule(static) reduction(+: log_pdf_total, num_snp_failed, num_infinite, func_evals)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      double tag_weight = static_cast<double>(weights_[tag_index]);

      const int causal_index = tag_index; // yes, causal==tag in this case -- see a large comment in calc_univariate_characteristic_function_times_cosinus function.
      ld_matrix_csr_.extract_row(causal_index, data.ld_matrix_row);
      data.tag_index = tag_index;
      data.func_evals = 0;

      double tag_pdf = 0, tag_pdf_err = 0;
      const double xmin = 0, xmax = 1;
      const int integrand_fdim = 1, ndim = 1;
      int cubature_result = hcubature(integrand_fdim, calc_univariate_characteristic_function_for_integration,
        &data, ndim, &xmin, &xmax, cubature_max_evals_, cubature_abs_error_, cubature_rel_error_, ERROR_INDIVIDUAL, &tag_pdf, &tag_pdf_err);
      func_evals += (weights_[tag_index] * (double)data.func_evals);
      if (cubature_result != 0) { num_snp_failed++; continue; }

      if (tag_pdf <= 0)
        tag_pdf = 1e-100;

      double increment = static_cast<double>(-std::log(tag_pdf) * weights_[tag_index]);
      if (!std::isfinite(increment)) num_infinite++;
      log_pdf_total += increment;
    }
  }    

  if (num_snp_failed > 0)
    LOG << " warning: hcubature failed for " << num_snp_failed << " tag snps";
  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  double total_weight = 0.0;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) total_weight += weights_[tag_index];
  func_evals /= total_weight;

  LOG << "<calc_univariate_cost_convolve(" << ss.str() << "), cost=" << log_pdf_total << ", evals=" << func_evals << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

class BivariateCharacteristicFunctionData {
 public:
  const float* pi_vec;
  const float* sig2_zero;
  const float* sig2_beta;
  float rho_beta;
  float rho_zero;

  int tag_index;  // for which SNP to calculate the characteristic function
                  // note that use_complete_tag_indices_ must be enabled for "convolve" calculator, 
                  // so "tag" and "snp" are in the same indexing
  LdMatrixRow* ld_matrix_row;

  const std::vector<float>* hvec;
  const std::vector<float>* zvec1;
  const std::vector<float>* nvec1;
  const std::vector<float>* zvec2;
  const std::vector<float>* nvec2;
  const std::vector<float>* ld_tag_sum_r2_below_r2min_adjust_for_hvec;
  int func_evals;
};

inline float exp_quad_form(float t1, float t2, float a11, float a12, float a22) {
  return std::exp(-0.5f*(t1*t1 * a11 + 2.0*t1*t2*a12 + t2*t2*a22));
}

int calc_bivariate_characteristic_function_times_cosinus(unsigned ndim, const double *x, void *raw_data, unsigned fdim, double* fval) {
  assert(ndim == 2);
  assert(fdim == 1);
  const float t1 = (float)x[0];
  const float t2 = (float)x[1];
  BivariateCharacteristicFunctionData* data = (BivariateCharacteristicFunctionData *)raw_data;
  const float pi1 = data->pi_vec[0];
  const float pi2 = data->pi_vec[1];
  const float pi12 = data->pi_vec[2];
  const float pi0 = 1.0 - pi1 - pi2 - pi12;
  const float eff_sig2_beta1 = (*data->nvec1)[data->tag_index] * data->sig2_beta[0];
  const float eff_sig2_beta2 = (*data->nvec2)[data->tag_index] * data->sig2_beta[1];
  const float eff_sig2_beta_cov = data->rho_beta * sqrt(eff_sig2_beta1 * eff_sig2_beta2);
  const float sig2_zero_cov = data->rho_zero * sqrt(data->sig2_zero[0] * data->sig2_zero[1]);
  const float zval1 = (*data->zvec1)[data->tag_index];
  const float zval2 = (*data->zvec2)[data->tag_index];
  const float inf_adj_r2 = (*data->ld_tag_sum_r2_below_r2min_adjust_for_hvec)[data->tag_index];
  const float scaling_factor = 0.5 * M_1_PI * M_1_PI;

  double result = scaling_factor * cos(t1 * zval1 + t2 * zval2);
  result *= exp_quad_form(t1, t2, data->sig2_zero[0], sig2_zero_cov, data->sig2_zero[1]);

  // apply infinitesimal model to adjust tag_r2sum for all r2 that are below r2min (and thus do not contribute via resampling)
  result *= exp_quad_form(t1, t2, pi1 * eff_sig2_beta1 * inf_adj_r2, 0, 0);
  result *= exp_quad_form(t1, t2, 0, 0, pi2 * eff_sig2_beta2 * inf_adj_r2);
  result *= exp_quad_form(t1, t2, pi12 * eff_sig2_beta1 * inf_adj_r2,
                                  pi12 * eff_sig2_beta_cov * inf_adj_r2,
                                  pi12 * eff_sig2_beta2 * inf_adj_r2);
  
  auto iter_end = data->ld_matrix_row->end();
  for (auto iter = data->ld_matrix_row->begin(); iter < iter_end; iter++) {
    int snp_index = iter.tag_index();  // yes, this is correct - snp_index on the LHS, tag_index on the RHS. 
                                       // See comment in the univariate counterpart, calc_univariate_characteristic_function_times_cosinus
    const float r2 = iter.r2();
    const float hval = (*data->hvec)[snp_index];
    const float r2_times_hval = r2 * hval;
    result *= (double) (pi0 + 
                        pi1 * exp_quad_form(t1, t2, eff_sig2_beta1 * r2_times_hval, 0, 0) +
                        pi2 * exp_quad_form(t1, t2, 0, 0, eff_sig2_beta2 * r2_times_hval) + 
                        pi12* exp_quad_form(t1, t2, eff_sig2_beta1 * r2_times_hval,
                                                    eff_sig2_beta_cov * r2_times_hval,
                                                    eff_sig2_beta2 * r2_times_hval));
  }

  data->func_evals++;
  *fval = result;
  return 0; 
}

int calc_bivariate_characteristic_function_for_integration(unsigned ndim, const double *x, void *raw_data, unsigned fdim, double* fval) {
  double x_transformed[2];

  const double t1 = x[0];
  const double inv_1_minus_t1 = 1.0 / (1.0 - t1);
  x_transformed[0] = t1 * inv_1_minus_t1;
  const double jacob1 = inv_1_minus_t1 * inv_1_minus_t1;

  const double t2 = x[1];
  const double inv_1_minus_t2sqr = 1.0 / (1.0 - t2*t2);
  x_transformed[1] = t2 * inv_1_minus_t2sqr;
  const double jacob2 = (1 + t2*t2) * inv_1_minus_t2sqr * inv_1_minus_t2sqr;
  
  int retval = calc_bivariate_characteristic_function_times_cosinus(ndim, x_transformed, raw_data, fdim, fval);
  (*fval) *= (jacob1 * jacob2);
  return retval;
}

double BgmgCalculator::calc_bivariate_cost_convolve(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  if (!use_complete_tag_indices_) BGMG_THROW_EXCEPTION(::std::runtime_error("Convolve calculator require 'use_complete_tag_indices' option"));
  if (!causalbetavec1_.empty() || !causalbetavec2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("Convolve calculator does not support causalbetavec"));

  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost_convolve(" << ss << ")";

  double log_pdf_total = 0.0;
  int num_snp_failed = 0;
  int num_infinite = 0;
  double func_evals = 0.0;
  SimpleTimer timer(-1);

  std::vector<float> hvec;
  find_hvec(*this, &hvec);

  // in use_complete_tag_indices_ many tag indices are undefined
  // we compute deftag to avoid unbalanced load across OpenMP threads (static scheduler is still the way to go here.)
  std::vector<int> deftag_indices;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec1_[tag_index]) || !std::isfinite(nvec1_[tag_index])) continue;
    if (!std::isfinite(zvec2_[tag_index]) || !std::isfinite(nvec2_[tag_index])) continue;
    deftag_indices.push_back(tag_index);
  }
  const int num_deftag = deftag_indices.size();

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    BivariateCharacteristicFunctionData data;
    data.pi_vec = pi_vec;
    data.sig2_zero = sig2_zero;
    data.sig2_beta = sig2_beta;
    data.rho_zero = rho_zero;
    data.rho_beta = rho_beta;
    data.ld_matrix_row = &ld_matrix_row;
    data.hvec = &hvec;
    data.zvec1 = &zvec1_;
    data.nvec1 = &nvec1_;
    data.zvec2 = &zvec2_;
    data.nvec2 = &nvec2_;
    data.ld_tag_sum_r2_below_r2min_adjust_for_hvec = &ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2(LD_TAG_COMPONENT_BELOW_R2MIN);

#pragma omp for schedule(static) reduction(+: log_pdf_total, num_snp_failed, num_infinite, func_evals)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      double tag_weight = static_cast<double>(weights_[tag_index]);

      const int causal_index = tag_index; // yes, causal==tag in this case -- see a large comment in calc_univariate_characteristic_function_times_cosinus function.
      ld_matrix_csr_.extract_row(causal_index, data.ld_matrix_row);
      data.tag_index = tag_index;
      data.func_evals = 0;

      double tag_pdf = 0, tag_pdf_err = 0;
      const double xmin[2] = {0.0, -1.0}, xmax[2] = {1.0, 1.0};
      const int integrand_fdim = 1, ndim = 2;
      int cubature_result = hcubature(integrand_fdim, calc_bivariate_characteristic_function_for_integration,
        &data, ndim, xmin, xmax, cubature_max_evals_, cubature_abs_error_, cubature_rel_error_, ERROR_INDIVIDUAL, &tag_pdf, &tag_pdf_err);
      func_evals += (weights_[tag_index] * (double)data.func_evals);
      if (cubature_result != 0) { num_snp_failed++; continue; }

      if (tag_pdf <= 0)
        tag_pdf = 1e-100;

      double increment = static_cast<double>(-std::log(tag_pdf) * weights_[tag_index]);
      if (!std::isfinite(increment)) num_infinite++;
      log_pdf_total += increment;
    }
  }    

  if (num_snp_failed > 0)
    LOG << " warning: hcubature failed for " << num_snp_failed << " tag snps";
  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  double total_weight = 0.0;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) total_weight += weights_[tag_index];
  func_evals /= total_weight;

  LOG << "<calc_bivariate_cost_convolve(" << ss << "), cost=" << log_pdf_total << ", evals=" << func_evals << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

void BgmgCalculator::clear_state() {
  LOG << " clear_state";

  ld_matrix_csr_.clear();

  // clear ordering of SNPs
  snp_order_.clear();
  k_pdf_.clear();
  tag_r2sum_.clear();
  last_num_causals_.clear();
}

void BgmgCalculator::clear_tag_r2sum(int component_id) {
  if (component_id < 0 || component_id >= num_components_) BGMG_THROW_EXCEPTION(::std::runtime_error("clear_tag_r2sum: component_id must be between 0 and num_components_"));
  LOG << " clear_tag_r2sum(component_id=" << component_id << ")";
  if (cache_tag_r2sum_) {
    if (tag_r2sum_.empty()) {
      // Initialize
      for (int i = 0; i < num_components_; i++) {
        last_num_causals_.push_back(0.0f);
        tag_r2sum_.push_back(std::make_shared<DenseMatrix<float>>(num_tag_, k_max_));
        tag_r2sum_[i]->InitializeZeros();
      }
    } else {
      // Clear just the requested component
      tag_r2sum_[component_id]->InitializeZeros();
      last_num_causals_[component_id] = 0;
    }
  } else {
    // Cache disabled, clear tag_r2sum to free up memory
    last_num_causals_.clear();
    tag_r2sum_.clear();
  }
}

int64_t BgmgCalculator::set_weights_randprune(int n, float r2_threshold) {
  if (!ld_matrix_csr_.is_ready()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call set_weights_randprune before set_ld_r2_csr"));
  LOG << ">set_weights_randprune(n=" << n << ", r2=" << r2_threshold << ")";
  if (r2_threshold < r2_min_) BGMG_THROW_EXCEPTION(::std::runtime_error("set_weights_randprune: r2 < r2_min_"));
  if (n <= 0) BGMG_THROW_EXCEPTION(::std::runtime_error("set_weights_randprune: n <= 0"));
  SimpleTimer timer(-1);

  std::valarray<int> passed_random_pruning(0, num_tag_);  // count how many times an index  has passed random pruning

#pragma omp parallel
  {
    std::valarray<int> passed_random_pruning_local(0, num_tag_);  // count how many times an index  has passed random pruning
    LdMatrixRow ld_matrix_row;

#pragma omp for schedule(static)
    for (int prune_i = 0; prune_i < n; prune_i++) {
      std::mt19937_64 random_engine;
      random_engine.seed(seed_ + prune_i);

      std::vector<int> candidate_tag_indices(num_tag_, 0);
      std::vector<char> processed_tag_indices(num_tag_, 0);
      for (int i = 0; i < num_tag_; i++) candidate_tag_indices[i] = i;
      std::set<int> non_processed_tag_indices(candidate_tag_indices.begin(), candidate_tag_indices.end());

      // random pruning should operate only on SNPs with defined zvec and nvec.
      int num_excluded_by_defvec = 0;
      for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
        bool exclude = false;
        if ((!zvec1_.empty()) && !std::isfinite(zvec1_[tag_index])) exclude = true;
        if ((!nvec1_.empty()) && !std::isfinite(nvec1_[tag_index])) exclude = true;
        if ((!zvec2_.empty()) && !std::isfinite(zvec2_[tag_index])) exclude = true;
        if ((!nvec2_.empty()) && !std::isfinite(nvec2_[tag_index])) exclude = true;
        if (exclude) {
          processed_tag_indices[tag_index] = 1;         // mark as processed, and
          non_processed_tag_indices.erase(tag_index);   // remove from the set
          num_excluded_by_defvec++;
        }
      }
      if ((prune_i == 0) && (num_excluded_by_defvec > 0))
        LOG << " set_weights_randprune excludes " << num_excluded_by_defvec << " variants due to undefined zvec or nvec";

      while (candidate_tag_indices.size() > 0) {
        // Here is the logic:
        // 1. select a random element X from the candidate_tag_indices
        // 2. if X is present in processed_tag_indices (collision case):
        //    - re-generate candidate_tag_indices from the set of non_processed_tag_indices
        //    - continue while loop.
        // 3. add X to passed_random_pruning
        // 4. query LD matrix for everything in LD with X (we asume that X will be part of that list). Then, for each Y in LD with X:
        //    - add Y to processed_tag_indices
        //    - remove Y from non_processed_tag_indices

        const int random_candidate_index = std::uniform_int_distribution<int>(0, candidate_tag_indices.size() - 1)(random_engine);
        const int random_tag_index = candidate_tag_indices[random_candidate_index];
        if (processed_tag_indices[random_tag_index]) {
          candidate_tag_indices.assign(non_processed_tag_indices.begin(), non_processed_tag_indices.end());
          // Validate that non_processed_tag_indices is consistent with processed_tag_indices.
          // for (int i = 0; i < num_tag_; i++) {
          //  const bool is_processed = (non_processed_tag_indices.find(i) == non_processed_tag_indices.end());
          //  if (processed_tag_indices[i] != is_processed) {
          //    LOG << " set_weights_randprune is stuck, processed_tag_indices inconsistent with non_processed_tag_indices. Cancel random pruning iteration " << prune_i;
          //    candidate_tag_indices.clear();
          //    break;
          //  }
          // }
          continue;
        }

        passed_random_pruning_local[random_tag_index] += 1;
        int causal_index = tag_to_snp_[random_tag_index];
        ld_matrix_csr_.extract_row(causal_index, &ld_matrix_row);
        auto iter_end = ld_matrix_row.end();
        int num_changes = 0;
        for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
          const int tag_index = iter.tag_index();
          const float r2_value = iter.r2();  // here we are interested in r2 (hvec is irrelevant)        
          if (r2_value < r2_threshold) continue;
          if (processed_tag_indices[tag_index]) continue;
          processed_tag_indices[tag_index] = 1;         // mark as processed, and
          non_processed_tag_indices.erase(tag_index);   // remove from the set
          num_changes++;
        }
        if (num_changes == 0) {
          LOG << " set_weights_randprune is stuck, num_changes=0. Cancel random pruning iteration " << prune_i;
          break;
        }
      }
    }

#pragma omp critical
    passed_random_pruning += passed_random_pruning_local;
  }

  weights_.clear(); weights_.resize(num_tag_, 0.0f);
  for (int i = 0; i < num_tag_; i++)
    weights_[i] = static_cast<float>(passed_random_pruning[i]) / static_cast<float>(n);

  LOG << "<set_weights_randprune(n=" << n << ", r2=" << r2_threshold << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::set_snp_order(int component_id, int64_t length, const int* buffer) {
  if ((component_id < 0) || (component_id >= num_components_)) BGMG_THROW_EXCEPTION(::std::runtime_error("component_id out of range"));
  if (length != (max_causals_ * k_max_)) BGMG_THROW_EXCEPTION(::std::runtime_error("buffer length must be max_causals_ * k_max_"));
  if ((snp_order_.size() != num_components_) || (snp_order_[component_id] == nullptr)) BGMG_THROW_EXCEPTION(::std::runtime_error("snp_order_.size() != num_components_, or is empty"));
  if (snp_order_[component_id]->size() != length) BGMG_THROW_EXCEPTION(::std::runtime_error("snp_order_[component_id] has a wrong size"));
  for (int64_t k = 0; k < k_max_; k++)
    for (int64_t j = 0; j < max_causals_; j++)
      (*snp_order_[component_id])(j, k) = buffer[k*max_causals_ + j];
  LOG << " set_snp_order(component_id" << component_id << ")";
  return 0;
}

int64_t BgmgCalculator::retrieve_snp_order(int component_id, int64_t length, int* buffer) {
  if ((component_id < 0) || (component_id >= num_components_)) BGMG_THROW_EXCEPTION(::std::runtime_error("component_id out of range"));
  if (length != (max_causals_ * k_max_)) BGMG_THROW_EXCEPTION(::std::runtime_error("buffer length must be max_causals_ * k_max_"));
  if ((snp_order_.size() != num_components_) || (snp_order_[component_id] == nullptr)) BGMG_THROW_EXCEPTION(::std::runtime_error("snp_order_.size() != num_components_, or is empty"));
  if (snp_order_[component_id]->size() != length) BGMG_THROW_EXCEPTION(::std::runtime_error("snp_order_[component_id] has a wrong size"));
  for (int64_t k = 0; k < k_max_; k++)
    for (int64_t j = 0; j < max_causals_; j++)
      buffer[k*max_causals_ + j] = (*snp_order_[component_id])(j, k);
  LOG << " retrieve_snp_order(component_id" << component_id << ")";
  return 0;
}

int64_t BgmgCalculator::retrieve_k_pdf(int length, double* buffer) {
  if (k_pdf_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("k_pdf_ is not generated; make sure 'calc_univariate_cost_nocache' is called."));
  if (length != k_pdf_.size()) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  LOG << " retrieve_k_pdf()";
  for (int i = 0; i < k_pdf_.size(); i++) buffer[i] = k_pdf_[i];
  return 0;
}

int64_t BgmgCalculator::retrieve_tag_indices(int num_tag, int* tag_indices) {
  if (num_tag != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (num_tag != tag_to_snp_.size()) BGMG_THROW_EXCEPTION(::std::runtime_error("num_tag != tag_to_snp_.size()"));
  LOG << " retrieve_tag_indices()";
  for (int i = 0; i < num_tag_; i++) tag_indices[i] = tag_to_snp_[i];
  return 0;
}

int64_t BgmgCalculator::retrieve_zvec(int trait, int length, float* buffer) {
  if (length != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  std::vector<float>& zvec(*get_zvec(trait));
  if (zvec.size() != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec.size() != num_tag_"));
  LOG << " retrieve_zvec()";
  for (int i = 0; i < num_tag_; i++) buffer[i] = zvec[i];
  return 0;
}

int64_t BgmgCalculator::retrieve_nvec(int trait, int length, float* buffer) {
  if (length != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  std::vector<float>& nvec(*get_nvec(trait));
  if (nvec.size() != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec.size() != num_tag_"));
  LOG << " retrieve_nvec()";
  for (int i = 0; i < num_tag_; i++) buffer[i] = nvec[i];
  return 0;
}

int64_t BgmgCalculator::retrieve_causalbetavec(int trait, int length, float* buffer) {
  if (length != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  const std::vector<float>& causalbetavec(*get_causalbetavec(trait));
  if (causalbetavec.size() != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("causalbetavec.size() != num_snp_"));
  LOG << " retrieve_causalbetavec()";
  for (int i = 0; i < num_snp_; i++) buffer[i] = causalbetavec[i];
  return 0;
}

int64_t BgmgCalculator::retrieve_mafvec(int length, float* buffer) {
  if (length != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (mafvec_.size() != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("mafvec_.size() != num_tag_"));
  LOG << " retrieve_mafvec()";
  for (int i = 0; i < num_snp_; i++) buffer[i] = mafvec_[i];
  return 0;
}

int64_t BgmgCalculator::retrieve_weights(int length, float* buffer) {
  if (length != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (weights_.size() != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("weights_.size() != num_tag_"));
  LOG << " retrieve_weights()";
  for (int i = 0; i < num_tag_; i++) buffer[i] = weights_[i];
  return 0;
}

void BgmgCalculator::find_tag_r2sum_no_cache(int component_id, float num_causal, int k_index, std::vector<float>* buffer) {
  assert(buffer->size() == num_tag_);
  if (snp_order_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("find_tag_r2sum_no_cache called before find_snp_order()"));

  std::vector<std::pair<int, float>> changeset;
  float floor_num_causals = floor(num_causal);
  for (int i = 0; i < (int)floor_num_causals; i++) changeset.push_back(std::make_pair(i, 1.0f));
  changeset.push_back(std::make_pair((int)floor_num_causals, num_causal - floor_num_causals));

  for (int i = 0; i < num_tag_; i++) buffer->at(i) = 0;

  std::vector<float> hvec;
  find_hvec(*this, &hvec);

  LdMatrixRow ld_matrix_row;
  for (auto change : changeset) {
    int scan_index = change.first;
    float scan_weight = change.second;
    int snp_index = (*snp_order_[component_id])(scan_index, k_index);
    ld_matrix_csr_.extract_row(snp_index, &ld_matrix_row);
    auto iter_end = ld_matrix_row.end();
    for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
      const float mafval = mafvec_[snp_index];
      const float hval = 2.0f * mafval * (1 - mafval);
      const int tag_index = iter.tag_index();
      const float r2 = iter.r2();
      buffer->at(tag_index) += (scan_weight * r2 * hval);
    }
  }

  // apply infinitesimal model to adjust tag_r2sum for all r2 that are below r2min (and thus do not contribute via resampling)
  const std::vector<float>& tag_sum_r2_below_r2min = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2(LD_TAG_COMPONENT_BELOW_R2MIN);
  const float pival = num_causal / static_cast<float>(num_snp_);
  for (int i = 0; i < num_tag_; i++) {
    buffer->at(i) += (pival * tag_sum_r2_below_r2min[i]);
  }
}

int64_t BgmgCalculator::retrieve_weighted_causal_r2(int length, float* buffer) {
  if (length != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_snp_: wrong buffer size"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
  if (!ld_matrix_csr_.is_ready()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call retrieve_weighted_causal_r2 before set_ld_r2_csr"));

  LOG << ">retrieve_weighted_causal_r2()";
  SimpleTimer timer(-1);

  for (int i = 0; i < num_snp_; i++) buffer[i] = 0.0f;

  LdMatrixRow ld_matrix_row;
  for (int causal_index = 0; causal_index < num_snp_; causal_index++) {
    ld_matrix_csr_.extract_row(causal_index, &ld_matrix_row);
    auto iter_end = ld_matrix_row.end();
    for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
      const int tag_index = iter.tag_index();
      const float r2 = iter.r2();  // here we are interested in r2 (hvec is irrelevant)          
      buffer[causal_index] += r2 * weights_[tag_index];
    }
  }

  LOG << "<retrieve_weighted_causal_r2(), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

LoglikeCacheElem::LoglikeCacheElem(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float  rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, double cost) {
  pi_vec1_ = pi_vec[0];
  pi_vec2_ = pi_vec[1];
  pi_vec3_ = pi_vec[2];
  sig2_beta1_ = sig2_beta[0];
  sig2_beta2_ = sig2_beta[1];
  sig2_zero1_ = sig2_zero[0];
  sig2_zero2_ = sig2_zero[1];
  rho_zero_ = rho_zero;
  rho_beta_ = rho_beta;
  cost_ = cost;
}

void LoglikeCacheElem::get(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float* rho_beta, int sig2_zero_len, float* sig2_zero, float* rho_zero, double* cost) {
  pi_vec[0] = pi_vec1_;
  pi_vec[1] = pi_vec2_;
  pi_vec[2] = pi_vec3_;
  sig2_beta[0] = sig2_beta1_;
  sig2_beta[1] = sig2_beta2_;
  sig2_zero[0] = sig2_zero1_;
  sig2_zero[1] = sig2_zero2_;
  *rho_zero = rho_zero_;
  *rho_beta = rho_beta_;
  *cost = cost_;
}

void    LoglikeCache::add_entry(float pi_vec, float sig2_zero, float sig2_beta, double cost) {
  float tmp_pi_vec[3] = { pi_vec, NAN, NAN };
  float tmp_sig2_zero[2] = { sig2_zero, NAN };
  float tmp_sig2_beta[2] = { sig2_beta, NAN };
  add_entry(3, tmp_pi_vec, 2, tmp_sig2_beta, NAN, 2, tmp_sig2_zero, NAN, cost);
}

int64_t LoglikeCache::get_entry(int entry_index, float* pi_vec, float* sig2_zero, float* sig2_beta, double* cost) {
  float tmp_pi_vec[3];
  float tmp_sig2_zero[2];
  float tmp_sig2_beta[2];
  float tmp_rho_beta;
  float tmp_rho_zero;
  double tmp_cost;
  int64_t retval = get_entry(entry_index, 3, tmp_pi_vec, 2, tmp_sig2_beta, &tmp_rho_beta, 2, tmp_sig2_zero, &tmp_rho_zero, &tmp_cost);
  *cost = tmp_cost;
  *pi_vec = tmp_pi_vec[0];
  *sig2_zero = tmp_sig2_zero[0];
  *sig2_beta = tmp_sig2_beta[0];
  return retval;
}
void    LoglikeCache::add_entry(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float  rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, double cost) {
  cache_.push_back(LoglikeCacheElem(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, cost));
}
int64_t LoglikeCache::get_entry(int entry_index, int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float* rho_beta, int sig2_zero_len, float* sig2_zero, float* rho_zero, double* cost) {
  if (entry_index < 0 || entry_index > cache_.size()) BGMG_THROW_EXCEPTION(::std::runtime_error("entry_index out of range"));
  cache_[entry_index].get(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, cost);
  return 0;
}

int64_t BgmgCalculator::set_chrnumvec(int num_snp, const int* chrlabel) {
  if (!chrnumvec_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can not set chrnumvec twice"));
  for (int i = 1; i < num_snp; i++) {
    if (chrlabel[i] < chrlabel[i - 1]) BGMG_THROW_EXCEPTION(::std::runtime_error("chrnumvec must be sorted"));
  }

  LOG << ">set_chrnumvec(" << num_snp << "); ";
  check_num_snp(num_snp);
  chrnumvec_.assign(chrlabel, chrlabel + num_snp);
  LOG << "<set_chrnumvec(" << num_snp << "); ";
  return 0;
}

int64_t BgmgCalculator::retrieve_chrnumvec(int length, int* buffer) {
  if (length != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (chrnumvec_.size() != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("chrnumvec_.size() != num_snp_"));
  LOG << " retrieve_chrnumvec()";
  for (int i = 0; i < num_snp_; i++) buffer[i] = chrnumvec_[i];
  return 0;
}

int64_t BgmgCalculator::num_ld_r2_snp(int snp_index) {
  if (!ld_matrix_csr_.is_ready()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call num_ld_r2_snp before set_ld_r2_csr"));
  CHECK_SNP_INDEX((*this), snp_index);
  return ld_matrix_csr_.num_ld_r2(snp_index);
}

int64_t BgmgCalculator::retrieve_ld_r2_snp(int snp_index, int length, int* tag_index, float* r2) {
  if (length != num_ld_r2_snp(snp_index)) BGMG_THROW_EXCEPTION(::std::runtime_error("length does not match num_ld_r2_snp"));
  LOG << " retrieve_ld_r2_snp(snp_index=" << snp_index << ")";
  
  LdMatrixRow ld_matrix_row;
  ld_matrix_csr_.extract_row(snp_index, &ld_matrix_row);
  auto iter_end = ld_matrix_row.end();
  int r2_index = 0;
  for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
    tag_index[r2_index] = iter.tag_index();
    r2[r2_index] = iter.r();
    r2_index++;
  }

  return 0;
}

int64_t BgmgCalculator::num_ld_r2_chr(int chr_label) {
  if (!ld_matrix_csr_.is_ready()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call num_ld_r2_chr before set_ld_r2_csr"));

  int64_t retval = 0;
  for (int snp_index = 0; snp_index < num_snp_; snp_index++) {
    if (chrnumvec_[snp_index] != chr_label) continue;
    retval += num_ld_r2_snp(snp_index);
  }

  return retval;
}

int64_t BgmgCalculator::retrieve_ld_r2_chr(int chr_label, int64_t length, int* snp_index, int* tag_index, float* r2) {
  if (length != num_ld_r2_chr(chr_label)) BGMG_THROW_EXCEPTION(::std::runtime_error("length does not match num_ld_r2_chr"));
  LOG << " retrieve_ld_r2_chr(chr_label=" << chr_label << ")";
  
  int r2_index = 0;
  LdMatrixRow ld_matrix_row;
  for (int causal_index = 0; causal_index < num_snp_; causal_index++) {
    if (chrnumvec_[causal_index] != chr_label) continue;
    ld_matrix_csr_.extract_row(causal_index, &ld_matrix_row);
    auto iter_end = ld_matrix_row.end();
    for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
      snp_index[r2_index] = causal_index;
      tag_index[r2_index] = iter.tag_index();
      r2[r2_index] = iter.r();
      r2_index++;
    }
  }

  return 0;
}

int64_t BgmgCalculator::init(std::string bim_file, std::string frq_file, std::string chr_labels, std::string trait1_file, std::string trait2_file, std::string exclude, std::string extract) {
  if (!trait2_file.empty() && trait1_file.empty())
    BGMG_THROW_EXCEPTION(std::runtime_error("trait2_file can be provided only together with trait1_file"));
  LOG << ">init(bim_file=" << bim_file << ", frq_file=" << frq_file << ", chr_labels=" << chr_labels << ", trait1_file=" << trait1_file << ", trait2_file=" << trait2_file << ", exclude=" << exclude << ", extract=" << extract << "); ";
  SimpleTimer timer(-1);

  std::vector<std::string> chr_labels_vector, bim_files, frq_files;

  if (chr_labels.empty()) {
    for (int i = 1; i <= 22; i++)
      chr_labels_vector.push_back(boost::lexical_cast<std::string>(i));
  } else {
    const std::string separators = " ,;\t\n\r";
    boost::trim_if(chr_labels, boost::is_any_of(separators));
    boost::split(chr_labels_vector, chr_labels, boost::is_any_of(separators), boost::token_compress_on);
  }

  if (bim_file.find("@") != std::string::npos) {
    for (auto chrlabel : chr_labels_vector) {
      bim_files.push_back(bim_file);
      boost::replace_all(bim_files.back(), "@", chrlabel);
    }
  }
  else {
    bim_files.push_back(bim_file);
  }

  if (frq_file.find("@") != std::string::npos) {
    for (auto chrlabel : chr_labels_vector) {
      frq_files.push_back(frq_file);
      boost::replace_all(frq_files.back(), "@", chrlabel);
    }
  }
  else if (!frq_file.empty()) {
    frq_files.push_back(frq_file);
  }

  for (auto& bim_file : bim_files) {
    if (!boost::filesystem::exists(bim_file)) {
      std::stringstream ss; ss << "ERROR: input file " << bim_file << " does not exist";
      BGMG_THROW_EXCEPTION(std::runtime_error(ss.str()));
    }
  }

  for (auto& frq_file : frq_files) {
    if (!boost::filesystem::exists(frq_file)) {
      std::stringstream ss; ss << "ERROR: input file " << frq_file << " does not exist";
      BGMG_THROW_EXCEPTION(std::runtime_error(ss.str()));
    }
  }

  if (!trait1_file.empty() && !boost::filesystem::exists(trait1_file))
    BGMG_THROW_EXCEPTION(std::runtime_error(trait1_file + " does not exist"));
  if (!trait2_file.empty() && !boost::filesystem::exists(trait2_file))
    BGMG_THROW_EXCEPTION(std::runtime_error(trait2_file + " does not exist"));

  bim_file_.clear();
  bim_file_.read(bim_files);
  bim_file_.find_snp_to_index_map();

  FrqFile frq_file_object;
  if (!frq_file.empty()) {
    frq_file_object.read(bim_file_, frq_files);
    frq_file_object.align_to_reference(bim_file_);  // raise an error if any of reference variants is not present in frq files.
  }

  std::vector<int> defvec(bim_file_.size(), 1);

  SumstatFile trait1_file_object;
  if (!trait1_file.empty()) {
    trait1_file_object.read(bim_file_, trait1_file);
    for (int i = 0; i < bim_file_.size(); i++) {
      if (!std::isfinite(trait1_file_object.zscore()[i])) defvec[i] = 0;
      if (!std::isfinite(trait1_file_object.sample_size()[i])) defvec[i] = 0;
    }

    LOG << " constrain analysis to " << std::accumulate(defvec.begin(), defvec.end(), 0) << " tag variants (due to trait1_file='" << trait1_file << "')";
  }

  SumstatFile trait2_file_object;
  if (!trait2_file.empty()) {
    trait2_file_object.read(bim_file_, trait2_file);
    for (int i = 0; i < bim_file_.size(); i++) {
      if (!std::isfinite(trait2_file_object.zscore()[i])) defvec[i] = 0;
      if (!std::isfinite(trait2_file_object.sample_size()[i])) defvec[i] = 0;
    }

    LOG << " constrain analysis to " << std::accumulate(defvec.begin(), defvec.end(), 0) << " tag variants (due to trait2_file='" << trait2_file << "')";
  }

  if (!extract.empty()) {
    SnpList extract_object;
    extract_object.read(extract);
    for (int i = 0; i < bim_file_.size(); i++) {
      if (!extract_object.contains(bim_file_.snp()[i]))
        defvec[i] = 0;
    }

    LOG << " constrain analysis to " << std::accumulate(defvec.begin(), defvec.end(), 0) << " tag variants (due to extract='" << extract << "')";
  }
  
  if (!exclude.empty()) {
    SnpList exclude_object;
    exclude_object.read(exclude);
    for (int i = 0; i < bim_file_.size(); i++) {
      if (exclude_object.contains(bim_file_.snp()[i]))
        defvec[i] = 0;
    }

    LOG << " constrain analysis to " << std::accumulate(defvec.begin(), defvec.end(), 0) << " tag variants (due to exclude='" << exclude << "')";
  }

  // Find tag indices
  std::vector<int> tag_indices;
  for (int i = 0; i < bim_file_.size(); i++)
    if (defvec[i] || use_complete_tag_indices_) tag_indices.push_back(i);

  // Initialize bgmg_calculator, e.i.
  // - set_tag_indices
  // - set_chrnumvec
  // - set_mafvec (if frq file is available)
  // - set_zvec, set_nvec (for each trait, if they are available)
  set_tag_indices(defvec.size(), tag_indices.size(), &tag_indices[0]);
  set_chrnumvec(bim_file_.size(), &bim_file_.chr_label()[0]);
  if (!frq_file.empty())
    set_mafvec(bim_file_.size(), &frq_file_object.frq()[0]);

  if (!trait1_file.empty()) {
    std::vector<float> zvec(tag_indices.size(), 0), nvec(tag_indices.size(), 0);
    for (int i = 0; i < tag_indices.size(); i++) {
      const int snp_index = tag_indices[i];
      zvec[i] = defvec[snp_index] ? trait1_file_object.zscore()[snp_index] : NAN;
      nvec[i] = defvec[snp_index] ? trait1_file_object.sample_size()[snp_index] : NAN;
    }
    set_zvec(1, tag_indices.size(), &zvec[0]);
    set_nvec(1, tag_indices.size(), &nvec[0]);
  }

  if (!trait2_file.empty()) {
    std::vector<float> zvec(tag_indices.size(), 0), nvec(tag_indices.size(), 0);
    for (int i = 0; i < tag_indices.size(); i++) {
      const int snp_index = tag_indices[i];
      zvec[i] = defvec[snp_index] ? trait2_file_object.zscore()[snp_index] : NAN;
      nvec[i] = defvec[snp_index] ? trait2_file_object.sample_size()[snp_index] : NAN;
    }
    set_zvec(2, tag_indices.size(), &zvec[0]);
    set_nvec(2, tag_indices.size(), &nvec[0]);
  }

  LOG << "<init(bim_file=" << bim_file << 
    ", frq_file=" << frq_file << 
    ", chr_labels=" << chr_labels << 
    ", trait1_file=" << trait1_file << 
    ", trait2_file=" << trait2_file << 
    ", exclude=" << exclude << 
    ", extract=" << extract <<
    ");  elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::convert_plink_ld(std::string plink_ld_gz, std::string plink_ld_bin) {
  PlinkLdFile plink_ld_file(bim_file_, plink_ld_gz);
  plink_ld_file.save_as_binary(plink_ld_bin);
  return 0;
}

int64_t BgmgCalculator::num_ld_r2_snp_range(int snp_index_from, int snp_index_to) {
  return retrieve_ld_r2_snp_range(snp_index_from, snp_index_to, -1, nullptr, nullptr, nullptr);
}

int64_t BgmgCalculator::retrieve_ld_r2_snp_range(int snp_index_from, int snp_index_to, int length, int* snp_index, int* tag_index, float* r2) {
  // length < 0 indicate that we just calculate num_ld_r2_snp_range (i.e. the buffer size needed to copy out the result)
  CHECK_SNP_INDEX((*this), snp_index_from);
  CHECK_SNP_INDEX((*this), snp_index_to);
  int64_t num_r2 = 0;
  LdMatrixRow ld_matrix_row;
  for (int iter_snp_index = snp_index_from; iter_snp_index < snp_index_to; iter_snp_index++) {
    ld_matrix_csr_.extract_row(iter_snp_index, &ld_matrix_row);
    auto iter_end = ld_matrix_row.end();
    for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
      const int iter_tag_index = iter.tag_index();
      if (tag_to_snp_[iter_tag_index] < snp_index_from || tag_to_snp_[iter_tag_index] >= snp_index_to) continue;
      num_r2++;
      if (length < 0) continue;
      if (num_r2 > length) BGMG_THROW_EXCEPTION(::std::runtime_error("insufficient length for retrieve_ld_r2_snp_range"));
      snp_index[num_r2 - 1] = iter_snp_index;
      tag_index[num_r2 - 1] = iter_tag_index;
      r2[num_r2 - 1] = iter.r(); // here we are interested in r2 (hvec is irrelevant)
    }
  }
  LOG << ((length < 0) ? " num_ld_r2_snp_range(" : " retrieve_ld_r2_snp_range(from=") << snp_index_from << ", to=" << snp_index_to << "), return " << num_r2;
  return (length < 0) ? num_r2 : 0;
}

int64_t BgmgCalculator::retrieve_fixed_effect_delta(int trait_index, int length, float* delta) {
  check_num_tag(length);
  std::valarray<float> fixed_effect_delta(0.0, num_tag_);
  calc_fixed_effect_delta_from_causalbetavec(trait_index, &fixed_effect_delta);
  for (int i = 0; i < num_tag_; i++) delta[i] = fixed_effect_delta[i];
}

void BgmgCalculator::calc_fixed_effect_delta_from_causalbetavec(int trait_index, std::valarray<float>* delta) {
  if (delta->size() != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_fixed_effect_delta_from_causalbetavec expect delta to be already initialized"));
  *delta = 0.0f;

  const std::vector<float>& causalbetavec(*get_causalbetavec(trait_index));
  if (causalbetavec.empty()) return;
  check_num_snp(causalbetavec.size());

  LOG << ">calc_fixed_effect_delta_from_causalbetavec(trait_index=" << trait_index << ")";
  SimpleTimer timer(-1);

  std::vector<float> sqrt_hvec;
  find_hvec(*this, &sqrt_hvec);
  for (int i = 0; i < sqrt_hvec.size(); i++) sqrt_hvec[i] = sqrt(sqrt_hvec[i]);

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    std::valarray<float> delta_local(0.0, num_tag_);

    // many entries in causalbetavec are expected to be zero, therefore static scheduler may give an unbalanced load
    // however it's fairly short operation anyway, so we don't bother too much.
#pragma omp for schedule(static)
    for (int causal_index = 0; causal_index < num_snp_; causal_index++) {
      if (causalbetavec[causal_index] == 0.0f) continue;
      ld_matrix_csr_.extract_row(causal_index, &ld_matrix_row);
      auto iter_end = ld_matrix_row.end();
      for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
        const int tag_index = iter.tag_index();
        const float r_value = iter.r();
        delta_local[tag_index] += r_value * sqrt_hvec[causal_index] * causalbetavec[causal_index];
      }
    }
#pragma omp critical
    (*delta) += delta_local;
  }

  const std::vector<float>& nvec(*get_nvec(trait_index));
  for (int i = 0; i < nvec.size(); i++) (*delta)[i] *= sqrt(nvec[i]);

  LOG << "<calc_fixed_effect_delta_from_causalbetavec(trait_index=" << trait_index  << "), elapsed time " << timer.elapsed_ms() << "ms";
}