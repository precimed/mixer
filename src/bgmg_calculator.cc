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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>

#if _OPENMP >= 200805
#include "parallel_stable_sort.h"
#endif

// Include namespace SEMT & global operators.
#define SEMT_DISABLE_PRINT 0
#include "semt/Semt.h"
// Include macros: INT, DINT, VAR, DVAR, PARAM, DPARAM
#include "semt/Shortcuts.h"

#include "bgmg_log.h"
#include "fmath.hpp"

#define OMP_CHUNK 1000

#define FLOAT_TYPE float

#define LD_TAG_COMPONENT_COUNT 2
#define LD_TAG_COMPONENT_BELOW_R2MIN 0
#define LD_TAG_COMPONENT_ABOVE_R2MIN 1

std::vector<float>* BgmgCalculator::get_zvec(int trait_index) {
  if ((trait_index != 1) && (trait_index != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  return (trait_index == 1) ? &zvec1_ : &zvec2_;
}

std::vector<float>* BgmgCalculator::get_nvec(int trait_index) {
  if ((trait_index != 1) && (trait_index != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  return (trait_index == 1) ? &nvec1_ : &nvec2_;
}

BgmgCalculator::BgmgCalculator() : num_snp_(-1), num_tag_(-1), k_max_(100), seed_(0), r2_min_(0.0), num_components_(1), max_causals_(100000), use_fast_cost_calc_(false), cache_tag_r2sum_(false) {
  boost::posix_time::ptime const time_epoch(boost::gregorian::date(1970, 1, 1));
  seed_ = (boost::posix_time::microsec_clock::local_time() - time_epoch).ticks();
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
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  LOG << " set_nvec(trait=" << trait << "); ";
  check_num_tag(length);
  get_nvec(trait)->assign(values, values + length);
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
    log_disgnostics(); return 0;
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
  } else if (!strcmp(option, "fast_cost")) {
    use_fast_cost_calc_ = (value != 0); return 0;
  } else if (!strcmp(option, "threads")) {
    omp_set_num_threads(static_cast<int>(value)); return 0;
  } else if (!strcmp(option, "cache_tag_r2sum")) {
    cache_tag_r2sum_ = (value != 0);
    for (int component_id = 0; component_id < num_components_; component_id++) clear_tag_r2sum(component_id);
    return 0;
  }

  BGMG_THROW_EXCEPTION(::std::runtime_error("unknown option"));
  return 0;
}

#define CHECK_SNP_INDEX(i) if (i < 0 || i >= num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("CHECK_SNP_INDEX failed"));
#define CHECK_TAG_INDEX(i) if (i < 0 || i >= num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("CHECK_TAG_INDEX failed"));

int64_t BgmgCalculator::set_tag_indices(int num_snp, int num_tag, int* tag_indices) {
  if (num_snp_ != -1 || num_tag_ != -1) BGMG_THROW_EXCEPTION(::std::runtime_error("can not call set_tag_indices twice"));

  LOG << " set_tag_indices(num_snp=" << num_snp << ", num_tag=" << num_tag << "); ";
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

  ld_tag_sum_adjust_for_hvec_ = std::make_shared<LdTagSum>(LD_TAG_COMPONENT_COUNT, num_tag_);
  ld_tag_sum_                 = std::make_shared<LdTagSum>(LD_TAG_COMPONENT_COUNT, num_tag_);
  return 0;
}

int64_t BgmgCalculator::set_ld_r2_coo(const std::string& filename) {
  std::ifstream is(filename, std::ifstream::binary);
  if (!is) BGMG_THROW_EXCEPTION(::std::runtime_error("can't open" + filename));
  if (sizeof(int) != 4) BGMG_THROW_EXCEPTION("sizeof(int) != 4, internal error in BGMG cpp"); // int -> int32_t
  
  int64_t numel;
  is.read(reinterpret_cast<char*>(&numel), sizeof(int64_t));
  LOG << " set_ld_r2_coo(filename=" << filename << "), reading " << numel << " elements...";

  std::vector<int> snp_index(numel, 0), tag_index(numel, 0);
  std::vector<float> r2(numel, 0.0f);

  is.read(reinterpret_cast<char*>(&snp_index[0]), numel * sizeof(int));
  is.read(reinterpret_cast<char*>(&tag_index[0]), numel * sizeof(int));
  is.read(reinterpret_cast<char*>(&r2[0]), numel * sizeof(float));

  if (!is) BGMG_THROW_EXCEPTION(::std::runtime_error("can't read from " + filename));
  is.close();

  return set_ld_r2_coo(numel, &snp_index[0], &tag_index[0], &r2[0]);
}

int64_t BgmgCalculator::set_ld_r2_coo(int64_t length, int* snp_index, int* tag_index, float* r2) {
  if (!csr_ld_r2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call set_ld_r2_coo after set_ld_r2_csr"));
  if (hvec_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call set_ld_r2_coo before set_hvec"));
  LOG << ">set_ld_r2_coo(length=" << length << "); ";

  for (int64_t i = 0; i < length; i++)
    if (snp_index[i] == tag_index[i])
      BGMG_THROW_EXCEPTION(::std::runtime_error("snp_index[i] == tag_index[i] --- unexpected for ld files created via plink"));

  if (snp_order_.empty()) find_snp_order();

  for (int64_t i = 0; i < length; i++) {
    if (!std::isfinite(r2[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  SimpleTimer timer(-1);

  int was = coo_ld_.size();
  for (int64_t i = 0; i < length; i++) {
    CHECK_SNP_INDEX(snp_index[i]); CHECK_SNP_INDEX(tag_index[i]);

    int ld_component = (r2[i] < r2_min_) ? LD_TAG_COMPONENT_BELOW_R2MIN : LD_TAG_COMPONENT_ABOVE_R2MIN;
    if (is_tag_[tag_index[i]]) ld_tag_sum_adjust_for_hvec_->store(ld_component, snp_to_tag_[tag_index[i]], r2[i] * hvec_[snp_index[i]]);
    if (is_tag_[snp_index[i]]) ld_tag_sum_adjust_for_hvec_->store(ld_component, snp_to_tag_[snp_index[i]], r2[i] * hvec_[tag_index[i]]);

    if (is_tag_[tag_index[i]]) ld_tag_sum_->store(ld_component, snp_to_tag_[tag_index[i]], r2[i]);
    if (is_tag_[snp_index[i]]) ld_tag_sum_->store(ld_component, snp_to_tag_[snp_index[i]], r2[i]);

    if (r2[i] < r2_min_) continue;
    // tricky part here is that we take into account snp_can_be_causal_
    // there is no reason to keep LD information about certain causal SNP if we never selecting it as causal
    // (see how snp_can_be_causal_ is created during find_snp_order() call)
    if (snp_can_be_causal_[snp_index[i]] && is_tag_[tag_index[i]]) coo_ld_.push_back(std::make_tuple(snp_index[i], snp_to_tag_[tag_index[i]], r2[i]));
    if (snp_can_be_causal_[tag_index[i]] && is_tag_[snp_index[i]]) coo_ld_.push_back(std::make_tuple(tag_index[i], snp_to_tag_[snp_index[i]], r2[i]));
  }
  LOG << "<set_ld_r2_coo: done; coo_ld_.size()=" << coo_ld_.size() << " (new: " << coo_ld_.size() - was << "), elapsed time " << timer.elapsed_ms() << " ms";
  return 0;
}

int64_t BgmgCalculator::set_ld_r2_csr() {
  if (coo_ld_.empty()) 
    BGMG_THROW_EXCEPTION(::std::runtime_error("coo_ld_ is empty"));

  LOG << ">set_ld_r2_csr (coo_ld_.size()==" << coo_ld_.size() << "); ";

  SimpleTimer timer(-1);

  LOG << " set_ld_r2_csr adds " << tag_to_snp_.size() << " elements with r2=1.0 to the diagonal of LD r2 matrix";
  for (int i = 0; i < tag_to_snp_.size(); i++) {
    coo_ld_.push_back(std::make_tuple(tag_to_snp_[i], i, 1.0f));
    ld_tag_sum_adjust_for_hvec_->store(LD_TAG_COMPONENT_ABOVE_R2MIN, i, 1.0f * hvec_[tag_to_snp_[i]]);
    ld_tag_sum_->store(LD_TAG_COMPONENT_ABOVE_R2MIN, i, 1.0f);
  }
  
  LOG << " sorting ld r2 elements... ";
  SimpleTimer timer2(-1);
  // Use parallel sort? https://software.intel.com/en-us/articles/a-parallel-stable-sort-using-c11-for-tbb-cilk-plus-and-openmp
#if _OPENMP >= 200805
  pss::parallel_stable_sort(coo_ld_.begin(), coo_ld_.end(), std::less<std::tuple<int, int, float>>());
  LOG << " pss::parallel_stable_sort took " << timer.elapsed_ms() << "ms.";
#else
  std::sort(coo_ld_.begin(), coo_ld_.end());
  LOG << " std::sort took " << timer.elapsed_ms() << "ms.";
  LOG << " (to enable parallel sort build bgmglib with compiler that supports OpenMP 3.0";
#endif

  csr_ld_tag_index_.reserve(coo_ld_.size());
  csr_ld_r2_.reserve(coo_ld_.size());

  for (int64_t i = 0; i < coo_ld_.size(); i++) {
    csr_ld_tag_index_.push_back(std::get<1>(coo_ld_[i]));
    csr_ld_r2_.push_back(std::get<2>(coo_ld_[i]));
  }

  // find starting position for each snp
  csr_ld_snp_index_.resize(snp_to_tag_.size() + 1, coo_ld_.size());
  for (int64_t i = (coo_ld_.size() - 1); i >= 0; i--) {
    int snp_index = std::get<0>(coo_ld_[i]);
    csr_ld_snp_index_[snp_index] = i;
  }

  for (int i = (csr_ld_snp_index_.size() - 2); i >= 0; i--)
    if (csr_ld_snp_index_[i] > csr_ld_snp_index_[i + 1])
      csr_ld_snp_index_[i] = csr_ld_snp_index_[i + 1];

  LOG << "<set_ld_r2_csr (coo_ld_.size()==" << coo_ld_.size() << "); elapsed time " << timer.elapsed_ms() << " ms"; 
  coo_ld_.clear();
  validate_ld_r2_csr();

  return 0;
}

int64_t BgmgCalculator::validate_ld_r2_csr() {
  LOG << ">validate_ld_r2_csr(); ";
  SimpleTimer timer(-1);

  // Test correctness of sparse representation
  if (csr_ld_snp_index_.size() != (num_snp_ + 1)) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_snp_index_.size() != (num_snp_ + 1))"));
  for (int i = 0; i < csr_ld_snp_index_.size(); i++) if (csr_ld_snp_index_[i] < 0 || csr_ld_snp_index_[i] > csr_ld_r2_.size()) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_snp_index_[i] < 0 || csr_ld_snp_index_[i] > csr_ld_r2_.size()"));
  for (int i = 1; i < csr_ld_snp_index_.size(); i++) if (csr_ld_snp_index_[i-1] > csr_ld_snp_index_[i]) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_snp_index_[i-1] > csr_ld_snp_index_[i]"));
  if (csr_ld_snp_index_.back() != csr_ld_r2_.size()) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_snp_index_.back() != csr_ld_r2_.size()"));
  if (csr_ld_tag_index_.size() != csr_ld_r2_.size()) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_tag_index_.size() != csr_ld_r2_.size()"));
  for (int64_t i = 0; i < csr_ld_tag_index_.size(); i++) if (csr_ld_tag_index_[i] < 0 || csr_ld_tag_index_[i] >= num_tag_) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_tag_index_ < 0 || csr_ld_tag_index_ >= num_tag_"));

  // Test that all values are between zero and r2min
  for (int64_t i = 0; i < csr_ld_r2_.size(); i++) if (csr_ld_r2_[i] < r2_min_ || csr_ld_r2_[i] > 1.0f) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_tag_index_ < 0 || csr_ld_tag_index_ >= num_tag_"));
  for (int64_t i = 0; i < csr_ld_r2_.size(); i++) if (!std::isfinite(csr_ld_r2_[i])) BGMG_THROW_EXCEPTION(std::runtime_error("!std::isfinite(csr_ld_r2_[i])"));

  // Test that LDr2 does not have duplicates
  for (int causal_index = 0; causal_index < num_snp_; causal_index++) {
    const int64_t r2_index_from = csr_ld_snp_index_[causal_index];
    const int64_t r2_index_to = csr_ld_snp_index_[causal_index + 1];
    for (int64_t r2_index = r2_index_from; r2_index < (r2_index_to - 1); r2_index++) {
      if (csr_ld_tag_index_[r2_index] == csr_ld_tag_index_[r2_index + 1])
        BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_tag_index_[r2_index] == csr_ld_tag_index_[r2_index + 1]"));
    }
  }

  // Test that LDr2 is symmetric (as long as both SNPs are tag)
  // Test that LDr2 contains the diagonal
  for (int causal_index = 0; causal_index < num_snp_; causal_index++) {
    if (!is_tag_[causal_index]) continue;
    const int tag_index_of_the_snp = snp_to_tag_[causal_index];

    const int64_t r2_index_from = csr_ld_snp_index_[causal_index];
    const int64_t r2_index_to = csr_ld_snp_index_[causal_index + 1];
    bool ld_r2_contains_diagonal = false;
    for (int64_t r2_index = r2_index_from; r2_index < r2_index_to; r2_index++) {
      const int tag_index = csr_ld_tag_index_[r2_index];
      const float r2 = csr_ld_r2_[r2_index];  // here we are interested in r2 (hvec is irrelevant)
      
      if (tag_index == tag_index_of_the_snp) ld_r2_contains_diagonal = true;
      float r2symm = find_and_retrieve_ld_r2(tag_to_snp_[tag_index], tag_index_of_the_snp);
      if (!std::isfinite(r2symm)) BGMG_THROW_EXCEPTION(std::runtime_error("!std::isfinite(r2symm)"));
      if (r2symm != r2) BGMG_THROW_EXCEPTION(std::runtime_error("r2symm != r2"));
    }

    if (!ld_r2_contains_diagonal) BGMG_THROW_EXCEPTION(std::runtime_error("!ld_r2_contains_diagonal"));
  }

  LOG << "<validate_ld_r2_csr (); elapsed time " << timer.elapsed_ms() << " ms";
  return 0;
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
  snp_can_be_causal_.resize(num_snp_, 1);
  const bool snp_can_be_causal_is_constant_1 = true;

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

    // Fill in snp_can_be_causal_
    if (!snp_can_be_causal_is_constant_1) {
      for (int k = 0; k < k_max_; k++) {
        for (int i = 0; i < max_causals_; i++) {
          snp_can_be_causal_[(*snp_order_[component_index])(i, k)] = 1;
        }
      }
    }
  }

  int num_can_be_causal = 0;
  for (int i = 0; i < num_snp_; i++) num_can_be_causal += snp_can_be_causal_[i];
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

  // if num_causal is more than twice lower than last_num_causals we should re-calculate tag_r2sum from scratch.
  if (num_causals < (last_num_causals / 2)) {
    clear_tag_r2sum(component_id);
    last_num_causals = 0.0f;
  }

  SimpleTimer timer(-1);

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
  const std::vector<float>& tag_sum_r2_below_r2min = ld_tag_sum_adjust_for_hvec_->ld_tag_sum_r2(LD_TAG_COMPONENT_BELOW_R2MIN);
  const float pival_delta = (num_causals_original - last_num_causals_original) / static_cast<float>(num_snp_);

  // it is OK to parallelize the following loop on k_index, because:
  // - all structures here are readonly, except tag_r2sum_ that we are accumulating
  // - two threads will never touch the same memory location (that's why we choose k_index as an outer loop)
#pragma omp parallel for schedule(static)
  for (int k_index = 0; k_index < k_max_; k_index++) {
    for (auto change : changeset) {
      int scan_index = change.first;
      float scan_weight = change.second;
      int snp_index = (*snp_order_[component_id])(scan_index, k_index);  // index of a causal snp
      int64_t r2_index_from = csr_ld_snp_index_[snp_index];
      int64_t r2_index_to = csr_ld_snp_index_[snp_index + 1];
      for (int64_t r2_index = r2_index_from; r2_index < r2_index_to; r2_index++) {
        int tag_index = csr_ld_tag_index_[r2_index];
        float r2 = csr_ld_r2_[r2_index];
        float hval = hvec_[snp_index];
        (*tag_r2sum_[component_id])(tag_index, k_index) += (scan_weight * r2 * hval);
      }
    }
    for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
      (*tag_r2sum_[component_id])(tag_index, k_index) += (pival_delta * tag_sum_r2_below_r2min[tag_index]);
    }
  }

  LOG << "<find_tag_r2sum(component_id=" << component_id << ", num_causals=" << num_causals_original << ", last_num_causals=" << last_num_causals << "), elapsed time " << timer.elapsed_ms() << "ms";

  last_num_causals_[component_id] = num_causals_original;
  return 0;
}

int64_t BgmgCalculator::set_hvec(int length, float* values) {
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  if (!hvec_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can not set hvec twice"));

  LOG << ">set_hvec(" << length << "); ";
  check_num_snp(length);
  hvec_.assign(values, values + length);
  LOG << "<set_hvec(" << length << "); ";
  return 0;
}

int64_t BgmgCalculator::retrieve_ld_tag_r2_sum(int length, float* buffer) {
  check_num_tag(length);
  LOG << " retrieve_ld_tag_r2_sum()";
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    buffer[tag_index] = ld_tag_sum_->ld_tag_sum_r2()[tag_index];
  }
  return 0;
}

int64_t BgmgCalculator::retrieve_ld_tag_r4_sum(int length, float* buffer) {
  check_num_tag(length);
  LOG << " retrieve_ld_tag_r4_sum()";
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    buffer[tag_index] = ld_tag_sum_->ld_tag_sum_r4()[tag_index];
  }
  return 0;
}

int64_t BgmgCalculator::retrieve_tag_r2_sum(int component_id, float num_causal, int length, float* buffer) {
  if (length != (k_max_ * num_tag_)) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (num_causal < 0 && !cache_tag_r2sum_) BGMG_THROW_EXCEPTION(::std::runtime_error("retrieve_tag_r2sum with num_causal<0 is meant for cache_tag_r2sum==true"));
  if (component_id < 0 || component_id >= num_components_ || tag_r2sum_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong component_id"));

  LOG << " retrieve_tag_r2_sum(component_id=" << component_id << ", num_causal=" << num_causal << ")";

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
  return pdf;
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
  return pdf;
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

  const float pdf = fmath::exp(log_pi + log_dt + log_exp);
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

  if (cache_tag_r2sum_) {
    find_tag_r2sum(component_id, num_causals);
  }

  SimpleTimer timer(-1);

  const float pi_k = 1. / static_cast<float>(k_max_);

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

        float tag_r2sum_value = tag_r2sum[tag_index];
        float sig2eff = tag_r2sum_value * nvec[tag_index] * sig2_beta + sig2_zero;
        float s = sqrt(sig2eff);

        for (int z_index = 0; z_index < length; z_index++) {
          FLOAT_TYPE pdf_tmp = pi_k * gaussian_pdf<FLOAT_TYPE>(zvec[z_index], s);
          pdf_double_local[z_index] += static_cast<double>(pdf_tmp * weights_[tag_index]);
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

double BgmgCalculator::calc_univariate_cost(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<float>& zvec(*get_zvec(trait_index));
  if (zvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec is not set"));
  if (nvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  double cost;
  if (use_fast_cost_calc_) cost = calc_univariate_cost_fast(trait_index,pi_vec, sig2_zero, sig2_beta);
  else if (!cache_tag_r2sum_) cost = calc_univariate_cost_nocache(trait_index, pi_vec, sig2_zero, sig2_beta);
  else cost = calc_univariate_cost_cache(trait_index, pi_vec, sig2_zero, sig2_beta);

  if (!use_fast_cost_calc_) loglike_cache_.add_entry(pi_vec, sig2_zero, sig2_beta, cost);
  return cost;
}

double BgmgCalculator::calc_univariate_cost_cache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<float>& zvec(*get_zvec(trait_index));

  float num_causals = pi_vec * static_cast<float>(num_snp_);
  if ((int)num_causals >= max_causals_) return 1e100; // too large pi_vec
  const int component_id = 0;   // univariate is always component 0.
    
  LOG << ">calc_univariate_cost(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";
  find_tag_r2sum(component_id, num_causals);

  SimpleTimer timer(-1);

  const float pi_k = 1. / static_cast<float>(k_max_);
  
  double log_pdf_total = 0.0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec[tag_index])) continue;

    double pdf_tag = 0.0f;
    for (int k_index = 0; k_index < k_max_; k_index++) {
      float tag_r2sum = (*tag_r2sum_[component_id])(tag_index, k_index);
      float sig2eff = tag_r2sum * nvec[tag_index] * sig2_beta + sig2_zero;

      float s = sqrt(sig2eff);
      FLOAT_TYPE pdf = pi_k * gaussian_pdf<FLOAT_TYPE>(zvec[tag_index], s);
      pdf_tag += static_cast<double>(pdf);
    }
    log_pdf_total += -std::log(pdf_tag) * static_cast<double>(weights_[tag_index]);
  }

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
  find_tag_r2sum(component_id, num_causals);

  SimpleTimer timer(-1);

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
      if (!std::isfinite(zvec[tag_index])) continue;

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
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec[tag_index])) continue;
    log_pdf_total += -std::log(pdf_double[tag_index]) * weights_[tag_index];
    pi_vec_io += (-pdf_deriv_pivec[tag_index] / pdf_double[tag_index]) * weights_[tag_index];
    sig2_zero_io += (-pdf_deriv_sig2zero[tag_index] / pdf_double[tag_index]) * weights_[tag_index];
    sig2_beta_io += (-pdf_deriv_sig2beta[tag_index] / pdf_double[tag_index]) * weights_[tag_index];
  }

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

  LOG << ">calc_univariate_cost_nocache(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";
  
  SimpleTimer timer(-1);

  const float pi_k = 1. / static_cast<float>(rhs.k_max_);

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
          if (!std::isfinite(zvec[tag_index])) continue;

          float tag_r2sum_value = tag_r2sum[tag_index];
          float sig2eff = tag_r2sum_value * nvec[tag_index] * sig2_beta + sig2_zero;

          float s = sqrt(sig2eff);
          T pdf = pi_k * gaussian_pdf<T>(zvec[tag_index], s);
          pdf_double_local[tag_index] += static_cast<double>(pdf);
        }
      }
#pragma omp critical
      pdf_double += pdf_double_local;
  }

  double log_pdf_total = 0.0;
  for (int tag_index = 0; tag_index < rhs.num_tag_; tag_index++) {
    if (rhs.weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec[tag_index])) continue;
    log_pdf_total += -std::log(pdf_double[tag_index]) * rhs.weights_[tag_index];
  }

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
  if (use_fast_cost_calc_) cost = calc_bivariate_cost_fast(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);
  else if (!cache_tag_r2sum_) cost = calc_bivariate_cost_nocache(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);
  else cost = calc_bivariate_cost_cache(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);

  if (!use_fast_cost_calc_) loglike_cache_.add_entry(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, cost);
  return cost;
}

double BgmgCalculator::calc_bivariate_cost_cache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {

  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost(" << ss << ")";

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) return 1e100; // too large pi_vec
  }

  for (int component_id = 0; component_id < 3; component_id++) {
    find_tag_r2sum(component_id, num_causals[component_id]);
  }

  SimpleTimer timer(-1);

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  // pi_k is mixture weight
  const float pi_k = 1. / static_cast<float>(k_max_);

  double log_pdf_total = 0.0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec1_[tag_index])) continue;
    if (!std::isfinite(zvec2_[tag_index])) continue;

    const float z1 = zvec1_[tag_index];
    const float z2 = zvec2_[tag_index];
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

      const FLOAT_TYPE pdf = pi_k * gaussian2_pdf<FLOAT_TYPE>(z1, z2, a11, a12, a22);
      pdf_tag += static_cast<double>(pdf);
    }

    log_pdf_total += static_cast<double>(-std::log(pdf_tag) * weights_[tag_index]);
  }

  LOG << "<calc_bivariate_cost(" << ss << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_bivariate_cost_nocache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost_nocache(" << ss << ")";

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) return 1e100; // too large pi_vec
  }

  SimpleTimer timer(-1);

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  // pi_k is mixture weight
  const float pi_k = 1. / static_cast<float>(k_max_);

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
        if (!std::isfinite(zvec1_[tag_index])) continue;
        if (!std::isfinite(zvec2_[tag_index])) continue;

        const float z1 = zvec1_[tag_index];
        const float z2 = zvec2_[tag_index];
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

        const FLOAT_TYPE pdf = pi_k * gaussian2_pdf<FLOAT_TYPE>(z1, z2, a11, a12, a22);
        pdf_double_local[tag_index] += static_cast<double>(pdf);
      }
    }
#pragma omp critical
    pdf_double += pdf_double_local;
  }

  double log_pdf_total = 0.0;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec1_[tag_index])) continue;
    if (!std::isfinite(zvec2_[tag_index])) continue;
    log_pdf_total += -std::log(pdf_double[tag_index]) * weights_[tag_index];
  }

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

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  }

  if (cache_tag_r2sum_) {
    for (int component_id = 0; component_id < 3; component_id++) {
      find_tag_r2sum(component_id, num_causals[component_id]);
    }
  }

  SimpleTimer timer(-1);

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  // pi_k is mixture weight
  const float pi_k = 1. / static_cast<float>(k_max_);

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
          FLOAT_TYPE pdf_tmp = pi_k * gaussian2_pdf<FLOAT_TYPE>(zvec1[z_index], zvec2[z_index], a11, a12, a22);
          pdf_double_local[z_index] += static_cast<double>(pdf_tmp * weights_[tag_index]);
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

void BgmgCalculator::log_disgnostics() {
  size_t mem_bytes = 0, mem_bytes_total = 0;
  LOG << " diag: num_snp_=" << num_snp_;
  LOG << " diag: num_tag_=" << num_tag_;
  LOG << " diag: csr_ld_snp_index_.size()=" << csr_ld_snp_index_.size();
  mem_bytes = csr_ld_tag_index_.size() * sizeof(int); mem_bytes_total += mem_bytes;
  LOG << " diag: csr_ld_tag_index_.size()=" << csr_ld_tag_index_.size() << " (mem usage = " << mem_bytes << " bytes)";
  mem_bytes = csr_ld_r2_.size() * sizeof(float); mem_bytes_total += mem_bytes;
  LOG << " diag: csr_ld_r2_.size()=" << csr_ld_r2_.size() << " (mem usage = " << mem_bytes << " bytes)";
  mem_bytes = coo_ld_.size() * (sizeof(float) + sizeof(int) + sizeof(int)); mem_bytes_total += mem_bytes;
  LOG << " diag: coo_ld_.size()=" << coo_ld_.size() << " (mem usage = " << mem_bytes << " bytes)";
  LOG << " diag: zvec1_.size()=" << zvec1_.size();
  LOG << " diag: zvec1_=" << std_vector_to_str(zvec1_);
  LOG << " diag: nvec1_.size()=" << nvec1_.size();
  LOG << " diag: nvec1_=" << std_vector_to_str(nvec1_);
  LOG << " diag: zvec2_.size()=" << zvec2_.size();
  LOG << " diag: zvec2_=" << std_vector_to_str(zvec2_);
  LOG << " diag: nvec2_.size()=" << nvec2_.size();
  LOG << " diag: nvec2_=" << std_vector_to_str(nvec2_);
  LOG << " diag: weights_.size()=" << weights_.size();
  LOG << " diag: weights_=" << std_vector_to_str(weights_);
  LOG << " diag: hvec_.size()=" << hvec_.size();
  LOG << " diag: hvec_=" << std_vector_to_str(hvec_);
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
  LOG << " diag: options.max_causals_=" << max_causals_;
  LOG << " diag: options.num_components_=" << num_components_;
  LOG << " diag: options.r2_min_=" << r2_min_;
  LOG << " diag: options.use_fast_cost_calc_=" << (use_fast_cost_calc_ ? "yes" : "no");
  LOG << " diag: options.cache_tag_r2sum_=" << (cache_tag_r2sum_ ? "yes" : "no");
  LOG << " diag: options.seed_=" << (seed_);
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

  int num_zero_tag_r2 = 0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec[tag_index])) continue;
    
    const float tag_r2 = ld_tag_sum_adjust_for_hvec_->ld_tag_sum_r2()[tag_index];
    const float tag_r4 = ld_tag_sum_adjust_for_hvec_->ld_tag_sum_r4()[tag_index];

    if (tag_r2 == 0 || tag_r4 == 0) {
      num_zero_tag_r2++; continue;
    }

    const float tag_chi = tag_r4 / tag_r2;

    const float tag_eta_factor = pi_vec * tag_r2 + (1.0f - pi_vec) * tag_chi;
    const float tag_pi1 = pi_vec * tag_r2 / tag_eta_factor;
    const float tag_pi0 = 1 - tag_pi1;
    const float tag_sig2beta = sig2_beta * tag_eta_factor;

    const float tag_z = zvec[tag_index];
    const float tag_n = nvec[tag_index];
    const FLOAT_TYPE tag_pdf0 = gaussian_pdf<FLOAT_TYPE>(tag_z, sqrt(sig2_zero));
    const FLOAT_TYPE tag_pdf1 = gaussian_pdf<FLOAT_TYPE>(tag_z, sqrt(sig2_zero + tag_n *tag_sig2beta));
    const FLOAT_TYPE tag_pdf = tag_pi0 * tag_pdf0 + tag_pi1 * tag_pdf1;
    log_pdf_total += static_cast<double>(-std::log(tag_pdf) * weights_[tag_index]);
  }

  if (num_zero_tag_r2 > 0)
    LOG << " warning: zero tag_r2 encountered " << num_zero_tag_r2 << " times";
  LOG << "<" << ss.str() << ", cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_bivariate_cost_fast(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost_fast(" << ss << ")";

  double log_pdf_total = 0.0;
  SimpleTimer timer(-1);

  int num_zero_tag_r2 = 0;

  const float s0_a11 = sig2_zero[0];
  const float s0_a22 = sig2_zero[1];
  const float s0_a12 = sqrt(sig2_zero[0] * sig2_zero[1]) * rho_zero;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec1_[tag_index])) continue;
    if (!std::isfinite(zvec2_[tag_index])) continue;

    const float z1 = zvec1_[tag_index];
    const float n1 = nvec1_[tag_index];
    const float z2 = zvec2_[tag_index];
    const float n2 = nvec2_[tag_index];

    const float tag_r2 = ld_tag_sum_adjust_for_hvec_->ld_tag_sum_r2()[tag_index];
    const float tag_r4 = ld_tag_sum_adjust_for_hvec_->ld_tag_sum_r4()[tag_index];

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

    FLOAT_TYPE tag_pdf = 0.0f;
    for (int i = 0; i < 8; i++) {
      const float pi1 = (f0[i] ? tag_pi1[0] : tag_pi0[0]);
      const float pi2 = (f1[i] ? tag_pi1[1] : tag_pi0[1]);
      const float pi3 = (f2[i] ? tag_pi1[2] : tag_pi0[2]);
      const float a11i = s0_a11 + f0[i] * a11[0] + f1[i] * a11[1] + f2[i] * a11[2];
      const float a22i = s0_a22 + f0[i] * a22[0] + f1[i] * a22[1] + f2[i] * a22[2];
      const float a12i = s0_a12 + f0[i] * a12[0] + f1[i] * a12[1] + f2[i] * a12[2];
      tag_pdf += static_cast<double>((pi1*pi2*pi3) * gaussian2_pdf<FLOAT_TYPE>(z1, z2, a11i, a12i, a22i));
    }

    if (tag_pdf <= 0)
      tag_pdf = 1e-100;

    log_pdf_total += static_cast<double>(-std::log(tag_pdf) * weights_[tag_index]);
  }

  if (num_zero_tag_r2 > 0)
    LOG << " warning: zero tag_r2 encountered " << num_zero_tag_r2 << " times";

  LOG << "<calc_bivariate_cost_fast(" << ss << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

void BgmgCalculator::clear_state() {
  LOG << " clear_state";

  // clear all info about LD structure
  csr_ld_snp_index_.clear();
  csr_ld_tag_index_.clear();
  csr_ld_r2_.clear();
  coo_ld_.clear();
  hvec_.clear();
  ld_tag_sum_adjust_for_hvec_->clear();
  ld_tag_sum_->clear();

  // clear ordering of SNPs
  snp_order_.clear();
  tag_r2sum_.clear();
  last_num_causals_.clear();
  snp_can_be_causal_.clear();
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
  if (csr_ld_r2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call set_weights_randprune before set_ld_r2_csr"));
  LOG << ">set_weights_randprune(n=" << n << ", r2=" << r2_threshold << ")";
  if (r2_threshold < r2_min_) BGMG_THROW_EXCEPTION(::std::runtime_error("set_weights_randprune: r2 < r2_min_"));
  if (n <= 0) BGMG_THROW_EXCEPTION(::std::runtime_error("set_weights_randprune: n <= 0"));
  SimpleTimer timer(-1);

  std::valarray<int> passed_random_pruning(0, num_tag_);  // count how many times an index  has passed random pruning

#pragma omp parallel
  {
    std::valarray<int> passed_random_pruning_local(0, num_tag_);  // count how many times an index  has passed random pruning

#pragma omp for schedule(static)
    for (int prune_i = 0; prune_i < n; prune_i++) {
      std::mt19937_64 random_engine;
      random_engine.seed(seed_ + prune_i);

      std::vector<int> candidate_tag_indices(num_tag_, 0);
      std::vector<char> processed_tag_indices(num_tag_, 0);
      for (int i = 0; i < num_tag_; i++) candidate_tag_indices[i] = i;
      std::set<int> non_processed_tag_indices(candidate_tag_indices.begin(), candidate_tag_indices.end());

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
        const int64_t r2_index_from = csr_ld_snp_index_[causal_index];
        const int64_t r2_index_to = csr_ld_snp_index_[causal_index + 1];
        int num_changes = 0;
        for (int64_t r2_index = r2_index_from; r2_index < r2_index_to; r2_index++) {
          const int tag_index = csr_ld_tag_index_[r2_index];
          const float r2_value = csr_ld_r2_[r2_index];  // here we are interested in r2 (hvec is irrelevant)
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

  LOG << ">set_weights_randprune(n=" << n << ", r2=" << r2_threshold << "), elapsed time " << timer.elapsed_ms() << "ms";
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

int64_t BgmgCalculator::retrieve_hvec(int length, float* buffer) {
  if (length != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (hvec_.size() != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("hvec_.size() != num_tag_"));
  LOG << " retrieve_hvec()";
  for (int i = 0; i < num_snp_; i++) buffer[i] = hvec_[i];
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

  std::vector<std::pair<int, float>> changeset;
  float floor_num_causals = floor(num_causal);
  for (int i = 0; i < (int)floor_num_causals; i++) changeset.push_back(std::make_pair(i, 1.0f));
  changeset.push_back(std::make_pair((int)floor_num_causals, num_causal - floor_num_causals));

  for (int i = 0; i < num_tag_; i++) buffer->at(i) = 0;

  for (auto change : changeset) {
    int scan_index = change.first;
    float scan_weight = change.second;
    int snp_index = (*snp_order_[component_id])(scan_index, k_index);
    int64_t r2_index_from = csr_ld_snp_index_[snp_index];
    int64_t r2_index_to = csr_ld_snp_index_[snp_index + 1];
    for (int64_t r2_index = r2_index_from; r2_index < r2_index_to; r2_index++) {
      int tag_index = csr_ld_tag_index_[r2_index];
      float r2 = csr_ld_r2_[r2_index];
      float hval = hvec_[snp_index];
      buffer->at(tag_index) += (scan_weight * r2 * hval);
    }
  }

  // apply infinitesimal model to adjust tag_r2sum for all r2 that are below r2min (and thus do not contribute via resampling)
  const std::vector<float>& tag_sum_r2_below_r2min = ld_tag_sum_adjust_for_hvec_->ld_tag_sum_r2(LD_TAG_COMPONENT_BELOW_R2MIN);
  const float pival = num_causal / static_cast<float>(num_snp_);
  for (int i = 0; i < num_tag_; i++) {
    buffer->at(i) += (pival * tag_sum_r2_below_r2min[i]);
  }
}

int64_t BgmgCalculator::retrieve_weighted_causal_r2(int length, float* buffer) {
  if (length != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_snp_: wrong buffer size"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
  if (csr_ld_r2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call retrieve_weighted_causal_r2 before set_ld_r2_csr"));

  LOG << ">retrieve_weighted_causal_r2()";
  SimpleTimer timer(-1);

  for (int i = 0; i < num_snp_; i++) buffer[i] = 0.0f;
  for (int causal_index = 0; causal_index < num_snp_; causal_index++) {
    const int64_t r2_index_from = csr_ld_snp_index_[causal_index];
    const int64_t r2_index_to = csr_ld_snp_index_[causal_index + 1];
    for (int64_t r2_index = r2_index_from; r2_index < r2_index_to; r2_index++) {
      const int tag_index = csr_ld_tag_index_[r2_index];
      const float r2 = csr_ld_r2_[r2_index];  // here we are interested in r2 (hvec is irrelevant)
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

int64_t BgmgCalculator::set_chrnumvec(int num_snp, int* chrlabel) {
  if (!chrnumvec_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can not set chrnumvec twice"));
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

float BgmgCalculator::find_and_retrieve_ld_r2(int snp_index, int tag_index) {
  auto r2_iter_from = csr_ld_tag_index_.begin() + csr_ld_snp_index_[snp_index];
  auto r2_iter_to = csr_ld_tag_index_.begin() + csr_ld_snp_index_[snp_index + 1];
  auto iter = std::lower_bound(r2_iter_from, r2_iter_to, tag_index);
  return (iter != r2_iter_to) ? csr_ld_r2_[iter - csr_ld_tag_index_.begin()] : NAN;
}