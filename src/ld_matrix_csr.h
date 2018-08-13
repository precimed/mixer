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

#include <memory>
#include <vector>
#include <tuple>

#include "bgmg_log.h"

#define LD_TAG_COMPONENT_COUNT 2
#define LD_TAG_COMPONENT_BELOW_R2MIN 0
#define LD_TAG_COMPONENT_ABOVE_R2MIN 1

// An interface to provide the following information:
// - how many snps are there in the reference, and which of them are available in GWAS ("tag snps")
// - chromosome label for each snp
// - mafvec
class TagToSnpMapping {
public:
  virtual ~TagToSnpMapping() {}
  virtual int num_snp() = 0;
  virtual int num_tag() = 0;
  virtual const std::vector<int>& tag_to_snp() = 0;
  virtual const std::vector<int>& snp_to_tag() = 0;
  virtual const std::vector<char>& is_tag() = 0;
  virtual const std::vector<int>& chrnumvec() = 0;
  virtual const std::vector<float>& mafvec() = 0;
  virtual const std::vector<char>& snp_can_be_causal() = 0;
};

#define CHECK_SNP_INDEX(mapping, i) if (i < 0 || i >= mapping.num_snp()) BGMG_THROW_EXCEPTION(::std::runtime_error("CHECK_SNP_INDEX failed"));
#define CHECK_TAG_INDEX(mapping, i) if (i < 0 || i >= mapping.num_tag)()) BGMG_THROW_EXCEPTION(::std::runtime_error("CHECK_TAG_INDEX failed"));

void find_hvec(TagToSnpMapping& mapping, std::vector<float>* hvec);

// pre-calculated sum of LD r2 and r4 for each tag snp
// This takes hvec into account, e.i. we store sum of r2*hvec and r4*hvec^2.
// The information is in several "components" (for example low and high r2).
class LdTagSum {
public:
  LdTagSum(int num_ld_components, int num_tag_snp) : total_ld_component(num_ld_components) {
    ld_tag_sum_r2_.resize(num_ld_components + 1);
    ld_tag_sum_r4_.resize(num_ld_components + 1);
    for (int ic = 0; ic < (num_ld_components + 1); ic++) {
      ld_tag_sum_r2_[ic].resize(num_tag_snp, 0.0f);
      ld_tag_sum_r4_[ic].resize(num_tag_snp, 0.0f);
    }
  }

  void store(int ld_component, int tag_index, float r2_times_hval) {
    ld_tag_sum_r2_[ld_component][tag_index] += r2_times_hval;
    ld_tag_sum_r4_[ld_component][tag_index] += (r2_times_hval * r2_times_hval);

    ld_tag_sum_r2_[total_ld_component][tag_index] += r2_times_hval;
    ld_tag_sum_r4_[total_ld_component][tag_index] += (r2_times_hval * r2_times_hval);
  }

  const std::vector<float>& ld_tag_sum_r2() const { return ld_tag_sum_r2_[total_ld_component]; }
  const std::vector<float>& ld_tag_sum_r4() const { return ld_tag_sum_r4_[total_ld_component]; }
  const std::vector<float>& ld_tag_sum_r2(int ld_component) const { return ld_tag_sum_r2_[ld_component]; }
  const std::vector<float>& ld_tag_sum_r4(int ld_component) const { return ld_tag_sum_r4_[ld_component]; }

  void clear() {
    for (int i = 0; i < ld_tag_sum_r2_.size(); i++) std::fill(ld_tag_sum_r2_[i].begin(), ld_tag_sum_r2_[i].end(), 0.0f);
    for (int i = 0; i < ld_tag_sum_r4_.size(); i++) std::fill(ld_tag_sum_r4_[i].begin(), ld_tag_sum_r4_[i].end(), 0.0f);
  }
private:
  std::vector<std::vector<float>> ld_tag_sum_r2_;
  std::vector<std::vector<float>> ld_tag_sum_r4_;
  const int total_ld_component;
};

// Class to store LD matrix for a given chromosome (or chunk) in CSR format
class LdMatrixCsrChunk {
 public:
  // csr_ld_snp_index_.size() == num_snp_ + 1; 
  // csr_ld_snp_index_[j]..csr_ld_snp_index_[j+1] is a range of values in CSR matrix corresponding to j-th variant
  // csr_ld_tag_index_.size() == csr_ld_r2_.size() == number of non-zero LD r2 values
  // csr_ld_tag_index_ contains values from 0 to num_tag_-1
  // csr_ld_r2_ contains values from 0 to 1, indicating LD r2 between snp and tag variants
  std::vector<int64_t> csr_ld_snp_index_;
  std::vector<int> csr_ld_tag_index_;  // NB! This array can be very long. Indeed more than 2e9 !
  std::vector<float> csr_ld_r2_;
  std::vector<std::tuple<int, int, float>> coo_ld_; // snp, tag, r2
  size_t log_diagnostics();
  void clear();
};

// Class for sparse LD matrix stored in CSR format (Compressed Sparse Row Format)
class LdMatrixCsr {
 public:
   LdMatrixCsr(TagToSnpMapping& mapping) : mapping_(mapping), chunks_(), combined_() {}

   int64_t set_ld_r2_coo(int64_t length, int* snp_index, int* tag_index, float* r2, float r2_min);
   int64_t set_ld_r2_coo(const std::string& filename, float r2_min);
   int64_t set_ld_r2_csr(float r2_min);  // finalize

   const std::vector<int64_t>& csr_ld_snp_index() { return combined_.csr_ld_snp_index_; }
   const std::vector<int>& csr_ld_tag_index() { return combined_.csr_ld_tag_index_; }
   const std::vector<float>& csr_ld_r2() { return combined_.csr_ld_r2_; }

   const LdTagSum* ld_tag_sum_adjust_for_hvec() { return ld_tag_sum_adjust_for_hvec_.get(); }
   const LdTagSum* ld_tag_sum() { return ld_tag_sum_.get(); }

   size_t log_diagnostics();
   void clear();
private:
  int64_t validate_ld_r2_csr(float r2_min);  // validate
  float find_and_retrieve_ld_r2(int snp_index, int tag_index);  // nan if doesn't exist.

  TagToSnpMapping& mapping_;
  std::vector<LdMatrixCsrChunk> chunks_;  // split per chromosomes (before aggregation)
  LdMatrixCsrChunk combined_;             // final LD matrix (after set_ld_r2_csr, e.i. aggregation across chromosomes)
  std::shared_ptr<LdTagSum> ld_tag_sum_adjust_for_hvec_;
  std::shared_ptr<LdTagSum> ld_tag_sum_;
};

