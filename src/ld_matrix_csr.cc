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

#include "ld_matrix_csr.h"

#include <assert.h>
#include <algorithm>
#include <numeric>

#include "TurboPFor/vsimple.h"

void find_hvec(TagToSnpMapping& mapping, std::vector<float>* hvec) {
  const std::vector<float>& mafvec = mapping.mafvec();
  hvec->resize(mafvec.size(), 0.0f);
  for (int snp_index = 0; snp_index < mafvec.size(); snp_index++) {
    hvec->at(snp_index) = 2.0f * mafvec[snp_index] * (1.0f - mafvec[snp_index]);
  }
}

int64_t LdMatrixCsr::set_ld_r2_coo(const std::string& filename, float r2_min) {
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

  return set_ld_r2_coo(numel, &snp_index[0], &tag_index[0], &r2[0], r2_min);
}

int64_t LdMatrixCsr::set_ld_r2_coo(int64_t length, int* snp_index, int* tag_index, float* r2, float r2_min) {
  if (mapping_.mafvec().empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call set_ld_r2_coo before set_mafvec"));
  if (mapping_.chrnumvec().empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call set_ld_r2_coo before set_chrnumvec"));
  LOG << ">set_ld_r2_coo(length=" << length << "); ";

  if (ld_tag_sum_adjust_for_hvec_ == nullptr) {
    ld_tag_sum_adjust_for_hvec_ = std::make_shared<LdTagSum>(LD_TAG_COMPONENT_COUNT, mapping_.num_tag());
    ld_tag_sum_ = std::make_shared<LdTagSum>(LD_TAG_COMPONENT_COUNT, mapping_.num_tag());
  }

  for (int64_t i = 0; i < length; i++)
    if (snp_index[i] == tag_index[i])
      BGMG_THROW_EXCEPTION(::std::runtime_error("snp_index[i] == tag_index[i] --- unexpected for ld files created via plink"));

  for (int64_t i = 0; i < length; i++) {
    if (!std::isfinite(r2[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  SimpleTimer timer(-1);

  std::vector<float> hvec;
  find_hvec(mapping_, &hvec);

  if (chunks_.empty()) {
    int max_chr_label = 0;
    for (int i = 0; i < mapping_.chrnumvec().size(); i++) {
      int chr_label = mapping_.chrnumvec()[i];
      if (chr_label >= max_chr_label) max_chr_label = chr_label;
    }
    chunks_.resize(max_chr_label + 1);  // here we most likely create one useless LD structure for chr_label==0. But that's fine, it'll just stay empty.
    LOG << " highest chr label: " << max_chr_label;

    // Find where each chunk starts and ends
    std::vector<int> chunk_snp_count(max_chr_label + 1, 0);
    for (int i = 0; i < mapping_.chrnumvec().size(); i++) {
      int chr_label = mapping_.chrnumvec()[i];
      chunk_snp_count[chr_label]++;
    }
    for (int chr_label = 0, snp_count_on_previous_chromosomes = 0; chr_label <= max_chr_label; chr_label++) {
      chunks_[chr_label].snp_index_from_inclusive_ = snp_count_on_previous_chromosomes;
      chunks_[chr_label].snp_index_to_exclusive_ = snp_count_on_previous_chromosomes + chunk_snp_count[chr_label];
      chunks_[chr_label].chr_label_ = chr_label;
      snp_count_on_previous_chromosomes += chunk_snp_count[chr_label];
    }

    LOG << " set_ld_r2_coo adds " << mapping_.tag_to_snp().size() << " elements with r2=1.0 to the diagonal of LD r2 matrix";
    for (int i = 0; i < mapping_.tag_to_snp().size(); i++) {
      int snp_index = mapping_.tag_to_snp()[i];
      chunks_[mapping_.chrnumvec()[snp_index]].coo_ld_.push_back(std::make_tuple(snp_index, i, 1.0f));
      ld_tag_sum_adjust_for_hvec_->store(LD_TAG_COMPONENT_ABOVE_R2MIN, i, 1.0f * hvec[mapping_.tag_to_snp()[i]]);
      ld_tag_sum_->store(LD_TAG_COMPONENT_ABOVE_R2MIN, i, 1.0f);
    }
  }

  int64_t new_elements = 0;
  int64_t elements_on_different_chromosomes = 0;
  for (int64_t i = 0; i < length; i++) {
    CHECK_SNP_INDEX(mapping_, snp_index[i]); CHECK_SNP_INDEX(mapping_, tag_index[i]);

    int chr_label = mapping_.chrnumvec()[tag_index[i]];
    if (chr_label != mapping_.chrnumvec()[snp_index[i]]) { elements_on_different_chromosomes++;  continue; }

    int ld_component = (r2[i] < r2_min) ? LD_TAG_COMPONENT_BELOW_R2MIN : LD_TAG_COMPONENT_ABOVE_R2MIN;
    if (mapping_.is_tag()[tag_index[i]]) ld_tag_sum_adjust_for_hvec_->store(ld_component, mapping_.snp_to_tag()[tag_index[i]], r2[i] * hvec[snp_index[i]]);
    if (mapping_.is_tag()[snp_index[i]]) ld_tag_sum_adjust_for_hvec_->store(ld_component, mapping_.snp_to_tag()[snp_index[i]], r2[i] * hvec[tag_index[i]]);

    if (mapping_.is_tag()[tag_index[i]]) ld_tag_sum_->store(ld_component, mapping_.snp_to_tag()[tag_index[i]], r2[i]);
    if (mapping_.is_tag()[snp_index[i]]) ld_tag_sum_->store(ld_component, mapping_.snp_to_tag()[snp_index[i]], r2[i]);

    if (r2[i] < r2_min) continue;
    if (mapping_.is_tag()[tag_index[i]]) { chunks_[chr_label].coo_ld_.push_back(std::make_tuple(snp_index[i], mapping_.snp_to_tag()[tag_index[i]], r2[i])); new_elements++; }
    if (mapping_.is_tag()[snp_index[i]]) { chunks_[chr_label].coo_ld_.push_back(std::make_tuple(tag_index[i], mapping_.snp_to_tag()[snp_index[i]], r2[i])); new_elements++; }
  }
  for (int i = 0; i < chunks_.size(); i++) chunks_[i].coo_ld_.shrink_to_fit();
  if (elements_on_different_chromosomes > 0) LOG << " ignore " << elements_on_different_chromosomes << " r2 elements on between snps located on different chromosomes";
  LOG << "<set_ld_r2_coo: done; (new_elements: " << new_elements << "), elapsed time " << timer.elapsed_ms() << " ms";
  return 0;
}

int64_t LdMatrixCsrChunk::set_ld_r2_csr() {
  if (!csr_ld_snp_index_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call set_ld_r2_csr twice"));
  csr_ld_snp_index_.resize(num_snps_in_chunk() + 1, 0);
  if (coo_ld_.empty()) {
    return 0;
  }

  LOG << ">set_ld_r2_csr(chr_label=" << chr_label_ << "); ";

  SimpleTimer timer(-1);

  LOG << " sorting ld r2 elements... ";
  // Use parallel sort? https://software.intel.com/en-us/articles/a-parallel-stable-sort-using-c11-for-tbb-cilk-plus-and-openmp
#if _OPENMP >= 200805
  {
    SimpleTimer timer2(-1);
    pss::parallel_stable_sort(coo_ld_.begin(), coo_ld_.end(), std::less<std::tuple<int, int, packed_r2_value>>());
    LOG << " pss::parallel_stable_sort took " << timer2.elapsed_ms() << "ms.";
  }
#else
  {
    SimpleTimer timer2(-1);
    std::sort(coo_ld_.begin(), coo_ld_.end());
    LOG << " std::sort took " << timer2.elapsed_ms() << "ms.";
  }
  static bool first_call = true;
  if (first_call) { first_call = false; LOG << " To enable parallel sort within each chr label build bgmglib with compiler that supports OpenMP 3.0";}
#endif

  csr_ld_tag_index_.reserve(coo_ld_.size());
  csr_ld_r2_.reserve(coo_ld_.size());

  for (int64_t i = 0; i < coo_ld_.size(); i++) {
    csr_ld_tag_index_.push_back(std::get<1>(coo_ld_[i]));
    csr_ld_r2_.push_back(std::get<2>(coo_ld_[i]));
  }

  // find starting position for each snp
  std::fill(csr_ld_snp_index_.begin(), csr_ld_snp_index_.end(), coo_ld_.size());
  for (int64_t ld_index = coo_ld_.size() - 1; ld_index >= 0; ld_index--) {
    int snp_index = std::get<0>(coo_ld_[ld_index]);
    if (snp_index < snp_index_from_inclusive_ || snp_index >= snp_index_to_exclusive_) BGMG_THROW_EXCEPTION(std::runtime_error("bgmglib internal error: snp_index < snp_index_from_inclusive_ || snp_index >= snp_index_to_exclusive_"));
    csr_ld_snp_index_[snp_index - snp_index_from_inclusive_] = ld_index;
  }

  for (int i = (csr_ld_snp_index_.size() - 2); i >= 0; i--)
    if (csr_ld_snp_index_[i] > csr_ld_snp_index_[i + 1])
      csr_ld_snp_index_[i] = csr_ld_snp_index_[i + 1];

  coo_ld_.clear();

  LOG << "<set_ld_r2_csr(chr_label=" << chr_label_ << "); elapsed time " << timer.elapsed_ms() << " ms";
  return 0;
}

int64_t LdMatrixCsr::set_ld_r2_csr(float r2_min, int chr_label) {
  if (chr_label < 0) {
    for (int i = 0; i < chunks_.size(); i++) set_ld_r2_csr(r2_min, i);
  } else {  
    chunks_[chr_label].set_ld_r2_csr();
    chunks_[chr_label].validate_ld_r2_csr(r2_min, chr_label, mapping_);
  }
  return 0;
}

int64_t LdMatrixCsrChunk::validate_ld_r2_csr(float r2_min, int chr_label, TagToSnpMapping& mapping_) {
  LOG << ">validate_ld_r2_csr(); ";
  SimpleTimer timer(-1);

  if (csr_ld_r2_.empty()) return 0;  // allow empty chunks
  
  // Test correctness of sparse representation
  if (csr_ld_snp_index_.size() != (num_snps_in_chunk() + 1)) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_snp_index_.size() != (num_snp_ + 1))"));
  for (int i = 0; i < csr_ld_snp_index_.size(); i++) if (csr_ld_snp_index_[i] < 0 || csr_ld_snp_index_[i] > csr_ld_r2_.size()) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_snp_index_[i] < 0 || csr_ld_snp_index_[i] > csr_ld_r2_.size()"));
  for (int i = 1; i < csr_ld_snp_index_.size(); i++) if (csr_ld_snp_index_[i - 1] > csr_ld_snp_index_[i]) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_snp_index_[i-1] > csr_ld_snp_index_[i]"));
  if (csr_ld_snp_index_.back() != csr_ld_r2_.size()) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_snp_index_.back() != csr_ld_r2_.size()"));
  if (csr_ld_tag_index_.size() != csr_ld_r2_.size()) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_tag_index_.size() != csr_ld_r2_.size()"));
  for (int64_t i = 0; i < csr_ld_tag_index_.size(); i++) if (csr_ld_tag_index_[i] < 0 || csr_ld_tag_index_[i] >= mapping_.num_tag()) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_tag_index_ < 0 || csr_ld_tag_index_ >= num_tag_"));

  // Test that all values are between zero and r2min
  for (int64_t i = 0; i < csr_ld_r2_.size(); i++) if (csr_ld_r2_[i].get() < r2_min || csr_ld_r2_[i].get() > 1.0f) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_tag_index_ < 0 || csr_ld_tag_index_ >= num_tag_"));
  for (int64_t i = 0; i < csr_ld_r2_.size(); i++) if (!std::isfinite(csr_ld_r2_[i].get())) BGMG_THROW_EXCEPTION(std::runtime_error("!std::isfinite(csr_ld_r2_[i])"));

  // Test that LDr2 does not have duplicates
  for (int snp_index_in_chunk = 0; snp_index_in_chunk < num_snps_in_chunk(); snp_index_in_chunk++) {
    const int64_t r2_index_from = csr_ld_snp_index_[snp_index_in_chunk];
    const int64_t r2_index_to = csr_ld_snp_index_[snp_index_in_chunk + 1];
    for (int64_t r2_index = r2_index_from; r2_index < (r2_index_to - 1); r2_index++) {
      if (csr_ld_tag_index_[r2_index] == csr_ld_tag_index_[r2_index + 1])
        BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_tag_index_[r2_index] == csr_ld_tag_index_[r2_index + 1]"));
    }
  }

  // Test that LDr2 is symmetric (as long as both SNPs are tag)
  // Test that LDr2 contains the diagonal (as long as we looking at correct chr_label; remember that this validation is for a given LD chunk, e.i. for a given LD chromosome)
  for (int causal_index = snp_index_from_inclusive_; causal_index < snp_index_to_exclusive_; causal_index++) {
    if (!mapping_.is_tag()[causal_index]) continue;
    if (mapping_.chrnumvec()[causal_index] != chr_label) continue;
    const int tag_index_of_the_snp = mapping_.snp_to_tag()[causal_index];

    const int64_t r2_index_from = csr_ld_snp_index_[causal_index - snp_index_from_inclusive_];
    const int64_t r2_index_to = csr_ld_snp_index_[causal_index - snp_index_from_inclusive_ + 1];
    bool ld_r2_contains_diagonal = false;
    for (int64_t r2_index = r2_index_from; r2_index < r2_index_to; r2_index++) {
      const int tag_index = csr_ld_tag_index_[r2_index];
      const float r2 = csr_ld_r2_[r2_index].get();  // here we are interested in r2 (hvec is irrelevant)

      if (tag_index == tag_index_of_the_snp) ld_r2_contains_diagonal = true;
      float r2symm = find_and_retrieve_ld_r2(mapping_.tag_to_snp()[tag_index], tag_index_of_the_snp);
      if (!std::isfinite(r2symm)) BGMG_THROW_EXCEPTION(std::runtime_error("!std::isfinite(r2symm)"));
      if (r2symm != r2) BGMG_THROW_EXCEPTION(std::runtime_error("r2symm != r2"));
    }

    if (!ld_r2_contains_diagonal) BGMG_THROW_EXCEPTION(std::runtime_error("!ld_r2_contains_diagonal"));
  }

  LOG << "<validate_ld_r2_csr (); elapsed time " << timer.elapsed_ms() << " ms";
  return 0;
}

float LdMatrixCsrChunk::find_and_retrieve_ld_r2(int snp_index, int tag_index) {
  auto r2_iter_from = csr_ld_tag_index_.begin() + csr_ld_snp_index_[snp_index - snp_index_from_inclusive_];
  auto r2_iter_to = csr_ld_tag_index_.begin() + csr_ld_snp_index_[snp_index - snp_index_from_inclusive_ + 1];
  auto iter = std::lower_bound(r2_iter_from, r2_iter_to, tag_index);
  return (iter != r2_iter_to) ? csr_ld_r2_[iter - csr_ld_tag_index_.begin()].get() : NAN;
}

size_t LdMatrixCsr::log_diagnostics() {
  size_t mem_bytes = 0, mem_bytes_total = 0;
  for (int i = 0; i < chunks_.size(); i++) {
    LOG << " diag: LdMatrixCsr chunk " << i << ", snp_index in ["<< chunks_[i].snp_index_from_inclusive_ << ", " << chunks_[i].snp_index_to_exclusive_ << ")";
    mem_bytes_total += chunks_[i].log_diagnostics();
  }

  return mem_bytes_total;
}

size_t LdMatrixCsrChunk::log_diagnostics() {
  size_t mem_bytes = 0, mem_bytes_total = 0;
  mem_bytes = coo_ld_.size() * (sizeof(packed_r2_value) + sizeof(int) + sizeof(int)); mem_bytes_total += mem_bytes;
  LOG << " diag: coo_ld_.size()=" << coo_ld_.size() << " (mem usage = " << mem_bytes << " bytes)";
  LOG << " diag: csr_ld_snp_index_.size()=" << csr_ld_snp_index_.size();
  mem_bytes = csr_ld_tag_index_.size() * sizeof(int); mem_bytes_total += mem_bytes;
  LOG << " diag: csr_ld_tag_index_.size()=" << csr_ld_tag_index_.size() << " (mem usage = " << mem_bytes << " bytes)";
  mem_bytes = csr_ld_r2_.size() * sizeof(packed_r2_value); mem_bytes_total += mem_bytes;
  LOG << " diag: csr_ld_r2_.size()=" << csr_ld_r2_.size() << " (mem usage = " << mem_bytes << " bytes)";
  return mem_bytes_total;
}

void LdMatrixCsr::clear() {
  chunks_.clear();
  if (ld_tag_sum_adjust_for_hvec_ != nullptr) ld_tag_sum_adjust_for_hvec_->clear();
  if (ld_tag_sum_ != nullptr) ld_tag_sum_->clear();
}

void LdMatrixCsrChunk::clear() {
  coo_ld_.clear();
}

void LdMatrixCsr::extract_row(int snp_index, LdMatrixRow* row) {
  const int chr_label = mapping_.chrnumvec()[snp_index];
  const LdMatrixCsrChunk& chunk = chunks_[chr_label];
  const int64_t ld_index_begin = chunk.ld_index_begin(snp_index);
  const int64_t ld_index_end = chunk.ld_index_end(snp_index);
  row->tag_index_.resize(ld_index_end - ld_index_begin);  // std::vector.resize() never reduce capacity
  row->r2_.resize(ld_index_end - ld_index_begin);
  for (int64_t ld_index = ld_index_begin; ld_index < ld_index_end; ld_index++) {
    row->tag_index_.at(ld_index - ld_index_begin) = chunk.csr_ld_tag_index_[ld_index];
    row->r2_.at(ld_index - ld_index_begin) = chunk.csr_ld_r2_[ld_index];
  }
}

int LdMatrixCsr::num_ld_r2(int snp_index) {
  const int chr_label = mapping_.chrnumvec()[snp_index];
  const LdMatrixCsrChunk& chunk = chunks_[chr_label];
  const int64_t ld_index_begin = chunk.ld_index_begin(snp_index);
  const int64_t ld_index_end = chunk.ld_index_end(snp_index);
  return ld_index_end - ld_index_begin;
}

int LdMatrixIterator::tag_index() const { return parent_->tag_index_[ld_index_]; }
float LdMatrixIterator::r2() const { return parent_->r2_[ld_index_].get(); }
