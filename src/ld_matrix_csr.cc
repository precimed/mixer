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
#include "FastDifferentialCoding/fastdelta.h"

#include "ld_matrix.h"

// Very important to use correct sizes for output buffers.
// There boundaries were suggested in https://github.com/powturbo/TurboPFor/issues/31
#define VSENC_BOUND(n, size) ((n + 32) * ((size)+1) )
#define VSDEC_BOUND(n, size) ((n + 32) * (size))
#define VSDEC_NUMEL(n      ) (n + 32)

void find_hvec_per_chunk(TagToSnpMapping& mapping, std::vector<float>* hvec, int index_from, int index_to) {
  if (mapping.mafvec().empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("set_mafvec() must be called first"));  
  const std::vector<float>& mafvec = mapping.mafvec();
  hvec->resize(index_to - index_from, 0.0f);
  for (int snp_index = index_from; snp_index < index_to; snp_index++) {
    hvec->at(snp_index - index_from) = 2.0f * mafvec[snp_index] * (1.0f - mafvec[snp_index]);
  }
}

void find_hvec(TagToSnpMapping& mapping, std::vector<float>* hvec) {
  find_hvec_per_chunk(mapping, hvec, 0, mapping.num_snp());
}

int64_t LdMatrixCsr::set_ld_r2_coo_version0(int chr_label_data, const std::string& filename, float r2_min) {
  std::vector<int> snp_index;
  std::vector<int> snp_other_index;
  std::vector<float> r;

  // the old-format file contains r2 values, so we take a square root in the next line
  load_ld_matrix_version0(filename, &snp_index, &snp_other_index, &r);
  for (int i = 0; i < r.size(); i++) r[i] = sqrt(r[i]);

  if (chr_label_data < 0 || chr_label_data >= chunks_forward_.size()) BGMG_THROW_EXCEPTION(::std::runtime_error("invalid value for chr_label argument"));
  const int index0 = chunks_forward_[chr_label_data].key_index_from_inclusive_;

  LOG << " set_ld_r2_coo_version0 takes square root of input r2 values, and offset all indices by " << index0;

  // the following operation expects all indices to be 0-based within chunk, but that's not the case in version0 format.
  // version0 all indices were with respect to the overal reference.
  for (int i = 0; i < snp_index.size(); i++) snp_index[i] -= index0;
  for (int i = 0; i < snp_other_index.size(); i++) snp_other_index[i] -= index0;

  return set_ld_r2_coo(chr_label_data, r.size(), &snp_index[0], &snp_other_index[0], &r[0], r2_min);
}

int64_t LdMatrixCsr::set_ld_r2_coo_version1plus(int chr_label, const std::string& filename, float r2_min) {
  LdMatrixCsrChunk chunk;
  std::vector<float> freqvec, ld_r2_sum, ld_r2_sum_adjust_for_hvec;  // these are ignored for now
  load_ld_matrix(filename, &chunk, &freqvec, &ld_r2_sum, &ld_r2_sum_adjust_for_hvec);

  int64_t numel = chunk.csr_ld_r_.size();
  LOG << " set_ld_r2_coo(filename=" << filename << ", numel=" << numel << ")...";

  const int index0 = (chr_label < 0) ? 0 : chunks_forward_[chr_label].key_index_from_inclusive_;

  std::vector<int> snp_index(numel, 0), snp_other_index(numel, 0);
  std::vector<float> r(numel, 0.0f);

  LdMatrixRow ld_matrix_row;
  int64_t elem_index = 0;
  for (int chunk_snp_index = chunk.key_index_from_inclusive_; chunk_snp_index < chunk.key_index_to_exclusive_; chunk_snp_index++) {
    chunk.extract_row(chunk_snp_index, &ld_matrix_row);
    auto iter_end = ld_matrix_row.end();
    for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
      snp_index[elem_index] = chunk_snp_index;
      snp_other_index[elem_index] = iter.index();
      r[elem_index] = iter.r();
      elem_index++;
    }

    // apply freqvec and ld_r2_sum/ld_r2_sum_adjust_for_hvec
    const int global_snp_index = index0 + chunk_snp_index;
    mapping_.mutable_mafvec()->at(global_snp_index) = freqvec[chunk_snp_index];
    ld_sum_->store_below_r2min(global_snp_index, ld_r2_sum[chunk_snp_index], 0);
    ld_sum_adjust_for_hvec_->store_below_r2min(global_snp_index, ld_r2_sum_adjust_for_hvec[chunk_snp_index], 0);
  }

  if (elem_index != numel) BGMG_THROW_EXCEPTION(::std::runtime_error("internal error, elem_index != numel"));

  chunk.clear();  // all information is extracted to snp_index, snp_other_index and r. 
  return set_ld_r2_coo(chr_label, numel, &snp_index[0], &snp_other_index[0], &r[0], r2_min);
}

void LdMatrixCsr::init_chunks() {
  if (mapping_.chrnumvec().empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call init_chunks() before set_chrnumvec"));
  if (!chunks_forward_.empty() || !chunks_reverse_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call init_chunks() twice"));
  LOG << ">LdMatrixCsr::init_chunks(); ";

  if (ld_sum_adjust_for_hvec_ == nullptr) ld_sum_adjust_for_hvec_ = std::make_shared<LdSum>(mapping_);
  if (ld_sum_ == nullptr) ld_sum_ = std::make_shared<LdSum>(mapping_);
  ld_sum_adjust_for_hvec_->clear();
  ld_sum_->clear();  

  int max_chr_label = 0;
  for (int i = 0; i < mapping_.chrnumvec().size(); i++) {
    int chr_label = mapping_.chrnumvec()[i];
    if (chr_label >= max_chr_label) max_chr_label = chr_label;
  }
  chunks_forward_.resize(max_chr_label + 1);  // here we most likely create one useless LD structure for chr_label==0. But that's fine, it'll just stay empty.
  chunks_reverse_.resize(max_chr_label + 1);
  LOG << " highest chr label: " << max_chr_label;

  // Find where each chunk starts and ends
  std::vector<int> chunk_snp_count(max_chr_label + 1, 0);
  std::vector<int> chunk_tag_count(max_chr_label + 1, 0);
  for (int snp_index = 0; snp_index < mapping_.chrnumvec().size(); snp_index++) {
    int chr_label = mapping_.chrnumvec()[snp_index];
    chunk_snp_count[chr_label]++;

    if (mapping_.is_tag()[snp_index])
      chunk_tag_count[chr_label]++;
  }

  for (int chr_label = 0, snp_count_on_previous_chromosomes = 0; chr_label <= max_chr_label; chr_label++) {
    chunks_forward_[chr_label].key_index_from_inclusive_ = snp_count_on_previous_chromosomes;
    chunks_forward_[chr_label].key_index_to_exclusive_ = snp_count_on_previous_chromosomes + chunk_snp_count[chr_label];
    chunks_forward_[chr_label].chr_label_ = chr_label;
    snp_count_on_previous_chromosomes += chunk_snp_count[chr_label];
  }

  for (int chr_label = 0, tag_count_on_previous_chromosomes = 0; chr_label <= max_chr_label; chr_label++) {
    chunks_reverse_[chr_label].key_index_from_inclusive_ = tag_count_on_previous_chromosomes;
    chunks_reverse_[chr_label].key_index_to_exclusive_ = tag_count_on_previous_chromosomes + chunk_tag_count[chr_label];
    chunks_reverse_[chr_label].chr_label_ = chr_label;
    tag_count_on_previous_chromosomes += chunk_tag_count[chr_label];
  }
  LOG << "<LdMatrixCsr::init_chunks(); ";
}

void LdMatrixCsr::init_diagonal(int chr_label) {
  LOG << ">LdMatrixCsr::init_diagonal(chr_label=" << chr_label << "); ";
  if (chr_label < 0 || chr_label >= chunks_forward_.size()) BGMG_THROW_EXCEPTION(::std::runtime_error("invalid value for chr_label argument"));
  const int snp_index_from = chunks_forward_[chr_label].key_index_from_inclusive_;
  const int snp_index_to =  chunks_forward_[chr_label].key_index_to_exclusive_;

  std::vector<float> hvec_per_chunk;
  find_hvec_per_chunk(mapping_, &hvec_per_chunk, snp_index_from, snp_index_to);

  int added = 0;
  for (int snp_index = snp_index_from; snp_index < snp_index_to; snp_index++) {
    ld_sum_adjust_for_hvec_->store_above_r2min(snp_index, 1.0f * hvec_per_chunk[snp_index - snp_index_from], 0);
    ld_sum_->store_above_r2min(snp_index, 1.0f, 0);
    if (!mapping_.is_tag()[snp_index]) continue;
    added++;
    int tag_index = mapping_.snp_to_tag()[snp_index];
    if (!mapping_.disable_snp_to_tag_map())
      chunks_forward_[chr_label].coo_ld_.push_back(std::make_tuple(snp_index, tag_index, 1.0f));
    chunks_reverse_[chr_label].coo_ld_.push_back(std::make_tuple(tag_index, snp_index, 1.0f));
  }
  LOG << " added " << added << " tag (out of " << (snp_index_to - snp_index_from)  << " snps) elements with r2=1.0 to the diagonal of LD r2 matrix";  
  LOG << "<LdMatrixCsr::init_diagonal(chr_label=" << chr_label << "); ";
}

int64_t LdMatrixCsr::set_ld_r2_coo(int chr_label_data, int64_t length, int* snp_index_data, int* snp_other_index_data, float* r, float r2_min) {
  if (mapping_.mafvec().empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call set_ld_r2_coo before set_mafvec"));
  if (mapping_.chrnumvec().empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call set_ld_r2_coo before set_chrnumvec"));
  LOG << ">set_ld_r2_coo(length=" << length << "); ";

  for (int64_t i = 0; i < length; i++)
    if (snp_index_data[i] == snp_other_index_data[i])
      BGMG_THROW_EXCEPTION(::std::runtime_error("snp_index[i] == snp_other_index[i] --- unexpected for ld files created via plink"));

  for (int64_t i = 0; i < length; i++) {
    if (!std::isfinite(r[i]) || r[i] < -1.0f || r[i] > 1.0f) BGMG_THROW_EXCEPTION(::std::runtime_error("encounter undefined, below -1.0 or above 1.0 values"));
  }

  SimpleTimer timer(-1);

  if (chr_label_data < 0) {
    for (int i = 0; i < chunks_forward_.size(); i++) 
      init_diagonal(i);
  } else {
    init_diagonal(chr_label_data);
  }

  int64_t new_elements = 0;
  int64_t elements_on_different_chromosomes = 0;

  if (chr_label_data >= chunks_forward_.size()) BGMG_THROW_EXCEPTION(::std::runtime_error("invalid value for chr_label argument"));
  const int snp_index_from = (chr_label_data < 0) ? 0 : chunks_forward_[chr_label_data].key_index_from_inclusive_;
  const int snp_index_to = (chr_label_data < 0) ? mapping_.num_snp() : chunks_forward_[chr_label_data].key_index_to_exclusive_;

  std::vector<float> hvec_per_chunk;
  find_hvec_per_chunk(mapping_, &hvec_per_chunk, snp_index_from, snp_index_to);

  for (int64_t i = 0; i < length; i++) {
    const int snp_index = snp_index_data[i] + snp_index_from;
    const int snp_other_index = snp_other_index_data[i] + snp_index_from;
    CHECK_SNP_INDEX(mapping_, snp_index); CHECK_SNP_INDEX(mapping_, snp_other_index);

    int chr_label = mapping_.chrnumvec()[snp_other_index];
    if (chr_label != mapping_.chrnumvec()[snp_index]) { elements_on_different_chromosomes++;  continue; }

    const float r2 = r[i] * r[i];

    ld_sum_adjust_for_hvec_->store(r2 < r2_min, snp_other_index, r2 * hvec_per_chunk[snp_index - snp_index_from], 0);
    ld_sum_adjust_for_hvec_->store(r2 < r2_min, snp_index, r2 * hvec_per_chunk[snp_other_index - snp_index_from], 0);
    ld_sum_->store(r2 < r2_min, snp_other_index, r2, 0);
    ld_sum_->store(r2 < r2_min, snp_index, r2, 0);

    if (r2 < r2_min) continue;

    if (mapping_.is_tag()[snp_other_index]) { 
      const int tag_other_index = mapping_.snp_to_tag()[snp_other_index];
      if (!mapping_.disable_snp_to_tag_map())
        chunks_forward_[chr_label].coo_ld_.push_back(std::make_tuple(snp_index, tag_other_index, r[i]));
      chunks_reverse_[chr_label].coo_ld_.push_back(std::make_tuple(tag_other_index, snp_index, r[i]));
      new_elements++;
    }

    if (mapping_.is_tag()[snp_index]) {
      const int tag_index = mapping_.snp_to_tag()[snp_index];
      if (!mapping_.disable_snp_to_tag_map())
        chunks_forward_[chr_label].coo_ld_.push_back(std::make_tuple(snp_other_index, tag_index, r[i]));
      chunks_reverse_[chr_label].coo_ld_.push_back(std::make_tuple(tag_index, snp_other_index, r[i]));
      new_elements++;
    }
  }

  for (int i = 0; i < chunks_forward_.size(); i++) chunks_forward_[i].coo_ld_.shrink_to_fit();
  for (int i = 0; i < chunks_reverse_.size(); i++) chunks_reverse_[i].coo_ld_.shrink_to_fit();
  if (elements_on_different_chromosomes > 0) LOG << " ignore " << elements_on_different_chromosomes << " r2 elements on between snps located on different chromosomes";
  LOG << "<set_ld_r2_coo: done; (new_elements: " << new_elements << "), elapsed time " << timer.elapsed_ms() << " ms";
  return 0;
}

int64_t LdMatrixCsrChunk::set_ld_r2_csr() {
  if (!csr_ld_key_index_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't call set_ld_r2_csr twice"));
  csr_ld_key_index_.resize(num_keys_in_chunk() + 1, 0);
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
    pss::parallel_stable_sort(coo_ld_.begin(), coo_ld_.end(), std::less<std::tuple<int, int, packed_r_value>>());
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

  std::vector<uint32_t> csr_ld_val_index_;
  csr_ld_val_index_.reserve(coo_ld_.size());
  csr_ld_r_.reserve(coo_ld_.size());

  for (int64_t i = 0; i < coo_ld_.size(); i++) {
    csr_ld_val_index_.push_back(std::get<1>(coo_ld_[i]));
    csr_ld_r_.push_back(std::get<2>(coo_ld_[i]));
  }

  // find starting position for each key
  std::fill(csr_ld_key_index_.begin(), csr_ld_key_index_.end(), coo_ld_.size());
  for (int64_t ld_index = coo_ld_.size() - 1; ld_index >= 0; ld_index--) {
    int key_index = std::get<0>(coo_ld_[ld_index]);
    if (key_index < key_index_from_inclusive_ || key_index >= key_index_to_exclusive_)
      BGMG_THROW_EXCEPTION(std::runtime_error("bgmglib internal error: key_index < key_index_from_inclusive_ || key_index >= key_index_to_exclusive_"));
    csr_ld_key_index_[key_index - key_index_from_inclusive_] = ld_index;
  }

  for (int i = (csr_ld_key_index_.size() - 2); i >= 0; i--)
    if (csr_ld_key_index_[i] > csr_ld_key_index_[i + 1])
      csr_ld_key_index_[i] = csr_ld_key_index_[i + 1];

  coo_ld_.clear();

  validate_ld_r2_csr(csr_ld_val_index_);

  {
    LOG << ">pack_ld_r2_csr(); ";
    SimpleTimer timer3(-1);
    // pack LD structure
    csr_ld_val_index_offset_.reserve(csr_ld_key_index_.size()); csr_ld_val_index_offset_.push_back(0);
    int64_t buffer_size = 0;
    for (int key_index_in_chunk = 0; key_index_in_chunk < num_keys_in_chunk(); key_index_in_chunk++) {
      int64_t ld_index_from = csr_ld_key_index_[key_index_in_chunk];
      int64_t ld_index_to = csr_ld_key_index_[key_index_in_chunk + 1];
      int64_t num_ld_indices = ld_index_to - ld_index_from;

      if (num_ld_indices > 0) {
        compute_deltas_inplace(&csr_ld_val_index_[ld_index_from], num_ld_indices, 0);

        const int growth_factor = 2;
        csr_ld_val_index_packed_.resize(growth_factor * buffer_size + VSENC_BOUND(num_ld_indices, sizeof(uint32_t)));
        unsigned char *inptr = &csr_ld_val_index_packed_[buffer_size];
        unsigned char *outptr = vsenc32(&csr_ld_val_index_[ld_index_from], num_ld_indices, inptr);
        buffer_size += (outptr - inptr);
      }

      csr_ld_val_index_offset_.push_back(buffer_size);
    }
    csr_ld_val_index_packed_.resize(buffer_size);
    csr_ld_val_index_packed_.shrink_to_fit();
    LOG << "<pack_ld_r2_csr(); elapsed time " << timer3.elapsed_ms() << " ms";
  }

  LOG << "<set_ld_r2_csr(chr_label=" << chr_label_ << "); elapsed time " << timer.elapsed_ms() << " ms";
  return 0;
}

int64_t LdMatrixCsr::set_ld_r2_csr(float r2_min, int chr_label) {
  if (chr_label < 0) {
    for (int i = 0; i < chunks_forward_.size(); i++) set_ld_r2_csr(r2_min, i);
  } else {  
    chunks_forward_[chr_label].set_ld_r2_csr();
    chunks_reverse_[chr_label].set_ld_r2_csr();
  }
  return 0;
}

int64_t LdMatrixCsrChunk::validate_ld_r2_csr(const std::vector<uint32_t>& csr_ld_val_index_) {
  LOG << ">validate_ld_r2_csr(); ";
  SimpleTimer timer(-1);

  if (csr_ld_r_.empty()) return 0;  // allow empty chunks
  
  // Test correctness of sparse representation
  if (csr_ld_key_index_.size() != (num_keys_in_chunk() + 1)) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_key_index_.size() != (num_keys_in_chunk() + 1))"));
  for (int i = 0; i < csr_ld_key_index_.size(); i++) if (csr_ld_key_index_[i] < 0 || csr_ld_key_index_[i] > csr_ld_r_.size()) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_key_index_[i] < 0 || csr_ld_key_index_[i] > csr_ld_r_.size()"));
  for (int i = 1; i < csr_ld_key_index_.size(); i++) if (csr_ld_key_index_[i - 1] > csr_ld_key_index_[i]) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_key_index_[i-1] > csr_ld_key_index_[i]"));
  if (csr_ld_key_index_.back() != csr_ld_r_.size()) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_key_index_.back() != csr_ld_r_.size()"));
  if (csr_ld_val_index_.size() != csr_ld_r_.size()) BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_val_index_.size() != csr_ld_r_.size()"));

  // Test that LDr2 does not have duplicates
  for (int key_index_in_chunk = 0; key_index_in_chunk < num_keys_in_chunk(); key_index_in_chunk++) {
    const int64_t r2_index_from = csr_ld_key_index_[key_index_in_chunk];
    const int64_t r2_index_to = csr_ld_key_index_[key_index_in_chunk + 1];
    for (int64_t r2_index = r2_index_from; r2_index < (r2_index_to - 1); r2_index++) {
      if (csr_ld_val_index_[r2_index] == csr_ld_val_index_[r2_index + 1])
        BGMG_THROW_EXCEPTION(std::runtime_error("csr_ld_val_index_[r2_index] == csr_ld_val_index_[r2_index + 1]"));
    }
  }

  LOG << "<validate_ld_r2_csr (); elapsed time " << timer.elapsed_ms() << " ms";
  return 0;
}

float LdMatrixCsrChunk::find_and_retrieve_ld_r2(int key_index, int val_index, const std::vector<uint32_t>& csr_ld_val_index_) {
  auto r2_iter_from = csr_ld_val_index_.begin() + csr_ld_key_index_[key_index - key_index_from_inclusive_];
  auto r2_iter_to = csr_ld_val_index_.begin() + csr_ld_key_index_[key_index - key_index_from_inclusive_ + 1];
  auto iter = std::lower_bound(r2_iter_from, r2_iter_to, val_index);
  return (iter != r2_iter_to) ? csr_ld_r_[iter - csr_ld_val_index_.begin()].get() : NAN;
}

size_t LdMatrixCsr::log_diagnostics() {
  size_t mem_bytes = 0, mem_bytes_total = 0;
  LOG << " diag: LdMatrixCsr " << chunks_forward_.size() <<  " chunks in total. Logging futher info for non-empty chunks only.";
  for (int i = 0; i < chunks_forward_.size(); i++) {
    if (chunks_forward_[i].is_empty()) continue;
    LOG << " diag: LdMatrixCsr forward chunk " << i << ", key_index in ["<< chunks_forward_[i].key_index_from_inclusive_ << ", " << chunks_forward_[i].key_index_to_exclusive_ << ")";
    mem_bytes_total += chunks_forward_[i].log_diagnostics();
  }
  for (int i = 0; i < chunks_reverse_.size(); i++) {
    if (chunks_reverse_[i].is_empty()) continue;
    LOG << " diag: LdMatrixCsr reverse chunk " << i << ", key_index in ["<< chunks_reverse_[i].key_index_from_inclusive_ << ", " << chunks_reverse_[i].key_index_to_exclusive_ << ")";
    mem_bytes_total += chunks_reverse_[i].log_diagnostics();
  }

  return mem_bytes_total;
}

size_t LdMatrixCsrChunk::log_diagnostics() {
  size_t mem_bytes = 0, mem_bytes_total = 0;
  mem_bytes = coo_ld_.size() * (sizeof(packed_r_value) + sizeof(int) + sizeof(int)); mem_bytes_total += mem_bytes;
  LOG << " diag: coo_ld_.size()=" << coo_ld_.size() << " (mem usage = " << mem_bytes << " bytes)";
  LOG << " diag: csr_ld_key_index_.size()=" << csr_ld_key_index_.size();
  LOG << " diag: csr_ld_val_index_offset_.size()=" << csr_ld_val_index_offset_.size();
  mem_bytes = csr_ld_val_index_packed_.size() * sizeof(unsigned char); mem_bytes_total += mem_bytes;
  LOG << " diag: csr_ld_val_index_packed_.size()=" << csr_ld_val_index_packed_.size() << " (mem usage = " << mem_bytes << " bytes)";
  mem_bytes = csr_ld_r_.size() * sizeof(packed_r_value); mem_bytes_total += mem_bytes;
  LOG << " diag: csr_ld_r_.size()=" << csr_ld_r_.size() << " (mem usage = " << mem_bytes << " bytes)";
  return mem_bytes_total;
}

void LdMatrixCsr::clear() {
  chunks_forward_.clear();
  chunks_reverse_.clear();
  if (ld_sum_adjust_for_hvec_ != nullptr) ld_sum_adjust_for_hvec_->clear();
  if (ld_sum_ != nullptr) ld_sum_->clear();
}

void LdMatrixCsrChunk::clear() {
  coo_ld_.clear();
}

void LdMatrixCsrChunk::extract_row(int key_index, LdMatrixRow* row) {
  const int64_t num_ld_r2 = this->num_ld_r2(key_index);
  row->index_vector_.reserve(VSDEC_NUMEL(num_ld_r2));
  row->index_vector_.resize(num_ld_r2);  // std::vector.resize() never reduce capacity
  row->r_.resize(num_ld_r2);

  // empty LD entry
  if (num_ld_r2 == 0) return;

  vsdec32(this->csr_ld_val_index_packed(key_index), num_ld_r2, reinterpret_cast<unsigned int*>(&row->index_vector_[0]));

  compute_prefix_sum_inplace(reinterpret_cast<uint32_t*>(&row->index_vector_[0]), num_ld_r2, 0);

  const int64_t ld_index_begin = this->ld_index_begin(key_index);
  const int64_t ld_index_end = this->ld_index_end(key_index);
  for (int64_t ld_index = ld_index_begin; ld_index < ld_index_end; ld_index++) {
    row->r_.at(ld_index - ld_index_begin) = this->csr_ld_r_[ld_index];
  }
}

void LdMatrixCsr::extract_snp_row(SnpIndex snp_index, LdMatrixRow* row) {
  if (mapping_.num_tag() == mapping_.num_snp()) {
    extract_tag_row(TagIndex(snp_index.index()), row);
    return;
  }

  const int chr_label = mapping_.chrnumvec()[snp_index.index()];
  LdMatrixCsrChunk& chunk = chunks_forward_[chr_label];
  if (chunk.is_empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("chunks_forward_ is empty()"));  
  chunk.extract_row(snp_index.index(), row);
}

void LdMatrixCsr::extract_tag_row(TagIndex tag_index, LdMatrixRow* row) {
  const int chr_label = mapping_.chrnumvec()[mapping_.tag_to_snp()[tag_index.index()]];
  LdMatrixCsrChunk& chunk = chunks_reverse_[chr_label];
  if (chunk.is_empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("chunks_reverse_ is empty()"));  
  chunk.extract_row(tag_index.index(), row);
}

int LdMatrixCsr::num_ld_r2_snp(int snp_index) {
  const int chr_label = mapping_.chrnumvec()[snp_index];
  const LdMatrixCsrChunk& chunk = chunks_forward_[chr_label];
  return chunk.num_ld_r2(snp_index);
}

int LdMatrixIterator::index() const { return parent_->index_vector_[ld_index_]; }
float LdMatrixIterator::r() const { return parent_->r_[ld_index_].get(); }
float LdMatrixIterator::r2() const { float r = parent_->r_[ld_index_].get(); return r*r; }

