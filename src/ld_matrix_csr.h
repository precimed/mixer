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
#include <numeric>

#include "bgmg_log.h"

#if _OPENMP >= 200805
#include "parallel_stable_sort.h"
#endif

// TagIndex and SnpIndex classes wrap an integer that represents either a tag or a snp,
// to enable compile-time checks (protection agains typos in the code).

class TagIndex {
 public:
  explicit TagIndex(int index) : index_(index) {}
  int operator()() const { return index_; }
  int index() const { return index_; }
 private:
  int index_;
};

class SnpIndex {
 public:
  explicit SnpIndex(int index) : index_(index) {}
  int operator()() const { return index_; }
  int index() const { return index_; }
 private:
  int index_;
};

class packed_r_value {
 public:
   packed_r_value() : value_(0) {}
   packed_r_value(float value) {
    assert((value >= -1.0f) && (value <= 1.0f));
    value_ = static_cast<uint16_t>(roundf((0.5f + 0.5f * value) * 65535.0f));
  }

  uint16_t raw_value() const { return value_; };
  float get() const {
    static float data_[65536];
    static bool initialized;
    if (!initialized) {
      initialized = true;
      for (int i = 0; i < 65536; i++) data_[i] = -1.0f + 2.0f * static_cast<float>(i) / 65535.0f;
    }

    return data_[value_];
  }

 private:
  uint16_t value_;
  friend inline bool operator <(const packed_r_value& lhs, const packed_r_value& rhs);
};

inline bool operator <(const packed_r_value& lhs, const packed_r_value& rhs)
{
  return lhs.value_ < rhs.value_;
}

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
  virtual std::vector<float>* mutable_mafvec() = 0;
};

#define CHECK_SNP_INDEX(mapping, i) if (i < 0 || i >= mapping.num_snp()) BGMG_THROW_EXCEPTION(::std::runtime_error("CHECK_SNP_INDEX failed"));
#define CHECK_TAG_INDEX(mapping, i) if (i < 0 || i >= mapping.num_tag()) BGMG_THROW_EXCEPTION(::std::runtime_error("CHECK_TAG_INDEX failed"));

void find_hvec_per_chunk(TagToSnpMapping& mapping, std::vector<float>* hvec, int index_from, int index_to);
void find_hvec(TagToSnpMapping& mapping, std::vector<float>* hvec);

// pre-calculated sum of LD r2 and r4 for each snp
// This takes hvec into account, e.i. we store sum of r2*hvec and r4*hvec^2.
class LdSum {
public:
  LdSum(TagToSnpMapping& mapping) : mapping_(mapping) {
    const int num_snp = mapping.num_snp();
    const int num_tag = mapping.num_tag();
    ld_sum_r2_below_r2min_.resize(num_snp, 0.0f);
    ld_sum_r2_above_r2min_.resize(num_snp, 0.0f);
    ld_sum_r4_above_r2min_.resize(num_snp, 0.0f);
    ld_tag_sum_r2_below_r2min_.resize(num_tag, 0.0f);
    ld_tag_sum_r2_above_r2min_.resize(num_tag, 0.0f);
    ld_tag_sum_r4_above_r2min_.resize(num_tag, 0.0f);
  }

  void store_below_r2min(int snp_index, float r2_times_hval, int to_be_removed) {
    ld_sum_r2_below_r2min_[snp_index] += r2_times_hval;
    if (mapping_.is_tag()[snp_index]) {
      ld_tag_sum_r2_below_r2min_[mapping_.snp_to_tag()[snp_index]] += r2_times_hval;
    }
  }

  void store_above_r2min(int snp_index, float r2_times_hval, int to_be_removed) {
    ld_sum_r2_above_r2min_[snp_index] += r2_times_hval;
    ld_sum_r4_above_r2min_[snp_index] += (r2_times_hval * r2_times_hval);
    if (mapping_.is_tag()[snp_index]) {
      const int tag_index = mapping_.snp_to_tag()[snp_index];
      ld_tag_sum_r2_above_r2min_[tag_index] += r2_times_hval;
      ld_tag_sum_r4_above_r2min_[tag_index] += (r2_times_hval * r2_times_hval);
    }
  }

  void store(bool below_r2min, int index, float r2_times_hval, int to_be_removed) {
    if (below_r2min) store_below_r2min(index, r2_times_hval, 0);
    else store_above_r2min(index, r2_times_hval, 0);
  }

  const std::vector<float>& ld_sum_r2_below_r2min() const { return ld_sum_r2_below_r2min_; }
  const std::vector<float>& ld_sum_r2_above_r2min() const { return ld_sum_r2_above_r2min_; }
  const std::vector<float>& ld_sum_r4_above_r2min() const { return ld_sum_r4_above_r2min_; }

  const std::vector<float>& ld_tag_sum_r2_below_r2min() const { return ld_tag_sum_r2_below_r2min_; }
  const std::vector<float>& ld_tag_sum_r2_above_r2min() const { return ld_tag_sum_r2_above_r2min_; }
  const std::vector<float>& ld_tag_sum_r4_above_r2min() const { return ld_tag_sum_r4_above_r2min_; }

  void clear() {
    std::fill(ld_sum_r2_below_r2min_.begin(), ld_sum_r2_below_r2min_.end(), 0.0f);
    std::fill(ld_sum_r2_above_r2min_.begin(), ld_sum_r2_above_r2min_.end(), 0.0f);
    std::fill(ld_sum_r4_above_r2min_.begin(), ld_sum_r4_above_r2min_.end(), 0.0f);
    std::fill(ld_tag_sum_r2_below_r2min_.begin(), ld_tag_sum_r2_below_r2min_.end(), 0.0f);
    std::fill(ld_tag_sum_r2_above_r2min_.begin(), ld_tag_sum_r2_above_r2min_.end(), 0.0f);
    std::fill(ld_tag_sum_r4_above_r2min_.begin(), ld_tag_sum_r4_above_r2min_.end(), 0.0f);
  }
 private:
  TagToSnpMapping& mapping_;

  std::vector<float> ld_sum_r2_below_r2min_;  // master data
  std::vector<float> ld_sum_r2_above_r2min_;
  std::vector<float> ld_sum_r4_above_r2min_;

  std::vector<float> ld_tag_sum_r2_below_r2min_;  // derived data
  std::vector<float> ld_tag_sum_r2_above_r2min_;
  std::vector<float> ld_tag_sum_r4_above_r2min_;
};

class LdMatrixRow; 

// Class to store LD matrix for a given chromosome (or chunk) in CSR format
// Iterpret this class as a mapping from key to val, where "key" is a snp, and "val" is a tag.
class LdMatrixCsrChunk {
 public:
  std::vector<std::tuple<int, int, packed_r_value>> coo_ld_; // key, val, r

  // csr_ld_key_index_.size() == num_keys_in_chunk() + 1; 
  // csr_ld_key_index_[j]..csr_ld_key_index_[j+1] is a range of values in CSR matrix corresponding to j-th variant
  // csr_ld_val_index_.size() == csr_ld_r_.size() == number of non-zero LD r values
  // csr_ld_val_index_ contains values from 0 to num_val_-1
  // csr_ld_r_ contains values from -1 to 1, indicating LD allelic correlation (r) between snp and tag variants
  std::vector<int64_t> csr_ld_key_index_;
  std::vector<uint64_t> csr_ld_val_index_offset_;      // pointers to csr_ld_val_index_packed_ (location where to decompress)
                                                       // number of elements to decompress can be deduced from csr_ld_key_index_
  std::vector<unsigned char> csr_ld_val_index_packed_;  // packed csr_ld_val_index (delta-encoded, then compressed with TurboPFor vsenc32 algorithm).
                                                        // The buffer has some extra capacity (as required by TurboPFor vsenc32/vdec32 algorithms).
  std::vector<packed_r_value> csr_ld_r_;
  
  unsigned char* csr_ld_val_index_packed(int key_index) {
    return &csr_ld_val_index_packed_[csr_ld_val_index_offset_[key_index - key_index_from_inclusive_]];
  }

  int64_t num_ld_r2(int key_index) const {
    return ld_index_end(key_index) - ld_index_begin(key_index);
  }

  int64_t ld_index_begin(int key_index) const {
    return csr_ld_key_index_[key_index - key_index_from_inclusive_];
  }

  int64_t ld_index_end(int key_index) const {
    return csr_ld_key_index_[key_index - key_index_from_inclusive_ + 1];
  }

  // indices where chunk starts (inclusive) and ends (exclusive)
  // [key_index_from_inclusive_, key_index_to_exclusive_)
  int key_index_from_inclusive_;
  int key_index_to_exclusive_;
  int chr_label_;
  int num_keys_in_chunk() const { return key_index_to_exclusive_ - key_index_from_inclusive_; }
  bool is_empty() const { return key_index_to_exclusive_ == key_index_from_inclusive_; }

  int64_t set_ld_r2_csr();
  int64_t validate_ld_r2_csr(const std::vector<uint32_t>& csr_ld_val_index);
  float find_and_retrieve_ld_r2(int key_index, int val_index, const std::vector<uint32_t>& csr_ld_val_index);  // nan if doesn't exist.
  void extract_row(int key_index, LdMatrixRow* row);

  size_t log_diagnostics();
  void clear();
};

class LdMatrixIterator {
public:
  LdMatrixIterator(int64_t ld_index, const LdMatrixRow* parent) : ld_index_(ld_index), parent_(parent) {}

  int index() const;  // this can be tag_index or snp_index, depending on how LdMatrix is stored (as snp->tag mapping or as tag->snp mapping)
  float r() const;
  float r2() const;

  LdMatrixIterator& operator++ () {
    ld_index_++;
    return *this;
  }
  LdMatrixIterator  operator++ (int)
  {
    LdMatrixIterator result(*this);
    ++(*this);
    return result;
  }

private:
  int64_t ld_index_;
  const LdMatrixRow* parent_;
  friend inline bool operator <(const LdMatrixIterator& lhs, const LdMatrixIterator& rhs);
  friend inline int operator -(const LdMatrixIterator& lhs, const LdMatrixIterator& rhs);
};

inline bool operator <(const LdMatrixIterator& lhs, const LdMatrixIterator& rhs)
{
  return (lhs.parent_ == rhs.parent_) && (lhs.ld_index_ < rhs.ld_index_);
}
inline int operator -(const LdMatrixIterator& lhs, const LdMatrixIterator& rhs)
{
  return lhs.ld_index_ - rhs.ld_index_;
}

class LdMatrixRow {
public:
  LdMatrixIterator begin() { return LdMatrixIterator(0, this); }
  LdMatrixIterator end() { return LdMatrixIterator(index_vector_.size(), this); }
private:
  std::vector<int> index_vector_;
  std::vector<packed_r_value> r_;
  friend class LdMatrixIterator;
  friend class LdMatrixCsr;
  friend class LdMatrixCsrChunk;
};

// Class for sparse LD matrix stored in CSR format (Compressed Sparse Row Format)
class LdMatrixCsr {
 public:
   LdMatrixCsr(TagToSnpMapping& mapping) : mapping_(mapping) {}

   int64_t set_ld_r2_coo(int chr_label, int64_t length, int* snp_index, int* snp_other_index, float* r, float r2_min);
   int64_t set_ld_r2_coo_version1plus(int chr_label, const std::string& filename, float r2_min);
   int64_t set_ld_r2_coo_version0(int chr_label, const std::string& filename, float r2_min);
   int64_t set_ld_r2_csr(float r2_min, int chr_label);  // finalize

   void extract_snp_row(SnpIndex snp_index, LdMatrixRow* row);  // retrieve all LD r2 entries for given snp_index or tag_index
   void extract_tag_row(TagIndex tag_index, LdMatrixRow* row);

   int num_ld_r2_snp(int snp_index);  // how many LD r2 entries is there for snp_index

   const LdSum* ld_sum_adjust_for_hvec() { return ld_sum_adjust_for_hvec_.get(); }
   const LdSum* ld_sum() { return ld_sum_.get(); }

   size_t log_diagnostics();
   void clear();
   void init_chunks();
   void init_diagonal(int chr_label);
private:

  TagToSnpMapping& mapping_;
  std::vector<LdMatrixCsrChunk> chunks_forward_;   // split per chromosomes, mapping from snp to tag
  std::vector<LdMatrixCsrChunk> chunks_reverse_;   // mapping from tag to snp
  
  std::shared_ptr<LdSum> ld_sum_adjust_for_hvec_;
  std::shared_ptr<LdSum> ld_sum_;
};
