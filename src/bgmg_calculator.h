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

#include <unordered_map>
#include <memory>

#include <iostream>
#include <sstream>

#include "boost/utility.hpp"

// Singleton class to manage a collection of objects, identifiable with some integer ID.
template<class Type>
class TemplateManager : boost::noncopyable {
 public:
  static TemplateManager<Type>& singleton() {
    static TemplateManager<Type> manager;
    return manager;
  }

  std::shared_ptr<Type> Get(int id) {
    if (map_.find(id) == map_.end()) {
      map_[id] = std::make_shared<Type>();
    }
    return map_[id];
  }

  void Erase(int id) {
    map_.erase(id);
  }

 private:
  TemplateManager() { }  // Singleton (make constructor private)
  std::unordered_map<int, std::shared_ptr<Type>> map_;
};

template<typename T>
class DenseMatrix {
 public:
  explicit DenseMatrix(int no_rows = 0, int no_columns = 0, bool store_by_rows = true)
    : no_rows_(no_rows),
    no_columns_(no_columns),
    store_by_rows_(store_by_rows),
    data_(nullptr) {
    if (no_rows > 0 && no_columns > 0) {
      data_ = new T[no_rows_ * no_columns_];
    }
  }

  DenseMatrix(const DenseMatrix<T>& src_matrix) {
    no_rows_ = src_matrix.no_rows();
    no_columns_ = src_matrix.no_columns();
    store_by_rows_ = src_matrix.store_by_rows_;
    if (no_columns_ >0 && no_rows_ > 0) {
      data_ = new T[no_rows_ * no_columns_];

      for (int i = 0; i < no_rows_ * no_columns_; ++i) {
        data_[i] = src_matrix.get_data()[i];
      }
    } else {
      data_ = nullptr;
    }
  }

  virtual ~DenseMatrix() {
    delete[] data_;
  }

  void InitializeZeros() {
    memset(data_, 0, sizeof(T)* no_rows_ * no_columns_);
  }

  T& operator() (int index_row, int index_col) {
    if (store_by_rows_) {
      return data_[index_row * no_columns_ + index_col];
    }
    return data_[index_col * no_rows_ + index_row];
  }

  const T& operator() (int index_row, int index_col) const {
    assert(index_row < no_rows_);
    assert(index_col < no_columns_);
    if (store_by_rows_) {
      return data_[index_row * no_columns_ + index_col];
    }
    return data_[index_col * no_rows_ + index_row];
  }

  DenseMatrix<T>& operator= (const DenseMatrix<T>& src_matrix) {
    no_rows_ = src_matrix.no_rows();
    no_columns_ = src_matrix.no_columns();
    store_by_rows_ = src_matrix.store_by_rows_;
    if (data_ != nullptr) {
      delete[] data_;
    }
    if (no_columns_ >0 && no_rows_ > 0) {
      data_ = new  T[no_rows_ * no_columns_];

      for (int i = 0; i < no_rows_ * no_columns_; ++i) {
        data_[i] = src_matrix.get_data()[i];
      }
    } else {
      data_ = nullptr;
    }

    return *this;
  }

  int no_rows() const { return no_rows_; }
  int no_columns() const { return no_columns_; }
  int size() const { return no_rows_ * no_columns_; }
  bool is_equal_size(const DenseMatrix<T>& rhs) const {
    return no_rows_ == rhs.no_rows_ && no_columns_ == rhs.no_columns_;
  }

  T* get_data() {
    return data_;
  }

  const T* get_data() const {
    return data_;
  }

 private:
  int no_rows_;
  int no_columns_;
  bool store_by_rows_;
  T* data_;
};

class BgmgCalculator {
 public:
  BgmgCalculator() : num_snp_(-1), num_tag_(-1), k_max_(100), r2_min_(0.0), max_causals_(100000) { status_.push_back(std::stringstream()); }
  
  // num_snp = total size of the reference (e.i. the total number of genotyped variants)
  // num_tag = number of tag variants to include in the inference (must be a subset of the reference)
  // indices = array of size num_tag, containing indices from 0 to num_snp-1
  // NB: all tag variants must have defined zvec, hvec, hvec and weights.
  int64_t set_tag_indices(int num_snp, int num_tag, int* tag_indices);

  // consume input in plink format, e.i. lower triangular LD r2 matrix
  // - snp_index is less than tag_index;
  // - does not contain r2=1.0 of the variant with itself
  // must be called after set_tag_indices
  // must be called one for each chromosome, sequentially, starting from lowest chromosome number
  // non-tag variants will be ignored
  int64_t set_ld_r2_coo(int length, int* snp_index, int* tag_index, float* r2);
  int64_t set_ld_r2_csr();  // finalize

  // must be called after set_ld_r2, as it adjusts r2 matrix
  // one value for each snp (tag and non-tag)
  int64_t set_hvec(int length, float* values);
  
  // zvec, nvec, weights for tag variants
  // all values must be defined
  int64_t set_zvec(int trait, int length, float* values);
  int64_t set_nvec(int trait, int length, float* values);
  int64_t set_weights(int length, float* values);

  int64_t find_snp_sample();

  int64_t set_option(char* option, double value);
  const char* status();

  double calc_univariate_cost(float pi_vec, float sig2_zero, float sig2_beta);
  double calc_bivariate_cost(int num_components, double* pi_vec, double* sig2_beta, double* rho_beta, double* sig2_zero, double rho_zero);
 private:

  int num_snp_;
  int num_tag_;
  std::vector<int> tag_to_snp_; // 0..num_snp_-1, size=num_tag_
  std::vector<int> snp_to_tag_; // 0..num_tag-1,  size=num_snp_, -1 indicate non-tag snp
  std::vector<char> is_tag_;    // true or false, size=num_snp_, is_tag_[i] == (snp_to_tag_[i] != -1)

  // csr_ld_snp_index_.size() == num_snp_ + 1; 
  // csr_ld_snp_index_[j]..csr_ld_snp_index_[j+1] is a range of values in CSR matrix corresponding to j-th variant
  // csr_ld_tag_index_.size() == csr_ld_r2_.size() == number of non-zero LD r2 values
  // csr_ld_tag_index_ contains values from 0 to num_tag_-1
  // csr_ld_r2_ contains values from 0 to 1, indicating LD r2 between snp and tag variants
  std::vector<int> csr_ld_snp_index_;
  std::vector<int> csr_ld_tag_index_;
  std::vector<float> csr_ld_r2_;
  std::vector<std::tuple<int, int, float>> coo_ld_; // snp, tag, r2

  // all stored for for tag variants (only)
  std::vector<float> zvec1_;
  std::vector<float> nvec1_;
  std::vector<float> weights_;

  std::shared_ptr<DenseMatrix<int>> snp_sample_;  // permutation matrix; #rows = pimax*num_snp; #cols=k_max_

  int k_max_;
  int max_causals_;
  float r2_min_;
  void set_num_snp(int length);
  void set_num_tag(int length);

  std::vector<std::stringstream> status_;
};

typedef TemplateManager<BgmgCalculator> BgmgCalculatorManager;