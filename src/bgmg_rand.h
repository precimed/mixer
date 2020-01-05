#pragma once

#include <random>

#if defined(__AVX2__) 
#include "SIMDxorshift/simdxorshift128plus.h"
class avx_xorshift128plus_engine {
 public:
  avx_xorshift128plus_engine(uint64_t key1, uint64_t key2) {
    index_ = 8;
    avx_xorshift128plus_init(key1, key2, &key_);
    avx_xorshift128plus_jump(&key_);  // without "jump" xorshift128+ seem to generate few zeros on first iteration (when keys are small numbers)
  }
  typedef uint32_t result_type;
  static constexpr result_type min() { return result_type(0UL); }
  static constexpr result_type max() { return result_type(~result_type(0UL)); }
  inline result_type operator()() { 
    if (index_==8) {
      __m256i R = avx_xorshift128plus(&key_);
      _mm256_storeu_si256((__m256i *) randomsource_, R);
      index_=0;
    }
    return randomsource_[index_++];
  }
  avx_xorshift128plus_key_t* key() { return &key_; }

 private:
  avx_xorshift128plus_key_t key_;
  uint32_t randomsource_[8];
  uint32_t index_;
};
class SubsetSampler {
 public:
  SubsetSampler(uint64_t key1, uint64_t key2, uint32_t size) : data_(size, 0), generator_(key1, key2) {
    for (uint32_t i = 0; i < size; i++) data_[i]=i;
  }

  // sample elements from [0, 1, ..., size) with probability p
  // return number of selected elements (count)
  // NB! stores output in the END of data_ array, i.e. data[size-count, .., size-1] 
  inline uint32_t sample_shuffle(double p) {
    uint32_t count = std::binomial_distribution<uint32_t>(data_.size(), p)(generator_);
    avx_xorshift128plus_shuffle32_partial(generator_.key(), &data_[0], data_.size(), data_.size() - count);
    return count;
  }

  const uint32_t* data() const { return &data_[0]; }

 private:
  std::vector<uint32_t> data_;
  avx_xorshift128plus_engine generator_;
};
class MultinomialSampler {
 public:
  MultinomialSampler(uint64_t key1, uint64_t key2, uint32_t n, uint32_t len) : data_(n, 0), p_(len, 0), counts_(len, 0), generator_(key1, key2) {
    for (uint32_t i = 0; i < n; i++) data_[i]=i;
  }

  inline uint32_t sample_shuffle() {
    uint32_t total_sampled = 0;
    const int len = p_.size();
    for (int i = 0; i < len; i++) {
      counts_[i] = std::binomial_distribution<uint32_t>(data_.size() - total_sampled, p_[i])(generator_);
      const double factor = 1.0 / (1.0-p_[i]);
      for (int j = (i + 1); j < len; j++) p_[j] *= factor;
      total_sampled += counts_[i];
    }

    avx_xorshift128plus_shuffle32_partial(generator_.key(), &data_[0], data_.size(), data_.size() - total_sampled);
    return total_sampled;
  }

  const uint32_t* data() const { return &data_[0]; }
  double* p() { return &p_[0]; }
  const uint32_t* counts() const { return &counts_[0]; }  

 private:
  std::vector<uint32_t> data_;
  std::vector<double> p_;
  std::vector<uint32_t> counts_;
  avx_xorshift128plus_engine generator_;
};
#else
#include "SIMDxorshift/xorshift128plus.h"
class xorshift128plus_engine {
 public:
  xorshift128plus_engine(uint64_t key1, uint64_t key2) {
    index_ = 2;
    xorshift128plus_init(key1, key2, &key_);
    xorshift128plus_jump(&key_);  // without "jump" xorshift128+ seem to generate few zeros on first iteration (when keys are small numbers)
  }
  typedef uint32_t result_type;
  static constexpr result_type min() { return result_type(0UL); }
  static constexpr result_type max() { return result_type(~result_type(0UL)); }
  inline result_type operator()() { 
    if (index_==2) {
      uint64_t R = xorshift128plus(&key_);
      randomsource_[0] = R & UINT64_C(0xFFFFFFFF);
      randomsource_[1] = R >> 32;
      index_=0;
    }
    return randomsource_[index_++];
  }
  xorshift128plus_key_t* key() { return &key_; }

 private:
  xorshift128plus_key_t key_;
  uint32_t randomsource_[2];
  uint32_t index_;
};
class SubsetSampler {
 public:
  SubsetSampler(uint64_t key1, uint64_t key2, uint32_t size) : data_(size, 0), generator_(key1, key2) {
    for (uint32_t i = 0; i < size; i++) data_[i]=i;
  }

  // sample elements from [0, 1, ..., size) with probability p
  // return number of selected elements (count)
  // NB! stores output in the END of data_ array, i.e. data[size-count, .., size-1] 
  inline uint32_t sample_shuffle(double p) {
    uint32_t count = std::binomial_distribution<uint32_t>(data_.size(), p)(generator_);
    xorshift128plus_shuffle32_partial(generator_.key(), &data_[0], data_.size(), data_.size() - count);
    return count;
  }

  const uint32_t* data() const { return &data_[0]; }

 private:
  std::vector<uint32_t> data_;
  xorshift128plus_engine generator_;
};

class MultinomialSampler {
 public:
  MultinomialSampler(uint64_t key1, uint64_t key2, uint32_t n, uint32_t len) : data_(n, 0), p_(len, 0), counts_(len, 0), generator_(key1, key2) {
    for (uint32_t i = 0; i < n; i++) data_[i]=i;
  }

  // multinomial sample elements from [0, 1, ..., size) with probabilities given by vector p_ of length len, where len = this->p_.size()
  // Effect:
  // - NB! Values in this->p_ array will be updated with some new values - make sure to restore original values if you call sample_shuffle repeatedly)
  // - this->counts_[i] is populated with the number of elements sampled with probability this->p_[i]
  // - the function returns total number of elements sampled, i.e. "total_count = sum_i this->counts_[i]"
  // - the elements are patrly shuffled so that the last total_count elements in this->data_ are now random
  // Note:
  // - p must be non-negative, and sum of elements in p must not exceed 1
  // - the API of MultinomialSampler class is ugly, partly because it's important to avoid array allocations both here and at caller's site.
  // Algorithm:
  // - https://math.stackexchange.com/questions/934941/conditional-probability-in-multinomial-distribution
  // - https://github.com/numba/numba/blob/master/numba/targets/randomimpl.py#L1433-L1508
  inline uint32_t sample_shuffle() {
    uint32_t total_sampled = 0;
    const int len = p_.size();
    for (int i = 0; i < len; i++) {
      counts_[i] = std::binomial_distribution<uint32_t>(data_.size() - total_sampled, p_[i])(generator_);
      const double factor = 1.0 / (1.0-p_[i]);
      for (int j = (i + 1); j < len; j++) p_[j] *= factor;
      total_sampled += counts_[i];
    }

    xorshift128plus_shuffle32_partial(generator_.key(), &data_[0], data_.size(), data_.size() - total_sampled);
    return total_sampled;
  }

  const uint32_t* data() const { return &data_[0]; }
  double* p() { return &p_[0]; }
  const uint32_t* counts() const { return &counts_[0]; }  

 private:
  std::vector<uint32_t> data_;
  std::vector<double> p_;
  std::vector<uint32_t> counts_;
  xorshift128plus_engine generator_;
};

#endif





