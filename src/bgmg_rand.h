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
#endif





