// This file is part of PLINK 1.90, copyright (C) 2005-2019 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// The original source code from PLINK was modified by Oleksandr Frei, April 2019

#include "plink_ld.h"
#include "snp_lookup.h"


#define MULTIPLEX_LD 1920

#ifdef __LP64__
static inline void ld_dot_prod_batch(__m128i* vec1, __m128i* vec2, __m128i* mask1, __m128i* mask2, int32_t* return_vals, uint32_t iters) {
  // Main routine for computation of \sum_i^M (x_i - \mu_x)(y_i - \mu_y), where
  // x_i, y_i \in \{-1, 0, 1\}, but there are missing values.
  //
  //
  // We decompose this sum into
  //   \sum_i x_iy_i - \mu_y\sum_i x_i - \mu_x\sum_i y_i +
  //   (M - # missing)\mu_x\mu_y.
  // *Without* missing values, this can be handled very cleanly.  The last
  // three terms can all be precomputed, and \sum_i x_iy_i can be handled in a
  // manner very similar to bitwise Hamming distance.  This is several times as
  // fast as the lookup tables used for relationship matrices.
  //
  // Unfortunately, when missing values are present,
  // \mu_y\sum_{i: nonmissing from y} x_i and
  // \mu_x\sum_{i: nonmissing from x} y_i must also be evaluated (and, in
  // practice, \mu_y\sum_{i: nonmissing from y} x_i^2 and
  // \mu_x\sum_{i: nonmissing from x} y_i^2 should be determined here as well);
  // this removes much of the speed advantage, and the best applications of the
  // underlying ternary dot product algorithm used here lie elsewhere.
  // Nevertheless, it is still faster, so we use it.
  // (possible todo: accelerated function when there really are no missing
  // values, similar to what is now done for --fast-epistasis)
  //
  //
  // Input:
  // * vec1 and vec2 are encoded -1 -> 00, 0/missing -> 01, 1 -> 10.
  // * mask1 and mask2 mask out missing values (i.e. 00 for missing, 11 for
  //   nonmissing).
  // * return_vals provides space for return values.
  // * iters is the number of 48-byte windows to process, anywhere from 1 to 10
  //   inclusive.
  //
  // This function performs the update
  //   return_vals[0] += (-N) + \sum_i x_iy_i
  //   return_vals[1] += N_y + \sum_{i: nonmissing from y} x_i
  //   return_vals[2] += N_x + \sum_{i: nonmissing from x} y_i
  //   return_vals[3] += N_y - \sum_{i: nonmissing from y} x_i^2
  //   return_vals[4] += N_x - \sum_{i: nonmissing from x} y_i^2
  // where N is the number of samples processed after applying the missingness
  // masks indicated by the subscripts.
  //
  // Computation of terms [1]-[4] is based on the identity
  //   N_y + \sum_{i: nonmissing from y} x_i = popcount2(vec1 & mask2)
  // where "popcount2" refers to starting with two-bit integers instead of
  // one-bit integers in our summing process (this allows us to skip a few
  // operations).  (Once we can assume the presence of hardware popcount, a
  // slightly different implementation may be better.)
  //
  // The trickier [0] computation currently proceeds as follows:
  //
  // 1. zcheck := (vec1 | vec2) & 0x5555...
  // Detects whether at least one member of the pair has a 0/missing value.
  //
  // 2. popcount2(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
  // Subtracting this *from* a bias will give us our desired \sum_i x_iy_i dot
  // product.
  //
  // MULTIPLEX_LD sets of values are usually handled per function call.  If
  // fewer values are present, the ends of all input vectors should be zeroed
  // out.

  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader1;
  __m128i loader2;
  __m128i sum1;
  __m128i sum2;
  __m128i sum11;
  __m128i sum22;
  __m128i sum12;
  __m128i tmp_sum1;
  __m128i tmp_sum2;
  __m128i tmp_sum12;
  __univec acc;
  __univec acc1;
  __univec acc2;
  __univec acc11;
  __univec acc22;
  acc.vi = _mm_setzero_si128();
  acc1.vi = _mm_setzero_si128();
  acc2.vi = _mm_setzero_si128();
  acc11.vi = _mm_setzero_si128();
  acc22.vi = _mm_setzero_si128();
  do {
    loader1 = *vec1++;
    loader2 = *vec2++;
    sum1 = *mask2++;
    sum2 = *mask1++;
    sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    // sum11 = _mm_and_si128(_mm_and_si128(_mm_xor_si128(sum1, m1), m1), loader1);
    // sum22 = _mm_and_si128(_mm_and_si128(_mm_xor_si128(sum2, m1), m1), loader2);
    sum1 = _mm_and_si128(sum1, loader1);
    sum2 = _mm_and_si128(sum2, loader2);
    sum11 = _mm_and_si128(sum1, m1);
    sum22 = _mm_and_si128(sum2, m1);
    // use andnot to eliminate need for 0xaaaa... to occupy an xmm register
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12), _mm_xor_si128(loader1, loader2));
    sum12 = _mm_or_si128(sum12, loader1);

    // sum1, sum2, and sum12 now store the (biased) two-bit sums of
    // interest; merge to 4 bits to prevent overflow.  this merge can be
    // postponed for sum11 and sum22 because the individual terms are 0/1
    // instead of 0/1/2.
    sum1 = _mm_add_epi64(_mm_and_si128(sum1, m2), _mm_and_si128(_mm_srli_epi64(sum1, 2), m2));
    sum2 = _mm_add_epi64(_mm_and_si128(sum2, m2), _mm_and_si128(_mm_srli_epi64(sum2, 2), m2));
    sum12 = _mm_add_epi64(_mm_and_si128(sum12, m2), _mm_and_si128(_mm_srli_epi64(sum12, 2), m2));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
    tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
    sum11 = _mm_add_epi64(sum11, _mm_and_si128(tmp_sum1, m1));
    sum22 = _mm_add_epi64(sum22, _mm_and_si128(tmp_sum2, m1));
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

    sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
    sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
    tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
    sum11 = _mm_add_epi64(sum11, _mm_and_si128(tmp_sum1, m1));
    sum22 = _mm_add_epi64(sum22, _mm_and_si128(tmp_sum2, m1));
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

    sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
    sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
    sum11 = _mm_add_epi64(_mm_and_si128(sum11, m2), _mm_and_si128(_mm_srli_epi64(sum11, 2), m2));
    sum22 = _mm_add_epi64(_mm_and_si128(sum22, m2), _mm_and_si128(_mm_srli_epi64(sum22, 2), m2));
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    acc1.vi = _mm_add_epi64(acc1.vi, _mm_add_epi64(_mm_and_si128(sum1, m4), _mm_and_si128(_mm_srli_epi64(sum1, 4), m4)));
    acc2.vi = _mm_add_epi64(acc2.vi, _mm_add_epi64(_mm_and_si128(sum2, m4), _mm_and_si128(_mm_srli_epi64(sum2, 4), m4)));
    acc11.vi = _mm_add_epi64(acc11.vi, _mm_add_epi64(_mm_and_si128(sum11, m4), _mm_and_si128(_mm_srli_epi64(sum11, 4), m4)));
    acc22.vi = _mm_add_epi64(acc22.vi, _mm_add_epi64(_mm_and_si128(sum22, m4), _mm_and_si128(_mm_srli_epi64(sum22, 4), m4)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sum12, m4), _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));
  } while (--iters);
  // moved down because we've almost certainly run out of xmm registers
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
#if MULTIPLEX_LD > 960
  acc1.vi = _mm_add_epi64(_mm_and_si128(acc1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
  acc2.vi = _mm_add_epi64(_mm_and_si128(acc2.vi, m8), _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  acc1.vi = _mm_and_si128(_mm_add_epi64(acc1.vi, _mm_srli_epi64(acc1.vi, 8)), m8);
  acc2.vi = _mm_and_si128(_mm_add_epi64(acc2.vi, _mm_srli_epi64(acc2.vi, 8)), m8);
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  acc11.vi = _mm_and_si128(_mm_add_epi64(acc11.vi, _mm_srli_epi64(acc11.vi, 8)), m8);
  acc22.vi = _mm_and_si128(_mm_add_epi64(acc22.vi, _mm_srli_epi64(acc22.vi, 8)), m8);

  return_vals[0] -= ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[1] += ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[2] += ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[3] += ((acc11.u8[0] + acc11.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[4] += ((acc22.u8[0] + acc22.u8[1]) * 0x1000100010001LLU) >> 48;
}

void ld_dot_prod(uintptr_t* vec1, uintptr_t* vec2, uintptr_t* mask1, uintptr_t* mask2, int32_t* return_vals, uint32_t batch_ct_m1, uint32_t last_batch_size) {
  while (batch_ct_m1--) {
    ld_dot_prod_batch((__m128i*)vec1, (__m128i*)vec2, (__m128i*)mask1, (__m128i*)mask2, return_vals, MULTIPLEX_LD / 192);
    vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
    vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
    mask1 = &(mask1[MULTIPLEX_LD / BITCT2]);
    mask2 = &(mask2[MULTIPLEX_LD / BITCT2]);
  }
  ld_dot_prod_batch((__m128i*)vec1, (__m128i*)vec2, (__m128i*)mask1, (__m128i*)mask2, return_vals, last_batch_size);
}

static inline int32_t ld_dot_prod_nm_batch(__m128i* vec1, __m128i* vec2, uint32_t iters) {
  // faster ld_dot_prod_batch() for no-missing-calls case.
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader1;
  __m128i loader2;
  __m128i sum12;
  __m128i tmp_sum12;
  __univec acc;
  acc.vi = _mm_setzero_si128();
  do {
    loader1 = *vec1++;
    loader2 = *vec2++;
    sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12), _mm_xor_si128(loader1, loader2));
    sum12 = _mm_or_si128(sum12, loader1);
    sum12 = _mm_add_epi64(_mm_and_si128(sum12, m2), _mm_and_si128(_mm_srli_epi64(sum12, 2), m2));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sum12, m4), _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));
  } while (--iters);
#if MULTIPLEX_LD > 960
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  return (uint32_t)(((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48);
}

int32_t ld_dot_prod_nm(uintptr_t* vec1, uintptr_t* vec2, uint32_t founder_ct, uint32_t batch_ct_m1, uint32_t last_batch_size) {
  // accelerated implementation for no-missing-loci case
  int32_t result = (int32_t)founder_ct;
  while (batch_ct_m1--) {
    result -= ld_dot_prod_nm_batch((__m128i*)vec1, (__m128i*)vec2, MULTIPLEX_LD / 192);
    vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
    vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
  }
  result -= ld_dot_prod_nm_batch((__m128i*)vec1, (__m128i*)vec2, last_batch_size);
  return result;
}
#else  // __LP64__
static inline void ld_dot_prod_batch(uintptr_t* vec1, uintptr_t* vec2, uintptr_t* mask1, uintptr_t* mask2, int32_t* return_vals, uint32_t iters) {
  uint32_t final_sum1 = 0;
  uint32_t final_sum2 = 0;
  uint32_t final_sum11 = 0;
  uint32_t final_sum22 = 0;
  uint32_t final_sum12 = 0;
  uintptr_t loader1;
  uintptr_t loader2;
  uintptr_t sum1;
  uintptr_t sum2;
  uintptr_t sum11;
  uintptr_t sum22;
  uintptr_t sum12;
  uintptr_t tmp_sum1;
  uintptr_t tmp_sum2;
  uintptr_t tmp_sum12;
  do {
    // (The important part of the header comment on the 64-bit version is
    // copied below.)
    //
    // Input:
    // * vec1 and vec2 are encoded -1 -> 00, 0/missing -> 01, 1 -> 10.
    // * mask1 and mask2 mask out missing values (i.e. 00 for missing, 11 for
    //   nonmissing).
    // * return_vals provides space for return values.
    // * iters is the number of 12-byte windows to process, anywhere from 1 to
    //   40 inclusive.  (No, this is not the interface you'd use for a
    //   general-purpose library.)  [32- and 64-bit differ here.]
    //
    // This function performs the update
    //   return_vals[0] += (-N) + \sum_i x_iy_i
    //   return_vals[1] += N_y + \sum_{i: nonmissing from y} x_i
    //   return_vals[2] += N_x + \sum_{i: nonmissing from x} y_i
    //   return_vals[3] += N_y - \sum_{i: nonmissing from y} x_i^2
    //   return_vals[4] += N_x - \sum_{i: nonmissing from x} y_i^2
    // where N is the number of samples processed after applying the
    // missingness masks indicated by the subscripts.
    //
    // Computation of terms [1]-[4] is based on the identity
    //   N_y + \sum_{i: nonmissing from y} x_i = popcount2(vec1 & mask2)
    // where "popcount2" refers to starting with two-bit integers instead of
    // one-bit integers in our summing process (this allows us to skip a few
    // operations).  (Once we can assume the presence of hardware popcount, a
    // slightly different implementation may be better.)
    //
    // The trickier [0] computation currently proceeds as follows:
    //
    // 1. zcheck := (vec1 | vec2) & 0x5555...
    // Detects whether at least one member of the pair has a 0/missing value.
    //
    // 2. popcount2(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
    // Subtracting this *from* a bias will give us our desired \sum_i x_iy_i
    // dot product.


    loader1 = *vec1++;
    loader2 = *vec2++;
    sum1 = *mask2++;
    sum2 = *mask1++;
    sum12 = (loader1 | loader2) & FIVEMASK;

    sum1 = sum1 & loader1;
    sum2 = sum2 & loader2;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
    sum12 = sum12 | loader1;
    sum11 = sum1 & FIVEMASK;
    sum22 = sum2 & FIVEMASK;

    sum1 = (sum1 & 0x33333333) + ((sum1 >> 2) & 0x33333333);
    sum2 = (sum2 & 0x33333333) + ((sum2 >> 2) & 0x33333333);
    sum12 = (sum12 & 0x33333333) + ((sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;
    tmp_sum1 = tmp_sum1 & loader1;
    tmp_sum2 = tmp_sum2 & loader2;

    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;
    sum11 += tmp_sum1 & FIVEMASK;
    sum22 += tmp_sum2 & FIVEMASK;

    sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
    sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;

    tmp_sum1 = tmp_sum1 & loader1;
    tmp_sum2 = tmp_sum2 & loader2;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;
    sum11 += tmp_sum1 & FIVEMASK;
    sum22 += tmp_sum2 & FIVEMASK;

    sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
    sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
    sum11 = (sum11 & 0x33333333) + ((sum11 >> 2) & 0x33333333);
    sum22 = (sum22 & 0x33333333) + ((sum22 >> 2) & 0x33333333);
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    sum1 = (sum1 & 0x0f0f0f0f) + ((sum1 >> 4) & 0x0f0f0f0f);
    sum2 = (sum2 & 0x0f0f0f0f) + ((sum2 >> 4) & 0x0f0f0f0f);
    sum11 = (sum11 & 0x0f0f0f0f) + ((sum11 >> 4) & 0x0f0f0f0f);
    sum22 = (sum22 & 0x0f0f0f0f) + ((sum22 >> 4) & 0x0f0f0f0f);
    sum12 = (sum12 & 0x0f0f0f0f) + ((sum12 >> 4) & 0x0f0f0f0f);

    // technically could do the multiply-and-shift only once every two rounds
    final_sum1 += (sum1 * 0x01010101) >> 24;
    final_sum2 += (sum2 * 0x01010101) >> 24;
    final_sum11 += (sum11 * 0x01010101) >> 24;
    final_sum22 += (sum22 * 0x01010101) >> 24;
    final_sum12 += (sum12 * 0x01010101) >> 24;
  } while (--iters);
  return_vals[0] -= final_sum12;
  return_vals[1] += final_sum1;
  return_vals[2] += final_sum2;
  return_vals[3] += final_sum11;
  return_vals[4] += final_sum22;
}

void ld_dot_prod(uintptr_t* vec1, uintptr_t* vec2, uintptr_t* mask1, uintptr_t* mask2, int32_t* return_vals, uint32_t batch_ct_m1, uint32_t last_batch_size) {
  while (batch_ct_m1--) {
    ld_dot_prod_batch(vec1, vec2, mask1, mask2, return_vals, MULTIPLEX_LD / 48);
    vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
    vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
    mask1 = &(mask1[MULTIPLEX_LD / BITCT2]);
    mask2 = &(mask2[MULTIPLEX_LD / BITCT2]);
  }
  ld_dot_prod_batch(vec1, vec2, mask1, mask2, return_vals, last_batch_size);
}

static inline int32_t ld_dot_prod_nm_batch(uintptr_t* vec1, uintptr_t* vec2, uint32_t iters) {
  uint32_t final_sum12 = 0;
  uintptr_t loader1;
  uintptr_t loader2;
  uintptr_t sum12;
  uintptr_t tmp_sum12;
  do {
    loader1 = *vec1++;
    loader2 = *vec2++;
    sum12 = (loader1 | loader2) & FIVEMASK;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
    sum12 = sum12 | loader1;
    sum12 = (sum12 & 0x33333333) + ((sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);
    sum12 = (sum12 & 0x0f0f0f0f) + ((sum12 >> 4) & 0x0f0f0f0f);

    final_sum12 += (sum12 * 0x01010101) >> 24;
  } while (--iters);
  return final_sum12;
}

int32_t ld_dot_prod_nm(uintptr_t* vec1, uintptr_t* vec2, uint32_t founder_ct, uint32_t batch_ct_m1, uint32_t last_batch_size) {
  int32_t result = (int32_t)founder_ct;
  while (batch_ct_m1--) {
    result -= ld_dot_prod_nm_batch(vec1, vec2, MULTIPLEX_LD / 48);
    vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
    vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
  }
  result -= ld_dot_prod_nm_batch(vec1, vec2, last_batch_size);
  return result;
}

#endif // __LP64__

void ld_process_load2(uintptr_t* geno_buf, uintptr_t* mask_buf, uint32_t* missing_ct_ptr, uint32_t founder_ct, uint32_t is_x, uintptr_t* founder_male_include2) {
  // ld_process_load(), except no missing_buf[] to conserve memory (and no
  // --ld-xchr 3 support yet), and no zero-variance check (we just want to
  // dump nans in that case)
  uintptr_t* geno_ptr = geno_buf;
  uintptr_t founder_ctl2 = QUATERCT_TO_WORDCT(founder_ct);
  uintptr_t* geno_end = &(geno_buf[founder_ctl2]);
  uintptr_t* mask_buf_ptr = mask_buf;
  uintptr_t cur_geno;
  uintptr_t shifted_masked_geno;
  uintptr_t new_geno;
  uintptr_t new_mask;
  do {
    // Desired encodings:
    // new_geno: nonset homozygote -> 00
    //           het/missing       -> 01
    //           set homozygote    -> 10
    // Given PLINK encoding xx, this is (xx - ((xx >> 1) & FIVEMASK)).
    //
    // new_mask: missing   -> 00
    //           otherwise -> 11
    // ...and this is (((xx >> 1) & FIVEMASK) | ((~xx) & FIVEMASK)) * 3.
    //
    // new_missing: missing   -> 1
    //              otherwise -> 0
    // This can be assembled via repeated CTZLU on ~new_mask.

    cur_geno = *geno_ptr;
    shifted_masked_geno = (cur_geno >> 1) & FIVEMASK;
    new_geno = cur_geno - shifted_masked_geno;
    *geno_ptr++ = new_geno;
    new_mask = (((~cur_geno) & FIVEMASK) | shifted_masked_geno) * 3;
    *mask_buf_ptr++ = new_mask;
  } while (geno_ptr < geno_end);
  if (is_x) {
    geno_ptr = geno_buf;
    do {
      new_geno = *geno_ptr;
      *geno_ptr++ = new_geno + ((~(new_geno | (new_geno >> 1))) & (*founder_male_include2++));
    } while (geno_ptr < geno_end);
  }
  if (founder_ct % BITCT2) {
    mask_buf[founder_ct / BITCT2] &= (ONELU << (2 * (founder_ct % BITCT2))) - ONELU;
  }
  *missing_ct_ptr = founder_ct - (popcount_longs(mask_buf, founder_ctl2) / 2);
}

uint32_t ld_missing_ct_intersect(uintptr_t* lptr1, uintptr_t* lptr2, uintptr_t word12_ct, uintptr_t word12_rem, uintptr_t lshift_last) {
  // variant of popcount_longs_intersect()
  uintptr_t tot = 0;
  uintptr_t* lptr1_end2;
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i* vptr1 = (__m128i*)lptr1;
  __m128i* vptr2 = (__m128i*)lptr2;
  __m128i* vend1;
  __m128i loader1;
  __m128i loader2;
  __univec acc;

  while (word12_ct >= 10) {
    word12_ct -= 10;
    vend1 = &(vptr1[60]);
  ld_missing_ct_intersect_main_loop:
    acc.vi = _mm_setzero_si128();
    do {
      loader1 = _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1);
      loader2 = _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1);
      loader1 = _mm_add_epi64(loader1, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
      loader2 = _mm_add_epi64(loader2, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
      loader1 = _mm_add_epi64(loader1, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
      loader2 = _mm_add_epi64(loader2, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
      loader1 = _mm_add_epi64(_mm_and_si128(loader1, m2), _mm_and_si128(_mm_srli_epi64(loader1, 2), m2));
      loader1 = _mm_add_epi64(loader1, _mm_add_epi64(_mm_and_si128(loader2, m2), _mm_and_si128(_mm_srli_epi64(loader2, 2), m2)));
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(loader1, m4), _mm_and_si128(_mm_srli_epi64(loader1, 4), m4)));
    } while (vptr1 < vend1);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (word12_ct) {
    vend1 = &(vptr1[word12_ct * 6]);
    word12_ct = 0;
    goto ld_missing_ct_intersect_main_loop;
  }
  lptr1 = (uintptr_t*)vptr1;
  lptr2 = (uintptr_t*)vptr2;
#else
  uintptr_t* lptr1_end = &(lptr1[word12_ct * 12]);
  uintptr_t tmp_stor;
  uintptr_t loader1;
  uintptr_t loader2;
  while (lptr1 < lptr1_end) {
    loader1 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
    loader1 += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
    tmp_stor = (loader1 & 0x0f0f0f0f) + ((loader1 >> 4) & 0x0f0f0f0f);

    loader1 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
    loader1 = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
    loader1 += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
    tmp_stor += (loader1 & 0x0f0f0f0f) + ((loader1 >> 4) & 0x0f0f0f0f);
    tot += (tmp_stor * 0x01010101) >> 24;
  }
#endif
  lptr1_end2 = &(lptr1[word12_rem]);
  while (lptr1 < lptr1_end2) {
    tot += popcount2_long((~((*lptr1++) | (*lptr2++))) & FIVEMASK);
  }
  if (lshift_last) {
    tot += popcount2_long(((~((*lptr1) | (*lptr2))) & FIVEMASK) << lshift_last);
  }
  return tot;
}

SampleCountInfo::SampleCountInfo(int num_subjects) {
  unfiltered_sample_ct = num_subjects;
  unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  unfiltered_sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(unfiltered_sample_ct);

  founder_ct = unfiltered_sample_ct;  // here we only work with simple plink files where everyone is founder
  final_mask = get_final_mask(founder_ct);

  founder_ct_mld = (founder_ct + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
  founder_ct_mld_m1 = ((uint32_t)founder_ct_mld) - 1;
#ifdef __LP64__
  founder_ct_mld_rem = (MULTIPLEX_LD / 192) - (founder_ct_mld * MULTIPLEX_LD - founder_ct) / 192;
#else
  founder_ct_mld_rem = (MULTIPLEX_LD / 48) - (founder_ct_mld * MULTIPLEX_LD - founder_ct) / 48;
#endif
  founder_ct_192_long = founder_ct_mld_m1 * (MULTIPLEX_LD / BITCT2) + founder_ct_mld_rem * (192 / BITCT2);

  founder_ctwd = founder_ct / BITCT2;
  founder_ctwd12 = founder_ctwd / 12;
  founder_ctwd12_rem = founder_ctwd - (12 * founder_ctwd12);
  lshift_last = 2 * ((0x7fffffc0 - founder_ct) % BITCT2);

  unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
}

void single_marker_3freqs(size_t sample_ctl2, uintptr_t* lptr, uintptr_t* include2, uint32_t* hom2p, uint32_t* hetp, uint32_t* missingp) {
  // adapted from the original plink_assoc.c::single_marker_cc_3freqs
  // Counts the number of heterozygotes, A2 homozygotes, and missing calls.  Assumes marker is diploid.  The caller is
  // expected to calculate the A1 allele count.
  // See single_marker_freqs_and_hwe() for discussion.
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & B
  //   popcount(C) = homozyg major ct
  //   popcount(B) = het ct + homozyg major ct
  //   popcount(A) = missing_ct + homozyg major ct
  // hom2: popcount(C)
  // het: popcount(B) - popcount(C)
  // missing: popcount(A) - popcount(C)

  uint32_t tot_a = 0;
  uint32_t tot_b = 0;
  uint32_t tot_c = 0;
  uintptr_t* lptr_end = &(lptr[sample_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
#ifdef __LP64__
  uintptr_t cur_decr = 120;
  uintptr_t* lptr_12x_end;
  sample_ctl2 -= sample_ctl2 % 12;
  while (sample_ctl2 >= 120) {
  single_marker_3freqs_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    count_3freq_1920b((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)include2, &tot_a, &tot_b, &tot_c);
    lptr = lptr_12x_end;
    include2 = &(include2[cur_decr]);
    sample_ctl2 -= cur_decr;
  }
  if (sample_ctl2) {
    cur_decr = sample_ctl2;
    goto single_marker_3freqs_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[sample_ctl2 - (sample_ctl2 % 12)]);
  while (lptr < lptr_twelve_end) {
    count_3freq_48b(lptr, include2, &tot_a, &tot_b, &tot_c);
    lptr = &(lptr[12]);
    include2 = &(include2[12]);
  }
#endif
  while (lptr < lptr_end) {
    //   A := genotype & 0x5555...
    //   B := (genotype >> 1) & 0x5555...
    //   C := A & B
    //   popcount(C) = homozyg major ct
    //   popcount(B) = het ct + homozyg major ct
    //   popcount(A) = missing_ct + homozyg major ct
    loader = *lptr++;
    loader2 = *include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_a += popcount2_long(loader2);
    tot_b += popcount2_long(loader3);
    tot_c += popcount2_long(loader2 & loader3);
  }
  *hom2p = tot_c;
  *hetp = tot_b - tot_c;
  *missingp = tot_a - tot_c;
}


uint32_t PlinkLdBedFileChunk::init(int num_subjects, int snp_start_index, int num_snps_in_chunk, FILE* bedfile) {
  num_subj_ = num_subjects;
  num_snps_in_chunk_ = num_snps_in_chunk;

  const SampleCountInfo sc(num_subjects);
  std::vector<uintptr_t> loadbuf_vec(sc.unfiltered_sample_ctv2, 0);
  std::vector<uintptr_t> founder_info_vec(sc.unfiltered_sample_ctl, 0);
  fill_all_bits(sc.unfiltered_sample_ct, &founder_info_vec[0]);

  std::vector<uintptr_t> quatervec(sc.unfiltered_sample_ctv2, 0);
  init_quaterarr_from_bitarr(&(founder_info_vec[0]), sc.unfiltered_sample_ct, &quatervec[0]);

  const bool is_x = false;  // no special processing for X chromosome
  uintptr_t* founder_male_include2 = nullptr;  // not used when is_x == false
  const bool is_marker_reverse = false;
  const int bed_offset = 3;
 
  geno_vec.resize(num_snps_in_chunk * sc.founder_ct_192_long, 0);
  geno_masks_vec.resize(num_snps_in_chunk * sc.founder_ct_192_long, 0);
  ld_missing_cts_vec.resize(num_snps_in_chunk, 0);
  freq_.resize(num_snps_in_chunk_, 0);

  if (fseeko(bedfile, bed_offset + (snp_start_index * ((uint64_t)sc.unfiltered_sample_ct4)), SEEK_SET)) {
    return RET_READ_FAIL;
  }

  for (int snp_index = 0; snp_index < num_snps_in_chunk; snp_index++) {
    uintptr_t* mainbuf = &(geno()[snp_index * sc.founder_ct_192_long]);
    uint32_t error_code = load_and_collapse_incl(
      sc.unfiltered_sample_ct, sc.founder_ct, &founder_info_vec[0], sc.final_mask, is_marker_reverse,
      bedfile, &loadbuf_vec[0], mainbuf);
    if (error_code != 0)
      return error_code;

    // single_marker_3freqs must happen before ld_process_load2, because the later will re-code data in the mainbuf.
    uint32_t hom2 = 0, het = 0, missing = 0;
    single_marker_3freqs(sc.unfiltered_sample_ctv2, mainbuf, &quatervec[0], &hom2, &het, &missing);
    uint32_t nonmissing = num_subjects-missing;
    freq_[snp_index] = (nonmissing > 0) ? (float)(het + 2*hom2) / (float)(2*num_subjects-2*missing) : 0.5f;

    ld_process_load2(&(geno_vec[snp_index * sc.founder_ct_192_long]), 
                     &(geno_masks_vec[snp_index * sc.founder_ct_192_long]),
                     &(ld_missing_cts_vec[snp_index]), 
                     sc.founder_ct, is_x, founder_male_include2);
  }

  return 0;
}

double PlinkLdBedFileChunk::calculate_ld_corr(PlinkLdBedFileChunk& fixed_chunk, PlinkLdBedFileChunk& var_chunk, int snp_fixed_index, int snp_var_index) {
  // The following routine is combined from plink's ld_block_thread() and ld_report_regular() in plink_ld.c
  const SampleCountInfo sc(fixed_chunk.num_subj());

  const bool is_r2 = false;
  const bool keep_sign = false;

  uintptr_t* mask_fixed_vec_ptr = &(fixed_chunk.geno_masks()[snp_fixed_index * sc.founder_ct_192_long]);
  uintptr_t* mask_var_vec_ptr = &(var_chunk.geno_masks()[snp_var_index * sc.founder_ct_192_long]);
  uintptr_t* geno_fixed_vec_ptr = &(fixed_chunk.geno()[snp_fixed_index * sc.founder_ct_192_long]);
  uintptr_t* geno_var_vec_ptr = &(var_chunk.geno()[snp_var_index * sc.founder_ct_192_long]);

  uint32_t fixed_missing_ct = fixed_chunk.ld_missing_cts()[snp_fixed_index];
  uint32_t var_missing_ct = var_chunk.ld_missing_cts()[snp_var_index];
  uint32_t fixed_non_missing_ct = sc.founder_ct - fixed_missing_ct;
  uint32_t var_non_missing_ct = sc.founder_ct - var_missing_ct;
  uint32_t non_missing_ct = fixed_non_missing_ct - var_missing_ct;
  if (fixed_missing_ct && var_missing_ct) {
    non_missing_ct += ld_missing_ct_intersect(mask_var_vec_ptr, mask_fixed_vec_ptr, sc.founder_ctwd12, sc.founder_ctwd12_rem, sc.lshift_last);
  }

  int32_t dp_result[5];
  dp_result[0] = sc.founder_ct;
  dp_result[1] = -fixed_non_missing_ct;
  dp_result[2] = -var_non_missing_ct;
  dp_result[3] = dp_result[1];
  dp_result[4] = dp_result[2];
  ld_dot_prod(geno_var_vec_ptr, geno_fixed_vec_ptr, mask_var_vec_ptr, mask_fixed_vec_ptr, dp_result, sc.founder_ct_mld_m1, sc.founder_ct_mld_rem);

  double non_missing_ctd = (double)((int32_t)non_missing_ct);
  double dxx = dp_result[1];
  double dyy = dp_result[2];
  double cov12 = dp_result[0] * non_missing_ctd - dxx * dyy;
  dxx = (dp_result[3] * non_missing_ctd + dxx * dxx) * (dp_result[4] * non_missing_ctd + dyy * dyy);
  if (!is_r2) {
    dxx = cov12 / sqrt(dxx);
  } else if (!keep_sign) {
    dxx = (cov12 * cov12) / dxx;
  } else {
    dxx = (fabs(cov12) * cov12) / dxx;
  }

  return dxx;
}
