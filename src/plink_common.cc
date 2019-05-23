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

#include "plink_common.h"

#ifdef __LP64__
// Basic SSE2 implementation of Lauradoux/Walisch popcount.
static inline uintptr_t popcount_vecs(const __m128i* vptr, uintptr_t ct) {
  // popcounts vptr[0..(ct-1)].  Assumes ct is a multiple of 3 (0 ok).
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  uintptr_t tot = 0;
  const __m128i* vend;
  __m128i count1;
  __m128i count2;
  __m128i half1;
  __m128i half2;
  __univec acc;

  while (ct >= 30) {
    ct -= 30;
    vend = &(vptr[30]);
  popcount_vecs_main_loop:
    acc.vi = _mm_setzero_si128();
    do {
      count1 = *vptr++;
      count2 = *vptr++;
      half1 = *vptr++;
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1
      // count2 store a partial bitcount covering themselves AND another bit
      // from elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (ct) {
    vend = &(vptr[ct]);
    ct = 0;
    goto popcount_vecs_main_loop;
  }
  return tot;
}

#endif

uintptr_t popcount_longs(const uintptr_t* lptr, uintptr_t word_ct) {
  // Efficiently popcounts lptr[0..(word_ct - 1)].  In the 64-bit case, lptr[]
  // must be 16-byte aligned.
  // The popcount_longs_nzbase() wrapper takes care of starting from a later
  // index.
  uintptr_t tot = 0;
  const uintptr_t* lptr_end = &(lptr[word_ct]);
#ifdef __LP64__
  uintptr_t six_ct;
  const __m128i* vptr;
  vptr = (const __m128i*)lptr;
  six_ct = word_ct / 6;
  tot += popcount_vecs(vptr, six_ct * 3);
  lptr = &(lptr[six_ct * 6]);
#else
  // The humble 16-bit lookup table actually beats
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
  // on my development machine by a hair.
  // However, if we take the hint from Lauradoux/Walisch and postpone the
  // multiply and right shift, this is no longer true.  Ah well.
  const uintptr_t* lptr_six_end;
  uintptr_t tmp_stor;
  uintptr_t loader;
  uintptr_t ulii;
  uintptr_t uljj;
  lptr_six_end = &(lptr[word_ct - (word_ct % 6)]);
  while (lptr < lptr_six_end) {
    loader = *lptr++;
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = *lptr++;
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101 is
    // equivalent to the left-shifts and adds we need to sum those four 8-bit
    // numbers in the high-order slot.
    tot += (tmp_stor * 0x01010101) >> 24;
  }
#endif
  while (lptr < lptr_end) {
    tot += popcount_long(*lptr++);
  }
  return tot;
}

void fill_all_bits(uintptr_t ct, uintptr_t* bitarr) {
  // leaves bits beyond the end unset
  // ok for ct == 0
  uintptr_t quotient = ct / BITCT;
  uintptr_t remainder = ct % BITCT;
  fill_ulong_one(quotient, bitarr);
  if (remainder) {
    bitarr[quotient] = (ONELU << remainder) - ONELU;
  }
}

void copy_quaterarr_nonempty_subset(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict output_quaterarr) {
  // in plink 2.0, we probably want (0-based) bit raw_quaterarr_size of
  // subset_mask to be always allocated and unset.  This removes a few special
  // cases re: iterating past the end of arrays.
  assert(subset_size);
  assert(raw_quaterarr_size >= subset_size);
  uintptr_t cur_output_word = 0;
  uintptr_t* output_quaterarr_last = &(output_quaterarr[subset_size / BITCT2]);
  const uint32_t word_write_halfshift_end = subset_size % BITCT2;
  uint32_t word_write_halfshift = 0;
  // if < 2/3-filled, use sparse copy algorithm
  if (subset_size * (3 * ONELU) < raw_quaterarr_size * (2 * ONELU)) {
    uint32_t subset_mask_widx = 0;
    while (1) {
      const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
      if (cur_include_word) {
	uint32_t wordhalf_idx = 0;
#ifdef __LP64__
	uint32_t cur_include_halfword = (uint32_t)cur_include_word;
#else
	uint32_t cur_include_halfword = (uint16_t)cur_include_word;
#endif
	while (1) {
	  if (cur_include_halfword) {
	    uintptr_t raw_quaterarr_word = raw_quaterarr[subset_mask_widx * 2 + wordhalf_idx];
	    do {
	      uint32_t rqa_idx_lowbits = __builtin_ctz(cur_include_halfword);
	      cur_output_word |= ((raw_quaterarr_word >> (rqa_idx_lowbits * 2)) & 3) << (word_write_halfshift * 2);
	      if (++word_write_halfshift == BITCT2) {
		*output_quaterarr++ = cur_output_word;
		word_write_halfshift = 0;
		cur_output_word = 0;
	      }
	      cur_include_halfword &= cur_include_halfword - 1;
	    } while (cur_include_halfword);
	  }
	  if (wordhalf_idx) {
	    break;
	  }
	  wordhalf_idx++;
#ifdef __LP64__
	  cur_include_halfword = cur_include_word >> 32;
#else
	  cur_include_halfword = cur_include_word >> 16;
#endif
	}
	if (output_quaterarr == output_quaterarr_last) {
	  if (word_write_halfshift == word_write_halfshift_end) {
            if (word_write_halfshift_end) {
	      *output_quaterarr_last = cur_output_word;
	    }
	    return;
	  }
	}
      }
      subset_mask_widx++;
    }
  }
  // blocked copy
  while (1) {
    const uintptr_t cur_include_word = *subset_mask++;
    uint32_t wordhalf_idx = 0;
#ifdef __LP64__
    uintptr_t cur_include_halfword = (uint32_t)cur_include_word;
#else
    uint32_t cur_include_halfword = (uint16_t)cur_include_word;
#endif
    while (1) {
      uintptr_t raw_quaterarr_word = *raw_quaterarr++;
      while (cur_include_halfword) {
	uint32_t rqa_idx_lowbits = CTZLU(cur_include_halfword);
	uintptr_t halfword_invshifted = (~cur_include_halfword) >> rqa_idx_lowbits;
	uintptr_t raw_quaterarr_curblock_unmasked = raw_quaterarr_word >> (rqa_idx_lowbits * 2);
	uint32_t rqa_block_len = CTZLU(halfword_invshifted);
	uint32_t block_len_limit = BITCT2 - word_write_halfshift;
	cur_output_word |= raw_quaterarr_curblock_unmasked << (2 * word_write_halfshift);
	if (rqa_block_len < block_len_limit) {
	  word_write_halfshift += rqa_block_len;
	  cur_output_word &= (ONELU << (word_write_halfshift * 2)) - ONELU;
	} else {
	  // no need to mask, extra bits vanish off the high end
	  *output_quaterarr++ = cur_output_word;
	  word_write_halfshift = rqa_block_len - block_len_limit;
	  if (word_write_halfshift) {
	    cur_output_word = (raw_quaterarr_curblock_unmasked >> (2 * block_len_limit)) & ((ONELU << (2 * word_write_halfshift)) - ONELU);
	  } else {
	    // avoid potential right-shift-64
	    cur_output_word = 0;
	  }
	}
	cur_include_halfword &= (~(ONELU << (rqa_block_len + rqa_idx_lowbits))) + ONELU;
      }
      if (wordhalf_idx) {
	break;
      }
      wordhalf_idx++;
#ifdef __LP64__
      cur_include_halfword = cur_include_word >> 32;
#else
      cur_include_halfword = cur_include_word >> 16;
#endif
    }
    if (output_quaterarr == output_quaterarr_last) {
      if (word_write_halfshift == word_write_halfshift_end) {
	if (word_write_halfshift_end) {
	  *output_quaterarr_last = cur_output_word;
	}
	return;
      }
    }
  }
}

void reverse_loadbuf(uintptr_t unfiltered_sample_ct, unsigned char* loadbuf) {
  // unfiltered_sample_ct can be zero
  uintptr_t sample_bidx = 0;
  unsigned char* loadbuf_end = &(loadbuf[(unfiltered_sample_ct + 3) / 4]);
  unsigned char ucc;
  unsigned char ucc2;
  uintptr_t unfiltered_sample_ctd;
  uint32_t* loadbuf_alias32;
  uint32_t uii;
  uint32_t ujj;
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  __m128i* loadbuf_alias;
  __m128i vii;
  __m128i vjj;
  // todo: use this vector loop even when loadbuf is unaligned, so stuff like
  // recode_load_to() is faster
  if (!(((uintptr_t)loadbuf) & 15)) {
    loadbuf_alias = (__m128i*)loadbuf;
    unfiltered_sample_ctd = unfiltered_sample_ct / 64;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      vii = *loadbuf_alias;
      // we want to exchange 00 and 11, and leave 01/10 untouched.  So make
      // vjj := 11 iff vii is 00/11, and vjj := 00 otherwise; then xor.
      vjj = _mm_andnot_si128(_mm_xor_si128(vii, _mm_srli_epi64(vii, 1)), m1);
      vjj = _mm_or_si128(vjj, _mm_slli_epi64(vjj, 1));
      *loadbuf_alias++ = _mm_xor_si128(vii, vjj);
    }
    loadbuf = (unsigned char*)loadbuf_alias;
  } else if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      uii = *loadbuf_alias32;
      ujj = 0x55555555 & (~(uii ^ (uii >> 1)));
      ujj *= 3;
      *loadbuf_alias32++ = uii ^ ujj;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
#else
  if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      uii = *loadbuf_alias32;
      ujj = 0x55555555 & (~(uii ^ (uii >> 1)));
      ujj *= 3;
      *loadbuf_alias32++ = uii ^ ujj;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
#endif
  for (; loadbuf < loadbuf_end;) {
    ucc = *loadbuf;
    ucc2 = 0x55 & (~(ucc ^ (ucc >> 1)));
    ucc2 *= 3;
    *loadbuf++ = ucc ^ ucc2;
  }
  uii = unfiltered_sample_ct & 3;
  if (uii) {
    loadbuf[-1] &= (0xff >> (8 - 2 * uii));
  }
}

uint32_t load_and_collapse_incl(uint32_t unfiltered_sample_ct, uint32_t sample_ct, const uintptr_t* __restrict sample_include, uintptr_t final_mask, uint32_t do_reverse, FILE* bedfile, uintptr_t* __restrict rawbuf, uintptr_t* __restrict mainbuf) {
  assert(unfiltered_sample_ct);
  uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  if (unfiltered_sample_ct == sample_ct) {
    rawbuf = mainbuf;
  }
  if (load_raw(unfiltered_sample_ct4, bedfile, rawbuf)) {
    return RET_READ_FAIL;
  }
  if (unfiltered_sample_ct != sample_ct) {
    copy_quaterarr_nonempty_subset(rawbuf, sample_include, unfiltered_sample_ct, sample_ct, mainbuf);
  } else {
    mainbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
  }
  if (do_reverse) {
    reverse_loadbuf(sample_ct, (unsigned char*)mainbuf);
  }
  return 0;
}
