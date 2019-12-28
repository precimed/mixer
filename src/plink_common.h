#ifndef __PLINK_COMMON_H__
#define __PLINK_COMMON_H__

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


// Resources needed across all plink modules.

// The original source code from PLINK was modified by Oleksandr Frei, April 2019

#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#ifndef __STDC_FORMAT_MACROS
#  define __STDC_FORMAT_MACROS 1
#endif
#include <inttypes.h>

// avoid compiler warning
#ifndef NDEBUG
  #define NDEBUG
#endif
#include <assert.h>

#ifdef _WIN32
  // needed for MEMORYSTATUSEX
  #ifndef _WIN64
    #define WINVER 0x0500
  #else
    #define __LP64__
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
#else // Unix
  #include <sys/stat.h>
#endif

#ifndef HAVE_NULLPTR
  #ifndef __cplusplus
    #define nullptr NULL
  #else
    #if __cplusplus <= 199711L
      #ifndef nullptr
        #define nullptr NULL
      #endif
    #endif
  #endif
#endif

#ifdef _WIN32
  #define fseeko fseeko64
  #define ftello ftello64
  #include <process.h>
#  undef PRId64
#  undef PRIu64
#  define PRId64 "I64d"
#  define PRIu64 "I64u"
  #define pthread_t HANDLE
  #define THREAD_RET_TYPE unsigned __stdcall
  #define THREAD_RETURN return 0
  #define EOLN_STR "\r\n"
  #define FOPEN_RB "rb"
  #define FOPEN_WB "wb"
  #ifdef _WIN64
    #define getc_unlocked _fgetc_nolock
    #define putc_unlocked _fputc_nolock
  #else
    #define getc_unlocked getc
    #define putc_unlocked putc
  #endif
  #if __cplusplus < 201103L
    #define uint64_t unsigned long long
    #define int64_t long long
  #endif
#else
  #include <pthread.h>
  #define THREAD_RET_TYPE void*
  #define THREAD_RETURN return nullptr
  #ifdef __cplusplus
    #ifndef PRId64
      #define PRId64 "lld"
    #endif
  #endif
  #define EOLN_STR "\n"
  #define FOPEN_RB "r"
  #define FOPEN_WB "w"
  #ifndef __APPLE__
    // argh
    // not sure what the right threshold actually is, but this works for now
    // (may break on gcc <3.0?  but that shouldn't matter anymore)
    // tried defining GCC_VERSION, but that didn't always work
    #if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8)
      #define uint64_t unsigned long long
      #define int64_t long long
    #endif
  #endif
#endif

#ifdef _WIN64
  #define __LP64__
  #define CTZLU __builtin_ctzll
  #define CLZLU __builtin_clzll
#else
  #define CTZLU __builtin_ctzl
  #define CLZLU __builtin_clzl
  #ifndef __LP64__
    // attempt to patch GCC 6 build failure
    #if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8)
      #ifndef uintptr_t
        #define uintptr_t unsigned long
      #endif
      #ifndef intptr_t
        #define intptr_t long
      #endif
    #endif
  #endif
#endif

#ifdef __cplusplus
  #include <algorithm>
#  ifdef _WIN32
// Windows C++11 <algorithm> resets these values :(
#    undef PRIu64
#    undef PRId64
#    define PRIu64 "I64u"
#    define PRId64 "I64d"
#    undef PRIuPTR
#    undef PRIdPTR
#    ifdef __LP64__
#      define PRIuPTR PRIu64
#      define PRIdPTR PRId64
#    else
#      if __cplusplus < 201103L
#        define PRIuPTR "lu"
#        define PRIdPTR "ld"
#      else
#        define PRIuPTR "u"
#        define PRIdPTR "d"
#      endif
#    endif
#  endif
  #define HEADER_INLINE inline
#else
  #define HEADER_INLINE static inline
#endif

// It would be useful to disable compilation on big-endian platforms, but I
// don't see a decent portable way to do this (see e.g. discussion at
// http://esr.ibiblio.org/?p=5095 ).

#ifdef __LP64__
  #ifndef __SSE2__
    // It's obviously possible to support this by writing 64-bit non-SSE2 code
    // shadowing each SSE2 intrinsic, but this almost certainly isn't worth the
    // development/testing effort until regular PLINK 2.0 development is
    // complete.  No researcher has ever asked me for this feature.
    #error "64-bit builds currently require SSE2.  Try producing a 32-bit build instead."
  #endif
  #include <emmintrin.h>

  #define VECFTYPE __m128
  #define VECITYPE __m128i
  #define VECDTYPE __m128d

  // useful because of its bitwise complement: ~ZEROLU is a word with all 1
  // bits, while ~0 is always 32 1 bits.
  #define ZEROLU 0LLU

  // mainly useful for bitshifts: (ONELU << 32) works in 64-bit builds, while
  // (1 << 32) is undefined.  also used to cast some numbers/expressions to
  // uintptr_t (e.g. multiplying an int constant by ONELU widens it to 64 bits
  // only in 64-bit builds; note that 1LU fails on Win64 while 1LLU doesn't do
  // the right thing for 32-bit builds).
  #define ONELU 1LLU

  #ifdef _WIN32 // i.e. Win64

    #ifndef PRIuPTR
      #define PRIuPTR PRIu64
    #endif
    #ifndef PRIdPTR
      #define PRIdPTR PRId64
    #endif
    #define PRIxPTR2 "016I64x"

  #else // not _WIN32

    #ifndef PRIuPTR
      #define PRIuPTR "lu"
    #endif
    #ifndef PRIdPTR
      #define PRIdPTR "ld"
    #endif
    #define PRIxPTR2 "016lx"

  #endif // Win64

  #define VEC_BYTES 16

#else // not __LP64__

  #define ZEROLU 0LU
  #define ONELU 1LU
#  if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8) && (__cplusplus < 201103L)
#    undef PRIuPTR
#    undef PRIdPTR
#    define PRIuPTR "lu"
#    define PRIdPTR "ld"
#  endif
  #define PRIxPTR2 "08lx"

  // todo: update code so this still works when reduced to 4
  #define VEC_BYTES 8

#endif // __LP64__

// use constexpr for these as soon as compiler support is available on all
// platforms
#define FIVEMASK ((~ZEROLU) / 3)
#define AAAAMASK (FIVEMASK * 2)

#define VEC_BYTES_M1 (VEC_BYTES - 1)
#define VEC_BITS (VEC_BYTES * 8)
#define VEC_BITS_M1 (VEC_BITS - 1)

// 64MB of non-workspace memory guaranteed for now.
// Currently also serves as the maximum allele length.
#define NON_BIGSTACK_MIN 67108864

#define PI 3.1415926535897932
#define RECIP_2_32 0.00000000023283064365386962890625
#define RECIP_2_53 0.00000000000000011102230246251565404236316680908203125
// floating point comparison-to-nonzero tolerance, currently 2^{-30}
#define EPSILON 0.000000000931322574615478515625
// less tolerant versions (2^{-35}, 2^{-44}) for some exact calculations
#define SMALLISH_EPSILON 0.00000000002910383045673370361328125
#define SMALL_EPSILON 0.00000000000005684341886080801486968994140625
// at least sqrt(SMALL_EPSILON)
#define BIG_EPSILON 0.000000476837158203125
// 53-bit double precision limit
#define DOUBLE_PREC_LIMIT 0.00000000000000011102230246251565404236316680908203125
#define TWO_63 9223372036854775808.0
#define SQRT_HALF 0.70710678118654746

// 2^{-83} bias to give exact tests maximum ability to determine tiny p-values.
// (~2^{-53} is necessary to take advantage of denormalized small numbers, then
// allow tail sum to be up to 2^30.)
#define EXACT_TEST_BIAS 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

// occasionally used as an infinity substitute that avoids the 32-bit Windows
// performance penalty
// can import from limits.h, we don't bother to include that for now
#ifndef DBL_MAX
  #define DBL_MAX 1.7976931348623157e308
#endif

// not quite the same as FLT_MAX since it's a double-precision constant
#define FLT_MAXD 3.4028234663852886e38

#define RET_SUCCESS 0
#define RET_NOMEM 1
#define RET_OPEN_FAIL 2
#define RET_INVALID_FORMAT 3
#define RET_CALC_NOT_YET_SUPPORTED 4
#define RET_INVALID_CMDLINE 5
#define RET_WRITE_FAIL 6
#define RET_READ_FAIL 7
#define RET_THREAD_CREATE_FAIL 8
#define RET_ALLELE_MISMATCH 9
#define RET_NULL_CALC 10
#define RET_ALL_SAMPLES_EXCLUDED 11
#define RET_ALL_MARKERS_EXCLUDED 12
#define RET_NETWORK 13
#define LOAD_PHENO_LAST_COL 127

#define BIGSTACK_MIN_MB 64
#define BIGSTACK_DEFAULT_MB 2048

#ifdef __LP64__
  #define BITCT 64

  // unions generally shouldn't be used for reinterpret_cast's job (memcpy is
  // the right C-compatible way), but vectors are an exception to this rule.
  typedef union {
    VECFTYPE vf;
    VECITYPE vi;
    VECDTYPE vd;
    uintptr_t u8[VEC_BITS / BITCT];
    double d8[VEC_BYTES / sizeof(double)];
    float f4[VEC_BYTES / sizeof(float)];
    uint32_t u4[VEC_BYTES / sizeof(int32_t)];
  } __univec;
#else
  #define BITCT 32
#endif

#define BITCT2 (BITCT / 2)
#define BYTECT (BITCT / 8)
#define BYTECT4 (BITCT / 32)
#define VEC_WORDS (VEC_BITS / BITCT)
#define VEC_INT32 (VEC_BYTES / 4)

// assumed number of bytes per cache line, for alignment
#define CACHELINE 64

#define CACHELINE_BIT (CACHELINE * 8)
#define CACHELINE_INT32 (CACHELINE / 4)
#define CACHELINE_INT64 (CACHELINE / 8)
#define CACHELINE_WORD (CACHELINE / BYTECT)
#define CACHELINE_DBL (CACHELINE / 8)

// alignment must be a power of 2
HEADER_INLINE uintptr_t round_up_pow2(uintptr_t val, uintptr_t alignment) {
  uintptr_t alignment_m1 = alignment - 1;
  assert(!(alignment & alignment_m1));
  return (val + alignment_m1) & (~alignment_m1);
}

#define BITCT_TO_VECCT(val) (((val) + (VEC_BITS - 1)) / VEC_BITS)
#define BITCT_TO_WORDCT(val) (((val) + (BITCT - 1)) / BITCT)
#define BITCT_TO_ALIGNED_WORDCT(val) (VEC_WORDS * BITCT_TO_VECCT(val))

#define QUATERCT_TO_VECCT(val) (((val) + ((VEC_BITS / 2) - 1)) / (VEC_BITS / 2))
#define QUATERCT_TO_WORDCT(val) (((val) + (BITCT2 - 1)) / BITCT2)
#define QUATERCT_TO_ALIGNED_WORDCT(val) (VEC_WORDS * QUATERCT_TO_VECCT(val))

// todo: get rid of (BITCT_TO_WORDCT(x) == QUATERCT_TO_VECCT(x)) and similar
// assumptions, in preparation for AVX2

#ifdef __LP64__
#define round_up_pow2_ull round_up_pow2
#else
HEADER_INLINE uint64_t round_up_pow2_ull(uint64_t val, uint64_t alignment) {
  uint64_t alignment_m1 = alignment - 1;
  assert(!(alignment & alignment_m1));
  return (val + alignment_m1) & (~alignment_m1);
}
#endif

// 32-bit instead of word-length bitwise not here, when val can be assumed to
// be 32-bit.
// (note that the sizeof operator "returns" an uintptr_t, not a uint32_t; hence
// the lack of sizeof in the CACHELINE_INT32, etc. definitions.)
HEADER_INLINE uint32_t round_up_pow2_ui(uint32_t val, uint32_t alignment) {
  uint32_t alignment_m1 = alignment - 1;
  assert(!(alignment & alignment_m1));
  return (val + alignment_m1) & (~alignment_m1);
}

#define MAXV(aa, bb) (((bb) > (aa))? (bb) : (aa))
#define MINV(aa, bb) (((aa) > (bb))? (bb) : (aa))

#ifdef _WIN32
// if MAX_THREADS > 65, single WaitForMultipleObjects calls must be converted
// into loops
  #define MAX_THREADS 64
  #define MAX_THREADS_P1 65
#else
// shouldn't be larger than MODEL_BLOCKSIZE for now
  #define MAX_THREADS 512
  #define MAX_THREADS_P1 513
#endif

// defined as a macro since type of idx can vary; might want a debug
// compilation mode which performs type-checking, though
#define EXTRACT_2BIT_GENO(ulptr, idx) (((ulptr)[(idx) / BITCT2] >> (2 * ((idx) % BITCT2))) & 3)

// generic maximum line length.  .ped/.vcf/etc. lines can of course be longer
#define MAXLINELEN 131072

// must be at least 2 * MAXLINELEN + 2 to support generic token loader.
#define TEXTBUF_SIZE (2 * MAXLINELEN + 256)

#ifdef __LP64__
  // number of snp-major .bed lines to read at once for distance calc if
  // exponent is nonzero.
  #define MULTIPLEX_DIST_EXP 64
  // number of snp-major .bed lines to read at once for relationship calc
  #define MULTIPLEX_REL 60
#else
  // N.B. 32-bit version not as carefully tested or optimized, but I'll try to
  // make sure it works properly
  #define MULTIPLEX_DIST_EXP 28
  #define MULTIPLEX_REL 30
#endif

// load markers in blocks to enable multithreading and, for quantitative
// phenotypes, PERMORY-style LD exploitation
#define MODEL_BLOCKSIZE 1024
#define MODEL_BLOCKKEEP 64

// string hash table constants, currently only relevant for merge operations
// and annotate()
// (dynamic sizing used for main marker name lookup)

// last prime before 2^19
// size chosen to be likely to fit in L3 cache
#define HASHSIZE 524287
#define HASHSIZE_S 524287

#ifdef __LP64__
#define HASHMEM 4194304
#else
#define HASHMEM 2097152
#endif

HEADER_INLINE uint32_t popcount2_long(uintptr_t val) {
#ifdef __LP64__
  val = (val & 0x3333333333333333LLU) + ((val >> 2) & 0x3333333333333333LLU);
  return (((val + (val >> 4)) & 0x0f0f0f0f0f0f0f0fLLU) * 0x0101010101010101LLU) >> 56;
#else
  val = (val & 0x33333333) + ((val >> 2) & 0x33333333);
  return (((val + (val >> 4)) & 0x0f0f0f0f) * 0x01010101) >> 24;
#endif
}

HEADER_INLINE uint32_t popcount_long(uintptr_t val) {
  // the simple version, good enough for all non-time-critical stuff
  return popcount2_long(val - ((val >> 1) & FIVEMASK));
}

uintptr_t popcount_longs(const uintptr_t* lptr, uintptr_t word_ct);

HEADER_INLINE void fill_ulong_one(size_t size, uintptr_t* ularr) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *ularr++ = ~ZEROLU;
  }
}

#ifdef __LP64__
void count_2freq_dbl_960b(const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end, const VECITYPE* __restrict mask1vp, const VECITYPE* __restrict mask2vp, uint32_t* __restrict ct1abp, uint32_t* __restrict ct1cp, uint32_t* __restrict ct2abp, uint32_t* __restrict ct2cp);

void count_3freq_1920b(const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end, const VECITYPE* __restrict maskvp, uint32_t* __restrict ctap, uint32_t* __restrict ctbp, uint32_t* __restrict ctcp);
#else
void count_2freq_dbl_24b(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict mask1p, const uintptr_t* __restrict mask2p, uint32_t* __restrict ct1abp, uint32_t* __restrict ct1cp, uint32_t* __restrict ct2abp, uint32_t* __restrict ct2cp);

void count_3freq_48b(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict maskp, uint32_t* __restrict ctap, uint32_t* __restrict ctbp, uint32_t* __restrict ctcp);
#endif

void genovec_3freq(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict include_quatervec, uintptr_t sample_ctl2, uint32_t* __restrict missing_ctp, uint32_t* __restrict het_ctp, uint32_t* __restrict homset_ctp);

void fill_all_bits(uintptr_t ct, uintptr_t* bitarr);

HEADER_INLINE uintptr_t get_final_mask(uint32_t sample_ct) {
  uint32_t uii = sample_ct % BITCT2;
  if (uii) {
    return (ONELU << (2 * uii)) - ONELU;
  } else {
    return ~ZEROLU;
  }
}

// "quaterarr" refers to a packed group of base-4 (2-bit) elements, analogous
// to "bitarr".  (Based on "quaternary", not "quarter".)  "quatervec"
// indicates that vector-alignment is also required.
void fill_quatervec_55(uint32_t ct, uintptr_t* quatervec);
void init_quaterarr_from_bitarr(const uintptr_t* __restrict bitarr, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict new_quaterarr);
void init_quaterarr_from_inverted_bitarr(const uintptr_t* __restrict inverted_bitarr, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict new_quaterarr);

HEADER_INLINE uint32_t load_raw(uintptr_t unfiltered_sample_ct4, FILE* bedfile, uintptr_t* rawbuf) {
  // only use this if all accesses to the data involve
  // 1. some sort of mask, or
  // 2. explicit iteration from 0..(unfiltered_sample_ct-1).
  // otherwise improper trailing bits might cause a segfault, when we should
  // be ignoring them or just issuing a warning.
  return (fread(rawbuf, 1, unfiltered_sample_ct4, bedfile) < unfiltered_sample_ct4);
}

void reverse_loadbuf(uintptr_t unfiltered_sample_ct, unsigned char* loadbuf);
void copy_quaterarr_nonempty_subset(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict output_quaterarr);
uint32_t load_and_collapse_incl(uint32_t unfiltered_sample_ct, uint32_t sample_ct, const uintptr_t* __restrict sample_include, uintptr_t final_mask, uint32_t do_reverse, FILE* bedfile, uintptr_t* __restrict rawbuf, uintptr_t* __restrict mainbuf);

#endif  // __PLINK_COMMON_H__