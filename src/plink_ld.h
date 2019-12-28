#ifndef __PLINK_LD_H__
#define __PLINK_LD_H__

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

#endif  // __PLINK_LD_H__

#include "plink_common.h"

#include <vector>

// Calculates several derived measures from the number of subjects.
// Has several simplifications compared to what is typically handled in plink:
// * Assumes that all individuals are founders.
// * Performs no filtering on individuals.
class SampleCountInfo {
public:
  uintptr_t unfiltered_sample_ct;
  uintptr_t unfiltered_sample_ctl;
  uintptr_t unfiltered_sample_ctl2;
  uintptr_t unfiltered_sample_ctv2;
  uint32_t unfiltered_sample_ct4;

  uintptr_t founder_ct;
  uintptr_t final_mask ;

  uintptr_t founder_ct_mld;
  uint32_t founder_ct_mld_m1;
  uintptr_t founder_ct_mld_rem;
  uintptr_t founder_ct_192_long;

  uintptr_t founder_ctwd;
  uintptr_t founder_ctwd12;
  uintptr_t founder_ctwd12_rem;
  uintptr_t lshift_last;
public:
  explicit SampleCountInfo(int num_subjects);
};

// A class that wraps a chunk of a plink BED file, and stores it into a format suitable for computing LD allelic correlation.
class PlinkLdBedFileChunk {
 public:
 PlinkLdBedFileChunk() {}
  explicit PlinkLdBedFileChunk(int num_subjects, int snp_start_index, int num_snps_in_chunk, FILE* bedfile) { init(num_subjects, snp_start_index, num_snps_in_chunk, bedfile); }
  uint32_t init(int num_subjects, int snp_start_index, int num_snps_in_chunk, FILE* bedfile);

  uintptr_t* geno() {return &geno_vec[0];}
  uintptr_t* geno_masks() {return &geno_masks_vec[0];}
  uint32_t* ld_missing_cts() {return &ld_missing_cts_vec[0];}
  float* freq() {return &freq_[0];}
  float hetval(int snp_index) { return 2.0f * freq_[snp_index] * (1.0f - freq_[snp_index]);} // heterozygosity
  int num_subj() { return num_subj_; }

  static double calculate_ld_corr(PlinkLdBedFileChunk& fixed_chunk, PlinkLdBedFileChunk& var_chunk, int snp_fixed_index, int snp_var_index);

 private:
  int num_subj_;
  int num_snps_in_chunk_;
  std::vector<uintptr_t> geno_vec;        // geno_vec and geno_masks_vec has special encoding for LD structure, see ld_process_load2()
  std::vector<uintptr_t> geno_masks_vec;
  std::vector<uint32_t> ld_missing_cts_vec;
  std::vector<float> freq_;
};
