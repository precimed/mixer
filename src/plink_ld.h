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

class PlinkLdBaseBedFile {
protected:
  uintptr_t unfiltered_sample_ct;
  uintptr_t unfiltered_sample_ctl;
  uintptr_t unfiltered_sample_ctl2;
  uintptr_t unfiltered_sample_ctv2;

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

  std::vector<uintptr_t> founder_info_vec;
  uintptr_t* founder_info;

  virtual ~PlinkLdBaseBedFile() {}

public:
  explicit PlinkLdBaseBedFile(int num_subjects);
};

// A class that wraps plink BED file.
// Assumes that all individuals are founders.
// Performs no filtering on individuals.
// Performs no filtering on variants.
class PlinkLdBedFile : public PlinkLdBaseBedFile {
 public:
  explicit PlinkLdBedFile(int num_subjects, int num_snps, FILE* bedfile);
  virtual ~PlinkLdBedFile() {}
  double ld_corr(int snp_fixed_index, int snp_var_index);

 private:
  FILE* bedfile;

  int num_snps_;
  std::vector<uintptr_t> geno_vec;
  uintptr_t* geno;

  std::vector<uintptr_t> loadbuf_vec;
  uintptr_t* loadbuf;

  std::vector<uintptr_t> geno_masks_vec;
  uintptr_t* geno_masks;

  std::vector<uint32_t> ld_missing_cts_vec;
  uint32_t* ld_missing_cts;
};

