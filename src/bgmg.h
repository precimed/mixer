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
#include "bgmg_export.h"

#define VERSION "v0.9.0"

extern "C" {
  DLL_PUBLIC const char* bgmg_get_last_error();
  DLL_PUBLIC int64_t bgmg_set_zvec(int context_id, int trait, int length, double* values);
  DLL_PUBLIC int64_t bgmg_set_nvec(int context_id, int trait, int length, double* values);
  DLL_PUBLIC int64_t bgmg_set_hvec(int context_id, int length, double* values);
  DLL_PUBLIC int64_t bgmg_set_w_ld(int context_id, int length, double* values);
  DLL_PUBLIC int64_t bgmg_set_ref_ld_sum_r2(int context_id, int length, int r2bins, double* values);
  DLL_PUBLIC int64_t bgmg_set_ref_ld_sum_r4(int context_id, int length, int r2bins, double* values);
  DLL_PUBLIC int64_t bgmg_set_option(int context_id, char* option, int value);
  DLL_PUBLIC double bgmg_calc_univariate_cost(int context_id, double pi_vec, double sig2_zero, double sig2_beta);
  DLL_PUBLIC double bgmg_calc_bivariate_cost(int context_id, int num_components, int num_traits, double* pi_vec, double* sig2_beta, double* rho_beta, double* sig2_zero, double rho_zero);
  DLL_PUBLIC int64_t dispose(int context_id);
}

