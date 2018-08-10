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
  // "System" functions - error diagnostics, logging debug info to a file, memory management, status.
  DLL_PUBLIC const char* bgmg_get_last_error();
  DLL_PUBLIC void bgmg_init_log(const char* file);
  DLL_PUBLIC void bgmg_log_message(const char* message);
  DLL_PUBLIC int64_t bgmg_dispose(int context_id);
  DLL_PUBLIC const char* bgmg_status(int context_id);

  // API to work with "defvec". Here 
  // - num_snp is how many SNPs there is in the reference (particularly, in LD files and hvec)
  // - num_tag is how many SNPs there is in the GWAS (zvec, nvec, weights)
  // Tag snps must be a subset of the reference. tag_indices give indices of tag SNPs in the reference.
  DLL_PUBLIC int64_t bgmg_set_tag_indices(int context_id, int num_snp, int num_tag, int* tag_indices);
  DLL_PUBLIC int64_t bgmg_get_num_tag(int context_id);
  DLL_PUBLIC int64_t bgmg_get_num_snp(int context_id);
  DLL_PUBLIC int64_t bgmg_retrieve_tag_indices(int context_id, int num_tag, int* tag_indices);

  DLL_PUBLIC int64_t bgmg_set_chrnumvec(int context_id, int length, int* values);
  DLL_PUBLIC int64_t bgmg_retrieve_chrnumvec(int context_id, int length, int* buffer);

  // Set variouns options:
  // diag, kmax, r2min, max_causals, num_components, seed, fast_cost, threads, cache_tag_r2sum; refer to BgmgCalculator::set_option for a full list.
  // NB. Most options reset LD structure, you'll have to bgmg_set_ld_r2_coo / bgmg_set_ld_r2_csr again.
  DLL_PUBLIC int64_t bgmg_set_option(int context_id, char* option, double value);

  // API to populate and retrieve zvec, nvec, hvec
  DLL_PUBLIC int64_t bgmg_set_zvec(int context_id, int trait, int length, float* values);
  DLL_PUBLIC int64_t bgmg_set_nvec(int context_id, int trait, int length, float* values);
  DLL_PUBLIC int64_t bgmg_set_hvec(int context_id, int length, float* values);

  DLL_PUBLIC int64_t bgmg_retrieve_zvec(int context_id, int trait, int length, float* buffer);
  DLL_PUBLIC int64_t bgmg_retrieve_nvec(int context_id, int trait, int length, float* buffer);
  DLL_PUBLIC int64_t bgmg_retrieve_hvec(int context_id, int length, float* buffer);

  // API to populate LD structure
  // "from_file" expect a binary file of the following structure:
  // First 8 bytes - an integer, numel, that says how many LD r2 is in the file.
  // size(int32)*numel bytes - vector of snp indices, snpA
  // size(int32)*numel bytes - vector of snp indices, snpB
  // size(float)*numel bytes - vector of LDr2 between snpA and snpB
  // NB. snp indices must be zero-based.
  // NB. we expect that this data originates from plink, where LD matrix is lower triangular, diagonal not included. So snpA must be always lower than snpB.
  DLL_PUBLIC int64_t bgmg_set_ld_r2_coo(int context_id, int64_t length, int* snp_index, int* tag_index, float* r2);
  DLL_PUBLIC int64_t bgmg_set_ld_r2_coo_from_file(int context_id, const char* filename);
  DLL_PUBLIC int64_t bgmg_set_ld_r2_csr(int context_id);

  // Set weights, either explicitly or based on random pruning
  DLL_PUBLIC int64_t bgmg_set_weights(int context_id, int length, float* values);
  DLL_PUBLIC int64_t bgmg_set_weights_randprune(int context_id, int n, float r2);
  DLL_PUBLIC int64_t bgmg_retrieve_weights(int context_id, int length, float* buffer);

  // Retrieve certain aspects of the DL structure. Mainly for debugging purpose.
  DLL_PUBLIC int64_t bgmg_retrieve_tag_r2_sum(int context_id, int component_id, float num_causal, int length, float* buffer);
  DLL_PUBLIC int64_t bgmg_retrieve_ld_tag_r2_sum(int context_id, int length, float* buffer);  // LD scores (r2 and r4)
  DLL_PUBLIC int64_t bgmg_retrieve_ld_tag_r4_sum(int context_id, int length, float* buffer);
  DLL_PUBLIC int64_t bgmg_retrieve_weighted_causal_r2(int context_id, int length, float* buffer);

  // Calc univariate cost function and pdf
  DLL_PUBLIC double bgmg_calc_univariate_cost(int context_id, int trait_index, double pi_vec, double sig2_zero, double sig2_beta);
  DLL_PUBLIC double bgmg_calc_univariate_cost_with_deriv(int context_id, int trait_index, double pi_vec, double sig2_zero, double sig2_beta, int deriv_length, double* deriv);
  DLL_PUBLIC double bgmg_calc_univariate_pdf(int context_id, int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* zvec, float* pdf);

  // Calc bivariate cost function and pdf
  DLL_PUBLIC double bgmg_calc_bivariate_cost(int context_id, int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);
  DLL_PUBLIC int64_t bgmg_calc_bivariate_pdf(int context_id, int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, int length, float* zvec1, float* zvec2, float* pdf);
  
  // functions to work with loglikelihood cache
  DLL_PUBLIC int64_t bgmg_clear_loglike_cache(int context_id);
  DLL_PUBLIC int64_t bgmg_get_loglike_cache_size(int context_id);
  DLL_PUBLIC int64_t bgmg_get_loglike_cache_univariate_entry(int context_id, int entry_index, float* pi_vec, float* sig2_zero, float* sig2_beta, double* cost);
  DLL_PUBLIC int64_t bgmg_get_loglike_cache_bivariate_entry(int context_id, int entry_index, int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float* rho_beta, int sig2_zero_len, float* sig2_zero, float* rho_zero, double* cost);
}

