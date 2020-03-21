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

#include "bgmg_calculator_impl.h"

/*
// This file contains the following functions, all of the are legacy and are used only in unit-tests.
int64_t BgmgCalculator::find_snp_order()
int64_t BgmgCalculator::find_tag_r2sum(int component_id, float num_causals)
int64_t BgmgCalculator::retrieve_tag_r2_sum(int component_id, float num_causal, int length, float* buffer)
int64_t BgmgCalculator::calc_univariate_pdf(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* zvec, float* pdf) {
int64_t BgmgCalculator::calc_univariate_power(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, float zthresh, int length, float* nvec, float* svec) {
int64_t BgmgCalculator::calc_univariate_delta_posterior(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* c0, float* c1, float* c2) {
double BgmgCalculator::calc_univariate_cost(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
double BgmgCalculator::calc_univariate_cost_cache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
template<typename T> double calc_univariate_cost_nocache_template(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, BgmgCalculator& rhs) {
double BgmgCalculator::calc_univariate_cost_nocache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
double BgmgCalculator::calc_univariate_cost_nocache_float(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
double BgmgCalculator::calc_univariate_cost_nocache_double(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
double BgmgCalculator::calc_bivariate_cost(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
double BgmgCalculator::calc_bivariate_cost_cache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
double BgmgCalculator::calc_bivariate_cost_nocache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
int64_t BgmgCalculator::calc_bivariate_pdf(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, int length, float* zvec1, float* zvec2, float* pdf) {
int64_t BgmgCalculator::calc_bivariate_delta_posterior(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, int length, float* c00, float* c10, float* c01, float* c20, float* c11, float* c02) {
double BgmgCalculator::calc_univariate_cost_fast(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
double BgmgCalculator::calc_bivariate_cost_fast(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
double BgmgCalculator::calc_univariate_cost_convolve(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
void BgmgCalculator::clear_tag_r2sum(int component_id)
double BgmgCalculator::calc_bivariate_cost_convolve(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
int64_t BgmgCalculator::set_snp_order(int component_id, int64_t length, const int* buffer) {
int64_t BgmgCalculator::retrieve_snp_order(int component_id, int64_t length, int* buffer) {
int64_t BgmgCalculator::retrieve_k_pdf(int length, double* buffer) {
void BgmgCalculator::find_tag_r2sum_no_cache(int component_id, float num_causal, int k_index, std::vector<float>* buffer) {
int BgmgCalculator::find_deftag_indices_znw(int trait_index, std::vector<int>* deftag_indices) {
int BgmgCalculator::find_deftag_indices_nw(int trait_index, std::vector<int>* deftag_indices) {
int BgmgCalculator::find_deftag_indices_w(std::vector<int>* deftag_indices) {
int BgmgCalculator::find_deftag_indices_znw(std::vector<int>* deftag_indices) {
int BgmgCalculator::find_deftag_indices_nw(std::vector<int>* deftag_indices) {
void BgmgCalculator::clear_state() {
*/

int64_t BgmgCalculator::find_snp_order() {
  if (max_causals_ <= 0 || max_causals_ > num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("find_snp_order: max_causals_ <= 0 || max_causals_ > num_snp_"));
  if (num_components_ <= 0 || num_components_ > 3) BGMG_THROW_EXCEPTION(::std::runtime_error("find_snp_order: num_components_ must be between 1 and 3"));
  if (snp_order_.size() > 0) BGMG_THROW_EXCEPTION(::std::runtime_error("find_snp_order: called twice"));

  LOG << ">find_snp_order(num_components_=" << num_components_ << ", k_max_=" << k_max_ << ", max_causals_=" << max_causals_ << ")";
  SimpleTimer timer(-1);

  // Right now all SNPs must be included in snp_can_be_causal_.
  // Remember that the only purpose of snp_can_be_causal_ is to limit the information that store about the LD matrix.
  // At some point we use LD structure only to calculate tag_r2, so we only store r2 for SNPs that are selected as causal by find_snp_order.
  // Later we've started to use LD structure to 
  //  - perform random pruning (=> LD must be stored for all tag variants)
  //    "for (int i = 0; i < num_tag_; i++) snp_can_be_causal_[tag_to_snp_[i]] = 1;"
  //  - calculate sum_r2, sum_r4 (=> LD must be stored for all variants, OR we need to change the logic and load hvec so that sum_r2 and sum_r4 are calculated on the fly in set_ld_r2_coo.
  // For now we simply set snp_can_be_causal_ to 1 and store LD structure for all variants.
  std::vector<char> snp_can_be_causal(num_snp_, 0);  // mask of SNPs that may be causal (e.i. included in snp_order array)

  SimpleTimer log_timer(10000); // log some message each 10 seconds
  for (int component_index = 0; component_index < num_components_; component_index++) {
    if (log_timer.fire())
      LOG << " find_snp_order still working, component_id=" << component_index;

    snp_order_.push_back(std::make_shared<DenseMatrix<int>>(max_causals_, k_max_));

    if (cache_tag_r2sum_) clear_tag_r2sum(component_index);
    
#pragma omp parallel
    {
      std::vector<int> perm(num_snp_, 0);

#pragma omp for schedule(static)
      for (int k = 0; k < k_max_; k++) {
        for (int i = 0; i < num_snp_; i++) perm[i] = i;

        std::mt19937_64 random_engine;
        random_engine.seed(seed_ + component_index * k_max_ + k);  // ensure each k in each component starts with its own seed.

        // perform partial Fisher Yates shuffle (must faster than full std::shuffle)
        // swap_offset is a random integer, with max of n-1, n-2, n-3, ..., n-max_causals
        for (int i = 0; i < max_causals_; i++) {
          const int swap_offset = std::uniform_int_distribution<int>(0, num_snp_ - i - 1)(random_engine);
          std::iter_swap(perm.begin() + i, perm.begin() + i + swap_offset);
        }

        for (int i = 0; i < max_causals_; i++) {
          (*snp_order_[component_index])(i, k) = perm[i];
        }
      }
    }

    // Fill in snp_can_be_causal
    for (int k = 0; k < k_max_; k++) {
      for (int i = 0; i < max_causals_; i++) {
        snp_can_be_causal[(*snp_order_[component_index])(i, k)] = 1;
      }
    }
  }

  int num_can_be_causal = 0;
  for (int i = 0; i < num_snp_; i++) num_can_be_causal += snp_can_be_causal[i];
  LOG << "<find_snp_order: num_can_be_causal = " << num_can_be_causal << ", elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::find_tag_r2sum(int component_id, float num_causals) {
  if (!cache_tag_r2sum_) BGMG_THROW_EXCEPTION(::std::runtime_error("find_tag_r2sum can be used only with cache_tag_r2sum==true"));
  if (num_causals < 0 || num_causals >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("find_tag_r2sum: num_causals < 0 || num_causals >= max_causals_"));
  if (component_id < 0 || component_id >= num_components_) BGMG_THROW_EXCEPTION(::std::runtime_error("find_tag_r2sum: component_id must be between 0 and num_components_"));

  const float num_causals_original = num_causals;
  if (snp_order_.empty()) find_snp_order();

  float last_num_causals = last_num_causals_[component_id]; 
  const float last_num_causals_original = last_num_causals;
  
  LOG << ">find_tag_r2sum(component_id=" << component_id << ", num_causals=" << num_causals << ", last_num_causals=" << last_num_causals << ")";
  SimpleTimer timer(-1);

  // if num_causal is more than twice lower than last_num_causals we should re-calculate tag_r2sum from scratch.
  if (num_causals < (last_num_causals / 2)) {
    clear_tag_r2sum(component_id);
    last_num_causals = 0.0f;
  }

  // changeset contains a list of indices with corresponding weight
  // indices apply to snp_order_[component_id] array.
  // weights are typicaly +1 (to increase by r2) or -1 (to decrease by r2).
  // First and last weights is float-point number between 1 and -1,
  // to handle cases when num_causals is float-point number (derived from pivec).
  // This is important for fminsearch which get's confused if cost is a stepwise of pivec.
  std::vector<std::pair<int, float>> changeset;
  
  // Decreasing number of causals from B to A has an opposite effect to increasing from A to B.
  // To handle decreasing case we just swap num_causals and last_num_causals, and set sign to -1.0f.
  float sign = 1.0f;
  if (num_causals < last_num_causals) {
    float tmp = num_causals; num_causals = last_num_causals; last_num_causals = tmp;
    sign = -1.0f;
  }

  // There are 3 cases
  // 1. floor(num_causals) == floor(last_num_causals)
  // 2. floor(num_causals) == floor(last_num_causals) + 1
  // 3. floor(num_causals) >= floor(last_num_causals) + 2

  float floor_num_causals = floor(num_causals);
  float floor_last_num_causals = floor(last_num_causals);
  if ((int)floor_num_causals == (int)floor_last_num_causals) {
    changeset.push_back(std::make_pair((int)floor_last_num_causals, sign * (num_causals - last_num_causals)));
  }
  else if ((int)floor_num_causals >= ((int)floor_last_num_causals + 1)) {
    // handle case 2 and case 3 - lower boundary
    changeset.push_back(std::make_pair((int)floor_last_num_causals, sign * (floor_last_num_causals + 1.0f - last_num_causals)));

    // happends for the case 3 - bulk change (empty loop in case 2)
    for (int i = ((int)floor_last_num_causals + 1); i < (int)floor_num_causals; i++) {
      changeset.push_back(std::make_pair(i, sign));
    }

    // handle case 2 and case 3 - upper boundary
    changeset.push_back(std::make_pair((int)floor_num_causals, sign * (num_causals - floor_num_causals)));
  }
  else {
    BGMG_THROW_EXCEPTION(::std::runtime_error("floor_num_causals < floor_last_num_causals"));
  }

  // apply infinitesimal model to adjust tag_r2sum for all r2 that are below r2min (and thus do not contribute via resampling)
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  const float pival_delta = (num_causals_original - last_num_causals_original) / static_cast<float>(num_snp_);

  std::vector<float> hvec; find_hvec(*this, &hvec);

  // it is OK to parallelize the following loop on k_index, because:
  // - all structures here are readonly, except tag_r2sum_ that we are accumulating
  // - two threads will never touch the same memory location (that's why we choose k_index as an outer loop)
#pragma omp parallel
{
  LdMatrixRow ld_matrix_row;

#pragma omp for schedule(static)
  for (int k_index = 0; k_index < k_max_; k_index++) {
    for (auto change : changeset) {
      int scan_index = change.first;
      float scan_weight = change.second;
      int snp_index = (*snp_order_[component_id])(scan_index, k_index);  // index of a causal snp
      ld_matrix_csr_.extract_snp_row(SnpIndex(snp_index), &ld_matrix_row);
      auto iter_end = ld_matrix_row.end();
      for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
        int tag_index = iter.index();
        float r2 = iter.r2();
        float hval = hvec[snp_index];
        (*tag_r2sum_[component_id])(tag_index, k_index) += (scan_weight * r2 * hval);
      }
    }
  }
}

{
  SimpleTimer timer2(-1);
  int num_tag_inf_adjusted = 0;  // number of snps adjusted according to infinitesimal model
#pragma omp parallel for schedule(static) reduction(+: num_tag_inf_adjusted)
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    float inf_adj = pival_delta * ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index];
    if (inf_adj != 0) {
      num_tag_inf_adjusted++;
      for (int k_index = 0; k_index < k_max_; k_index++) {
        (*tag_r2sum_[component_id])(tag_index, k_index) += inf_adj;
      }
    }
  }
  if (num_tag_inf_adjusted > 0)
    LOG << " apply infinitesimal model to " << num_tag_inf_adjusted << " tag SNPs, to adjust tag_r2sum for all r2 that are below r2min, elapsed time " << timer2.elapsed_ms() << "ms";
}

  LOG << "<find_tag_r2sum(component_id=" << component_id << ", num_causals=" << num_causals_original << ", last_num_causals=" << last_num_causals << "), elapsed time " << timer.elapsed_ms() << "ms";

  last_num_causals_[component_id] = num_causals_original;
  return 0;
}

int64_t BgmgCalculator::retrieve_tag_r2_sum(int component_id, float num_causal, int length, float* buffer) {
  if (length != (k_max_ * num_tag_)) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (num_causal < 0 && !cache_tag_r2sum_) BGMG_THROW_EXCEPTION(::std::runtime_error("retrieve_tag_r2sum with num_causal<0 is meant for cache_tag_r2sum==true"));
  if (component_id < 0 || component_id >= num_components_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong component_id"));

  LOG << " retrieve_tag_r2_sum(component_id=" << component_id << ", num_causal=" << num_causal << ")";

  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) {
    // use negative to retrieve tag_r2_sum for last_num_causal (for debugging purpose)
    if (num_causal >= 0) {
      find_tag_r2sum(component_id, num_causal);
    }

    for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
      for (int k_index = 0; k_index < k_max_; k_index++) {
        float tag_r2sum = (*tag_r2sum_[component_id])(tag_index, k_index);
        buffer[tag_index * k_max_ + k_index] = tag_r2sum;
      }
    }
  } else {
#pragma omp parallel
    {
      std::vector<float> tag_r2sum(num_tag_, 0.0f);

#pragma omp for schedule(static)
      for (int k_index = 0; k_index < k_max_; k_index++) {
        find_tag_r2sum_no_cache(0, num_causal, k_index, &tag_r2sum);
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          buffer[tag_index * k_max_ + k_index] = tag_r2sum[tag_index];
        }
      }
    }
  }
  return 0;
}

int64_t BgmgCalculator::calc_univariate_pdf(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* zvec, float* pdf) {
  // input buffer "zvec" contains z scores (presumably an equally spaced grid)
  // output buffer contains pdf(z), aggregated across all SNPs with corresponding weights
  
  std::vector<float>& nvec(*get_nvec(trait_index));
  if (nvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  float num_causals = pi_vec * static_cast<float>(num_snp_);
  if ((int)num_causals >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  const int component_id = 0;   // univariate is always component 0.

  LOG << ">calc_univariate_pdf(trait_index="<< trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ", length(zvec)=" << length << ")";
  SimpleTimer timer(-1);

  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) find_tag_r2sum(component_id, num_causals);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices_nw(trait_index, &deftag_indices);

  const double pi_k = 1.0 / static_cast<double>(k_max_);

// omp reduction on std::vector ( https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector ) - did not work for microsoft compiler
// #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
//                               std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
//                     initializer(omp_priv = omp_orig)
// Final solution is to do a the reduction with omp critical (see here http://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-reduction.html )

  // we accumulate crazy many small values - each of them is OK as float; the sum is also OK as float;  
  // but accumulation must be done with double precision.
  // std::vector<double> pdf_double(length, 0.0);
  std::valarray<double> pdf_double(0.0, length);

#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, length);
    std::vector<float> tag_r2sum(num_tag_, 0.0f);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {

      if (cache_tag_r2sum_) {
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          tag_r2sum[tag_index] = (*tag_r2sum_[0])(tag_index, k_index);
        }
      } else {
        find_tag_r2sum_no_cache(0, num_causals, k_index, &tag_r2sum);
      }

      for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
        int tag_index = deftag_indices[deftag_index];

        const double tag_weight = static_cast<double>(weights_[tag_index]);

        float tag_r2sum_value = tag_r2sum[tag_index];
        float sig2eff = tag_r2sum_value * nvec[tag_index] * sig2_beta + sig2_zero;
        float s = sqrt(sig2eff);

        for (int z_index = 0; z_index < length; z_index++) {
          double pdf_tmp = static_cast<double>(gaussian_pdf<FLOAT_TYPE>(zvec[z_index], s));
          pdf_double_local[z_index] += pi_k * pdf_tmp * tag_weight;
        }
      }
    }
#pragma omp critical
    pdf_double += pdf_double_local;
  }

  for (int i = 0; i < length; i++) pdf[i] = static_cast<float>(pdf_double[i]);
  LOG << "<calc_univariate_pdf(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::calc_univariate_power(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, float zthresh, int length, float* nvec, float* svec) {
  // input buffer "nvec" contains a set of sample sizes (N) to calculate power
  // output buffer svec(n), e.i. a fraction of heritability explained by genome-wide significant SNPs, aggregated across all SNPs with corresponding weights.
  // note a special case: when  nvec contains just one element, the svec buffer is expected to length num_tag, so that S(N) can be saved for each SNP independently.
  //
  // NB! This function uses analytical formula for the following double integral:
  // C(z) = E(\delta^2 | z) * P(z) = \int_{z : |z| > zt} \int_{delta} p(z|delta) p(delta) delta^2 d[delta] d[z]
  // As long as we fix the set of causal variants p(delta) is just a normal distribution with zero mean variance "delta2eff" (see code below).
  // In this case the integral can be taken analytically, and benefit from the fact that both numerator and denominator in S(N) formula are additive across
  // tag SNPs and across resampling iterations (1..kmax).

  float num_causals = pi_vec * static_cast<float>(num_snp_);
  if ((int)num_causals >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  const int component_id = 0;   // univariate is always component 0.

  LOG << ">calc_univariate_power(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ", zthresh=" << zthresh << ", length(nvec)=" << length << ")";
  SimpleTimer timer(-1);

  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) find_tag_r2sum(component_id, num_causals);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices_w(&deftag_indices);
  const double pi_k = 1.0 / static_cast<double>(k_max_);

  const int out_length = (length > 1) ? length : num_tag_;

  std::valarray<double> s_numerator_global(0.0, out_length);
  std::valarray<double> s_denominator_global(0.0, out_length);

#pragma omp parallel
  {
    std::valarray<double> s_numerator_local(0.0, out_length);
    std::valarray<double> s_denominator_local(0.0, out_length);
    std::vector<float> tag_r2sum(num_tag_, 0.0f);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {

      if (cache_tag_r2sum_) {
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          tag_r2sum[tag_index] = (*tag_r2sum_[0])(tag_index, k_index);
        }
      }
      else {
        find_tag_r2sum_no_cache(0, num_causals, k_index, &tag_r2sum);
      }

      for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
        int tag_index = deftag_indices[deftag_index];
        const double tag_weight = static_cast<double>(weights_[tag_index]);

        float tag_r2sum_value = tag_r2sum[tag_index];
        for (int n_index = 0; n_index < length; n_index++) {
          const int out_index = (length > 1) ? n_index : tag_index;
          float delta2eff = tag_r2sum_value * nvec[n_index] * sig2_beta;
          float sig2eff = delta2eff + sig2_zero;
          float sqrt_sig2eff = sqrt(sig2eff);
          static const float sqrt_2 = sqrtf(2.0);
          float numerator1 = gaussian_pdf<FLOAT_TYPE>(zthresh, sqrt_sig2eff) * 2 * delta2eff * delta2eff * zthresh / sig2eff;
          float numerator2 = std::erfcf(zthresh / (sqrt_2 * sqrt_sig2eff)) * delta2eff;
          s_numerator_local[out_index] += tag_weight*(numerator1 + numerator2);
          s_denominator_local[out_index] += tag_weight*delta2eff;
        }
      }
    }
#pragma omp critical
    {
      s_numerator_global += s_numerator_local;
      s_denominator_global += s_denominator_local;
    }
  }

  for (int i = 0; i < out_length; i++) svec[i] = static_cast<float>(s_numerator_global[i] / s_denominator_global[i]);
  LOG << "<calc_univariate_power(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ", zthresh=" << zthresh << ", length(nvec)=" << length << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::calc_univariate_delta_posterior(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* c0, float* c1, float* c2) {
  // c0 = c(0), c1=c(1), c2=c(2), where c(q) = \int_\delta \delta^q P(z|delta) P(delta)
  // c(q) is define so that:
  //  E(\delta^2|z_j) = c2[j]/c0[j];
  //  E(\delta  |z_j) = c1[j]/c0[j];

  if ((length == 0) || (length != num_tag_)) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_tag_"));

  float num_causals = pi_vec * static_cast<float>(num_snp_);
  if ((int)num_causals >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  const int component_id = 0;   // univariate is always component 0.

  LOG << ">calc_univariate_delta_posterior(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ", length(nvec)=" << length << ")";
  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices_znw(trait_index, &deftag_indices);
  
  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) find_tag_r2sum(component_id, num_causals);

  std::valarray<double> c0_global(0.0f, num_tag_);
  std::valarray<double> c1_global(0.0f, num_tag_);
  std::valarray<double> c2_global(0.0f, num_tag_);

#pragma omp parallel
  {
    std::vector<float> tag_r2sum(num_tag_, 0.0f);
    std::valarray<double> c0_local(0.0f, num_tag_);
    std::valarray<double> c1_local(0.0f, num_tag_);
    std::valarray<double> c2_local(0.0f, num_tag_);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {
      if (cache_tag_r2sum_) {
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          tag_r2sum[tag_index] = (*tag_r2sum_[0])(tag_index, k_index);
        }
      }
      else {
        find_tag_r2sum_no_cache(0, num_causals, k_index, &tag_r2sum);
      }

      for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
        int tag_index = deftag_indices[deftag_index];

        const float tag_r2sum_value = tag_r2sum[tag_index];
        const float delta2eff = tag_r2sum_value * nvec[tag_index] * sig2_beta;  // S^2_kj
        const float sig2eff = delta2eff + sig2_zero;
        const float sig2eff_1_2 = sqrt(sig2eff);
        const float sig2eff_3_2 = sig2eff_1_2 * sig2eff;
        const float sig2eff_5_2 = sig2eff_3_2 * sig2eff;

        const float z = z_minus_fixed_effect_delta[tag_index];
        const float exp_common = std::exp(-0.5f*z*z / sig2eff);

        c0_local[tag_index] += (exp_common / sig2eff_1_2);
        c1_local[tag_index] += (exp_common / sig2eff_3_2) * z * delta2eff;
        c2_local[tag_index] += (exp_common / sig2eff_5_2) *     delta2eff * (sig2_zero*sig2_zero + sig2_zero*delta2eff + z*z*delta2eff);
      }
    }

#pragma omp critical
    {
      c0_global += c0_local;
      c1_global += c1_local;
      c2_global += c2_local;
    }
  }

  // save results to output buffers
  const double pi_k = 1.0 / static_cast<double>(k_max_);
  static const double inv_sqrt_2pi = 0.3989422804014327;
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];
    c0[tag_index] = pi_k * inv_sqrt_2pi * c0_global[tag_index];
    c1[tag_index] = pi_k * inv_sqrt_2pi * c1_global[tag_index];
    c2[tag_index] = pi_k * inv_sqrt_2pi * c2_global[tag_index];
  }

  LOG << "<calc_univariate_delta_posterior(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ", length(nvec)=" << length << "), elapsed time " << timer.elapsed_ms() << "ms";
}

double BgmgCalculator::calc_univariate_cost(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  double cost;
  if (cost_calculator_ == CostCalculator_Gaussian) cost = calc_univariate_cost_fast(trait_index, pi_vec, sig2_zero, sig2_beta);
  else if (cost_calculator_ == CostCalculator_Convolve) cost = calc_univariate_cost_convolve(trait_index, pi_vec, sig2_zero, sig2_beta);
  else if (!cache_tag_r2sum_) cost = calc_univariate_cost_nocache(trait_index, pi_vec, sig2_zero, sig2_beta);
  else cost = calc_univariate_cost_cache(trait_index, pi_vec, sig2_zero, sig2_beta);
  return cost;
}

double BgmgCalculator::calc_univariate_cost_cache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  float num_causals = pi_vec * static_cast<float>(num_snp_);
  if ((int)num_causals >= max_causals_) return 1e100; // too large pi_vec
  const int component_id = 0;   // univariate is always component 0.
    
  LOG << ">calc_univariate_cost(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";
  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices_znw(trait_index, &deftag_indices);
  const double zmax = (trait_index==1) ? z1max_ : z2max_;
  const double pi_k = 1.0 / static_cast<double>(k_max_);

  find_tag_r2sum(component_id, num_causals);  // tag_r2sum_ is adjusted for ld_tag_sum_r2_below_r2min_adjust_for_hvec
  
  double log_pdf_total = 0.0;
  int num_infinite = 0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_infinite)
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];

    double pdf_tag = 0.0f;
    for (int k_index = 0; k_index < k_max_; k_index++) {
      float tag_r2sum = (*tag_r2sum_[component_id])(tag_index, k_index);
      float sig2eff = tag_r2sum * nvec[tag_index] * sig2_beta + sig2_zero;

      const float tag_z = z_minus_fixed_effect_delta[tag_index];
      float s = sqrt(sig2eff);
      const bool censoring = std::abs(tag_z) > zmax;

      double pdf = static_cast<double>(censoring ? censored_cdf<FLOAT_TYPE>(zmax, s) : gaussian_pdf<FLOAT_TYPE>(tag_z, s));
      pdf_tag += pi_k * pdf;
    }
    double increment = -std::log(pdf_tag) * static_cast<double>(weights_[tag_index]);
    if (!std::isfinite(increment)) num_infinite++;
    else log_pdf_total += increment;
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_univariate_cost(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

template<typename T>
double calc_univariate_cost_nocache_template(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, BgmgCalculator& rhs) {
  float num_causals = pi_vec * static_cast<float>(rhs.num_snp_);
  if ((int)num_causals >= rhs.max_causals_) return 1e100; // too large pi_vec
  const int component_id = 0;   // univariate is always component 0.
  if (rhs.snp_order_.empty()) rhs.find_snp_order();

  LOG << ">calc_univariate_cost_nocache(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";
  SimpleTimer timer(-1);

  const double pi_k = 1.0 / static_cast<double>(rhs.k_max_);

  rhs.k_pdf_.assign(rhs.k_max_, 0.0f);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; rhs.find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*rhs.get_nvec(trait_index));
  std::vector<int> deftag_indices; const int num_deftag = rhs.find_deftag_indices_znw(trait_index, &deftag_indices);
  const double zmax = (trait_index==1) ? rhs.z1max_ : rhs.z2max_;

  std::valarray<double> pdf_double(0.0, rhs.num_tag_);
#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, rhs.num_tag_);
    std::vector<float> tag_r2sum(rhs.num_tag_, 0.0f);

#pragma omp for schedule(static)
      for (int k_index = 0; k_index < rhs.k_max_; k_index++) {

        rhs.find_tag_r2sum_no_cache(component_id, num_causals, k_index, &tag_r2sum);  // tag_r2sum are adjusted for ld_tag_sum_r2_below_r2min_adjust_for_hvec
        for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
          int tag_index = deftag_indices[deftag_index];

          float tag_r2sum_value = tag_r2sum[tag_index];
          float sig2eff = tag_r2sum_value * nvec[tag_index] * sig2_beta + sig2_zero;

          const float tag_z = z_minus_fixed_effect_delta[tag_index];
          float s = sqrt(sig2eff);
          const bool censoring = std::abs(tag_z) > zmax;
          double pdf = static_cast<double>(censoring ? censored_cdf<T>(zmax, s) : gaussian_pdf<T>(tag_z, s));
          pdf_double_local[tag_index] += pdf * pi_k;

          if (rhs.calc_k_pdf_) rhs.k_pdf_[k_index] += (-std::log(pdf) * rhs.weights_[tag_index]);
        }
      }
#pragma omp critical
      pdf_double += pdf_double_local;
  }

  double log_pdf_total = 0.0;
  int num_infinite = 0;
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];
    double increment = -std::log(pdf_double[tag_index]) * rhs.weights_[tag_index];
    if (!std::isfinite(increment)) num_infinite++;
    else log_pdf_total += increment;
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_univariate_cost_nocache(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_univariate_cost_nocache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  return calc_univariate_cost_nocache_template<FLOAT_TYPE>(trait_index, pi_vec, sig2_zero, sig2_beta, *this);
}
double BgmgCalculator::calc_univariate_cost_nocache_float(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  return calc_univariate_cost_nocache_template<float>(trait_index, pi_vec, sig2_zero, sig2_beta, *this);
}
double BgmgCalculator::calc_univariate_cost_nocache_double(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  return calc_univariate_cost_nocache_template<double>(trait_index, pi_vec, sig2_zero, sig2_beta, *this);
}

std::string calc_bivariate_params_to_str(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, int length) {
  std::stringstream ss;
  ss << "pi_vec=[" << pi_vec[0] << ", " << pi_vec[1] << ", " << pi_vec[2] << "], "
     << "sig2_beta=[" << sig2_beta[0] << ", " << sig2_beta[1] << "], "
     << "rho_beta=" << rho_beta << ", "
     << "sig2_zero=[" << sig2_zero[0] << ", " << sig2_zero[1] << "], "
     << "rho_zero=" << rho_zero;
  if (length >= 0) ss << ", length(zvec)=" << length;
  return ss.str();
}

double BgmgCalculator::calc_bivariate_cost(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  if (num_components_ != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: require num_components == 3. Remember to call set_option('num_components', 3)."));
  if (sig2_beta_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_beta_len != 2"));
  if (sig2_zero_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_zero_len != 2"));
  if (pi_vec_len != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: pi_vec_len != 3"));

  double cost;
  if (cost_calculator_==CostCalculator_Gaussian) cost = calc_bivariate_cost_fast(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);
  else if (cost_calculator_==CostCalculator_Convolve) cost = calc_bivariate_cost_convolve(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);
  else if (!cache_tag_r2sum_) cost = calc_bivariate_cost_nocache(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);
  else cost = calc_bivariate_cost_cache(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero);
  return cost;
}

double BgmgCalculator::calc_bivariate_cost_cache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {

  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost(" << ss << ")";
  SimpleTimer timer(-1);

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) return 1e100; // too large pi_vec
  }

  for (int component_id = 0; component_id < 3; component_id++) {
    find_tag_r2sum(component_id, num_causals[component_id]);
  }

  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices_znw(&deftag_indices);
  const double pi_k = 1.0 / static_cast<double>(k_max_);

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  double log_pdf_total = 0.0;
  int num_infinite = 0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_infinite)
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];

    const float z1 = z1_minus_fixed_effect_delta[tag_index];
    const float z2 = z2_minus_fixed_effect_delta[tag_index];
    const float n1 = nvec1_[tag_index];
    const float n2 = nvec2_[tag_index];

    double pdf_tag = 0.0f;
    for (int k_index = 0; k_index < k_max_; k_index++) {
      const float tag_r2sum_c1 = (*tag_r2sum_[0])(tag_index, k_index);
      const float tag_r2sum_c2 = (*tag_r2sum_[1])(tag_index, k_index);
      const float tag_r2sum_c3 = (*tag_r2sum_[2])(tag_index, k_index);

      // Sigma  = [A1+A3  B3;  B3  C2+C3] + Sigma0 = ...
      //        = [a11    a12; a12   a22]
      const float A1 = tag_r2sum_c1 * n1 * sig2_beta[0];
      const float C2 = tag_r2sum_c2 * n2 * sig2_beta[1];
      const float A3 = tag_r2sum_c3 * n1 * sig2_beta[0];
      const float C3 = tag_r2sum_c3 * n2 * sig2_beta[1];
      const float B3 = sqrt(A3*C3) * rho_beta;

      const float a11 = A1 + A3 + a0;
      const float a22 = C2 + C3 + c0;
      const float a12 =      B3 + b0;

      const bool censoring = (std::abs(z1) > z1max_) || (std::abs(z2) > z2max_);
      const double pdf = static_cast<double>(censoring ? censored2_cdf<FLOAT_TYPE>(z1max_, z2max_, a11, a12, a22) : gaussian2_pdf<FLOAT_TYPE>(z1, z2, a11, a12, a22));
      pdf_tag += pi_k * pdf;
    }

    double increment = static_cast<double>(-std::log(pdf_tag) * weights_[tag_index]);
    if (!std::isfinite(increment)) num_infinite++;
    else log_pdf_total += increment;
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_bivariate_cost(" << ss << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_bivariate_cost_nocache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost_nocache(" << ss << ")";
  SimpleTimer timer(-1);

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) return 1e100; // too large pi_vec
  }

  if (snp_order_.empty()) find_snp_order();

  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices_znw(&deftag_indices);
  const double pi_k = 1.0 / static_cast<double>(k_max_);

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  std::valarray<double> pdf_double(0.0, num_tag_);
#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, num_tag_);
    std::vector<float> tag_r2sum0(num_tag_, 0.0f);
    std::vector<float> tag_r2sum1(num_tag_, 0.0f);
    std::vector<float> tag_r2sum2(num_tag_, 0.0f);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {

      find_tag_r2sum_no_cache(0, num_causals[0], k_index, &tag_r2sum0);
      find_tag_r2sum_no_cache(1, num_causals[1], k_index, &tag_r2sum1);
      find_tag_r2sum_no_cache(2, num_causals[2], k_index, &tag_r2sum2);

      for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
        int tag_index = deftag_indices[deftag_index];

        const float z1 = z1_minus_fixed_effect_delta[tag_index];
        const float z2 = z2_minus_fixed_effect_delta[tag_index];
        const float n1 = nvec1_[tag_index];
        const float n2 = nvec2_[tag_index];

        float pdf_tag = 0.0f;

        const float tag_r2sum_c1 = tag_r2sum0[tag_index];
        const float tag_r2sum_c2 = tag_r2sum1[tag_index];
        const float tag_r2sum_c3 = tag_r2sum2[tag_index];

        // Sigma  = [A1+A3  B3;  B3  C2+C3] + Sigma0 = ...
        //        = [a11    a12; a12   a22]
        const float A1 = tag_r2sum_c1 * n1 * sig2_beta[0];
        const float C2 = tag_r2sum_c2 * n2 * sig2_beta[1];
        const float A3 = tag_r2sum_c3 * n1 * sig2_beta[0];
        const float C3 = tag_r2sum_c3 * n2 * sig2_beta[1];
        const float B3 = sqrt(A3*C3) * rho_beta;

        const float a11 = A1 + A3 + a0;
        const float a22 = C2 + C3 + c0;
        const float a12 = B3 + b0;

        const bool censoring = (std::abs(z1) > z1max_) || (std::abs(z2) > z2max_);
        const double pdf = static_cast<double>(censoring ? censored2_cdf<FLOAT_TYPE>(z1max_, z2max_, a11, a12, a22) : gaussian2_pdf<FLOAT_TYPE>(z1, z2, a11, a12, a22));
        pdf_double_local[tag_index] += pi_k * pdf;
      }
    }
#pragma omp critical
    pdf_double += pdf_double_local;
  }

  double log_pdf_total = 0.0;
  int num_infinite = 0;
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];
    double increment = -std::log(pdf_double[tag_index]) * weights_[tag_index];
    if (!std::isfinite(increment)) num_infinite++;
    else log_pdf_total += increment;
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_bivariate_cost_nocache(" << ss << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

int64_t BgmgCalculator::calc_bivariate_pdf(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, int length, float* zvec1, float* zvec2, float* pdf) {
  // input buffer "zvec1" and "zvec2" contains z scores (presumably an equally spaced grid)
  // output buffer contains pdf(z), aggregated across all SNPs with corresponding weights
  if (nvec1_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec1 is not set"));
  if (nvec2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec2 is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
  if (num_components_ != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_pdf: require num_components == 3. Remember to call set_option('num_components', 3)."));
  if (sig2_beta_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_beta_len != 2"));
  if (sig2_zero_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_zero_len != 2"));
  if (pi_vec_len != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: pi_vec_len != 3"));

  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, length);
  LOG << ">calc_bivariate_pdf(" << ss << ")";
  SimpleTimer timer(-1);

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  }

  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices_nw(&deftag_indices);
  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) {
    for (int component_id = 0; component_id < 3; component_id++) {
      find_tag_r2sum(component_id, num_causals[component_id]);
    }
  }

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  // pi_k is mixture weight
  const double pi_k = 1.0 / static_cast<double > (k_max_);

  std::valarray<double> pdf_double(0.0, length);

#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, length);
    std::vector<float> tag_r2sum0(num_tag_, 0.0f);
    std::vector<float> tag_r2sum1(num_tag_, 0.0f);
    std::vector<float> tag_r2sum2(num_tag_, 0.0f);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {
      if (cache_tag_r2sum_) {
        for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
          tag_r2sum0[tag_index] = (*tag_r2sum_[0])(tag_index, k_index);
          tag_r2sum1[tag_index] = (*tag_r2sum_[1])(tag_index, k_index);
          tag_r2sum2[tag_index] = (*tag_r2sum_[2])(tag_index, k_index);
        }
      } else {
        find_tag_r2sum_no_cache(0, num_causals[0], k_index, &tag_r2sum0);
        find_tag_r2sum_no_cache(1, num_causals[1], k_index, &tag_r2sum1);
        find_tag_r2sum_no_cache(2, num_causals[2], k_index, &tag_r2sum2);
      }

      for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
        int tag_index = deftag_indices[deftag_index];

        double tag_weight = static_cast<double>(weights_[tag_index]);

        const float n1 = nvec1_[tag_index];
        const float n2 = nvec2_[tag_index];

        const float tag_r2sum_c1 = tag_r2sum0[tag_index];
        const float tag_r2sum_c2 = tag_r2sum1[tag_index];
        const float tag_r2sum_c3 = tag_r2sum2[tag_index];

        // Sigma  = [A1+A3  B3;  B3  C2+C3] + Sigma0 = ...
        //        = [a11    a12; a12   a22]
        const float A1 = tag_r2sum_c1 * n1 * sig2_beta[0];
        const float C2 = tag_r2sum_c2 * n2 * sig2_beta[1];
        const float A3 = tag_r2sum_c3 * n1 * sig2_beta[0];
        const float C3 = tag_r2sum_c3 * n2 * sig2_beta[1];
        const float B3 = sqrt(A3*C3) * rho_beta;

        const float a11 = A1 + A3 + a0;
        const float a22 = C2 + C3 + c0;
        const float a12 = B3 + b0;
        
        for (int z_index = 0; z_index < length; z_index++) {
          double pdf_tmp = static_cast<double>(gaussian2_pdf<FLOAT_TYPE>(zvec1[z_index], zvec2[z_index], a11, a12, a22));
          pdf_double_local[z_index] += pi_k * pdf_tmp * tag_weight;
        }
      }
    }
#pragma omp critical
    pdf_double += pdf_double_local;
  }

  for (int i = 0; i < length; i++) pdf[i] = static_cast<float>(pdf_double[i]);
  LOG << "<calc_bivariate_pdf(" << ss << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}


int64_t BgmgCalculator::calc_bivariate_delta_posterior(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero,
                                                       int length, float* c00, float* c10, float* c01, float* c20, float* c11, float* c02) {
  // where c(i,j) = \int_{\delta1, \delta2} \delta1^i \delta2^j P(z1, z2 | delta1, delta2) P(delta1, delta2)
  if (num_components_ != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: require num_components == 3. Remember to call set_option('num_components', 3)."));
  if (sig2_beta_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_beta_len != 2"));
  if (sig2_zero_len != 2) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: sig2_zero_len != 2"));
  if (pi_vec_len != 3) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_bivariate_cost: pi_vec_len != 3"));
  if ((length == 0) || (length != num_tag_)) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_tag_"));

  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, length);
  LOG << ">calc_bivariate_delta_posterior(" << ss << ")";
  SimpleTimer timer(-1);

  float num_causals[3];
  for (int component_id = 0; component_id < 3; component_id++) {
    num_causals[component_id] = pi_vec[component_id] * static_cast<float>(num_snp_);
    if ((int)num_causals[component_id] >= max_causals_) BGMG_THROW_EXCEPTION(::std::runtime_error("too large values in pi_vec"));
  }

  if (snp_order_.empty()) find_snp_order();
  if (cache_tag_r2sum_) {
    for (int component_id = 0; component_id < 3; component_id++) {
      find_tag_r2sum(component_id, num_causals[component_id]);
    }
  }

  std::vector<int> deftag_indices;
  const int num_deftag = find_deftag_indices_znw(&deftag_indices);

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zero[0];
  const float c0 = sig2_zero[1];
  const float b0 = sqrt(a0 * c0) * rho_zero;

  std::valarray<double> c00_global(0.0f, num_tag_);
  std::valarray<double> c10_global(0.0f, num_tag_);
  std::valarray<double> c01_global(0.0f, num_tag_);
  std::valarray<double> c20_global(0.0f, num_tag_);
  std::valarray<double> c11_global(0.0f, num_tag_);
  std::valarray<double> c02_global(0.0f, num_tag_);

#pragma omp parallel
  {
    std::vector<float> tag_r2sum0(num_tag_, 0.0f);
    std::vector<float> tag_r2sum1(num_tag_, 0.0f);
    std::vector<float> tag_r2sum2(num_tag_, 0.0f);
    std::valarray<double> c00_local(0.0f, num_tag_);
    std::valarray<double> c10_local(0.0f, num_tag_);
    std::valarray<double> c01_local(0.0f, num_tag_);
    std::valarray<double> c20_local(0.0f, num_tag_);
    std::valarray<double> c11_local(0.0f, num_tag_);
    std::valarray<double> c02_local(0.0f, num_tag_);

#pragma omp for schedule(static)
    for (int k_index = 0; k_index < k_max_; k_index++) {
      if (cache_tag_r2sum_) {
        for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
          int tag_index = deftag_indices[deftag_index];
          tag_r2sum0[tag_index] = (*tag_r2sum_[0])(tag_index, k_index);
          tag_r2sum1[tag_index] = (*tag_r2sum_[1])(tag_index, k_index);
          tag_r2sum2[tag_index] = (*tag_r2sum_[2])(tag_index, k_index);
        }
      } else {
        find_tag_r2sum_no_cache(0, num_causals[0], k_index, &tag_r2sum0);
        find_tag_r2sum_no_cache(1, num_causals[1], k_index, &tag_r2sum1);
        find_tag_r2sum_no_cache(2, num_causals[2], k_index, &tag_r2sum2);
      }

      for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
        int tag_index = deftag_indices[deftag_index];

        const float n1 = nvec1_[tag_index];
        const float n2 = nvec2_[tag_index];
        const float z1 = zvec1_[tag_index];
        const float z2 = zvec2_[tag_index];

        const float tag_r2sum_c1 = tag_r2sum0[tag_index];
        const float tag_r2sum_c2 = tag_r2sum1[tag_index];
        const float tag_r2sum_c3 = tag_r2sum2[tag_index];

        // Sigma  = [A1+A3  B3;  B3  C2+C3] + Sigma0 = ...
        //        = [a11    a12; a12   a22]
        const float A1 = tag_r2sum_c1 * n1 * sig2_beta[0];
        const float C2 = tag_r2sum_c2 * n2 * sig2_beta[1];
        const float A3 = tag_r2sum_c3 * n1 * sig2_beta[0];
        const float C3 = tag_r2sum_c3 * n2 * sig2_beta[1];
        const float B3 = sqrt(A3*C3) * rho_beta;

        const float A = A1 + A3;
        const float C = C2 + C3;
        const float B = B3;

        float c00buf, c10buf, c01buf, c20buf, c11buf, c02buf;
        BgmgCalculator::calc_bivariate_delta_posterior_integrals(a0, b0, c0, A, B, C, z1, z2, &c00buf, &c10buf, &c01buf, &c20buf, &c11buf, &c02buf);
        c00_local[tag_index] += static_cast<double>(c00buf);
        c10_local[tag_index] += static_cast<double>(c10buf);
        c01_local[tag_index] += static_cast<double>(c01buf);
        c20_local[tag_index] += static_cast<double>(c20buf);
        c11_local[tag_index] += static_cast<double>(c11buf);
        c02_local[tag_index] += static_cast<double>(c02buf);
      }
    }
#pragma omp critical
    {
      c00_global += c00_local;
      c10_global += c10_local;
      c01_global += c01_local;
      c20_global += c20_local;
      c11_global += c11_local;
      c02_global += c02_local;
    }
  }

  // save results to output buffers
  const double pi_k = 1.0 / static_cast<double>(k_max_);
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];
    c00[tag_index] = pi_k * c00_global[tag_index];
    c10[tag_index] = pi_k * c10_global[tag_index];
    c01[tag_index] = pi_k * c01_global[tag_index];    
    c20[tag_index] = pi_k * c20_global[tag_index];
    c11[tag_index] = pi_k * c11_global[tag_index];
    c02[tag_index] = pi_k * c02_global[tag_index];    
  }

  LOG << "<calc_bivariate_delta_posterior(" << ss << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

// Use an approximation that preserves variance and kurtosis.
// This gives a robust cost function that scales up to a very high pivec, including infinitesimal model pi==1.
double BgmgCalculator::calc_univariate_cost_fast(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  std::stringstream ss;
  ss << "calc_univariate_cost_fast(trait_index=" << trait_index << ", pi_vec=" << pi_vec << ", sig2_zero=" << sig2_zero << ", sig2_beta=" << sig2_beta << ")";
  LOG << ">" << ss.str();

  double log_pdf_total = 0.0;
  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);
  std::vector<float>& nvec(*get_nvec(trait_index));
  const double zmax = (trait_index==1) ? z1max_ : z2max_;

  int num_zero_tag_r2 = 0;
  int num_infinite = 0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_zero_tag_r2, num_infinite)
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];
    double tag_weight = static_cast<double>(weights_[tag_index]);
    
    const float tag_r2inf = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min()[tag_index];
    const float tag_r2 = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_above_r2min()[tag_index];
    const float tag_r4 = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r4_above_r2min()[tag_index];

    if (tag_r2 == 0 || tag_r4 == 0) {
      num_zero_tag_r2++; continue;
    }

    const float tag_chi = tag_r4 / tag_r2;

    const float tag_eta_factor = pi_vec * tag_r2 + (1.0f - pi_vec) * tag_chi;
    const double tag_pi1 = static_cast<double>(pi_vec * tag_r2 / tag_eta_factor);
    const double tag_pi0 = 1.0 - tag_pi1;
    const float tag_sig2beta = sig2_beta * tag_eta_factor;

    const float tag_z = z_minus_fixed_effect_delta[tag_index];
    const float tag_n = nvec[tag_index];

    const double tag_sig2inf = tag_r2inf * pi_vec * sig2_beta * tag_n;

    const bool censoring = std::abs(tag_z) > zmax;
    const float s1 = sqrt(sig2_zero + tag_sig2inf);
    const float s2 = sqrt(sig2_zero + tag_sig2inf + tag_n *tag_sig2beta);

    // printf("%.4f %.4f %.4f %.4f %.4f\n", tag_pi0, s1, tag_pi1, s2, tag_pi0*s1*s1+tag_pi1*s2*s2);
    const double tag_pdf0 = static_cast<double>(censoring ? censored_cdf<FLOAT_TYPE>(zmax, s1) : gaussian_pdf<FLOAT_TYPE>(tag_z, s1));
    const double tag_pdf1 = static_cast<double>(censoring ? censored_cdf<FLOAT_TYPE>(zmax, s2) : gaussian_pdf<FLOAT_TYPE>(tag_z, s2));
    const double tag_pdf = tag_pi0 * tag_pdf0 + tag_pi1 * tag_pdf1;
    const double increment = (-std::log(tag_pdf) * tag_weight);
    if (!std::isfinite(increment)) num_infinite++;
    else log_pdf_total += increment;
  }

  if (num_zero_tag_r2 > 0)
    LOG << " warning: zero tag_r2 encountered " << num_zero_tag_r2 << " times";
  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<" << ss.str() << ", cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_bivariate_cost_fast(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  std::string ss = calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost_fast(" << ss << ")";

  double log_pdf_total = 0.0;
  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices_znw(&deftag_indices);

  int num_zero_tag_r2 = 0;
  int num_infinite = 0;

  const float s0_a11 = sig2_zero[0];
  const float s0_a22 = sig2_zero[1];
  const float s0_a12 = sqrt(sig2_zero[0] * sig2_zero[1]) * rho_zero;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_zero_tag_r2, num_infinite)
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];

    const float z1 = z1_minus_fixed_effect_delta[tag_index];
    const float n1 = nvec1_[tag_index];
    const float z2 = z2_minus_fixed_effect_delta[tag_index];
    const float n2 = nvec2_[tag_index];

    const float tag_r2 = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_above_r2min()[tag_index]; // TBD: apply ld_tag_sum_r2_below_r2min as an infinitesimal model
    const float tag_r4 = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r4_above_r2min()[tag_index];

    if (tag_r2 == 0 || tag_r4 == 0) {
      num_zero_tag_r2++; continue;
    }

    const float tag_chi = tag_r4 / tag_r2;

    const float tag_eta_factor[3] = {
      pi_vec[0] * tag_r2 + (1.0f - pi_vec[0]) * tag_chi,
      pi_vec[1] * tag_r2 + (1.0f - pi_vec[1]) * tag_chi,
      pi_vec[2] * tag_r2 + (1.0f - pi_vec[2]) * tag_chi
    };

    const float tag_pi1[3] = {
      pi_vec[0] * tag_r2 / tag_eta_factor[0],
      pi_vec[1] * tag_r2 / tag_eta_factor[1],
      pi_vec[2] * tag_r2 / tag_eta_factor[2]
    };

    const float tag_pi0[3] = {
      1.0f - tag_pi1[0],
      1.0f - tag_pi1[1],
      1.0f - tag_pi1[2]
    };

    const float a11[3] = { tag_eta_factor[0] * n1 * sig2_beta[0], 0,                                     tag_eta_factor[2] * n1 * sig2_beta[0] };
    const float a22[3] = { 0,                                     tag_eta_factor[1] * n2 * sig2_beta[1], tag_eta_factor[2] * n2 * sig2_beta[1] };
    const float a12[3] = { 0,                                     0,                                     rho_beta * sqrt(a11[2] * a22[2]) };

    const float f0[8] = { 0,0,0,0,1,1,1,1 };
    const float f1[8] = { 0,0,1,1,0,0,1,1 };
    const float f2[8] = { 0,1,0,1,0,1,0,1 };

    const bool censoring = (std::abs(z1) > z1max_) || (std::abs(z2) > z2max_);

    FLOAT_TYPE tag_pdf = 0.0f;
    for (int i = 0; i < 8; i++) {
      const float pi1 = (f0[i] ? tag_pi1[0] : tag_pi0[0]);
      const float pi2 = (f1[i] ? tag_pi1[1] : tag_pi0[1]);
      const float pi3 = (f2[i] ? tag_pi1[2] : tag_pi0[2]);
      const float a11i = s0_a11 + f0[i] * a11[0] + f1[i] * a11[1] + f2[i] * a11[2];
      const float a22i = s0_a22 + f0[i] * a22[0] + f1[i] * a22[1] + f2[i] * a22[2];
      const float a12i = s0_a12 + f0[i] * a12[0] + f1[i] * a12[1] + f2[i] * a12[2];
      tag_pdf += static_cast<double>(pi1*pi2*pi3) * static_cast<double>(censoring ? censored2_cdf<FLOAT_TYPE>(z1max_, z2max_, a11i, a12i, a22i) : gaussian2_pdf<FLOAT_TYPE>(z1, z2, a11i, a12i, a22i));
    }

    if (tag_pdf <= 0)
      tag_pdf = 1e-100;

    double increment = static_cast<double>(-std::log(tag_pdf) * weights_[tag_index]);
    if (!std::isfinite(increment)) num_infinite++;
    else log_pdf_total += increment;
  }

  if (num_zero_tag_r2 > 0)
    LOG << " warning: zero tag_r2 encountered " << num_zero_tag_r2 << " times";
  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<calc_bivariate_cost_fast(" << ss << "), cost=" << log_pdf_total << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_univariate_cost_convolve(int trait_index, float pi_vec, float sig2_zero, float sig2_beta) {
  std::vector<float> pi_vec_constant_vector(num_snp_, pi_vec);
  std::vector<float> sig2_vec_constant_vector(num_snp_, sig2_beta);
  const float sig2_zeroA = sig2_zero;
  const float sig2_zeroC = 1.0f;
  const float sig2_zeroL = pi_vec * sig2_beta;
  return calc_unified_univariate_cost_convolve(trait_index, 1, num_snp_, &pi_vec_constant_vector[0], &sig2_vec_constant_vector[0], sig2_zeroA, sig2_zeroC, sig2_zeroL, nullptr);
}

void BgmgCalculator::clear_tag_r2sum(int component_id) {
  if (component_id < 0 || component_id >= num_components_) BGMG_THROW_EXCEPTION(::std::runtime_error("clear_tag_r2sum: component_id must be between 0 and num_components_"));
  LOG << " clear_tag_r2sum(component_id=" << component_id << ")";
  if (cache_tag_r2sum_) {
    if (tag_r2sum_.empty()) {
      // Initialize
      for (int i = 0; i < num_components_; i++) {
        last_num_causals_.push_back(0.0f);
        tag_r2sum_.push_back(std::make_shared<DenseMatrix<float>>(num_tag_, k_max_));
        tag_r2sum_[i]->InitializeZeros();
      }
    } else {
      // Clear just the requested component
      tag_r2sum_[component_id]->InitializeZeros();
      last_num_causals_[component_id] = 0;
    }
  } else {
    // Cache disabled, clear tag_r2sum to free up memory
    last_num_causals_.clear();
    tag_r2sum_.clear();
  }
}

double BgmgCalculator::calc_bivariate_cost_convolve(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  const int num_components = 3;
  const int num_traits = 2;
  std::vector<float> pi_vec_constant_vector(num_components * num_snp_);
  for (int comp_index = 0; comp_index < num_components; comp_index++)
    for (int snp_index = 0; snp_index < num_snp_; snp_index++)
      pi_vec_constant_vector[comp_index * num_snp_ + snp_index] = pi_vec[comp_index];

  std::vector<float> sig2_vec_constant_vector(num_traits * num_snp_);
  for (int trait_index = 0; trait_index < num_traits; trait_index++)
    for (int snp_index = 0; snp_index < num_snp_; snp_index++)
      sig2_vec_constant_vector[trait_index * num_snp_ + snp_index] = sig2_beta[trait_index];
    
  std::vector<float> rho_vec_constant_vector(num_snp_, rho_beta);
  float sig2_zeroC[2] = {1.0, 1.0};
  float sig2_zeroL[2] = {(pi_vec[0]+pi_vec[2]) * sig2_beta[0], (pi_vec[1] + pi_vec[2]) * sig2_beta[1]};
  float rho_zeroL = rho_beta * pi_vec[2] / sqrt((pi_vec[0]+pi_vec[2]) * (pi_vec[1]+pi_vec[2]));

  return calc_unified_bivariate_cost_convolve(num_snp_, &pi_vec_constant_vector[0], &sig2_vec_constant_vector[0], &rho_vec_constant_vector[0],
                                              sig2_zero, sig2_zeroC, sig2_zeroL, rho_zero, rho_zeroL, nullptr);
}

int64_t BgmgCalculator::set_snp_order(int component_id, int64_t length, const int* buffer) {
  if ((component_id < 0) || (component_id >= num_components_)) BGMG_THROW_EXCEPTION(::std::runtime_error("component_id out of range"));
  if (length != (max_causals_ * k_max_)) BGMG_THROW_EXCEPTION(::std::runtime_error("buffer length must be max_causals_ * k_max_"));
  if ((snp_order_.size() != num_components_) || (snp_order_[component_id] == nullptr)) BGMG_THROW_EXCEPTION(::std::runtime_error("snp_order_.size() != num_components_, or is empty"));
  if (snp_order_[component_id]->size() != length) BGMG_THROW_EXCEPTION(::std::runtime_error("snp_order_[component_id] has a wrong size"));
  for (int64_t k = 0; k < k_max_; k++)
    for (int64_t j = 0; j < max_causals_; j++)
      (*snp_order_[component_id])(j, k) = buffer[k*max_causals_ + j];
  LOG << " set_snp_order(component_id" << component_id << ")";
  return 0;
}

int64_t BgmgCalculator::retrieve_snp_order(int component_id, int64_t length, int* buffer) {
  if ((component_id < 0) || (component_id >= num_components_)) BGMG_THROW_EXCEPTION(::std::runtime_error("component_id out of range"));
  if (length != (max_causals_ * k_max_)) BGMG_THROW_EXCEPTION(::std::runtime_error("buffer length must be max_causals_ * k_max_"));
  if ((snp_order_.size() != num_components_) || (snp_order_[component_id] == nullptr)) BGMG_THROW_EXCEPTION(::std::runtime_error("snp_order_.size() != num_components_, or is empty"));
  if (snp_order_[component_id]->size() != length) BGMG_THROW_EXCEPTION(::std::runtime_error("snp_order_[component_id] has a wrong size"));
  for (int64_t k = 0; k < k_max_; k++)
    for (int64_t j = 0; j < max_causals_; j++)
      buffer[k*max_causals_ + j] = (*snp_order_[component_id])(j, k);
  LOG << " retrieve_snp_order(component_id" << component_id << ")";
  return 0;
}

int64_t BgmgCalculator::retrieve_k_pdf(int length, double* buffer) {
  if (k_pdf_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("k_pdf_ is not generated; make sure 'calc_univariate_cost_nocache' is called."));
  if (length != k_pdf_.size()) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  LOG << " retrieve_k_pdf()";
  for (int i = 0; i < k_pdf_.size(); i++) buffer[i] = k_pdf_[i];
  return 0;
}

void BgmgCalculator::find_tag_r2sum_no_cache(int component_id, float num_causal, int k_index, std::vector<float>* buffer) {
  assert(buffer->size() == num_tag_);
  if (snp_order_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("find_tag_r2sum_no_cache called before find_snp_order()"));

  std::vector<std::pair<int, float>> changeset;
  float floor_num_causals = floor(num_causal);
  for (int i = 0; i < (int)floor_num_causals; i++) changeset.push_back(std::make_pair(i, 1.0f));
  changeset.push_back(std::make_pair((int)floor_num_causals, num_causal - floor_num_causals));

  for (int i = 0; i < num_tag_; i++) buffer->at(i) = 0;

  std::vector<float> hvec;
  find_hvec(*this, &hvec);

  LdMatrixRow ld_matrix_row;
  for (auto change : changeset) {
    int scan_index = change.first;
    float scan_weight = change.second;
    int snp_index = (*snp_order_[component_id])(scan_index, k_index);
    ld_matrix_csr_.extract_snp_row(SnpIndex(snp_index), &ld_matrix_row);
    auto iter_end = ld_matrix_row.end();
    for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
      const float mafval = mafvec_[snp_index];
      const float hval = 2.0f * mafval * (1 - mafval);
      const int tag_index = iter.index();
      const float r2 = iter.r2();
      buffer->at(tag_index) += (scan_weight * r2 * hval);
    }
  }

  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();

  const float pival = num_causal / static_cast<float>(num_snp_);
  for (int i = 0; i < num_tag_; i++) {
    buffer->at(i) += (pival * ld_tag_sum_r2_below_r2min_adjust_for_hvec[i]);
  }
}

int BgmgCalculator::find_deftag_indices_znw(int trait_index, std::vector<int>* deftag_indices) {
  std::vector<float>& zvec(*get_zvec(trait_index));
  std::vector<float>& nvec(*get_nvec(trait_index));
  if (zvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec is not set"));
  if (nvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec[tag_index]) || !std::isfinite(nvec[tag_index])) continue;
    deftag_indices->push_back(tag_index);
  }
  return deftag_indices->size();
}

int BgmgCalculator::find_deftag_indices_nw(int trait_index, std::vector<int>* deftag_indices) {
  std::vector<float>& nvec(*get_nvec(trait_index));
  if (nvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(nvec[tag_index])) continue;
    deftag_indices->push_back(tag_index);
  }
  return deftag_indices->size();
}

int BgmgCalculator::find_deftag_indices_w(std::vector<int>* deftag_indices) {
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    deftag_indices->push_back(tag_index);
  }
  return deftag_indices->size();
}

int BgmgCalculator::find_deftag_indices_znw(std::vector<int>* deftag_indices) {
  if (zvec1_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec1 is not set"));
  if (nvec1_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec1 is not set"));
  if (zvec2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec2 is not set"));
  if (nvec2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec2 is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(zvec1_[tag_index]) || !std::isfinite(nvec1_[tag_index])) continue;
    if (!std::isfinite(zvec2_[tag_index]) || !std::isfinite(nvec2_[tag_index])) continue;
    deftag_indices->push_back(tag_index);
  }
  return deftag_indices->size();
}

int BgmgCalculator::find_deftag_indices_nw(std::vector<int>* deftag_indices) {
  if (nvec1_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec1 is not set"));
  if (nvec2_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec2 is not set"));
  if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));

  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights_[tag_index] == 0) continue;
    if (!std::isfinite(nvec1_[tag_index])) continue;
    if (!std::isfinite(nvec2_[tag_index])) continue;
    deftag_indices->push_back(tag_index);
  }
  return deftag_indices->size();
}

void BgmgCalculator::clear_state() {
  LOG << " clear_state";

  // clear ordering of SNPs
  snp_order_.clear();
  k_pdf_.clear();
  tag_r2sum_.clear();
  last_num_causals_.clear();
}