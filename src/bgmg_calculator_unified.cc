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

static const double kMinTagPdf = 1e-100;

double BgmgCalculator::calc_unified_univariate_cost(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux) {
  check_num_snp(num_snp);

  if (cost_calculator_ == CostCalculator_Gaussian) return calc_unified_univariate_cost_gaussian(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, aux);
  else if (cost_calculator_ == CostCalculator_Convolve) return calc_unified_univariate_cost_convolve(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, aux);
  else if (cost_calculator_ == CostCalculator_Sampling) return calc_unified_univariate_cost_sampling(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, aux, nullptr);
  else BGMG_THROW_EXCEPTION(::std::runtime_error("unsupported cost calculator in calc_unified_univariate_cost"));
}

class UnivariateCharacteristicFunctionData {
 public:
  int num_components;
  int num_snp;
  float* pi_vec;
  float* sig2_vec;
  float sig2_zeroA;
  float sig2_zeroC;
  float sig2_zeroL;
  int tag_index;  // for which SNP to calculate the characteristic function
                  // note that use_complete_tag_indices_ must be enabled for "convolve" calculator, 
                  // so "tag" and "snp" are in the same indexing
  LdMatrixRow* ld_matrix_row;
  const std::vector<float>* hvec;
  const std::vector<float>* nvec;
  const std::vector<float>* z_minus_fixed_effect_delta;
  const std::vector<float>* ld_tag_sum_r2_below_r2min_adjust_for_hvec;
  int func_evals;
};

int calc_univariate_characteristic_function_times_cosinus(unsigned ndim, const double *x, void *raw_data, unsigned fdim, double* fval) {
  assert(ndim == 1);
  assert(fdim == 1);
  const double t = x[0];
  const double minus_tsqr_half = -t*t/2.0;
  UnivariateCharacteristicFunctionData* data = (UnivariateCharacteristicFunctionData *)raw_data;
  const double nval = (*data->nvec)[data->tag_index];
  const double zval = (*data->z_minus_fixed_effect_delta)[data->tag_index];
  const double sig2_zero = data->sig2_zeroA + (*data->ld_tag_sum_r2_below_r2min_adjust_for_hvec)[data->tag_index] * nval * data->sig2_zeroL;
  const double m_1_pi = M_1_PI;

  double result = m_1_pi * cos(t * zval) * std::exp(minus_tsqr_half * sig2_zero);
  
  auto iter_end = data->ld_matrix_row->end();
  for (auto iter = data->ld_matrix_row->begin(); iter < iter_end; iter++) {
    int snp_index = iter.tag_index();  // yes, this is correct - snp_index on the LHS, tag_index on the RHS.
                                       // ld_matrix was designed to work with "sampling" calculator which require 
                                       // a mapping from a causal SNP "i" to all tag SNPs "j" that are in LD with "i".
                                       // however, for for "convolve" we are interested in mapping from a tag SNP "j"
                                       // to all causal SNPs "i" in LD with "j". To make this work we store a complete
                                       // LD matrix (e.g. _num_snp == _num_tag), and explore symmetry of the matrix. 
    const double r2 = iter.r2();
    const double hval = (*data->hvec)[snp_index];
    const double minus_tsqr_half_r2_hval_nval = minus_tsqr_half * r2 * hval * nval * (data->sig2_zeroC);
    double factor = 0.0;
    double pi_complement = 1.0;        // handle a situation where pi0 N(0, 0) is not specified as a column in pi_vec and sig2_vec.
    for (int comp_index = 0; comp_index < data->num_components; comp_index++) {
      const int index = (comp_index*data->num_snp + snp_index);
      const double pi_val = data->pi_vec[index];
      const double sig2_val = data->sig2_vec[index];
      factor += pi_val * std::exp(minus_tsqr_half_r2_hval_nval * sig2_val);
      pi_complement -= pi_val;
    }
    factor += pi_complement;

    result *= (double)(factor);
  }

  data->func_evals++;
  *fval = result;
  return 0; 
}

int calc_univariate_characteristic_function_for_integration(unsigned ndim, const double *x, void *raw_data, unsigned fdim, double* fval) {
  const double t = x[0];
  const double inv_1_minus_t = 1.0 / (1.0 - t);
  const double x_transformed = t * inv_1_minus_t;
  const double jacob = inv_1_minus_t * inv_1_minus_t;
  int retval = calc_univariate_characteristic_function_times_cosinus(ndim, &x_transformed, raw_data, fdim, fval);
  (*fval) *= jacob;
  return retval;
}

double BgmgCalculator::calc_unified_univariate_cost_convolve(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux) {
  if (!use_complete_tag_indices_) BGMG_THROW_EXCEPTION(::std::runtime_error("Convolve calculator require 'use_complete_tag_indices' option"));

  std::stringstream ss;
  ss << "trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL;
  LOG << ">calc_unified_univariate_cost_convolve(" << ss.str() << ")";

  double log_pdf_total = 0.0;
  int num_snp_failed = 0;
  int num_infinite = 0;
  double func_evals = 0.0;
  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*get_nvec(trait_index));
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  const double zmax = (trait_index==1) ? z1max_ : z2max_;

  std::vector<float> weights_convolve(weights_.begin(), weights_.end());
  std::vector<float> weights_sampling(weights_.begin(), weights_.end()); int num_deftag_sampling = 0;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    const float tag_z = z_minus_fixed_effect_delta[tag_index];
    const bool censoring = std::abs(tag_z) > zmax;
    if (censoring) {
      weights_convolve[tag_index] = 0;
      num_deftag_sampling++;
    } else {
      weights_sampling[tag_index] = 0;
    }
  }

  if (num_deftag_sampling > 0) {  // fall back to sampling approach for censored z-scores
    log_pdf_total += calc_unified_univariate_cost_sampling(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, aux, &weights_sampling[0]);
  }

  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(&weights_convolve[0], &deftag_indices);

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    UnivariateCharacteristicFunctionData data;
    data.num_components = num_components;
    data.num_snp = num_snp_;
    data.pi_vec = pi_vec;
    data.sig2_vec = sig2_vec;
    data.sig2_zeroA = sig2_zeroA;
    data.sig2_zeroC = sig2_zeroC;
    data.sig2_zeroL = sig2_zeroL;
    data.ld_matrix_row = &ld_matrix_row;
    data.hvec = &hvec;
    data.z_minus_fixed_effect_delta = &z_minus_fixed_effect_delta;
    data.nvec = &nvec;
    data.ld_tag_sum_r2_below_r2min_adjust_for_hvec = &ld_tag_sum_r2_below_r2min_adjust_for_hvec;

#pragma omp for schedule(static) reduction(+: log_pdf_total, num_snp_failed, num_infinite, func_evals)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      double tag_weight = static_cast<double>(weights_convolve[tag_index]);

      const int causal_index = tag_index; // yes,snp==tag in this case -- see a large comment in calc_univariate_characteristic_function_times_cosinus function.
      ld_matrix_csr_.extract_row(causal_index, data.ld_matrix_row);
      data.tag_index = tag_index;
      data.func_evals = 0;

      double tag_pdf = 0, tag_pdf_err = 0;
      const double xmin = 0, xmax = 1;
      const int integrand_fdim = 1, ndim = 1;
      int cubature_result = hcubature(integrand_fdim, calc_univariate_characteristic_function_for_integration,
        &data, ndim, &xmin, &xmax, cubature_max_evals_, cubature_abs_error_, cubature_rel_error_, ERROR_INDIVIDUAL, &tag_pdf, &tag_pdf_err);
      func_evals += (weights_convolve[tag_index] * (double)data.func_evals);
      if (cubature_result != 0) { num_snp_failed++; continue; }

      if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdf)) aux[tag_index] = tag_pdf;
      if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdfErr)) aux[tag_index] = tag_pdf_err;

      double increment = static_cast<double>(-std::log(tag_pdf) * weights_convolve[tag_index]);
      if (!std::isfinite(increment)) {
        increment = static_cast<double>(-std::log(kMinTagPdf) * weights_convolve[tag_index]);
        num_infinite++;
      }

      log_pdf_total += increment;
    }
  }    

  if (num_snp_failed > 0)
    LOG << " warning: hcubature failed for " << num_snp_failed << " tag snps";
  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  double total_weight = 0.0;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) total_weight += weights_convolve[tag_index];
  func_evals /= total_weight;

  LOG << "<calc_unified_univariate_cost_convolve(" << ss.str() << "), cost=" << log_pdf_total << ", evals=" << func_evals << ", num_deftag=" << num_deftag << "+" << num_deftag_sampling << ", elapsed time " << timer.elapsed_ms() << "ms";

  return log_pdf_total;
}

// Use an approximation that preserves variance and kurtosis.
// This gives a robust cost function that scales up to a very high pivec, including infinitesimal model pi==1.
double BgmgCalculator::calc_unified_univariate_cost_gaussian(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux) {
  std::stringstream ss;
  ss << "calc_unified_univariate_cost_gaussian(trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*get_nvec(trait_index));
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);
  const double zmax = (trait_index==1) ? z1max_ : z2max_;

  // Step 1. Calculate Ebeta2 and Ebeta4
  std::valarray<float> Ebeta2(0.0, num_snp_);
  std::valarray<float> Ebeta4(0.0, num_snp_);// Ebeta4 is a simplified name - in fact, this variable contains E(\beta^4) - 3 (E \beta^2)^2.
  for (int comp_index = 0; comp_index < num_components; comp_index++) {
    for (int snp_index = 0; snp_index < num_snp_; snp_index++) {
      const float p = pi_vec[comp_index*num_snp_ + snp_index];
      const float s2 = sig2_vec[comp_index*num_snp_ + snp_index];
      const float s4 = s2*s2;
      Ebeta2[snp_index] += p * s2;
      Ebeta4[snp_index] += 3.0f * p * s4;
    }
  }
  for (int snp_index = 0; snp_index < num_snp_; snp_index++) {
    Ebeta4[snp_index] -= (3.0f * Ebeta2[snp_index] * Ebeta2[snp_index]);
  }

  // Step 2. Calculate Edelta2 and Edelta4
  std::valarray<float> Edelta2(0.0, num_tag_);
  std::valarray<float> Edelta4(0.0, num_tag_);  // Edelta4 is a simplified name - see comment for Ebeta4

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    std::valarray<float> Edelta2_local(0.0, num_tag_);
    std::valarray<float> Edelta4_local(0.0, num_tag_);

    if (use_complete_tag_indices_) {
      // Optimize implemnetation when LD matrix is rectangular
      // In this case we may have SNPs with undefined zvec and nvec,
      // and for those we don't need to compute Edelta2 and Edelta4.
      // The alternative implementation (i.e. the "else" branch with !use_complete_tag_indices_)
      // has to compute Edelta2 and Edelta4 for all tag snps, because 
      // our LD matrix is stored as a mapping from causal to tag variants.

#pragma omp for schedule(static)    
      for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
        int tag_index = deftag_indices[deftag_index];

        const int snp_index = tag_index; // yes, snp==tag in this case -- same story here as in calc_univariate_characteristic_function_times_cosinus function.
        ld_matrix_csr_.extract_row(snp_index, &ld_matrix_row);
        auto iter_end = ld_matrix_row.end();
        for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
          const int causal_index = iter.tag_index();   // tag_index on RHS can be misleading here
          const float r2_value = iter.r2();
          const float a2ij = sig2_zeroC * nvec[tag_index] * hvec[causal_index] * r2_value;
          Edelta2_local[tag_index] += a2ij *        Ebeta2[causal_index];
          Edelta4_local[tag_index] += a2ij * a2ij * Ebeta4[causal_index];
        }
      }
    } else {  // use_complete_tag_indices_
#pragma omp for schedule(static)
      for (int causal_index = 0; causal_index < num_snp_; causal_index++) {
        ld_matrix_csr_.extract_row(causal_index, &ld_matrix_row);
        auto iter_end = ld_matrix_row.end();
        for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
          const int tag_index = iter.tag_index();
          const float r2_value = iter.r2();
          const float a2ij = sig2_zeroC * nvec[tag_index] * hvec[causal_index] * r2_value;
          Edelta2_local[tag_index] += a2ij *        Ebeta2[causal_index];
          Edelta4_local[tag_index] += a2ij * a2ij * Ebeta4[causal_index];
        }
      }
    }  // use_complete_tag_indices_

#pragma omp critical
    {
      Edelta2 += Edelta2_local;
      Edelta4 += Edelta4_local;
    }
  }  // parallel

  double log_pdf_total = 0.0;
  int num_zero_tag_r2 = 0;
  int num_infinite = 0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_zero_tag_r2, num_infinite)
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];
    if (Edelta2[tag_index] == 0) { num_zero_tag_r2++; continue;}
    double tag_weight = static_cast<double>(weights_[tag_index]);

    const float A = Edelta2[tag_index];
    const float B = Edelta4[tag_index];
    const float Ax3 = 3.0f*A;
    const float A2x3= A * Ax3;
    const float BplusA2x3 = B + A2x3;
    const float tag_pi0 = B / BplusA2x3;
    const float tag_pi1 = A2x3 / BplusA2x3;
    const float sig2_tag = BplusA2x3 / Ax3;

    // additive inflation, plus contribution from small LD r2 (those below r2min)
    const float sig2_zero = sig2_zeroA + ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index] * nvec[tag_index] * sig2_zeroL;

    // export the expected values of z^2 distribution
    if ((aux != nullptr) && (aux_option_ == AuxOption_Ezvec2)) aux[tag_index] = A + sig2_zero;

    const float tag_z = z_minus_fixed_effect_delta[tag_index];
    const float tag_n = nvec[tag_index];

    const bool censoring = (std::abs(tag_z) > zmax);
    const float s1 = sqrt(sig2_zero);
    const float s2 = sqrt(sig2_zero + sig2_tag);  // sqrt(tag_n) * 

    // printf("%.4f %.4f %.4f %.4f %.4f\n", tag_pi0, s1, tag_pi1, s2, tag_pi0*s1*s1+tag_pi1*s2*s2);
    const double tag_pdf0 = static_cast<double>(censoring ? censored_cdf<FLOAT_TYPE>(zmax, s1) : gaussian_pdf<FLOAT_TYPE>(tag_z, s1));
    const double tag_pdf1 = static_cast<double>(censoring ? censored_cdf<FLOAT_TYPE>(zmax, s2) : gaussian_pdf<FLOAT_TYPE>(tag_z, s2));
    const double tag_pdf = tag_pi0 * tag_pdf0 + tag_pi1 * tag_pdf1;
    if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdf)) aux[tag_index] = tag_pdf;
    double increment = (-std::log(tag_pdf) * tag_weight);
    if (!std::isfinite(increment)) {
      increment = static_cast<double>(-std::log(kMinTagPdf) * tag_weight);
      num_infinite++;
    }

    log_pdf_total += increment;
  }

  if (num_zero_tag_r2 > 0)
    LOG << " warning: zero tag_r2 encountered " << num_zero_tag_r2 << " times";
  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<" << ss.str() << ", cost=" << log_pdf_total << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_unified_univariate_cost_sampling(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux, const float* weights) {
  if (!use_complete_tag_indices_) BGMG_THROW_EXCEPTION(::std::runtime_error("Unified Sampling calculator require 'use_complete_tag_indices' option"));
  if (weights == nullptr) {
    if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
    weights = &weights_[0];
  }

  std::stringstream ss;
  ss << "calc_unified_univariate_cost_sampling(trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*get_nvec(trait_index));
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(weights, &deftag_indices);

  const double z_max = (trait_index==1) ? z1max_ : z2max_;
  const double pi_k = 1.0 / static_cast<double>(k_max_);
  double log_pdf_total = 0.0;
  int num_infinite = 0;

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + omp_get_thread_num(), k_max_, num_components);
    std::vector<float> tag_delta2(k_max_, 0.0f);

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_infinite)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      const int tag_index = deftag_indices[deftag_index];
      const float sig2_zero = sig2_zeroA + ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index] * nvec[tag_index] * sig2_zeroL;
      find_unified_univariate_tag_delta_sampling(num_components, pi_vec, sig2_vec, sig2_zeroC, tag_index, &nvec[0], &hvec[0], &tag_delta2, &subset_sampler, &ld_matrix_row);

      double pdf_tag = 0.0;
      double average_tag_delta2 = 0.0;
      for (int k = 0; k < k_max_; k++) {
        float s = sqrt(tag_delta2[k] + sig2_zero);
        const float tag_z = z_minus_fixed_effect_delta[tag_index];
        const bool censoring = std::abs(tag_z) > z_max;
        double pdf = static_cast<double>(censoring ? censored_cdf<FLOAT_TYPE>(z_max, s) : gaussian_pdf<FLOAT_TYPE>(tag_z, s));
        pdf_tag += pdf * pi_k;
        average_tag_delta2 += tag_delta2[k] * pi_k;
      }

      // export the expected values of z^2 distribution
      if ((aux != nullptr) && (aux_option_ == AuxOption_Ezvec2)) aux[tag_index] = average_tag_delta2 + sig2_zero;
      if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdf)) aux[tag_index] = pdf_tag;

      double increment = -std::log(pdf_tag) * static_cast<double>(weights[tag_index]);
      if (!std::isfinite(increment)) {
        increment = static_cast<double>(-std::log(kMinTagPdf) * static_cast<double>(weights[tag_index]));
        num_infinite++;
      }

      log_pdf_total += increment;
    }
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<" << ss.str() << ", cost=" << log_pdf_total << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

void BgmgCalculator::find_unified_univariate_tag_delta_sampling(int num_components, float* pi_vec, float* sig2_vec, float sig2_zeroC, int tag_index, const float* nvec, const float* hvec, std::vector<float>* tag_delta2, MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row) {
  if (!use_complete_tag_indices_) BGMG_THROW_EXCEPTION(::std::runtime_error("Unified Sampling calculator require 'use_complete_tag_indices' option"));
  tag_delta2->assign(k_max_, 0.0f);
  const int snp_index = tag_index; // yes, snp==tag in this case -- same story here as in calc_univariate_characteristic_function_times_cosinus function.
  ld_matrix_csr_.extract_row(snp_index, ld_matrix_row);
  auto iter_end = ld_matrix_row->end();

  float delta2_inf = 0.0;
  for (auto iter = ld_matrix_row->begin(); iter < iter_end; iter++) {
    const int causal_index = iter.tag_index();
    const float nval = nvec[tag_index];
    const float r2 = iter.r2();
    const float hval = hvec[causal_index];
    const float r2_hval_nval_sig2zeroC = (r2 * hval * nval * sig2_zeroC);

    for (int comp_index = 0; comp_index < num_components; comp_index++) {
      const int index = (comp_index*num_snp_ + causal_index);
      float pi_val = pi_vec[index];
      if (pi_val > 0.5f) pi_val = 1.0f - pi_val;
      subset_sampler->p()[comp_index] = static_cast<double>(pi_val);
    }

    int sample_global_index = k_max_ - subset_sampler->sample_shuffle();
    const uint32_t* indices = subset_sampler->data();
    for (int comp_index = 0; comp_index < num_components; comp_index++) {
      const int index = (comp_index*num_snp_ + causal_index);

      float pi_val = pi_vec[index];
      float delta2_val = r2_hval_nval_sig2zeroC * sig2_vec[index];
      if (pi_val > 0.5f) { // for pi_val close to 1.0 it'll be faster to compute total, and deduct selected (1-pi_val) samples at random
        pi_val = 1.0f - pi_val;
        delta2_inf += delta2_val;
        delta2_val *= -1;
      }

      const int num_samples=subset_sampler->counts()[comp_index];
      for (int sample_index=0; sample_index < num_samples; sample_index++, sample_global_index++) {
        tag_delta2->at(indices[sample_global_index]) += delta2_val;
      }
    }
    assert(sample_global_index==k_max_);
  }

  for (int k_index = 0; k_index < k_max_; k_index++) {
    float val = tag_delta2->at(k_index);
    val += delta2_inf;
    if (val < 0.0f) val=0.0f;
    tag_delta2->at(k_index) = val;
  }
}

void BgmgCalculator::find_unified_bivariate_tag_delta_sampling(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, int tag_index, const float* nvec1, const float* nvec2, const float* hvec, std::vector<float>* tag_delta20, std::vector<float>* tag_delta02, std::vector<float>* tag_delta11, MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row) {
  if (!use_complete_tag_indices_) BGMG_THROW_EXCEPTION(::std::runtime_error("Unified Sampling calculator require 'use_complete_tag_indices' option"));
  tag_delta20->assign(k_max_, 0.0f); tag_delta02->assign(k_max_, 0.0f); tag_delta11->assign(k_max_, 0.0f);
  const int snp_index = tag_index; // yes, snp==tag in this case -- same story here as in calc_univariate_characteristic_function_times_cosinus function.
  ld_matrix_csr_.extract_row(snp_index, ld_matrix_row);
  auto iter_end = ld_matrix_row->end();

  float delta20_inf = 0.0, delta02_inf = 0.0, delta11_inf = 0.0;
  for (auto iter = ld_matrix_row->begin(); iter < iter_end; iter++) {
    const int causal_index = iter.tag_index();
    const float nval1 = nvec1[tag_index];
    const float nval2 = nvec2[tag_index];
    const float r2 = iter.r2();
    const float hval = hvec[causal_index];

    const float r2_hval_nval1_sig2_zeroC = (r2 * hval * nval1 * sig2_zeroC[0]);
    const float r2_hval_nval2_sig2_zeroC = (r2 * hval * nval2 * sig2_zeroC[1]);

    const float sig2_beta1_val = sig2_vec[0*num_snp_ + causal_index];
    const float sig2_beta2_val = sig2_vec[1*num_snp_ + causal_index];
    const float sig2_beta1[3] = {sig2_beta1_val, 0, sig2_beta1_val};
    const float sig2_beta2[3] = {0, sig2_beta2_val, sig2_beta2_val};
    const float rho[3] = {0, 0, rho_vec[causal_index]};

    const int num_components = 3;
    for (int comp_index = 0; comp_index < num_components; comp_index++) {
      const int index = (comp_index*num_snp_ + causal_index);
      float pi_val = pi_vec[index];
      if (pi_val > 0.5f) pi_val = 1.0f - pi_val;
      subset_sampler->p()[comp_index] = static_cast<double>(pi_val);
    }

    int sample_global_index = k_max_ - subset_sampler->sample_shuffle();
    const uint32_t* indices = subset_sampler->data();
    for (int comp_index = 0; comp_index < num_components; comp_index++) {
      const int index = (comp_index*num_snp_ + causal_index);

      float pi_val = pi_vec[index];
      float delta20_val = r2_hval_nval1_sig2_zeroC * sig2_beta1[comp_index];
      float delta02_val = r2_hval_nval2_sig2_zeroC * sig2_beta2[comp_index];
      float delta11_val = rho[comp_index] * sqrt(delta20_val * delta02_val);

      if (pi_val > 0.5f) { // for pi_val close to 1.0 it'll be faster to compute total, and deduct selected (1-pi_val) samples at random
        pi_val = 1.0f - pi_val;
        delta20_inf += delta20_val; delta20_val *= -1;
        delta02_inf += delta02_val; delta02_val *= -1;
        delta11_inf += delta11_val; delta11_val *= -1;
      }

      const int num_samples=subset_sampler->counts()[comp_index];
      for (int sample_index=0; sample_index < num_samples; sample_index++, sample_global_index++) {
        tag_delta20->at(indices[sample_global_index]) += delta20_val;
        tag_delta02->at(indices[sample_global_index]) += delta02_val;
        tag_delta11->at(indices[sample_global_index]) += delta11_val;
      }
    }
    assert(sample_global_index==k_max_);    
  } 

  for (int k_index = 0; k_index < k_max_; k_index++) {
    tag_delta20->at(k_index) = std::max(0.0f, tag_delta20->at(k_index) + delta20_inf);
    tag_delta02->at(k_index) = std::max(0.0f, tag_delta02->at(k_index) + delta02_inf);
    tag_delta11->at(k_index) = std::max(0.0f, tag_delta11->at(k_index) + delta11_inf);
  }
}

int64_t BgmgCalculator::calc_unified_univariate_pdf(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, int length, float* zvec, float* pdf) {
  check_num_snp(num_snp);

  std::stringstream ss;
  ss << "calc_unified_univariate_pdf(trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL << ", length=" << length << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);

  const double pi_k = 1.0 / static_cast<double>(k_max_);
  std::vector<float>& nvec(*get_nvec(trait_index));
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);

  std::valarray<double> pdf_double(0.0, length);
#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, length);
    LdMatrixRow ld_matrix_row;
    MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + omp_get_thread_num(), k_max_, num_components);
    std::vector<float> tag_delta2(k_max_, 0.0f);

#pragma omp for schedule(static)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      find_unified_univariate_tag_delta_sampling(num_components, pi_vec, sig2_vec, sig2_zeroC, tag_index, &nvec[0], &hvec[0], &tag_delta2, &subset_sampler, &ld_matrix_row);
      const float sig2_zero = sig2_zeroA + ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index] * nvec[tag_index] * sig2_zeroL;
      const double tag_weight = static_cast<double>(weights_[tag_index]);

      for (int k_index = 0; k_index < k_max_; k_index++) {
        const float tag_delta2_value = tag_delta2[k_index];
        float s = sqrt(tag_delta2_value + sig2_zero);
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
  LOG << "<" << ss.str() << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::calc_unified_univariate_power(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float zthresh, int length, float* nvec, float* svec) {
  std::stringstream ss;
  ss << "calc_unified_univariate_power(trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL << ", zthresh=" << zthresh << ", length=" << length << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);
  const double pi_k = 1.0 / static_cast<double>(k_max_);
  std::vector<float> hvec; find_hvec(*this, &hvec);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> nvec_dummy(num_tag_, 1.0f);

  std::valarray<double> s_numerator_global(0.0, length);
  std::valarray<double> s_denominator_global(0.0, length);

#pragma omp parallel
  {
    std::valarray<double> s_numerator_local(0.0, length);
    std::valarray<double> s_denominator_local(0.0, length);
    LdMatrixRow ld_matrix_row;
    MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + omp_get_thread_num(), k_max_, num_components);
    std::vector<float> tag_delta2(k_max_, 0.0f);

#pragma omp for schedule(static)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      const double tag_weight = static_cast<double>(weights_[tag_index]);
      find_unified_univariate_tag_delta_sampling(num_components, pi_vec, sig2_vec, sig2_zeroC, tag_index, &nvec_dummy[0], &hvec[0], &tag_delta2, &subset_sampler, &ld_matrix_row);

      for (int k_index = 0; k_index < k_max_; k_index++) {
        for (int n_index = 0; n_index < length; n_index++) {
          float delta2eff = tag_delta2[k_index] * nvec[n_index] + ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index] * nvec[n_index] * sig2_zeroL;
          float sig2eff = delta2eff + sig2_zeroA;
          float sqrt_sig2eff = sqrt(sig2eff);
          static const float sqrt_2 = sqrtf(2.0);
          float numerator1 = gaussian_pdf<FLOAT_TYPE>(zthresh, sqrt_sig2eff) * 2 * delta2eff * delta2eff * zthresh / sig2eff;
          float numerator2 = std::erfcf(zthresh / (sqrt_2 * sqrt_sig2eff)) * delta2eff;
          s_numerator_local[n_index] += tag_weight*(numerator1 + numerator2);
          s_denominator_local[n_index] += tag_weight*delta2eff;
        }
      }
    }

#pragma omp critical
    {
      s_numerator_global += s_numerator_local;
      s_denominator_global += s_denominator_local;
    }
  }

  for (int i = 0; i < length; i++) svec[i] = static_cast<float>(s_numerator_global[i] / s_denominator_global[i]);
  LOG << "<" << ss.str() << ", elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

// c0 = c(0), c1=c(1), c2=c(2), where c(q) = \int_\delta \delta^q P(z|delta) P(delta)
// c(q) is define so that:
//  E(\delta^2|z_j) = c2[j]/c0[j];
//  E(\delta  |z_j) = c1[j]/c0[j];
int64_t BgmgCalculator::calc_unified_univariate_delta_posterior(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, int length, float* c0, float* c1, float* c2) {
  std::stringstream ss;
  ss << "calc_unified_univariate_delta_posterior(trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL << ", length=" << length << ")";
  LOG << ">" << ss.str();

  if ((length == 0) || (length != num_tag_)) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_tag_"));

  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);
  std::vector<float> hvec; find_hvec(*this, &hvec);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();

  std::valarray<double> c0_global(0.0f, num_tag_);
  std::valarray<double> c1_global(0.0f, num_tag_);
  std::valarray<double> c2_global(0.0f, num_tag_);

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + omp_get_thread_num(), k_max_, num_components);
    std::vector<float> tag_delta2(k_max_, 0.0f);
    std::valarray<double> c0_local(0.0f, num_tag_);
    std::valarray<double> c1_local(0.0f, num_tag_);
    std::valarray<double> c2_local(0.0f, num_tag_);

#pragma omp for schedule(static)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      find_unified_univariate_tag_delta_sampling(num_components, pi_vec, sig2_vec, sig2_zeroC, tag_index, &nvec[0], &hvec[0], &tag_delta2, &subset_sampler, &ld_matrix_row);
    
      for (int k_index = 0; k_index < k_max_; k_index++) {

        const float delta2eff = tag_delta2[k_index] + ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index] * nvec[tag_index] * sig2_zeroL;  // S^2_kj
        const float sig2eff = delta2eff + sig2_zeroA;
        const float sig2eff_1_2 = sqrt(sig2eff);
        const float sig2eff_3_2 = sig2eff_1_2 * sig2eff;
        const float sig2eff_5_2 = sig2eff_3_2 * sig2eff;

        const float z = z_minus_fixed_effect_delta[tag_index];
        const float exp_common = std::exp(-0.5f*z*z / sig2eff);

        c0_local[tag_index] += (exp_common / sig2eff_1_2);
        c1_local[tag_index] += (exp_common / sig2eff_3_2) * z * delta2eff;
        c2_local[tag_index] += (exp_common / sig2eff_5_2) *     delta2eff * (sig2_zeroA*sig2_zeroA + sig2_zeroA*delta2eff + z*z*delta2eff);
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

  LOG << "<" << ss.str() << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

// pi_vec     : num_components X num_snp,  - weights of the three mixture components (num_components = 3)
// sig2_vec   : num_traits X num_snp,      - variance of cauasal effects for the two traits (num_traits = 2)
// rho_vec    : 1 x num_snps,              - correlation of genetic effects
// sig2_zeroX : 1 x num_traits,            - inflation parameters (additive, multiplicative, truncation of the LD structure)
// rho_zeroA, sig2_zeroL                   - correlation of sig2_zeroA and sig2_zeroL across the two traits
double BgmgCalculator::calc_unified_bivariate_cost(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux) { 
  check_num_snp(num_snp);
  if (cost_calculator_ == CostCalculator_Gaussian) return calc_unified_bivariate_cost_gaussian(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, aux);
  else if (cost_calculator_ == CostCalculator_Sampling) return calc_unified_bivariate_cost_sampling(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, aux, nullptr);
  else BGMG_THROW_EXCEPTION(::std::runtime_error("unsupported cost calculator in calc_unified_bivariate_cost"));
  return 0;
}

std::string find_bivariate_params_description(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL) {
  std::stringstream ss;
  ss << "num_snp=" << num_snp
     << ", sig2_zeroA=[" << sig2_zeroA[0] << "," << sig2_zeroA[1]  << "]"
     << ", sig2_zeroC=[" << sig2_zeroC[0] << "," << sig2_zeroC[1]  << "]"
     << ", sig2_zeroL=[" << sig2_zeroL[0] << "," << sig2_zeroL[1]  << "]"
     << ", rho_zeroA=" << rho_zeroA << ", rho_zeroL=" << rho_zeroL;
  return ss.str();
}

double BgmgCalculator::calc_unified_bivariate_cost_gaussian(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux) {
  std::stringstream ss;
  ss << "calc_unified_bivariate_cost_gaussian(" << find_bivariate_params_description(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL) << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);

  // Step 1. Calculate Ebeta20, Ebeta02, Ebeta11
  std::valarray<float> Ebeta20(0.0, num_snp_);
  std::valarray<float> Ebeta02(0.0, num_snp_);
  std::valarray<float> Ebeta11(0.0, num_snp_);
  for (int snp_index = 0; snp_index < num_snp_; snp_index++) {
    const float p1 = pi_vec[0*num_snp_ + snp_index];
    const float p2 = pi_vec[1*num_snp_ + snp_index];
    const float p12 = pi_vec[2*num_snp_ + snp_index];

    const float s1 = sig2_vec[0*num_snp_ + snp_index];
    const float s2 = sig2_vec[1*num_snp_ + snp_index];

    const float rho = rho_vec[snp_index];
    
    Ebeta20[snp_index] = (p1 + p12) * s1;
    Ebeta02[snp_index] = (p2 + p12) * s2;
    Ebeta11[snp_index] = p12 * rho * sqrt(s1*s2);
  }

  // Step 2. Calculate Edelta20, Edelta02, Edelta11
  std::valarray<float> Edelta20(0.0, num_snp_);
  std::valarray<float> Edelta02(0.0, num_snp_);
  std::valarray<float> Edelta11(0.0, num_snp_);
  
#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    std::valarray<float> Edelta20_local(0.0, num_snp_);
    std::valarray<float> Edelta02_local(0.0, num_snp_);
    std::valarray<float> Edelta11_local(0.0, num_snp_);

    if (use_complete_tag_indices_) {  // see discussion in calc_unified_univariate_cost_gaussian
#pragma omp for schedule(static)    
      for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
        int tag_index = deftag_indices[deftag_index];

        const int snp_index = tag_index; // yes, snp==tag in this case -- same story here as in calc_univariate_characteristic_function_times_cosinus function.
        ld_matrix_csr_.extract_row(snp_index, &ld_matrix_row);
        auto iter_end = ld_matrix_row.end();
        for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
          const int causal_index = iter.tag_index();   // tag_index on RHS can be misleading here
          const float r2_value = iter.r2();
          const float a2ij1 = sig2_zeroC[0] * nvec1_[tag_index] * hvec[causal_index] * r2_value;
          const float a2ij2 = sig2_zeroC[1] * nvec2_[tag_index] * hvec[causal_index] * r2_value;
          Edelta20_local[tag_index] += a2ij1 * Ebeta20[causal_index];
          Edelta02_local[tag_index] += a2ij2 * Ebeta02[causal_index];
          Edelta11_local[tag_index] += sqrt(a2ij1 * a2ij2) * Ebeta11[causal_index];
        }
      }
    } else {  // use_complete_tag_indices_
#pragma omp for schedule(static)        
      for (int causal_index = 0; causal_index < num_snp_; causal_index++) {
        ld_matrix_csr_.extract_row(causal_index, &ld_matrix_row);
        auto iter_end = ld_matrix_row.end();
        for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
          const int tag_index = iter.tag_index();
          const float r2_value = iter.r2();
          const float a2ij1 = sig2_zeroC[0] * nvec1_[tag_index] * hvec[causal_index] * r2_value;
          const float a2ij2 = sig2_zeroC[1] * nvec2_[tag_index] * hvec[causal_index] * r2_value;
          Edelta20_local[tag_index] += a2ij1 * Ebeta20[causal_index];
          Edelta02_local[tag_index] += a2ij2 * Ebeta02[causal_index];
          Edelta11_local[tag_index] += sqrt(a2ij1 * a2ij2) * Ebeta11[causal_index];
        }
      }
    }  // use_complete_tag_indices_
#pragma omp critical
    {
      Edelta20 += Edelta20_local;
      Edelta02 += Edelta02_local;
      Edelta11 += Edelta11_local;
    }
  }

  double log_pdf_total;
  int num_zero_tag_r2 = 0;
  int num_infinite = 0;

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_zero_tag_r2, num_infinite)
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];
    if (Edelta20[tag_index] == 0 && Edelta02[tag_index] == 0) { num_zero_tag_r2++; continue;}
    double tag_weight = static_cast<double>(weights_[tag_index]);

    // the indices (20, 02, 11) denote powers "p, q" in raw moment E[x^p y^q]
    const float A20 = Edelta20[tag_index];
    const float A02 = Edelta02[tag_index];
    const float A11 = Edelta11[tag_index];

    // additive inflation, plus contribution from small LD r2 (those below r2min)
    // the indices (11, 12, 22) denote indices in a matrix [a11 a12; a21 a22]. Usually a21=a12 so we don't compute it.
    const float adj_hval = ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index];
    const float sig2_zero_11 = sig2_zeroA[0] + adj_hval * nvec1_[tag_index] * sig2_zeroL[0];
    const float sig2_zero_22 = sig2_zeroA[1] + adj_hval * nvec2_[tag_index] * sig2_zeroL[1];
    const float sig2_zero_12 =            rho_zeroA * sqrt(sig2_zeroA[0] * sig2_zeroA[1]) + 
                               adj_hval * rho_zeroL * sqrt(nvec1_[tag_index] * nvec2_[tag_index] * sig2_zeroL[0] * sig2_zeroL[1]);

    const float tag_z1 = z1_minus_fixed_effect_delta[tag_index];
    const float tag_z2 = z2_minus_fixed_effect_delta[tag_index];
    const float tag_n1 = nvec1_[tag_index];
    const float tag_n2 = nvec2_[tag_index];

    const bool censoring = (std::abs(tag_z1) > z1max_) || (std::abs(tag_z2) > z2max_);

    const float a11 = Edelta20[tag_index] + sig2_zero_11;
    const float a12 = Edelta11[tag_index] + sig2_zero_12;
    const float a22 = Edelta02[tag_index] + sig2_zero_22;

    // export the expected values of z^2 distribution
    if ((aux != nullptr) && (aux_option_ == AuxOption_Ezvec2)) {
      aux[0 * num_tag_ + tag_index] = a11;
      aux[1 * num_tag_ + tag_index] = a12;
      aux[2 * num_tag_ + tag_index] = a22;
    }

    const double tag_pdf = static_cast<double>(censoring ? censored2_cdf<FLOAT_TYPE>(z1max_, z2max_, a11, a12, a22) : gaussian2_pdf<FLOAT_TYPE>(tag_z1, tag_z2, a11, a12, a22));

    if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdf)) aux[tag_index] = tag_pdf;

    double increment = (-std::log(tag_pdf) * tag_weight);
    if (!std::isfinite(increment)) {
      increment = static_cast<double>(-std::log(kMinTagPdf) * tag_weight);
      num_infinite++;
    }

    log_pdf_total += increment;
  }

  if (num_zero_tag_r2 > 0)
    LOG << " warning: zero tag_r2 encountered " << num_zero_tag_r2 << " times";
  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";


  LOG << "<" << ss.str() << ", cost=" << log_pdf_total << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

class BivariateCharacteristicFunctionData {
 public:
  const float* pi_vec;
  const float* sig2_zero;
  const float* sig2_beta;
  float rho_beta;
  float rho_zero;

  int tag_index;  // for which SNP to calculate the characteristic function
                  // note that use_complete_tag_indices_ must be enabled for "convolve" calculator, 
                  // so "tag" and "snp" are in the same indexing
  LdMatrixRow* ld_matrix_row;

  const std::vector<float>* hvec;
  const std::vector<float>* z1_minus_fixed_effect_delta;
  const std::vector<float>* nvec1;
  const std::vector<float>* z2_minus_fixed_effect_delta;
  const std::vector<float>* nvec2;
  const std::vector<float>* ld_tag_sum_r2_below_r2min_adjust_for_hvec;
  int func_evals;
};

inline float exp_quad_form(float t1, float t2, float a11, float a12, float a22) {
  return std::exp(-0.5f*(t1*t1 * a11 + 2.0*t1*t2*a12 + t2*t2*a22));
}

int calc_bivariate_characteristic_function_times_cosinus(unsigned ndim, const double *x, void *raw_data, unsigned fdim, double* fval) {
  assert(ndim == 2);
  assert(fdim == 1);
  const float t1 = (float)x[0];
  const float t2 = (float)x[1];
  BivariateCharacteristicFunctionData* data = (BivariateCharacteristicFunctionData *)raw_data;
  const float pi1 = data->pi_vec[0];
  const float pi2 = data->pi_vec[1];
  const float pi12 = data->pi_vec[2];
  const float pi0 = 1.0 - pi1 - pi2 - pi12;
  const float eff_sig2_beta1 = (*data->nvec1)[data->tag_index] * data->sig2_beta[0];
  const float eff_sig2_beta2 = (*data->nvec2)[data->tag_index] * data->sig2_beta[1];
  const float eff_sig2_beta_cov = data->rho_beta * sqrt(eff_sig2_beta1 * eff_sig2_beta2);
  const float sig2_zero_cov = data->rho_zero * sqrt(data->sig2_zero[0] * data->sig2_zero[1]);
  const float zval1 = (*data->z1_minus_fixed_effect_delta)[data->tag_index];
  const float zval2 = (*data->z2_minus_fixed_effect_delta)[data->tag_index];
  const float inf_adj_r2 = (*data->ld_tag_sum_r2_below_r2min_adjust_for_hvec)[data->tag_index];
  const float scaling_factor = 0.5 * M_1_PI * M_1_PI;

  double result = scaling_factor * cos(t1 * zval1 + t2 * zval2);
  result *= exp_quad_form(t1, t2, data->sig2_zero[0], sig2_zero_cov, data->sig2_zero[1]);

  // apply infinitesimal model to adjust tag_r2sum for all r2 that are below r2min (and thus do not contribute via resampling)
  result *= exp_quad_form(t1, t2, pi1 * eff_sig2_beta1 * inf_adj_r2, 0, 0);
  result *= exp_quad_form(t1, t2, 0, 0, pi2 * eff_sig2_beta2 * inf_adj_r2);
  result *= exp_quad_form(t1, t2, pi12 * eff_sig2_beta1 * inf_adj_r2,
                                  pi12 * eff_sig2_beta_cov * inf_adj_r2,
                                  pi12 * eff_sig2_beta2 * inf_adj_r2);
  
  auto iter_end = data->ld_matrix_row->end();
  for (auto iter = data->ld_matrix_row->begin(); iter < iter_end; iter++) {
    int snp_index = iter.tag_index();  // yes, causal==tag this is correct - snp_index on the LHS, tag_index on the RHS. 
                                       // See comment in the univariate counterpart, calc_univariate_characteristic_function_times_cosinus
    const float r2 = iter.r2();
    const float hval = (*data->hvec)[snp_index];
    const float r2_times_hval = r2 * hval;
    result *= (double) (pi0 + 
                        pi1 * exp_quad_form(t1, t2, eff_sig2_beta1 * r2_times_hval, 0, 0) +
                        pi2 * exp_quad_form(t1, t2, 0, 0, eff_sig2_beta2 * r2_times_hval) + 
                        pi12* exp_quad_form(t1, t2, eff_sig2_beta1 * r2_times_hval,
                                                    eff_sig2_beta_cov * r2_times_hval,
                                                    eff_sig2_beta2 * r2_times_hval));
  }

  data->func_evals++;
  *fval = result;
  return 0; 
}

int calc_bivariate_characteristic_function_for_integration(unsigned ndim, const double *x, void *raw_data, unsigned fdim, double* fval) {
  double x_transformed[2];

  const double t1 = x[0];
  const double inv_1_minus_t1 = 1.0 / (1.0 - t1);
  x_transformed[0] = t1 * inv_1_minus_t1;
  const double jacob1 = inv_1_minus_t1 * inv_1_minus_t1;

  const double t2 = x[1];
  const double inv_1_minus_t2sqr = 1.0 / (1.0 - t2*t2);
  x_transformed[1] = t2 * inv_1_minus_t2sqr;
  const double jacob2 = (1 + t2*t2) * inv_1_minus_t2sqr * inv_1_minus_t2sqr;
  
  int retval = calc_bivariate_characteristic_function_times_cosinus(ndim, x_transformed, raw_data, fdim, fval);
  (*fval) *= (jacob1 * jacob2);
  return retval;
}

double BgmgCalculator::calc_bivariate_cost_convolve(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero) {
  if (!use_complete_tag_indices_) BGMG_THROW_EXCEPTION(::std::runtime_error("Convolve calculator require 'use_complete_tag_indices' option"));

  std::string ss = ""; // calc_bivariate_params_to_str(pi_vec_len, pi_vec, sig2_beta_len, sig2_beta, rho_beta, sig2_zero_len, sig2_zero, rho_zero, -1);
  LOG << ">calc_bivariate_cost_convolve(" << ss << ")";

  double log_pdf_total = 0.0;
  int num_snp_failed = 0;
  int num_infinite = 0;
  double func_evals = 0.0;
  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    BivariateCharacteristicFunctionData data;
    data.pi_vec = pi_vec;
    data.sig2_zero = sig2_zero;
    data.sig2_beta = sig2_beta;
    data.rho_zero = rho_zero;
    data.rho_beta = rho_beta;
    data.ld_matrix_row = &ld_matrix_row;
    data.hvec = &hvec;
    data.z1_minus_fixed_effect_delta = &z1_minus_fixed_effect_delta;
    data.nvec1 = &nvec1_;
    data.z2_minus_fixed_effect_delta = &z2_minus_fixed_effect_delta;
    data.nvec2 = &nvec2_;
    data.ld_tag_sum_r2_below_r2min_adjust_for_hvec = &ld_tag_sum_r2_below_r2min_adjust_for_hvec;

#pragma omp for schedule(static) reduction(+: log_pdf_total, num_snp_failed, num_infinite, func_evals)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      double tag_weight = static_cast<double>(weights_[tag_index]);

      const int causal_index = tag_index; // yes, causal==tag in this case -- see a large comment in calc_univariate_characteristic_function_times_cosinus function.
      ld_matrix_csr_.extract_row(causal_index, data.ld_matrix_row);
      data.tag_index = tag_index;
      data.func_evals = 0;

      double tag_pdf = 0, tag_pdf_err = 0;
      const double xmin[2] = {0.0, -1.0}, xmax[2] = {1.0, 1.0};
      const int integrand_fdim = 1, ndim = 2;
      int cubature_result = hcubature(integrand_fdim, calc_bivariate_characteristic_function_for_integration,
        &data, ndim, xmin, xmax, cubature_max_evals_, cubature_abs_error_, cubature_rel_error_, ERROR_INDIVIDUAL, &tag_pdf, &tag_pdf_err);
      func_evals += (weights_[tag_index] * (double)data.func_evals);
      if (cubature_result != 0) { num_snp_failed++; continue; }

      if (tag_pdf <= 0)
        tag_pdf = 1e-100;

      double increment = static_cast<double>(-std::log(tag_pdf) * weights_[tag_index]);
      if (!std::isfinite(increment)) num_infinite++;
      else log_pdf_total += increment;
    }
  }    

  if (num_snp_failed > 0)
    LOG << " warning: hcubature failed for " << num_snp_failed << " tag snps";
  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  double total_weight = 0.0;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) total_weight += weights_[tag_index];
  func_evals /= total_weight;

  LOG << "<calc_bivariate_cost_convolve(" << ss << "), cost=" << log_pdf_total << ", evals=" << func_evals << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_unified_bivariate_cost_sampling(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux, const float* weights) {
  if (!use_complete_tag_indices_) BGMG_THROW_EXCEPTION(::std::runtime_error("Unified Sampling calculator require 'use_complete_tag_indices' option"));
  if (weights == nullptr) {
    if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
    weights = &weights_[0];
  }

  std::stringstream ss;
  ss << "calc_unified_bivariate_cost_sampling(" << find_bivariate_params_description(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL) << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);
  
  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_tag_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);

  const double pi_k = 1.0 / static_cast<double>(k_max_);
  double log_pdf_total = 0.0;
  int num_infinite = 0;
  const int num_components = 3;

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + omp_get_thread_num(), k_max_, num_components);
    std::vector<float> tag_delta20(k_max_, 0.0f);
    std::vector<float> tag_delta02(k_max_, 0.0f);
    std::vector<float> tag_delta11(k_max_, 0.0f);

#pragma omp parallel for schedule(static) reduction(+: log_pdf_total, num_infinite)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      const int tag_index = deftag_indices[deftag_index];
      const float adj_hval = ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index];
      const float sig2_zero_11 = sig2_zeroA[0] + adj_hval * nvec1_[tag_index] * sig2_zeroL[0];
      const float sig2_zero_22 = sig2_zeroA[1] + adj_hval * nvec2_[tag_index] * sig2_zeroL[1];
      const float sig2_zero_12 =            rho_zeroA * sqrt(sig2_zeroA[0] * sig2_zeroA[1]) + 
                                 adj_hval * rho_zeroL * sqrt(nvec1_[tag_index] * nvec2_[tag_index] * sig2_zeroL[0] * sig2_zeroL[1]);

      const float tag_z1 = z1_minus_fixed_effect_delta[tag_index];
      const float tag_z2 = z2_minus_fixed_effect_delta[tag_index];
      const float tag_n1 = nvec1_[tag_index];
      const float tag_n2 = nvec2_[tag_index];

      const bool censoring = (std::abs(tag_z1) > z1max_) || (std::abs(tag_z2) > z2max_);

      find_unified_bivariate_tag_delta_sampling(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, tag_index, &nvec1_[0], &nvec2_[0], &hvec[0], &tag_delta20, &tag_delta02, &tag_delta11, &subset_sampler, &ld_matrix_row);

      double pdf_tag = 0.0;
      double average_tag_delta20 = 0.0, average_tag_delta02 = 0.0, average_tag_delta11 = 0.0;
      for (int k = 0; k < k_max_; k++) {
        const float a11 = tag_delta20[k] + sig2_zero_11;
        const float a12 = tag_delta11[k] + sig2_zero_12;
        const float a22 = tag_delta02[k] + sig2_zero_22;
        const double pdf = static_cast<double>(censoring ? censored2_cdf<FLOAT_TYPE>(z1max_, z2max_, a11, a12, a22) : gaussian2_pdf<FLOAT_TYPE>(tag_z1, tag_z2, a11, a12, a22));
        pdf_tag += pdf * pi_k;
        average_tag_delta20 += a11 * pi_k;
        average_tag_delta11 += a12 * pi_k;
        average_tag_delta02 += a22 * pi_k;
      }

      // export the expected values of z^2 distribution
      if ((aux != nullptr) && (aux_option_ == AuxOption_Ezvec2)) {
        aux[0 * num_tag_ + tag_index] = average_tag_delta20;
        aux[1 * num_tag_ + tag_index] = average_tag_delta11;
        aux[2 * num_tag_ + tag_index] = average_tag_delta02;
      }
      if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdf)) aux[tag_index] = pdf_tag;

      double increment = -std::log(pdf_tag) * static_cast<double>(weights[tag_index]);
      if (!std::isfinite(increment)) {
        increment = static_cast<double>(-std::log(kMinTagPdf) * static_cast<double>(weights[tag_index]));
        num_infinite++;
      }

      log_pdf_total += increment;
    }
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<" << ss.str() << ", cost=" << log_pdf_total << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

int64_t BgmgCalculator::calc_unified_bivariate_pdf(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, int length, float* zvec1, float* zvec2, float* pdf) {
  return 0;
}

int64_t BgmgCalculator::calc_unified_bivariate_delta_posterior(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL,
                                                               int length, float* c00, float* c10, float* c01, float* c20, float* c11, float* c02) {
  return 0;
}
