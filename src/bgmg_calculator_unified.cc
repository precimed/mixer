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
static const int kOmpDynamicChunk = 512;

double BgmgCalculator::calc_unified_univariate_cost(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux) {
  check_num_snp(num_snp);

  if (cost_calculator_ == CostCalculator_Gaussian) return calc_unified_univariate_cost_gaussian(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, aux);
  else if (cost_calculator_ == CostCalculator_Convolve) return calc_unified_univariate_cost_convolve(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, aux);
  else if (cost_calculator_ == CostCalculator_Sampling) return calc_unified_univariate_cost_sampling(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, aux, nullptr);
  else if (cost_calculator_ == CostCalculator_Smplfast) return calc_unified_univariate_cost_smplfast(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, aux, nullptr);
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
    int snp_index = iter.index();

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
  std::stringstream ss;
  ss << "trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL;
  LOG << ">calc_unified_univariate_cost_convolve(" << ss.str() << ")";

  double log_pdf_total = 0.0;
  int num_snp_failed = 0;
  int num_infinite = 0;
  double func_evals = 0.0;
  double total_weight = 0.0;
  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*get_nvec(trait_index));
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
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

#pragma omp for schedule(dynamic, kOmpDynamicChunk) reduction(+: log_pdf_total, num_snp_failed, num_infinite, func_evals, total_weight)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      double tag_weight = static_cast<double>(weights_convolve[tag_index]);

      ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), data.ld_matrix_row);
      data.tag_index = tag_index;
      data.func_evals = 0;

      double tag_pdf = 0, tag_pdf_err = 0;
      const double xmin = 0, xmax = 1;
      const int integrand_fdim = 1, ndim = 1;
      int cubature_result = hcubature(integrand_fdim, calc_univariate_characteristic_function_for_integration,
        &data, ndim, &xmin, &xmax, cubature_max_evals_, cubature_abs_error_, cubature_rel_error_, ERROR_INDIVIDUAL, &tag_pdf, &tag_pdf_err);
      func_evals += (weights_convolve[tag_index] * (double)data.func_evals);
      total_weight += weights_convolve[tag_index];
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

  LOG << "<calc_unified_univariate_cost_convolve(" << ss.str() << "), cost=" << log_pdf_total << ", evals=" << func_evals / total_weight << ", num_deftag=" << num_deftag << "+" << num_deftag_sampling << ", elapsed time " << timer.elapsed_ms() << "ms";

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
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
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

#pragma omp for schedule(dynamic, kOmpDynamicChunk)    
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];

      ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), &ld_matrix_row);
      auto iter_end = ld_matrix_row.end();
      for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
        const int snp_index = iter.index();
        const float r2_value = iter.r2();
        const float a2ij = sig2_zeroC * nvec[tag_index] * hvec[snp_index] * r2_value;
        Edelta2_local[tag_index] += a2ij *        Ebeta2[snp_index];
        Edelta4_local[tag_index] += a2ij * a2ij * Ebeta4[snp_index];
      }
    }

#pragma omp critical
    {
      Edelta2 += Edelta2_local;
      Edelta4 += Edelta4_local;
    }
  }  // parallel

  double log_pdf_total = 0.0;
  int num_zero_tag_r2 = 0;
  int num_infinite = 0;

#pragma omp parallel for schedule(dynamic, kOmpDynamicChunk) reduction(+: log_pdf_total, num_zero_tag_r2, num_infinite)
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
  if (weights == nullptr) {
    if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
    weights = &weights_[0];
  }

  std::stringstream ss;
  ss << "calc_unified_univariate_cost_sampling(trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL << ", k_max_=" << k_max_ << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*get_nvec(trait_index));
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(weights, &deftag_indices);

  const double z_max = (trait_index==1) ? z1max_ : z2max_;
  const double pi_k = 1.0 / static_cast<double>(k_max_);
  double log_pdf_total = 0.0;
  int num_infinite = 0;

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    std::vector<float> tag_delta2(k_max_, 0.0f);
    std::valarray<double> tag_kpdf(0.0f, k_max_);

#pragma omp for schedule(dynamic, kOmpDynamicChunk) reduction(+: log_pdf_total, num_infinite)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      const int tag_index = deftag_indices[deftag_index];
      
      // for those who are woundering what's the point of tag_to_snp_[tag_index]...
      // each tag SNPs should have its own sequence of random values, and we control this by second seed
      // however, in unit-tests we validate that use_complete_tag_indices doesn't change anything => it's best to parametrize each tag variant by its snp index
      MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + tag_to_snp_[tag_index], k_max_, num_components);
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
        tag_kpdf[k] = pdf;
      }

      // export the expected values of z^2 distribution
      if ((aux != nullptr) && (aux_option_ == AuxOption_Ezvec2)) aux[tag_index] = average_tag_delta2 + sig2_zero;
      if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdf)) aux[tag_index] = pdf_tag;
      if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdfErr)) aux[tag_index] = sqrt(((tag_kpdf - pdf_tag) * (tag_kpdf - pdf_tag)).sum() / (static_cast<double>(k_max_) - 1.0));

      double increment = -std::log(pdf_tag) * static_cast<double>(weights[tag_index]);
      if (!std::isfinite(increment)) {
        increment = static_cast<double>(-std::log(kMinTagPdf) * static_cast<double>(weights[tag_index]));
        num_infinite++;
      }

#pragma omp critical
      {
        log_pdf_total += increment;
      }
    }
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<" << ss.str() << ", cost=" << log_pdf_total << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_unified_univariate_cost_smplfast(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux, const float* weights) {
  if (weights == nullptr) {
    if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
    weights = &weights_[0];
  }

  std::stringstream ss;
  ss << "calc_unified_univariate_cost_smplfast(trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL << ", k_max_=" << k_max_ << ")";
  LOG << ">" << ss.str();

  for (int comp_index = 0; comp_index < num_components; comp_index++) {
    const float ref_pi = pi_vec[comp_index*num_snp_];
    if (ref_pi > 0.5) LOG << " warning: pi_vec has values above 0.5, which is not optimized in calc_unified_univariate_cost_smplfast; consider using calc_unified_univariate_cost_sampling instead";
    for (int snp_index = 1; snp_index < num_snp; snp_index++) {
      if (ref_pi != pi_vec[comp_index*num_snp + snp_index])
        BGMG_THROW_EXCEPTION(::std::runtime_error("pi_vec doesn't appear to be constant across SNPs; use calc_unified_univariate_cost_sampling instead"));
    }
  }

  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*get_nvec(trait_index));
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(weights, &deftag_indices);

  const double z_max = (trait_index==1) ? z1max_ : z2max_;
  const double pi_k = 1.0 / static_cast<double>(k_max_);

  std::valarray<double> pdf_double(0.0, num_tag_);
  std::valarray<double> aux_Ezvec2(0.0, num_tag_);

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    std::vector<float> tag_delta2(num_tag_, 0.0f);
    std::valarray<double> pdf_double_local(0.0, num_tag_);
    std::valarray<double> aux_Ezvec2_local(0.0, num_tag_);

#pragma omp for schedule(dynamic)
    for (int k_index = 0; k_index < k_max_; k_index++) {
      MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + k_index, num_snp_, num_components);
      find_unified_univariate_tag_delta_smplfast(num_components, pi_vec, sig2_vec, sig2_zeroC, k_index, &nvec[0], &hvec[0], &tag_delta2, &subset_sampler, &ld_matrix_row);
      for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
        const int tag_index = deftag_indices[deftag_index];
        const float sig2_zero = sig2_zeroA + ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index] * nvec[tag_index] * sig2_zeroL;

        float s = sqrt(tag_delta2[tag_index] + sig2_zero);
        const float tag_z = z_minus_fixed_effect_delta[tag_index];
        const bool censoring = std::abs(tag_z) > z_max;
        double pdf = static_cast<double>(censoring ? censored_cdf<FLOAT_TYPE>(z_max, s) : gaussian_pdf<FLOAT_TYPE>(tag_z, s));
        pdf_double_local[tag_index] += pdf * pi_k;
        aux_Ezvec2_local[tag_index] += (tag_delta2[tag_index] + sig2_zero) * pi_k;
      }
    }

#pragma omp critical
    {
      pdf_double += pdf_double_local;
      aux_Ezvec2 += aux_Ezvec2_local;
    }
  }

  double log_pdf_total = 0.0;
  int num_infinite = 0;
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];
    double increment = -std::log(pdf_double[tag_index]) * static_cast<double>(weights[tag_index]);
    if (!std::isfinite(increment)) {
      increment = static_cast<double>(-std::log(kMinTagPdf) * static_cast<double>(weights[tag_index]));
      num_infinite++;
    }
    log_pdf_total += increment;

    // export the expected values of z^2 distribution
    if ((aux != nullptr) && (aux_option_ == AuxOption_Ezvec2)) aux[tag_index] = aux_Ezvec2[tag_index];
    if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdf)) aux[tag_index] = pdf_double[tag_index];
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<" << ss.str() << ", cost=" << log_pdf_total << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

void BgmgCalculator::find_unified_univariate_tag_delta_sampling(int num_components, float* pi_vec, float* sig2_vec, float sig2_zeroC, int tag_index, const float* nvec, const float* hvec, std::vector<float>* tag_delta2, MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row) {
  tag_delta2->assign(k_max_, 0.0f);
  ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), ld_matrix_row);
  auto iter_end = ld_matrix_row->end();

  float delta2_inf = 0.0;
  for (auto iter = ld_matrix_row->begin(); iter < iter_end; iter++) {
    const int snp_index = iter.index();
    const float nval = nvec[tag_index];
    const float r2 = iter.r2();
    const float hval = hvec[snp_index];
    const float r2_hval_nval_sig2zeroC = (r2 * hval * nval * sig2_zeroC);

    for (int comp_index = 0; comp_index < num_components; comp_index++) {
      const int index = (comp_index*num_snp_ + snp_index);
      float pi_val = pi_vec[index];
      if (pi_val > 0.5f) pi_val = 1.0f - pi_val;
      subset_sampler->p()[comp_index] = static_cast<double>(pi_val);
    }

    int sample_global_index = k_max_ - subset_sampler->sample_shuffle();
    const uint32_t* indices = subset_sampler->data();
    for (int comp_index = 0; comp_index < num_components; comp_index++) {
      const int index = (comp_index*num_snp_ + snp_index);

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

void BgmgCalculator::find_unified_univariate_tag_delta_smplfast(int num_components, float* pi_vec, float* sig2_vec, float sig2_zeroC, int k_index, const float* nvec, const float* hvec, std::vector<float>* tag_delta2, MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row) {
  tag_delta2->assign(num_tag_, 0.0f);

  const int global_snp_pi_index = 0; // here in smplfast we assume that all SNPs share the same pi
  for (int comp_index = 0; comp_index < num_components; comp_index++) {
    const int index = (comp_index*num_snp_ + global_snp_pi_index);
    const float pi_val = pi_vec[index];
    subset_sampler->p()[comp_index] = static_cast<double>(pi_val);
  }

  int sample_global_index = num_snp_ - subset_sampler->sample_shuffle();
  const uint32_t* indices = subset_sampler->data();
  for (int comp_index = 0; comp_index < num_components; comp_index++) {
    const int num_samples=subset_sampler->counts()[comp_index];
    for (int sample_index=0; sample_index < num_samples; sample_index++, sample_global_index++) {
      const int snp_index = indices[sample_global_index];
      ld_matrix_csr_.extract_snp_row(SnpIndex(snp_index), ld_matrix_row);
      auto iter_end = ld_matrix_row->end();

      for (auto iter = ld_matrix_row->begin(); iter < iter_end; iter++) {
        const int tag_index = iter.index();
        const float nval_times_sig2_zeroC = nvec[tag_index] * sig2_zeroC; // can be moved to the caller of find_unified_univariate_tag_delta_smplfast
        const float hval = hvec[snp_index];
        const float r2 = iter.r2();
        const float delta2_val = r2 * hval * nval_times_sig2_zeroC * sig2_vec[comp_index*num_snp_ + snp_index];
        tag_delta2->at(tag_index) += delta2_val;
      }
    }
  }
  assert(sample_global_index==num_snp_);
}

void BgmgCalculator::find_unified_bivariate_tag_delta_sampling(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, int tag_index, const float* nvec1, const float* nvec2, const float* hvec, std::vector<float>* tag_delta20, std::vector<float>* tag_delta02, std::vector<float>* tag_delta11, MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row) {
  tag_delta20->assign(k_max_, 0.0f); tag_delta02->assign(k_max_, 0.0f); tag_delta11->assign(k_max_, 0.0f);
  ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), ld_matrix_row);
  auto iter_end = ld_matrix_row->end();

  float delta20_inf = 0.0, delta02_inf = 0.0, delta11_inf = 0.0;
  for (auto iter = ld_matrix_row->begin(); iter < iter_end; iter++) {
    const int snp_index = iter.index();
    const float nval1 = nvec1[tag_index];
    const float nval2 = nvec2[tag_index];
    const float r2 = iter.r2();
    const float hval = hvec[snp_index];

    const float r2_hval_nval1_sig2_zeroC = (r2 * hval * nval1 * sig2_zeroC[0]);
    const float r2_hval_nval2_sig2_zeroC = (r2 * hval * nval2 * sig2_zeroC[1]);

    const float sig2_beta1_val = sig2_vec[0*num_snp_ + snp_index];
    const float sig2_beta2_val = sig2_vec[1*num_snp_ + snp_index];
    const float sig2_beta1[3] = {sig2_beta1_val, 0, sig2_beta1_val};
    const float sig2_beta2[3] = {0, sig2_beta2_val, sig2_beta2_val};
    const float rho[3] = {0, 0, rho_vec[snp_index]};

    const int num_components = 3;
    for (int comp_index = 0; comp_index < num_components; comp_index++) {
      const int index = (comp_index*num_snp_ + snp_index);
      float pi_val = pi_vec[index];
      if (pi_val > 0.5f) pi_val = 1.0f - pi_val;
      subset_sampler->p()[comp_index] = static_cast<double>(pi_val);
    }

    int sample_global_index = k_max_ - subset_sampler->sample_shuffle();
    const uint32_t* indices = subset_sampler->data();
    for (int comp_index = 0; comp_index < num_components; comp_index++) {
      const int index = (comp_index*num_snp_ + snp_index);

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

void BgmgCalculator::find_unified_bivariate_tag_delta_smplfast(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, int k_index, const float* nvec1, const float* nvec2, const float* hvec, std::vector<float>* tag_delta20, std::vector<float>* tag_delta02, std::vector<float>* tag_delta11, MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row) {
  tag_delta20->assign(num_tag_, 0.0f); tag_delta02->assign(num_tag_, 0.0f); tag_delta11->assign(num_tag_, 0.0f);

  const int global_snp_pi_index = 0; // here in smplfast we assume that all SNPs share the same pi
  const int num_components = 3;
  for (int comp_index = 0; comp_index < num_components; comp_index++) {
    const int index = (comp_index*num_snp_ + global_snp_pi_index);
    const float pi_val = pi_vec[index];
    subset_sampler->p()[comp_index] = static_cast<double>(pi_val);
  }

  int sample_global_index = num_snp_ - subset_sampler->sample_shuffle();
  const uint32_t* indices = subset_sampler->data();

  // int comp_index = 0;
  const int num_samples0=subset_sampler->counts()[0];
  for (int sample_index=0; sample_index < num_samples0; sample_index++, sample_global_index++) {
    const int snp_index = indices[sample_global_index];

    const float hval = hvec[snp_index];
    const float sig2_beta1 = sig2_vec[0*num_snp_ + snp_index] * hval;

    ld_matrix_csr_.extract_snp_row(SnpIndex(snp_index), ld_matrix_row);
    auto iter_end = ld_matrix_row->end();
    for (auto iter = ld_matrix_row->begin(); iter < iter_end; iter++) {
      const int tag_index = iter.index();
      const float r2 = iter.r2();
      tag_delta20->at(tag_index) += r2 * sig2_beta1;
    }
  }

  // int comp_index = 1;
  const int num_samples1=subset_sampler->counts()[1];
  for (int sample_index=0; sample_index < num_samples1; sample_index++, sample_global_index++) {
    const int snp_index = indices[sample_global_index];

    const float hval = hvec[snp_index];
    const float sig2_beta2 = sig2_vec[1*num_snp_ + snp_index] * hval;

    ld_matrix_csr_.extract_snp_row(SnpIndex(snp_index), ld_matrix_row);
    auto iter_end = ld_matrix_row->end();
    for (auto iter = ld_matrix_row->begin(); iter < iter_end; iter++) {
      const int tag_index = iter.index();
      const float r2 = iter.r2();
      tag_delta02->at(tag_index) += r2 * sig2_beta2;
    }
  }

  // int comp_index = 2;
  const int num_samples2=subset_sampler->counts()[2];
  for (int sample_index=0; sample_index < num_samples2; sample_index++, sample_global_index++) {
    const int snp_index = indices[sample_global_index];

    const float hval = hvec[snp_index];
    const float rho_val = rho_vec[snp_index];
    const float sig2_beta1 = sig2_vec[0*num_snp_ + snp_index] * hval;
    const float sig2_beta2 = sig2_vec[1*num_snp_ + snp_index] * hval;
    const float rho_beta12 = rho_vec[snp_index] * sqrt(sig2_beta1 * sig2_beta2);

    ld_matrix_csr_.extract_snp_row(SnpIndex(snp_index), ld_matrix_row);
    auto iter_end = ld_matrix_row->end();
    for (auto iter = ld_matrix_row->begin(); iter < iter_end; iter++) {
      const int tag_index = iter.index();
      const float r2 = iter.r2();
      tag_delta20->at(tag_index) += r2 * sig2_beta1;
      tag_delta02->at(tag_index) += r2 * sig2_beta2;
      tag_delta11->at(tag_index) += r2 * rho_beta12;
    }
  }

  assert(sample_global_index==num_snp_);
}

int64_t BgmgCalculator::calc_unified_univariate_pdf(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, int length, float* zvec, float* pdf) {
  check_num_snp(num_snp);

  std::stringstream ss;
  ss << "calc_unified_univariate_pdf(trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL << ", length=" << length << ", k_max_=" << k_max_ << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);

  const double pi_k = 1.0 / static_cast<double>(k_max_);
  std::vector<float>& nvec(*get_nvec(trait_index));
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);

  std::valarray<double> pdf_double(0.0, length);
#pragma omp parallel
  {
    std::valarray<double> pdf_double_local(0.0, length);
    LdMatrixRow ld_matrix_row;
    std::vector<float> tag_delta2(k_max_, 0.0f);

#pragma omp for schedule(dynamic, kOmpDynamicChunk)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + tag_to_snp_[tag_index], k_max_, num_components);
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
    {
      pdf_double += pdf_double_local;
    }
  }

  for (int i = 0; i < length; i++) pdf[i] = static_cast<float>(pdf_double[i]);
  LOG << "<" << ss.str() << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::calc_unified_univariate_power(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float zthresh, int length, float* nvec, float* svec) {
  std::stringstream ss;
  ss << "calc_unified_univariate_power(trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL << ", zthresh=" << zthresh << ", length=" << length << ", k_max_=" << k_max_ << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);
  const double pi_k = 1.0 / static_cast<double>(k_max_);
  std::vector<float> hvec; find_hvec(*this, &hvec);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> nvec_dummy(num_tag_, 1.0f);

  std::valarray<double> s_numerator_global(0.0, length);
  std::valarray<double> s_denominator_global(0.0, length);

#pragma omp parallel
  {
    std::valarray<double> s_numerator_local(0.0, length);
    std::valarray<double> s_denominator_local(0.0, length);
    LdMatrixRow ld_matrix_row;
    std::vector<float> tag_delta2(k_max_, 0.0f);

#pragma omp for schedule(dynamic, kOmpDynamicChunk)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + tag_to_snp_[tag_index], k_max_, num_components);
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
  ss << "calc_unified_univariate_delta_posterior(trait_index=" << trait_index << ", num_components=" << num_components << ", num_snp=" << num_snp << ", sig2_zeroA=" << sig2_zeroA << ", sig2_zeroC=" << sig2_zeroC << ", sig2_zeroL=" << sig2_zeroL << ", length=" << length << ", k_max_=" << k_max_ << ")";
  LOG << ">" << ss.str();

  if ((length == 0) || (length != num_tag_)) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_tag_"));

  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  std::vector<float>& nvec(*get_nvec(trait_index));
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);
  std::vector<float> hvec; find_hvec(*this, &hvec);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();

  std::valarray<double> c0_global(0.0f, num_tag_);
  std::valarray<double> c1_global(0.0f, num_tag_);
  std::valarray<double> c2_global(0.0f, num_tag_);

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    std::vector<float> tag_delta2(k_max_, 0.0f);
    std::valarray<double> c0_local(0.0f, num_tag_);
    std::valarray<double> c1_local(0.0f, num_tag_);
    std::valarray<double> c2_local(0.0f, num_tag_);

#pragma omp for schedule(dynamic, kOmpDynamicChunk)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      const float z = z_minus_fixed_effect_delta[tag_index];
      if (!std::isfinite(z)) continue;

      MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + tag_to_snp_[tag_index], k_max_, num_components);
      find_unified_univariate_tag_delta_sampling(num_components, pi_vec, sig2_vec, sig2_zeroC, tag_index, &nvec[0], &hvec[0], &tag_delta2, &subset_sampler, &ld_matrix_row);
    
      for (int k_index = 0; k_index < k_max_; k_index++) {

        const float delta2eff = tag_delta2[k_index] + ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index] * nvec[tag_index] * sig2_zeroL;  // S^2_kj
        const float sig2eff = delta2eff + sig2_zeroA;
        const float sig2eff_1_2 = sqrt(sig2eff);
        const float sig2eff_3_2 = sig2eff_1_2 * sig2eff;
        const float sig2eff_5_2 = sig2eff_3_2 * sig2eff;

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
  else if (cost_calculator_ == CostCalculator_Convolve) return calc_unified_bivariate_cost_convolve(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, aux);
  else if (cost_calculator_ == CostCalculator_Sampling) return calc_unified_bivariate_cost_sampling(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, aux, nullptr);
  else if (cost_calculator_ == CostCalculator_Smplfast) return calc_unified_bivariate_cost_smplfast(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, aux, nullptr);  
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
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
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

#pragma omp for schedule(dynamic, kOmpDynamicChunk)    
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];

      ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), &ld_matrix_row);
      auto iter_end = ld_matrix_row.end();
      for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
        const int snp_index = iter.index();
        const float r2_value = iter.r2();
        const float a2ij1 = sig2_zeroC[0] * nvec1_[tag_index] * hvec[snp_index] * r2_value;
        const float a2ij2 = sig2_zeroC[1] * nvec2_[tag_index] * hvec[snp_index] * r2_value;
        Edelta20_local[tag_index] += a2ij1 * Ebeta20[snp_index];
        Edelta02_local[tag_index] += a2ij2 * Ebeta02[snp_index];
        Edelta11_local[tag_index] += sqrt(a2ij1 * a2ij2) * Ebeta11[snp_index];
      }
    }

#pragma omp critical
    {
      Edelta20 += Edelta20_local;
      Edelta02 += Edelta02_local;
      Edelta11 += Edelta11_local;
    }
  }

  double log_pdf_total = 0;
  int num_zero_tag_r2 = 0;
  int num_infinite = 0;

#pragma omp parallel for schedule(dynamic, kOmpDynamicChunk) reduction(+: log_pdf_total, num_zero_tag_r2, num_infinite)
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
  int num_snp;
  const float* pi_vec;
  const float* sig2_vec;
  const float* rho_vec;
  const float* sig2_zeroA;
  const float* sig2_zeroC;
  const float* sig2_zeroL;
  float rho_zeroA;
  float rho_zeroL;
  int tag_index;  // for which SNP to calculate the characteristic function
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
  const int num_snp = data->num_snp;
  const float nval1 = (*data->nvec1)[data->tag_index];
  const float nval2 = (*data->nvec2)[data->tag_index];
  const float zval1 = (*data->z1_minus_fixed_effect_delta)[data->tag_index];
  const float zval2 = (*data->z2_minus_fixed_effect_delta)[data->tag_index];
  const float inf_adj_r2 = (*data->ld_tag_sum_r2_below_r2min_adjust_for_hvec)[data->tag_index];
  const float scaling_factor = 0.5 * M_1_PI * M_1_PI;

  double result = scaling_factor * cos(t1 * zval1 + t2 * zval2);

  result *= exp_quad_form(t1, t2,                        data->sig2_zeroA[0],
                                  data->rho_zeroA * sqrt(data->sig2_zeroA[0] * data->sig2_zeroA[1]),
                                                                               data->sig2_zeroA[1]);

  result *= exp_quad_form(t1, t2, inf_adj_r2 *                        nval1 * data->sig2_zeroL[0],
                                  inf_adj_r2 * data->rho_zeroL * sqrt(nval1 * data->sig2_zeroL[0] * nval2 * data->sig2_zeroL[1]),
                                  inf_adj_r2 *                                                      nval2 * data->sig2_zeroL[1]);

  auto iter_end = data->ld_matrix_row->end();
  for (auto iter = data->ld_matrix_row->begin(); iter < iter_end; iter++) {
    int snp_index = iter.index();

    const float r2 = iter.r2();
    const float hval = (*data->hvec)[snp_index];
    const float r2_times_hval = r2 * hval;

    const float pi1 = data->pi_vec[0*num_snp + snp_index];
    const float pi2 = data->pi_vec[1*num_snp + snp_index];
    const float pi12= data->pi_vec[2*num_snp + snp_index];
    const float pi0 = 1.0f - pi1 - pi2 - pi12;

    const float eff_sig2_beta1 = nval1 * r2_times_hval * data->sig2_zeroC[0] * data->sig2_vec[0*num_snp + snp_index];
    const float eff_sig2_beta2 = nval2 * r2_times_hval * data->sig2_zeroC[1] * data->sig2_vec[1*num_snp + snp_index];
    const float eff_sig2_beta_cov = data->rho_vec[snp_index] * sqrt(eff_sig2_beta1 * eff_sig2_beta2);

    result *= (double) (pi0 + 
                        pi1 * exp_quad_form(t1, t2, eff_sig2_beta1, 0, 0) +
                        pi2 * exp_quad_form(t1, t2, 0, 0, eff_sig2_beta2) + 
                        pi12* exp_quad_form(t1, t2, eff_sig2_beta1, eff_sig2_beta_cov, eff_sig2_beta2));
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

double BgmgCalculator::calc_unified_bivariate_cost_convolve(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux) {
  std::string ss = find_bivariate_params_description(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL);
  LOG << ">calc_unified_bivariate_cost_convolve(" << ss << ")";

  double log_pdf_total = 0.0;
  int num_snp_failed = 0;
  int num_infinite = 0;
  double func_evals = 0.0;
  double total_weight = 0.0;
  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);

  std::vector<float> weights_convolve(weights_.begin(), weights_.end());
  std::vector<float> weights_sampling(weights_.begin(), weights_.end()); int num_deftag_sampling = 0;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    const float tag_z1 = z1_minus_fixed_effect_delta[tag_index];
    const float tag_z2 = z2_minus_fixed_effect_delta[tag_index];
    const bool censoring = (std::abs(tag_z1) > z1max_) || (std::abs(tag_z2) > z2max_);
    if (censoring) {
      weights_convolve[tag_index] = 0;
      num_deftag_sampling++;
    } else {
      weights_sampling[tag_index] = 0;
    }
  }

  if (num_deftag_sampling > 0) {  // fall back to sampling approach for censored z-scores
    log_pdf_total += calc_unified_bivariate_cost_sampling(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, aux, &weights_sampling[0]);
  }

  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(&weights_convolve[0], &deftag_indices);

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    BivariateCharacteristicFunctionData data;
    data.num_snp = num_snp;
    data.pi_vec = pi_vec;
    data.sig2_vec = sig2_vec;
    data.rho_vec = rho_vec;
    data.sig2_zeroA = sig2_zeroA;
    data.sig2_zeroC = sig2_zeroC;
    data.sig2_zeroL = sig2_zeroL;
    data.rho_zeroA = rho_zeroA;
    data.rho_zeroL = rho_zeroL;
    data.ld_matrix_row = &ld_matrix_row;
    data.hvec = &hvec;
    data.z1_minus_fixed_effect_delta = &z1_minus_fixed_effect_delta;
    data.nvec1 = &nvec1_;
    data.z2_minus_fixed_effect_delta = &z2_minus_fixed_effect_delta;
    data.nvec2 = &nvec2_;
    data.ld_tag_sum_r2_below_r2min_adjust_for_hvec = &ld_tag_sum_r2_below_r2min_adjust_for_hvec;

#pragma omp for schedule(dynamic, kOmpDynamicChunk) reduction(+: log_pdf_total, num_snp_failed, num_infinite, func_evals, total_weight)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      int tag_index = deftag_indices[deftag_index];
      double tag_weight = static_cast<double>(weights_convolve[tag_index]);

      ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), data.ld_matrix_row);
      data.tag_index = tag_index;
      data.func_evals = 0;

      double tag_pdf = 0, tag_pdf_err = 0;
      const double xmin[2] = {0.0, -1.0}, xmax[2] = {1.0, 1.0};
      const int integrand_fdim = 1, ndim = 2;
      int cubature_result = hcubature(integrand_fdim, calc_bivariate_characteristic_function_for_integration,
        &data, ndim, xmin, xmax, cubature_max_evals_, cubature_abs_error_, cubature_rel_error_, ERROR_INDIVIDUAL, &tag_pdf, &tag_pdf_err);
      func_evals += (weights_convolve[tag_index] * (double)data.func_evals);
      total_weight += weights_convolve[tag_index];
      if (cubature_result != 0) { num_snp_failed++; continue; }

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

  LOG << "<calc_unified_bivariate_cost_convolve(" << ss << "), cost=" << log_pdf_total << ", evals=" << func_evals / total_weight << ", num_deftag=" << num_deftag << "+" << num_deftag_sampling << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

double BgmgCalculator::calc_unified_bivariate_cost_sampling(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux, const float* weights) {
  if (weights == nullptr) {
    if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
    weights = &weights_[0];
  }

  std::stringstream ss;
  ss << "calc_unified_bivariate_cost_sampling(" << find_bivariate_params_description(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL) << ", k_max_=" << k_max_ << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);
  
  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(weights, &deftag_indices);

  const double pi_k = 1.0 / static_cast<double>(k_max_);
  double log_pdf_total = 0.0;
  int num_infinite = 0;
  const int num_components = 3;

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    std::vector<float> tag_delta20(k_max_, 0.0f);
    std::vector<float> tag_delta02(k_max_, 0.0f);
    std::vector<float> tag_delta11(k_max_, 0.0f);
    std::valarray<double> tag_kpdf(0.0f, k_max_);

#pragma omp for schedule(dynamic, kOmpDynamicChunk) reduction(+: log_pdf_total, num_infinite)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      const int tag_index = deftag_indices[deftag_index];
      MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + tag_to_snp_[tag_index], k_max_, num_components);
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
        tag_kpdf[k] = pdf;        
      }

      // export the expected values of z^2 distribution
      if ((aux != nullptr) && (aux_option_ == AuxOption_Ezvec2)) {
        aux[0 * num_tag_ + tag_index] = average_tag_delta20;
        aux[1 * num_tag_ + tag_index] = average_tag_delta11;
        aux[2 * num_tag_ + tag_index] = average_tag_delta02;
      }
      if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdf)) aux[tag_index] = pdf_tag;
      if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdfErr)) aux[tag_index] = sqrt(((tag_kpdf - pdf_tag) * (tag_kpdf - pdf_tag)).sum() / (static_cast<double>(k_max_) - 1.0));

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

double BgmgCalculator::calc_unified_bivariate_cost_smplfast(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux, const float* weights) {
  if (weights == nullptr) {
    if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
    weights = &weights_[0];
  }

  std::stringstream ss;
  ss << "calc_unified_bivariate_cost_smplfast(" << find_bivariate_params_description(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL) << ", k_max_=" << k_max_ << ")";
  LOG << ">" << ss.str();

  const int num_components = 3;
  for (int comp_index = 0; comp_index < num_components; comp_index++) {
    const float ref_pi = pi_vec[comp_index*num_snp_];
    if (ref_pi > 0.5) LOG << " warning: pi_vec has values above 0.5, which is not optimized in calc_unified_bivariate_cost_smplfast; consider using calc_unified_bivariate_cost_sampling instead";
    for (int snp_index = 1; snp_index < num_snp; snp_index++) {
      if (ref_pi != pi_vec[comp_index*num_snp + snp_index])
        BGMG_THROW_EXCEPTION(::std::runtime_error("pi_vec doesn't appear to be constant across SNPs; use calc_unified_bivariate_cost_sampling instead"));
    }
  }

  SimpleTimer timer(-1);
  
  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(weights, &deftag_indices);

  const double pi_k = 1.0 / static_cast<double>(k_max_);
  std::valarray<double> pdf_double(0.0, num_tag_);
  std::valarray<double> average_tag_delta20(0.0, num_tag_);
  std::valarray<double> average_tag_delta11(0.0, num_tag_);
  std::valarray<double> average_tag_delta02(0.0, num_tag_);

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    std::vector<float> tag_delta20(num_tag_, 0.0f);
    std::vector<float> tag_delta02(num_tag_, 0.0f);
    std::vector<float> tag_delta11(num_tag_, 0.0f);
    std::valarray<double> pdf_double_local(0.0, num_tag_);
    std::valarray<double> average_tag_delta20_local(0.0, num_tag_);
    std::valarray<double> average_tag_delta11_local(0.0, num_tag_);
    std::valarray<double> average_tag_delta02_local(0.0, num_tag_);

#pragma omp for schedule(dynamic, 1) 
    for (int k_index = 0; k_index < k_max_; k_index++) {
      MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + k_index, num_snp_, num_components);
      find_unified_bivariate_tag_delta_smplfast(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, k_index, &nvec1_[0], &nvec2_[0], &hvec[0], &tag_delta20, &tag_delta02, &tag_delta11, &subset_sampler, &ld_matrix_row);
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

        const float a11 = tag_delta20[tag_index] * tag_n1 * sig2_zeroC[0] + sig2_zero_11;
        const float a12 = tag_delta11[tag_index] * sqrt(tag_n1 * tag_n2 * sig2_zeroC[0] * sig2_zeroC[1]) + sig2_zero_12;
        const float a22 = tag_delta02[tag_index] * tag_n2 * sig2_zeroC[1] + sig2_zero_22;
        const double pdf = static_cast<double>(censoring ? censored2_cdf<FLOAT_TYPE>(z1max_, z2max_, a11, a12, a22) : gaussian2_pdf<FLOAT_TYPE>(tag_z1, tag_z2, a11, a12, a22));
        pdf_double_local[tag_index] += pdf * pi_k;

        average_tag_delta20_local[tag_index] += a11 * pi_k;
        average_tag_delta11_local[tag_index] += a12 * pi_k;
        average_tag_delta02_local[tag_index] += a22 * pi_k;
      }
    }

#pragma omp critical
    {
      pdf_double += pdf_double_local;
      average_tag_delta20 += average_tag_delta20_local;
      average_tag_delta11 += average_tag_delta11_local;
      average_tag_delta02 += average_tag_delta02_local;
    }
  }

  double log_pdf_total = 0.0;
  int num_infinite = 0;
  for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
    int tag_index = deftag_indices[deftag_index];
    double increment = -std::log(pdf_double[tag_index]) * static_cast<double>(weights[tag_index]);
    if (!std::isfinite(increment)) {
      increment = static_cast<double>(-std::log(kMinTagPdf) * static_cast<double>(weights[tag_index]));
      num_infinite++;
    }
    log_pdf_total += increment;

    // export the expected values of z^2 distribution
    if ((aux != nullptr) && (aux_option_ == AuxOption_Ezvec2)) {
      aux[0 * num_tag_ + tag_index] = average_tag_delta20[tag_index];
      aux[1 * num_tag_ + tag_index] = average_tag_delta11[tag_index];
      aux[2 * num_tag_ + tag_index] = average_tag_delta02[tag_index];
    }
    if ((aux != nullptr) && (aux_option_ == AuxOption_TagPdf)) aux[tag_index] = pdf_double[tag_index];
  }

  if (num_infinite > 0)
    LOG << " warning: infinite increments encountered " << num_infinite << " times";

  LOG << "<" << ss.str() << ", cost=" << log_pdf_total << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return log_pdf_total;
}

void BgmgCalculator::calc_bivariate_delta_posterior_integrals(float a, float b, float c, float i, float j, float k, float z1, float z2,
                                                              float* c00, float* c10, float* c01, float* c20, float* c11, float* c02) {
  // [a b; b c] is variance-covariance matrix of eps in "z=delta+eps"
  // [i j; j k] is variance-covariance matrix of delta for a specific choice of causal variants in sampling
  // here we are calculating moments cPQ = E[delta1^P delta2^Q | z1, z2].
  // the formulas are pretty tricky to derive - loads of integrals involved.

  static const float inv_sqrt2_pi = 0.2250790790392765f;
  const float z12 = z1*z1;
  const float z22 = z2*z2;
  const float z1z2 = z1*z2;

  const float ci = c*i;
  const float bj = b*j;
  const float j2 = j*j;
  const float ik = i*k;
  const float bi = b*i;
  const float aj = a*j;
  const float cj = c*j;
  const float bk = b*k;
  const float ak = a*k;
  const float ij = i*j;
  const float i2 = i*i;
  const float k2 = k*k;
  const float jk = j*k;

  const float ci2 = ci*i;
  const float bij = bi*j;
  const float aj2 = aj*j;
  const float ij2 = ij*j;
  const float i2k = i2*k;
  const float cj2 = cj*j;
  const float bjk = bj*k;
  const float j2k = j2*k;
  const float ak2 = ak*k;
  const float ik2 = ik*k;
  const float cij = ci*j;
  const float bj2 = bj*j;
  const float bik = bi*k;
  const float ajk = aj*k;
  const float ijk = ij*k;

  const float j3 = j2*j;

  const float eN = (c + k) * z12 - 2 * (b + j) * z1z2 + (a + i) * z22;
  const float eD = -2.0f * (b + j)*(b + j) + 2.0f * (a + i) * (c + k);

  const float eD2 = eD*eD;
  const float inv_eD = 1.0f / eD;
  const float inv_eD2 = inv_eD*inv_eD;
  const float x2eD = 2.0f*eD;
  const float x2eD_x4eN = x2eD - 4.0f*eN;

  const float c00_tmp = inv_sqrt2_pi * exp(-eN * inv_eD) * sqrt(inv_eD);

  const float c10_pol = ci * z1 - bj * z1 - j2 * z1 + ik * z1 - bi * z2 + aj * z2;
  const float c01_pol = cj * z1 - bk * z1 - bj * z2 - j2 * z2 + ak * z2 + ik * z2;
  const float c20_pol = eD2*i + x2eD_x4eN*(bij + bij - ci2 - aj2 + ij2 - i2k) - x2eD*(j2*z12 - 2.0f*ij*z1z2           + i2*z22);
  const float c02_pol = eD2*k + x2eD_x4eN*(bjk + bjk - cj2 - ak2 + j2k - ik2) - x2eD*(k2*z12 - 2.0f*jk*z1z2           + j2*z22);
  const float c11_pol = eD2*j + x2eD_x4eN*(bj2 + j3  - cij - ajk + bik - ijk) - x2eD*(jk*z12 -      j2*z1z2 - ik*z1z2 + ij*z22);

  (*c00) = c00_tmp;
  (*c10) = c00_tmp * 2.0f * inv_eD  * c10_pol;
  (*c01) = c00_tmp * 2.0f * inv_eD  * c01_pol;
  (*c20) = c00_tmp *        inv_eD2 * c20_pol;
  (*c11) = c00_tmp *        inv_eD2 * c11_pol;
  (*c02) = c00_tmp *        inv_eD2 * c02_pol;  
}

int64_t BgmgCalculator::calc_unified_bivariate_pdf(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, int length, float* zvec1, float* zvec2, float* pdf) {
  std::stringstream ss;
  ss << "calc_unified_bivariate_pdf(" << find_bivariate_params_description(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL) << ", k_max_=" << k_max_ << ", length=" << length << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);

  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);
  const double pi_k = 1.0 / static_cast<double>(k_max_);
  const int num_components = 3;

  std::valarray<double> pdf_double(0.0, length);

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    std::valarray<double> pdf_double_local(0.0, length);
    std::vector<float> tag_delta20(k_max_, 0.0f);
    std::vector<float> tag_delta02(k_max_, 0.0f);
    std::vector<float> tag_delta11(k_max_, 0.0f);

#pragma omp for schedule(dynamic, kOmpDynamicChunk)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      const int tag_index = deftag_indices[deftag_index];
      MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + tag_to_snp_[tag_index], k_max_, num_components);
      double tag_weight = static_cast<double>(weights_[tag_index]);

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

      for (int k = 0; k < k_max_; k++) {
        const float a11 = tag_delta20[k] + sig2_zero_11;
        const float a12 = tag_delta11[k] + sig2_zero_12;
        const float a22 = tag_delta02[k] + sig2_zero_22;

        for (int z_index = 0; z_index < length; z_index++) {
          double pdf_tmp = static_cast<double>(gaussian2_pdf<FLOAT_TYPE>(zvec1[z_index], zvec2[z_index], a11, a12, a22));
          pdf_double_local[z_index] += pi_k * pdf_tmp * tag_weight;
        }
      }
    }
#pragma omp critical
    {
      pdf_double += pdf_double_local;
    }
  }

  for (int i = 0; i < length; i++) pdf[i] = static_cast<float>(pdf_double[i]);

  LOG << "<" << ss.str() << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::calc_unified_bivariate_delta_posterior(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL,
                                                               int length, float* c00, float* c10, float* c01, float* c20, float* c11, float* c02) {
  std::stringstream ss;
  ss << "calc_unified_bivariate_delta_posterior(" << find_bivariate_params_description(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL) << ", k_max_=" << k_max_ << ")";
  LOG << ">" << ss.str();

  SimpleTimer timer(-1);
  
  // standard variables
  std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
  std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);
  const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
  std::vector<float> hvec; find_hvec(*this, &hvec);
  std::vector<int> deftag_indices; const int num_deftag = find_deftag_indices(nullptr, &deftag_indices);

  const double pi_k = 1.0 / static_cast<double>(k_max_);
  const int num_components = 3;

  // Sigma0  = [a0 b0; b0 c0];
  const float a0 = sig2_zeroA[0];
  const float c0 = sig2_zeroA[1];
  const float b0 = sqrt(a0 * c0) * rho_zeroA;

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    double c00_local, c10_local, c01_local, c20_local, c11_local, c02_local;
    std::vector<float> tag_delta20(k_max_, 0.0f);
    std::vector<float> tag_delta02(k_max_, 0.0f);
    std::vector<float> tag_delta11(k_max_, 0.0f);

#pragma omp for schedule(dynamic, kOmpDynamicChunk)
    for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
      const int tag_index = deftag_indices[deftag_index];
      MultinomialSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1 + tag_to_snp_[tag_index], k_max_, num_components);
      const float adj_hval = ld_tag_sum_r2_below_r2min_adjust_for_hvec[tag_index];
      const float sig2_zeroL_11 = adj_hval * nvec1_[tag_index] * sig2_zeroL[0];
      const float sig2_zeroL_22 = adj_hval * nvec2_[tag_index] * sig2_zeroL[1];
      const float sig2_zeroL_12 = adj_hval * rho_zeroL * sqrt(nvec1_[tag_index] * nvec2_[tag_index] * sig2_zeroL[0] * sig2_zeroL[1]);

      const float tag_z1 = z1_minus_fixed_effect_delta[tag_index];
      const float tag_z2 = z2_minus_fixed_effect_delta[tag_index];
      const float tag_n1 = nvec1_[tag_index];
      const float tag_n2 = nvec2_[tag_index];

      find_unified_bivariate_tag_delta_sampling(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, tag_index, &nvec1_[0], &nvec2_[0], &hvec[0], &tag_delta20, &tag_delta02, &tag_delta11, &subset_sampler, &ld_matrix_row);

      c00_local = 0; c10_local = 0; c01_local = 0; c20_local = 0; c11_local = 0; c02_local = 0;
      for (int k = 0; k < k_max_; k++) {
        const float A = tag_delta20[k] + sig2_zeroL_11;
        const float B = tag_delta11[k] + sig2_zeroL_12;
        const float C = tag_delta02[k] + sig2_zeroL_22;

        float c00buf, c10buf, c01buf, c20buf, c11buf, c02buf;
        BgmgCalculator::calc_bivariate_delta_posterior_integrals(a0, b0, c0, A, B, C, tag_z1, tag_z2, &c00buf, &c10buf, &c01buf, &c20buf, &c11buf, &c02buf);
        c00_local += static_cast<double>(c00buf);
        c10_local += static_cast<double>(c10buf);
        c01_local += static_cast<double>(c01buf);
        c20_local += static_cast<double>(c20buf);
        c11_local += static_cast<double>(c11buf);
        c02_local += static_cast<double>(c02buf);
      }
      
      c00[tag_index] = pi_k * c00_local;
      c10[tag_index] = pi_k * c10_local;
      c01[tag_index] = pi_k * c01_local;
      c20[tag_index] = pi_k * c20_local;
      c11[tag_index] = pi_k * c11_local;
      c02[tag_index] = pi_k * c02_local;
    }
  }

  LOG << "<" << ss.str() << ", num_deftag=" << num_deftag << ", elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

