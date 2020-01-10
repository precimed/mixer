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

/*
  Implementation notes:

  1. Loops across tag indices should be implemented as follows:
  
        std::vector<int> deftag_indices;
        const int num_deftag = find_deftag_indices(trait_index, &deftag_indices);   # univariate
        const int num_deftag = find_deftag_indices(             &deftag_indices);   # bivariate
        for (int deftag_index = 0; deftag_index < num_deftag; deftag_index++) {
          int tag_index = deftag_indices[deftag_index];
          ...
        }
  
     Note that find_deftag_indices() checks and rises an exception is zvec/nvec are undefined.
     Also note that only SNPs with positive weight[tag_index] are selected by find_deftag_indices().

  2. Access to zvec is, typically, via the following construct:

      std::vector<float> z_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(trait_index, &z_minus_fixed_effect_delta);
  
     or
  
       std::vector<float> z1_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(1, &z1_minus_fixed_effect_delta);
       std::vector<float> z2_minus_fixed_effect_delta; find_z_minus_fixed_effect_delta(2, &z2_minus_fixed_effect_delta);

  3. Access to LD matrix is as follows:

      LdMatrixRow ld_matrix_row;
      ld_matrix_csr_.extract_snp_row(snp_index, &ld_matrix_row);  // case A: mapping from snp to tag
      ld_matrix_csr_.extract_tag_row(tag_index, &ld_matrix_row);  // case B: mapping from tag to snp

      auto iter_end = ld_matrix_row.end();
      for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
        int snp_index = iter.index();   // for case A
        int tag_index = iter.index();   // for case B
        float r2 = iter.r2();
        ...
      }

  4. Other frequent constructs:

     const double zmax = (trait_index==1) ? z1max_ : z2max_;
     const double pi_k = 1.0 / static_cast<double>(k_max_);
     const std::vector<float>& ld_tag_sum_r2_below_r2min_adjust_for_hvec = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_tag_sum_r2_below_r2min();
     std::vector<float> hvec; find_hvec(*this, &hvec);
     if (snp_order_.empty()) find_snp_order();
     if (cache_tag_r2sum_) find_tag_r2sum(component_id, num_causals);

*/

#include "bgmg_calculator.h"

#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_set_num_threads(i)
#define omp_get_thread_num() 0
#endif

#include <chrono>
#include <random>
#include <limits>
#include <algorithm>
#include <vector>
#include <valarray>
#include <cmath>
#include <sstream>
#include <fstream>
#include <set>
#include <numeric>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include "cubature/cubature.h"

#include "bgmg_log.h"
#include "bgmg_parse.h"
#include "bgmg_math.h"
#include "bgmg_rand.h"
#include "fmath.hpp"

#define FLOAT_TYPE float

// pdf of gaussian (normal) distribution with 0 mean and std error of s
// z is the point where to calculate pdf
template<typename T>
inline T gaussian_pdf(const T z, const T s) {
  static const T inv_sqrt_2pi = static_cast<T>(0.3989422804014327);
  const T a = z / s;
  const T pdf = inv_sqrt_2pi / s * std::exp(static_cast<T>(-0.5) * a * a);
  return pdf + std::numeric_limits<T>::min();
}

// censored_cdf gives cumulated probability of |z|>zmax, where z is normal random variable
// http://www.public.iastate.edu/~stat415/meeker/ml_estimation_chapter.pdf
// Principles of Maximum Likelihood Estimation and The Analysis of Censored Data
template<typename T>
inline T censored_cdf(const T z, const T s) {
  assert(z >= 0);
  static const T inv_sqrt_2 = static_cast<T>(0.7071067811865475);
  return std::erfc((z / s) * inv_sqrt_2) + std::numeric_limits<T>::min();
}

/*
// partial specification for float, to use fmath::exp instead of std::exp
template<>
inline float gaussian_pdf<float>(const float z, const float s) {
  static const float inv_sqrt_2pi = static_cast<float>(0.3989422804014327);
  const float a = z / s;
  const float pdf = inv_sqrt_2pi / s * fmath::exp(static_cast<float>(-0.5) * a * a);
  return pdf;
}
*/

template<typename T>
inline T gaussian2_pdf(const T z1, const T z2, const T a11, const T a12, const T a22) {
  static const T log_pi = static_cast<T>(-1.0 * log(2.0 * 3.14159265358979323846));

  // Calculation of log - likelihood and pdf, specific to bivariate normal
  // distribution with zero mean.It takes into account an explicit formula
  // for inverse 2x2 matrix, S = [a b; c d], => S^-1 = [d - b; -c a] . / det(S)
  const T dt = a11 * a22 - a12 * a12;  // det(S)

  const T log_exp = -0.5 * (a22*z1*z1 + a11*z2*z2 - 2.0*a12*z1*z2) / dt;
  const T log_dt = -0.5 * std::log(dt);

  const T pdf = std::exp(log_pi + log_dt + log_exp);
  return pdf + std::numeric_limits<T>::min();
}

/*
template<>
inline float gaussian2_pdf(const float z1, const float z2, const float a11, const float a12, const float a22) {
  static const float log_pi = static_cast<float>(-1.0 * log(2.0 * 3.14159265358979323846));

  // Calculation of log - likelihood and pdf, specific to bivariate normal
  // distribution with zero mean.It takes into account an explicit formula
  // for inverse 2x2 matrix, S = [a b; c d], => S^-1 = [d - b; -c a] . / det(S)
  const float dt = a11 * a22 - a12 * a12;  // det(S)

  const float log_exp = -0.5 * (a22*z1*z1 + a11*z2*z2 - 2.0*a12*z1*z2) / dt;
  const float log_dt = -0.5 * fmath::log(dt);

  const float pdf = fmath::exp(log_pi + log_dt + log_exp) + std::numeric_limits<float>::min();
  return pdf;
}
*/
