#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/owens_t.hpp>

// from STAN: http://www.stat.columbia.edu/~gelman/bda.course/_book/custom-probability-functions-chapter.html
// P(Z1>=z1 && Z2 >=z2)
template<typename T>
inline T binormal_cdf(T z1, T z2, T rho) {
  static const T inv_sqrt_2 = static_cast<T>(0.7071067811865475);
  if (z1 != 0 || z2 != 0) {
    T denom = std::abs(rho) < 1.0 ? sqrt((1 + rho) * (1 - rho)) : NAN;
    T a1 = (z2 / z1 - rho) / denom;
    T a2 = (z1 / z2 - rho) / denom;
    T product = z1 * z2;
    T delta = product < 0 || (product == 0 && (z1 + z2) < 0);
    T Phi_z1 = 0.5 * (1 + std::erf(z1 * inv_sqrt_2));
    T Phi_z2 = 0.5 * (1 + std::erf(z2 * inv_sqrt_2));
    return 0.5 * (Phi_z1 + Phi_z2 - delta) - boost::math::owens_t(z1, a1) - boost::math::owens_t(z2, a2);
  }
  return 0.25 + std::asin(rho) / (2.0 * boost::math::constants::pi<T>());
}

template<typename T>
inline T censored2_cdf(T z1max, T z2max, T a11, T a12, T a22) {
  assert((z1max > 0) && (z2max > 0));
  const T sqrt_a11 = sqrt(a11);
  const T sqrt_a22 = sqrt(a22);
  const T rho = a12 / (sqrt_a11 * sqrt_a22);
  const T z1norm = z1max / sqrt_a11;
  const T z2norm = z2max / sqrt_a22;

  // return area OUTSIDE the rectangle [-z1norm, -z2norm] x [z1norm, z2norm]
  return 1.0
    + binormal_cdf(-z1norm, z2norm, rho)
    + binormal_cdf(z1norm, -z2norm, rho)
    - binormal_cdf(z1norm, z2norm, rho)
    - binormal_cdf(-z1norm, -z2norm, rho);
}
