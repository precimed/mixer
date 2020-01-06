#include <cmath>
#include <numeric>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
// #include <boost/math/special_functions/owens_t.hpp>

// This file contains two implementations of bivariate cumulative density function
// One is copied from STAN library, another from BVNcdf package for matlab (see links below).
// BVNcdf is Based on Matlab code from Alan Getz, e.i. http://www.math.wsu.edu/faculty/genz/software/matlab/bvnl.m

#define Inf     1e+20
#define pi      3.14159265358979323846
#define twopi   2*pi
#define MaxQuad 10

// from STAN: http://www.stat.columbia.edu/~gelman/bda.course/_book/custom-probability-functions-chapter.html
// P(Z1>=z1 && Z2 >=z2)
template<typename T>
inline T binormal_cdf_stan(T z1, T z2, T rho) {
  return 0;  // disable because boost 1.49.0 does not support owens_t
  /*
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
  return 0.25 + std::asin(rho) / (2.0 * pi);
  */
}


// phid and BVNcdf are copied from http://www.tjeconomics.com/code/ 
inline double phid(double x) {
  // Based on "BETTER APPROXIMATIONS TO CUMULATIVE NORMAL FUNCTIONS" by Craeme West, 2004--> double precision
  double XAbs = std::abs(x), build, Cumnorm;
  if (XAbs > 37) {
    return (x<0) ? 0 : 1;       // <------ here it was a bug in BVNcdf library
  }
  else {
    if (XAbs < 7.07106781186547) {
      build = 3.52624965998911E-02 * XAbs + 0.700383064443688;
      build = build * XAbs + 6.37396220353165;
      build = build * XAbs + 33.912866078383;
      build = build * XAbs + 112.079291497871;
      build = build * XAbs + 221.213596169931;
      build = build * XAbs + 220.206867912376;
      Cumnorm = std::exp(-XAbs*XAbs / 2) * build;
      build = 8.83883476483184E-02 * XAbs + 1.75566716318264;
      build = build * XAbs + 16.064177579207;
      build = build * XAbs + 86.7807322029461;
      build = build * XAbs + 296.564248779674;
      build = build * XAbs + 637.333633378831;
      build = build * XAbs + 793.826512519948;
      build = build * XAbs + 440.413735824752;
      Cumnorm = Cumnorm / build;
    }
    else {
      build = XAbs + 0.65;
      build = XAbs + 4 / build;
      build = XAbs + 3 / build;
      build = XAbs + 2 / build;
      build = XAbs + 1 / build;
      Cumnorm = std::exp(-XAbs*XAbs / 2) / build / 2.506628274631;
    }
  }
  if (x > 0) {
    Cumnorm = 1 - Cumnorm;
  }
  return Cumnorm;
}

// phid and BVNcdf are copied from http://www.tjeconomics.com/code/ 
// P(x < dh & y < dk)
inline double BVNcdf(double dh, double dk, double r) {
  // Based on Matlab code from Alan Getz
  double w[MaxQuad], x[MaxQuad], hk, bvn, hs, asr, sn, as, a, bs, c, d, b, sp, xs, rs, ep;
  int lg, i, j, is;

  if (dh >= Inf && dk >= Inf)
    return 1;
  else if (dh <= -Inf || dk <= -Inf)
    return 0;
  else if (dk >= Inf)
    return phid(dh);
  else if (dh >= Inf)
    return phid(dk);
  else if (r == 0) {
    return phid(dh)*phid(dk);          // <------ this step was missing in BVNcdf library, but present in the original Alan Getz code 
  }
  else {
    dh = -dh; dk = -dk;
    if (std::abs(r) < 0.3) {
      lg = 3;
      //       Gauss Legendre points and weights, n =  6
      w[0] = 0.1713244923791705; w[1] = 0.3607615730481384; w[2] = 0.4679139345726904;
      x[0] = 0.9324695142031522; x[1] = 0.6612093864662647; x[2] = 0.2386191860831970;
    }
    else if (std::abs(r) < 0.75) {
      lg = 6;
      //       Gauss Legendre points and weights, n = 12
      w[0] = .04717533638651177; w[1] = 0.1069393259953183; w[2] = 0.1600783285433464;
      w[3] = 0.2031674267230659; w[4] = 0.2334925365383547; w[5] = 0.2491470458134029;
      x[0] = 0.9815606342467191; x[1] = 0.9041172563704750; x[2] = 0.7699026741943050;
      x[3] = 0.5873179542866171; x[4] = 0.3678314989981802; x[5] = 0.1252334085114692;
    }
    else {
      lg = 10;
      //       Gauss Legendre points and weights, n = 20
      w[0] = .01761400713915212; w[1] = .04060142980038694; w[2] = .06267204833410906;
      w[3] = .08327674157670475; w[4] = 0.1019301198172404; w[5] = 0.1181945319615184;
      w[6] = 0.1316886384491766; w[7] = 0.1420961093183821; w[8] = 0.1491729864726037; w[9] = 0.1527533871307259;
      x[0] = 0.9931285991850949; x[1] = 0.9639719272779138; x[2] = 0.9122344282513259;
      x[3] = 0.8391169718222188; x[4] = 0.7463319064601508; x[5] = 0.6360536807265150;
      x[6] = 0.5108670019508271; x[7] = 0.3737060887154196; x[8] = 0.2277858511416451; x[9] = 0.07652652113349733;
    }
    hk = dh*dk; bvn = 0;
    if (std::abs(r) < 0.925) {
      hs = (dh*dh + dk*dk) / 2; asr = std::asin(r); // findes asin i math.h? Det ser det ud til
      for (i = 0; i<lg; i++) {
        sn = std::sin(asr*(1 - x[i]) / 2);
        bvn = bvn + w[i] * std::exp((sn*hk - hs) / (1 - sn*sn));
        sn = std::sin(asr*(1 + x[i]) / 2);
        bvn = bvn + w[i] * std::exp((sn*hk - hs) / (1 - sn*sn));
      }
      bvn = bvn*asr / (4 * pi);
      bvn = bvn + phid(-dh)*phid(-dk);
    }
    else {
      if (r < 0) { dk = -dk; hk = -hk; }
      if (std::abs(r) < 1) {
        as = (1.0 - r)*(1.0 + r); a = std::sqrt(as); bs = (dh - dk)*(dh - dk);
        c = (4.0 - hk) / 8.0; d = (12.0 - hk) / 16.0; asr = -(bs / as + hk) / 2.0;
        if (asr > -100) { bvn = a*std::exp(asr)*(1 - c*(bs - as)*(1 - d*bs / 5) / 3 + c*d*as*as / 5); }
        if (hk > -100) {
          b = std::sqrt(bs); sp = std::sqrt(twopi)*phid(-b / a);
          bvn = bvn - std::exp(-hk / 2.0)*sp*b*(1.0 - c*bs*(1.0 - d*bs / 5.0) / 3.0);
        }
        a = a / 2.0;
        for (i = 0; i<lg; i++) {
          for (j = 0; j <= 1; j++) {
            is = -1 * (j == 0) + 1 * (j == 1);
            xs = (a + a*is*x[i])*(a + a*is*x[i]);
            rs = std::sqrt(1.0 - xs); asr = -(bs / xs + hk) / 2.0;
            if (asr > -100) {
              sp = (1.0 + c*xs*(1.0 + d*xs));
              ep = std::exp(-hk*(1.0 - rs) / (2.0*(1.0 + rs))) / rs;
              bvn = bvn + a*w[i] * std::exp(asr)*(ep - sp);
            }
          }
        }
        bvn = -bvn / twopi;
      }
      if (r > 0) { bvn = bvn + phid(-std::max<double>(dh, dk)); }
      else if (r < 0) { bvn = -bvn + std::max<double>(0, phid(-dh) - phid(-dk)); }
    }
    return std::max<double>(0, std::min<double>(1, bvn));
  }
}

template<typename T>
inline T censored2_cdf_stan(T z1max, T z2max, T a11, T a12, T a22) {
  assert((z1max > 0) && (z2max > 0));
  const T sqrt_a11 = sqrt(a11);
  const T sqrt_a22 = sqrt(a22);
  const T rho = a12 / (sqrt_a11 * sqrt_a22);
  const T z1norm = z1max / sqrt_a11;
  const T z2norm = z2max / sqrt_a22;

  // return area OUTSIDE the rectangle [-z1norm, -z2norm] x [z1norm, z2norm]
  return std::max<T>(std::numeric_limits<T>::min(), 1.0
    + binormal_cdf_stan(-z1norm, z2norm, rho)
    + binormal_cdf_stan(z1norm, -z2norm, rho)
    - binormal_cdf_stan(z1norm, z2norm, rho)
    - binormal_cdf_stan(-z1norm, -z2norm, rho));
}

template<typename T>
inline T censored2_cdf_BVN(T z1max, T z2max, T a11, T a12, T a22) {
  assert((z1max > 0) && (z2max > 0));
  const T sqrt_a11 = sqrt(a11);
  const T sqrt_a22 = sqrt(a22);
  const T rho = a12 / (sqrt_a11 * sqrt_a22);
  const T z1norm = z1max / sqrt_a11;
  const T z2norm = z2max / sqrt_a22;

  // return area OUTSIDE the rectangle [-z1norm, -z2norm] x [z1norm, z2norm]
  // in makes sure that all numbers are close to 0 (not close to 1), thus ensure we do not lose precision 
  T x1 = BVNcdf(0, -z2norm, rho);
  T x2 = BVNcdf(-z1norm, -z2norm, rho);
  T x3 = BVNcdf(-z1norm, z2norm, rho);
  T x4 = BVNcdf(0, -z2norm, -rho);
  //if (x1>0.3 || x2>0.3 || x3 > 0.3 || x4>0.3) printf("%.2f: %.3e %.3e %.3e %.3e\n", x1, x2, x3, x3);  // this is not supposed to happen often -
  return std::max<T>(std::numeric_limits<T>::min(), 2 * (x1 - x2 + x3 + x4));

  // Here we utilize some of the "symmetric reflection" properties, e.i. (7), (8) or (9):
  // https://www.jstor.org/stable/pdf/2006276.pdf
  // (1)  F(h, k, r) = Pr(x1<h & x2<k)
  // (7)  F(h, k, r) = f(h) + f(k) - 1 + F(-h, -k, r)
  // (8)  F(h, k, r) =        f(k)     - F(-h,  k, -r)
  // (9)  F(h, k, r) = f(h)            - F(h,  -k, -r)
  // (10) f(h) = 1/(2 pi) \int_{-inf}^h exp(-x^2 / 2) dx
}

// bgmg_math_test does not reveal any significant differences between STAN and BVNcdf implementations.
// going for STAN because it looks simpler.
template<typename T>
inline T censored2_cdf(T z1max, T z2max, T a11, T a12, T a22) {
  return censored2_cdf_BVN<T>(z1max, z2max, a11, a12, a22);
}
