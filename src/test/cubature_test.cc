#include "cubature/cubature.h"

#include "gtest/gtest.h"
#include "assert.h"
#include "math.h"

#include <vector>

int f0 (unsigned ndim, const double *x, void *, unsigned fdim, double* fval) {
    assert(fdim == 1);
    double prod = 1.0;
    for (unsigned int i = 0; i < ndim; ++i)
        prod *= cos(x[i]);
    fval[0] = prod;
    return 0;
}

void TestHcubature() {
    const int integrand_fdim = 1;
    const int ndim = 2;
    const int maxEval = 0;
    const double reqAbsError = 0;
    const double reqRelError = 1e-6;
    std::vector<double> xmin = {0.0, 0.0}, xmax={1.0, 1.0};
    std::vector<double> val = {0.0}, err={0.0};
    ASSERT_TRUE(0 == hcubature(integrand_fdim, f0, nullptr, ndim, &xmin[0], &xmax[0], 
                               maxEval, reqAbsError, reqRelError, ERROR_INDIVIDUAL, &val[0], &err[0]));
    ASSERT_FLOAT_EQ(val[0], 0.70807342596362365938);
}

// bgmg-test.exe --gtest_filter=Cubature.Hcubature
TEST(Cubature, Hcubature) {
  TestHcubature();
}
