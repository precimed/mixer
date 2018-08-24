#include "gtest/gtest.h"

#include "nlopt/neldermead.h"

#include <vector>
#include <numeric>
#include <assert.h>

namespace {

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
  static int nevals = 0;
  assert(grad == nullptr);
  auto sqr = [](double x) {return x*x; };
  double cost = 100.0 * sqr(x[1] - sqr(x[0])) + sqr(1 - x[0]);
  // printf("%i: f(%.3f, %.3f)=%.3f\n", ++nevals, x[0], x[1], cost);
  return cost;
}

void TestRosenbrockMin() {
  // tricks from matlab's fminsearch
  double usual_delta = 0.05;             // 5 percent deltas for non-zero terms
  double zero_term_delta = 0.00025;      // Even smaller delta for zero elements of x

  int nargs = 2;
  double x[2] = { -1.2, 2 };
  double xtol_abs[2] = { 1e-4, 1e-4 };
  double xstep[2] = { usual_delta*x[0], usual_delta*x[1] };
  double ub[2] = { 10, 10 };
  double lb[2] = { -10, -10 };
  int numevals = 0;

  nlopt_stopping stop; 
  stop.n = nargs;
  stop.minf_max = -std::numeric_limits<double>::infinity();
  stop.ftol_rel = 0;
  stop.ftol_abs = 0;
  stop.xtol_rel = 0;
  stop.xtol_abs = xtol_abs;
  stop.nevals_p = &numevals;
  stop.maxeval = std::numeric_limits<int>::max();
  stop.maxtime = std::numeric_limits<double>::max();
  stop.start = nlopt_seconds();
  stop.force_stop = nullptr;
  stop.stop_msg = nullptr;

  double minf;
  nlopt_result nldrmd_result = nldrmd_minimize(nargs, myfunc,
    /*void *f_data*/ nullptr,
    lb, /* lower bound */
    ub, /* upper bound */
    x, /* in: initial guess, out: minimizer */
    &minf,
    xstep, /* initial step sizes */
    &stop);

  printf("solution: (1,1)+(%.3e, %.3e), neval=%i\n", x[0]-1, x[1]-1, *stop.nevals_p);

  ASSERT_TRUE(true);
}

// bgmg-test.exe --gtest_filter=Nlopt.NeldermeadRosenbrockMin
TEST(Nlopt, NeldermeadRosenbrockMin) {
  TestRosenbrockMin();
}

}  // namespace
