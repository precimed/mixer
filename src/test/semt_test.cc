#include "gtest/gtest.h"
#include <iostream>
using namespace std;

// Include namespace SEMT & global operators.
#define SEMT_DISABLE_PRINT 0
#include "semt/Semt.h"
// Include macros: INT, DINT, VAR, DVAR, PARAM, DPARAM
#include "semt/Shortcuts.h"
using namespace SEMT;

namespace {
  
int introduction()
{
    // Define an instance of variable with index 0 named x.
    DVAR(x, 0);

    // f = x^3 + sin(x)
    auto f = pow(x, INT(3)) + sin(x);
    cout << f << endl;

    // We have one variable and want four derivatives.
    DifferentiableVectorExpr<1, 4>DVE;
    // Insert the expression, the actual differentiation has taken place
    // long before the whole program is even executed ^^
    // But this line issues the compiler to instantiate the derivatives.
    DVE.push_back(f);

    // Print derivatives and evaluate.
    vector<SEMT_PRECISION>x0(1, 3.1415);
    for (int i = 0; i < 5; ++i)
        cout << DVE.get_derivative(i) << " == "
                << (DVE.get_derivative(i)(x0))
                << " for x0 = " << (x0[0]) << endl;
    return 0;
// output:
// ((x0)^(3) + sin(x0))
// [    ((x0)^(3) + sin(x0))    ] == [ 31.0036 ] for x0 = 3.1415
// [    ((3 * (x0)^(2)) + cos(x0))  ] == [ 28.6071 ] for x0 = 3.1415
// [    ((3 * (2 * x0)) + (-1 * sin(x0)))   ] == [ 18.8489 ] for x0 = 3.1415
// [    ((3 * 2) + (-1 * cos(x0)))  ] == [    7 ] for x0 = 3.1415
// [    (-1 * (-1 * sin(x0)))   ] == [ 9.26536e-05 ] for x0 = 3.1415
}


int macros()
{
  DVAR(x2, 2);
  DINT(Two, 2);
  DPARAM(t1, 1)
    cout << (VAR(0) * VAR(1) - t1 + pow(x2, Two)) << endl;
  cout << deriv_t(pow(VAR(0) * x2, PARAM(0)), x2) << endl;
  return 0;
  // output:
  // (((x0 * x1) - t1) + (x2)^(2))
  // (((x0 * x2))^(t0) * ((t0 * x0) / (x0 * x2)))
}

int binary_operators()
{
  DVAR(x0, 0);
  DVAR(x1, 1);
  auto f1 = INT(5) * x0 + x1;
  // instead of: Expr< Plus< Times <Integer<5>, Variable<0> >, Variable<0> > f1;
  auto f2 = INT(1) / x0;
  // instead of: Expr< Divide< Integer<1>, Variable<0> > f2;
  auto f3 = pow(f1, f2);
  // instead of: Too< Long< Type<...> > > f3;
  cout << "f1 = " << f1 << "\nf2 = " << f2 << "\nf3 = " << f3 << endl;
  return 0;
  // output:
  // f1 = ((5 * x0) + x1)
  // f2 = (1 / x0)
  // f3 = (((5 * x0) + x1))^((1 / x0))
}

int unary_operators()
{
  DVAR(x, 0);
  auto g1 = sin(x);
  // instead of: Expr< Sin_t< Variable<0> > > g1;
  auto g2 = abs(x) * exp(x);
  // instead of: Expr< Times< Abs_t< Variable<0> >, Exp_t< Variable<0> > > > g2;
  cout << "g1 = " << g1 << "\ng2 = " << g2 << endl;
  return 0;
  // output:
  // g1 = sin(x0)
  // g2 = (abs(x0) * exp(x0))
}

int simple_deriv()
{
  auto h = pow(VAR(0), INT(5));
  auto dh = deriv_t(h, VAR(0));
  Array x(1, 2.3);
  cout << "d/dx " << h << " = " << dh << ", which is " << dh.apply(x)
    << " for x0 = " << x[0] << endl;
  return 0;
  // output:
  // d/dx ((x0)^(5)) = (5 * (x0)^(4)), which is 139.92 for x0 = 2.3
}

int full_example()
{
  // Instances can be global, const, static, you don't even need them.
  DVAR(x2, 2);
  // Prepare a function: f1 = x0^2 + sin(t0 / x1) - t1 * cos(x2^2)
  auto f1 = pow(VAR(0), INT(2)) +
    sin(PARAM(0) / VAR(1)) - PARAM(1) * (cos(pow(x2, INT(2))));
  cout << "f1 = " << f1 << endl;

  // Obtain the derivative using a global method, print via operator<<.
  cout << "dx2(f) = " << deriv_t(f1, x2) << endl;

  // Set some values for the unknowns.
  std::vector<SEMT_PRECISION> x_0(3);
  x_0[0] = 1.4141;
  x_0[1] = 2.0021;
  x_0[2] = 2.3;

  // Evaluate our function at this argument.
  cout << "f(" << x_0 << ") = " << f1.apply(x_0) << endl;

  // Evaluate the derivative at the given point.
  cout << "dx0(f)(" << x_0 << ") = " << diff_at(f1, x2, x_0) << endl;

  // The parameter's value is zero, we'll change that now.
  set_parameters<2>::to(x_0);

  // Let's see the changes:
  cout << "set t=3.1415\n"
    << "dx0(f)(" << x_0 << ") = " << diff_at(f1, x2, x_0) << endl;

  // We can store the derivative in an expression:
  auto df1 = deriv_t(f1, x2);

  // And use it like any other expression:
  auto ddf1 = deriv_t(df1, VAR(0));
  cout << "Haha, there's " << abs(ddf1) << " left.\n";

  // output:
  // f1 = (((x0)^(2) + sin((t0 / x1))) - (t1 * cos((x2)^(2))))
  // dx0(f) = (-1 * (t1 * (-1 * (sin((x2)^(2)) * (2 * x2)))))
  // f([ 1.4141, 2.0021,  2.3 ]) = 1.99968
  // dx0(f)([ 1.4141, 2.0021,  2.3 ]) = -0
  // set t=3.1415
  // dx0(f)([ 1.4141, 2.0021,  2.3 ]) = -7.71557
  // Haha, there's abs(0) left.

  // Make some more functions.
  auto f2 = VAR(1) * PARAM(1) - pow(VAR(1), INT(3));
  auto f3 = foldr<Plus, Variable, 0, 2>();

  // We have 3 Variables and want to differentiate once.
  // The expressions are passed via operator=(Expr) and operator,(Expr).
  DifferentiableVectorExpr<3, 1>DVE;
  DVE = f1, f2, f3;

  // Get the derivative as a function object. This is no template,
  // no more differentiation by SEMT means possible.
  const VectorExpr& dF = DVE.get_derivative(1);
  cout << "dF has " << dF.vars() << " vars and "
    << dF.size() << " components.\n";

  // We can print and evaluate VectorExpr, just like the function itself:
  cout << "F = " << DVE.get_function() << endl;
  cout << "D[F] = \n" << dF << endl;
  cout << "D[F](" << x_0 << ") = \n\t" << dF(x_0) << endl;

  // Exercise your CPU.
  SEMT_PRECISION step = 0.001;
  int steps = 100000;
  cout << "Calculating grad[f] at " << steps << " points...";
  cout.flush();
  std::vector<SEMT_PRECISION> res(9);
  for (int i = 0; i < steps; ++i)
  {
    x_0[0] += step;
    x_0[1] += step;
    x_0[2] += step;

    dF.eval(x_0, res);
    // or dF(x_0); if you don't have a vector allocated already.
  }
  cout << " done. \n";

  return 0;

  // output:
  // dF has 3 vars and 9 components.
  // F = [    (((x0)^(2) + sin((t0 / x1))) - (t1 * cos((x2)^(2)))),
  //          ((x1 * t1) - (x1)^(3)),
  //          (x0 + (x1 + x2))    ]
  // D[F] =
  // [    (2 * x0),
  //      (cos((t0 / x1)) * ((-1 * t0) / (x1 * x1))),
  //      (-1 * (t1 * (-1 * (sin((x2)^(2)) * (2 * x2))))),
  //      0,
  //      (t1 - (3 * (x1)^(2))),
  //      0,
  //      1,
  //      1,
  //      1   ]
  // D[F]([ 1.4141, 2.0021,  2.3 ]) =
  //        [ 2.8282, -0.268385, -7.71557,    0, -10.0231,    0,    1,    1,    1 ]
  // Calculating grad[f] at 100000 points... done.
} // end

int fold()
{
  DVAR(x, 0);
  DINT(Five, 5);
  auto r = foldr<Plus>(x, Five);
  auto l = foldl<Minus>(x, Five);
  cout << "r = " << r << "\nl = " << l << endl;
  return 0;
  // output:
  // r = (x0 + (x1 + (x2 + (x3 + (x4 + x5)))))
  // l = (((((x0 - x1) - x2) - x3) - x4) - x5)
}

// xn / (n+1)
template<int n>struct xn_by_npp
{
  typedef Divide<Variable<n>, Integer<n + 1>>simple_type;
};

int iterators()
{
  auto s = foldr<Plus, xn_by_npp, 0, 5>();
  cout << "s = " << s << endl;
  return 0;
  // output:
  // s = (x0 + ((x1 / 2) + ((x2 / 3) + ((x3 / 4) + ((x4 / 5) + (x5 / 6))))))
} // end

  // xn^2
template<int n>struct xn_squared
{
  typedef Power<Variable<n>, 2>simple_type;
};
/*
// x0^2 + ... + x3^2
struct abs_x_squared
{
  typedef Fold_t<Plus,
    IntIterator<xn_squared, 0, 3>,
    true    // Post-order paranthesis.
  >::Result simple_type;
};

// expr / sum(xk^2)
template<class expr>struct normed
{
  typedef Divide<expr, abs_x_squared>simple_type;
};
*/
auto f1 = foldr<Plus, Variable, 0, 3>();
auto f2 = sin(INT(2) * VAR(0)) * VAR(1) + VAR(2) * pow(VAR(0), INT(3));
auto f3 = pow(VAR(0), INT(2)) + INT(2) * VAR(0) * VAR(1) + pow(VAR(2), INT(2));

auto list = (start_tl(f1), f2, f3);
//auto list_normed = for_each<normed>(list);
/*
int for_each()
{
  cout << fold_tl_r<Plus>(list_normed) << endl;
  return 0;
  // output: (it's actually only one line...)
  // (((x0 + (x1 + (x2 + x3))) / ((x0)^(2) + ((x1)^(2) + (x2)^(2)))) +
  // ((((sin((2 * x0)) * x1) + (x2 * (x0)^(3))) / ((x0)^(2) + ((x1)^(2) + (x2)^(2)))) +
  // ((((x0)^(2) + ((2 * x0) * x1)) + (x2)^(2)) / ((x0)^(2) + ((x1)^(2) + (x2)^(2))))))
} // end
*/

using namespace SEMT;

struct pi_struct
{
    constexpr static SEMT_PRECISION value = 3.14159265358979323846;
};
struct inv_sqrt_2pi_struct
{
  constexpr static SEMT_PRECISION value = 0.3989422804014327;
};


Expr<Literal<pi_struct>> pi_value;
Expr<Literal<inv_sqrt_2pi_struct>> inv_sqrt_2pi_value;


void SemtUgmgTest() {
  DVAR(sig2zero, 0);
  DVAR(sig2beta, 1);
  DVAR(pivec, 2);
  DPARAM(z2, 0)
    //RAT
    // auto f = inv_sqrt_2pi_value * POW(sig2beta, -0.5)
    // EPower(sig2beta, Rational()
  auto s = sig2zero + sig2beta * pivec;
  auto f = inv_sqrt_2pi_value * pow(s, RAT(-1, 2)) * exp(RAT(-1, 2) * z2 / s);

  cout << f << endl;
  cout << deriv_t(f, sig2zero) << endl;
  cout << deriv_t(f, sig2beta) << endl;
  cout << deriv_t(f, pivec) << endl;

}

// bgmg-test.exe --gtest_filter=Semt.Test
TEST(Semt, Test) {
  SemtUgmgTest();
  return;
  cout << "These are the examples you find in the documentation.\n";
  cout << "\nintroduction:\n";
  introduction();
  cout << "\nfull example:\n";
  full_example();

  cout << "\nvarious examples:\n";
  macros();
  binary_operators();
  unary_operators();
  simple_deriv();

  cout << "\nadvanced examples:\n";
  fold();
  iterators();
  //for_each();
}

}  // namespace
