/*!
 * @file
 * @todo annotations
 */

#include "semt_speed_func.h"
#include "semt/Semt.h"
#include "semt/Common.h"

using namespace SEMT;

void my_semt_vector::set_parameters(CAR new_values) const
{
    SEMT::set_parameters<10>::to(new_values);
}

template<int i> struct iter1
{
    typedef Times<Parameter<i> , Divide< Power<Variable<i>, i + 1> , Integer<i + 1> > > simple_type;
};

template<int i> struct iter2
{
    typedef Minus< Power<Variable<i>, i + 1> , Integer<i + 1> > simple_type;
};

auto f1 = foldr<Plus, iter1, 0, 5>();
auto f2 = PARAM(0) * VAR(9) + foldr<Times, iter2, 0, 5>();
auto f3 = f1 * (abs(x0 + PARAM(0)) + sgn(x1 + PARAM(1)) + exp(x2) + ln(abs(x3) + PARAM(3)));
auto f4 = pow(sin(x4) * (tan(PARAM(4)) + sinh(x5)), PARAM(5));
auto f5 = (acos(x0) + atan(x1) - PARAM(7) * acosh(x2)) + f2;

my_semt_vector::my_semt_vector()
{
    DifferentiableVectorExpr<6, 2> Fn;
    Fn = f1, f2, f3, f4, f5;
    functions_.push_back(Fn.get_function());
    functions_.push_back(Fn.get_derivative(1));
    functions_.push_back(Fn.get_derivative(2));
}


/*
auto f1 = PARAM(0) * x0 + x1 + x2;
auto f2 = (f1 * f1);
auto f3 = pow(f1, Ten);
auto f4 = Quarter * x0 - PARAM(1) * (pow(x1, Two) + pow(x2, mOne));
auto f5 = pow(sin(x0 / x1), f4);

my_semt_vector::my_semt_vector()
{
    DifferentiableVectorExpr<3, 2> Fn = f1, f2, f3, f4, f5;
    functions_.push_back(Fn.get_function());
    functions_.push_back(Fn.get_derivative(1));
    functions_.push_back(Fn.get_derivative(2));
}

*/
