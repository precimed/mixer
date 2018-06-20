#include "semt_func_impl.h"
#include "semt/Semt.h"
#include "semt/Common.h"

double global_var = 7.444;

double alphaf()
{
    return global_var;
}
struct alpha
{
    constexpr static SEMT_PRECISION value = 5.234;
};

using namespace SEMT;

Expr<Literal<alpha>>a;
Expr<Func<alphaf>>af;

auto g1 = Quarter * x0 + a * x1;
auto g2 = pow(x2, Two) - Quarter;
auto g3 = (x0 + x2);

auto g4 = sin(x3);
auto g5 = cos(x4);
auto g6 = exp(x5) - Ten;

auto g7 = sinh((x6)) - Ten;
auto g8 = cosh((x7)) - Ten;
auto g9 = pow(x8, Ten) - Quarter;
auto g10 = x9 + ln(x9);

my_semt_func::my_semt_func()
{
    DifferentiableVectorExpr<10, 1>DVE;
    DVE = g1, g2, g3, g4, g5, g6, g7, g8, g9, g10;
    bunch_of_functions.push_back(DVE.get_function());
    bunch_of_functions.push_back(DVE.get_derivative(1));
}

