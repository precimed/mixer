/*!
 * @file
 * @todo annotations
 */

#include <iostream>

using namespace std;

#include "semt_speed_func.h"

using namespace SEMT;


//#define LOG(msg)
#ifndef LOG
    #define LOG(msg) std::cout << msg << std::endl;
#endif

//#define DEBUG_LOG(msg)
#ifndef DEBUG_LOG
    #define DEBUG_LOG(msg) std::cout << msg << std::endl;
#endif

#define ERROR_LOG(msg) std::cerr << msg << std::endl;

int main()
{
    const my_semt_vector F1;
    const VectorExpr& f1 = F1.get_fn(0);
    const VectorExpr& G = F1.get_fn(1);
    const VectorExpr& G2 = F1.get_fn(2);

    const size_t vars = f1.vars();

    LOG("F = " << f1 );
    DEBUG_LOG("D F = " << G );
    DEBUG_LOG( "D^2 F = " << G2 );

    Array param(10, 0.123456789);
    Array x_0(vars, 0.123456789);
    Array res(f1.size());
    Array grad(G.size());
    Array grad2(G2.size());

    LOG("\nF(" << x_0 << ") = " << f1(x_0));
    DEBUG_LOG("D F(" << x_0 << ") = \n" << G(x_0));
    DEBUG_LOG("D^2 F(" << x_0 << ") = \n" << G2(x_0));

    int steps = 100000;
    double step = 0.5/steps;

    LOG("\nCalculating F, DF, D^2 F at " << steps << " points... ");

#if not SEMT_DISABLE_PRINT
    cout.flush();
#endif

    F1.set_parameters(param);
    size_t k = 0; //used in the inner loop.
    for (int i = 0; i < steps; ++i)
    {
        for(k = 0; k < vars; ++k)
            x_0[k] += step;

        f1.eval(x_0,res);
        G.eval(x_0, grad);
        G2.eval(x_0, grad2);
    }
    LOG(" done.");

    DEBUG_LOG( "F(" << x_0 << ") = \n" << f1(x_0) );
    DEBUG_LOG( "D F(" << x_0 << ") = \n" << G(x_0) );
    DEBUG_LOG( "D^2 F(" << x_0 << ") = \n" << G2(x_0) );
}
