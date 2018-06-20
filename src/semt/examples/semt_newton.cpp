#include <cassert>
#include <cmath>
#include <limits>

#include <iostream>
#include <sstream>

#include "semt_func_impl.h"

using namespace std;
using namespace SEMT;

#define DEBUG_LOG(msg)
#ifndef DEBUG_LOG
#define DEBUG_LOG(msg) std::cout << msg << std::endl;
#endif

#define ERROR_LOG(msg) std::cerr << msg << std::endl;

template<size_t dim>
size_t search_pivot_column(CAR A, size_t column)
{
    numeric_limits<SEMT_PRECISION>real_info;
    size_t j = 0;
    size_t p = column;

    for (j = column + 1; j < dim; ++j)
    {
        if (fabs(A[p * dim + column]) < fabs(A[j * dim + column]))
            p = j;
    }
    if (fabs(A[p * dim + column]) < (real_info.epsilon()))
        throw "singular matrix";
    return p;
}

template<size_t dim>
int LRDecomp(Array& A, size_t * ipiv)
{
    size_t i, j, p;
    // Initialize pivot vector.
    for (i = 0; i < dim; ++i)
        ipiv[i] = i;

    for (i = 0; i < dim; ++i)
    {
        try
        {
            p = search_pivot_column<dim>(A, i);
        }
        catch (char*)
        {
            ERROR_LOG("\tMatrix numerically singular.");
            return 1;
        }

        //swap rows
        if (p != i)
        {
            swap_ranges(A.begin() + (i * dim), A.begin() + (i * dim + dim), A.begin() + (p * dim));
            swap(ipiv[i], ipiv[p]);
        }

        // eliminate column
        for (p = i + 1; p < dim; ++p)
        {
            A[p * dim + i] /= A[i * dim + i];
            for (j = i + 1; j < dim; ++j)
            {
                A[p * dim + j] -= A[p * dim + i] * A[i * dim + j];
            }
        }
    }
    return 0;
}

template<size_t dim>
void LRSolve(CAR LR, CAR a, Array& y, const size_t * const ipiv)
{
    //forward substitution: solve Ly = a
    for (size_t i = 0; i < dim; ++i)
    {
        y[i] = a[ipiv[i]];
        for (size_t j = 0; j < i; ++j)
        {
            y[i] -= (LR[i * dim + j] * y[j]);
        }
    }
    //backward substitution: solve Ry = y
    for (int i = dim - 1; i >= 0; --i)
    {
        for (int j = dim - 1; j > i; --j)
            y[i] -= (LR[i * dim + j] * y[j]);
        y[i] /= LR[i * dim + i];
    }
}

template<size_t dim>inline long double norm2(CAR x)
{
    return norm2<dim - 1>(x) + x[dim - 1] * x[dim - 1];
}
template<>inline long double norm2<0>(CAR x)
{
    return 0;
}

template<size_t dim>
int newtons_method(const VectorExpr& f, const VectorExpr& df, Array& x)
{
    assert(x.size() == dim);
    assert(f.size() == dim);
    assert(df.size() == dim * dim);

    double residual = 0.0;

    Array fv(f(x)); //initialize f at x
    Array dfv(df(x)); //initialize df at x
    Array tmp(dim, 0);

    cout << "Starting newton iteration with initial guess: " << x << endl;
    size_t ipiv[dim];
    for (int it = 0;; ++it)
    {
        //solve Df * y = f
        if (LRDecomp<dim>(dfv, ipiv))
        {
            cerr << "encountered singular jacobian in iteration " << it << "\nfor x = " << x
                    << endl;
            return 1;
        }
        LRSolve<dim>(dfv, fv, tmp, ipiv);
        DEBUG_LOG( "iteration " << it << "\tincrement = " << tmp );

        //apply Newton-correction
        for (size_t i = 0; i < dim; ++i)
            x[i] -= tmp[i];
        DEBUG_LOG( "\tcurrent 'solution' = " << x );

        //calculate new value
        f.eval(x, fv);
        residual = norm2<dim>(fv);
        DEBUG_LOG( "\tresidual^2 = " << residual );

        //exit loop if it overflows or residual is very small
        if ((it == 1000))
        {
            ERROR_LOG("Leaving after " << it << " iterations. \n\n");
            return 1;
        }
        if (residual < SEMT_EPS)
        {
            cout << "Leaving with residual^2 = " << residual << " after " << it
                    << " iterations.\n\n";
            return 0;
        }

        df.eval(x, dfv);
        if (!(it % 5))
            DEBUG_LOG( "Jacobian:\n "<< (matrix_str<dim , dim>(dfv)) << "\n");

    }
    return 1;
}

int main()
{
    cout.setf(ios::scientific, ios::floatfield);
    cout.precision(9);

    const size_t dim = 10;
    std::vector<double> x_0(dim);
    std::vector<double> res(dim);
    std::vector<double> grad(dim * dim);

    for (size_t i = 0; i < dim; ++i)
        x_0[i] = (rand() % 10000) * 0.001;

    cout << "\n\n Newtons-method:\n";
    const my_semt_func F;
    const VectorExpr& f = F.get_func(0);
    const VectorExpr& df = F.get_func(1);
    cout << "f = " << f << endl;
    cout << "Df = " << (matrix_str<dim, dim>(df)) << endl;
    if (!newtons_method<dim>(f, df, x_0))
    {
        cout << "Root of f found at \n  x  = " << x_0 << endl;
        cout << "f(x) = " << f(x_0) << endl;
    }
    else
        cout << "No solution found." << endl;
    return 0;
}
