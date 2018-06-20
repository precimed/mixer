#ifndef __SEMT_EXPRESSION__H__
#define __SEMT_EXPRESSION__H__

#include "Forwards.h"

namespace SEMT
{

/*!
 * A proxy for complex %SEMT types, example for the static interface that all unary, binary and node types conform to.
 * This type is wrapped around all other SEMT types,
 * everything is just forwarded to it's template parameter.
 * @struct Expr
 * @ingroup advapi
 * @see nodes
 * @see untypes
 * @see bintypes
 */
template<typename expr>
struct Expr
{
    /// First (and Second) are the argument types of a unary (or binary) operator.
    typedef typename RecursiveSimplifier<expr>::Result First;

    /// Used to get bounds for evaluation.
    static const int LastVar = First::LastVar;

    /// This is the heart, \c var is Variable\<i\> for an integer i, we just forward here.
    template<typename var> struct partial
    {
        /// Indicates if var has any influence in this expression.
        static const bool dependent = First::template partial<var>::dependent;
        /// Contains the derivative in respect to var.
        typedef typename First::template partial<var>::deriv deriv;
    };
    /// Calculation is done by recursion.
    static SEMT_INLINE_STATIC SEMT_PRECISION apply(CAR x)
    {
        return First::apply(x);
    }

    /// Used by operator<<.
    SEMT_DEFINE_TOSTRING(First::toString())
};

/*!
 * @addtogroup api
 * @{
 */

/*!
 * Evaluate a partial derivative at a given point.
 * @tparam  expr    The function to be differentiated.
 * @tparam  V       Differentiate for the Variable with this index.
 * @param   x       The value vector for the variables present in expr.
 * @return  The value of the partial derivative of expr with respect to Variable<V> at x.
 */
template<typename expr, int V>
SEMT_PRECISION diff_at(Expr<expr>, Expr<Variable<V> >, CAR x)
{
    return Expr<typename expr::template partial<Variable<V> >::deriv>::apply(x);
}

/*!
 * Give the expression of a partial derivative.
 * @tparam  expr    The function to be differentiated.
 * @tparam  V       Differentiate for the Variable with this index.
 * @return  The expression of the derivative.
 * @par Usage
 *      @dontinclude semt_examples.cpp
 *      @skip int simple_deriv()
 *      @until }
 */
template<typename expr, int V>
Expr<typename expr::template partial<Variable<V> >::deriv>
deriv_t(Expr<expr>, Expr<Variable<V> > = Expr<Variable<V> >())
{
    return Expr<typename expr::template partial<Variable<V> >::deriv>();
}

/// @} // group api

}// namesapce SEMT

/*!
 * Pass an expression as a std::string to std::ostream.
 * @fn std::ostream& operator<<(std::ostream &, SEMT::Expr)
 * @tparam  expr    The function.
 * @param   os      The std::ostream reference that will receive the string.
 * @return  The modified std::ostream.
 * @ingroup api
 */
template<typename expr>
std::ostream& operator<<(std::ostream &os, SEMT::Expr<expr> ex)
{
    os << ex.toString();
    return os;
}

#endif // __SEMT_EXPRESSION__H__
