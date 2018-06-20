#ifndef __SEMT_BINARY_OPERATORS__H__
#define __SEMT_BINARY_OPERATORS__H__

#include "Forwards.h"

/*!
 * @addtogroup binops
 * @{
 */

SEMT_DEFINE_BINARYOP(Plus, operator+);
SEMT_DEFINE_BINARYOP(Minus, operator-);
SEMT_DEFINE_BINARYOP(Times, operator*);
SEMT_DEFINE_BINARYOP(Divide, operator/);
SEMT_DEFINE_BINARYOP(EPower, pow);

///@bug This one is actually not intended, but it's necessary for deriv_t() to work with Power.
template<typename LHS, int exp>
inline SEMT::Expr<SEMT::Power<LHS, exp> >
pow(SEMT::Expr<LHS>, SEMT::Expr<SEMT::Integer<exp> >)
{
    return SEMT::Expr<SEMT::Power<LHS, exp> >();
}

/// @}

#endif // __SEMT_BINARY_OPERATORS__H__
