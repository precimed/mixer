#ifndef __SEMT_CONSTEXPR__H__
#define __SEMT_CONSTEXPR__H__

#include "Forwards.h"

namespace SEMT
{
/*!
 * @addtogroup  constexpr
 * @todo Write some examples.
 * @{
 */

/*!
 * Compile-time constant integer value.
 * Used extensively across this library, so use it anywhere you can.
 */
SEMT_DEFINE_CONSTEXPR(Integer, int n, n, n*1.0, n);

/// Not very interesting rational type, may be useful for rational exponents.
SEMT_DEFINE_CONSTEXPR(Rational,
        int nom COMMA size_t den,
        nom COMMA den,
        (nom*1.0)/den,
        nom << "/" << den);

/// Bind a type with some static const member \c value to a SEMT-compatible type.
SEMT_DEFINE_CONSTEXPR(Literal, typename T, T, T::value, T::value);

/// Bind a global function to a SEMT-compatible type.
SEMT_DEFINE_CONSTEXPR(Func, SEMT_PRECISION (f)(), f, f(), f());

/// Bind a global function expecting the data for the variables to a SEMT-compatible type.
SEMT_DEFINE_CONSTEXPR(VFunc, SEMT_PRECISION (f)(CAR), f, f(x), " f(x) ");

/// @} // group constexpr

} // namespace SEMT

#endif // __SEMT_CONSTEXPR__H__
