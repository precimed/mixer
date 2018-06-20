#ifndef __SEMT_BINARY_TYPES__H__
#define __SEMT_BINARY_TYPES__H__

#include "Forwards.h"

namespace SEMT
{

/*!
 * @addtogroup bintypes
 * @{
 */

/// Add two expressions.
SEMT_DEFINE_BINARY_TYPE(Plus,
        (First::apply(x) + Second::apply(x)),
        Plus<First_d COMMA Second_d>,
        "(" << First::toString() << " + " << Second::toString() << ")");

/// Substract two expressions.
SEMT_DEFINE_BINARY_TYPE(Minus,
        (First::apply(x) - Second::apply(x)),
        Minus<First_d COMMA Second_d>,
        "(" << First::toString() << " - " << Second::toString() << ")");

/// Multiply two expressions.
SEMT_DEFINE_BINARY_TYPE(Times,
        (First::apply(x) * Second::apply(x)),
        Plus<Times<First_d COMMA Second> COMMA Times<First COMMA Second_d> >,
        "(" << First::toString() << " * " << Second::toString() << ")");

/// The quotient of two expressions.
SEMT_DEFINE_BINARY_TYPE(Divide,
        (First::apply(x) / Second::apply(x)),
        Divide<Minus<Times<First_d COMMA Second> COMMA Times<First COMMA Second_d> > COMMA Times<Second COMMA Second> >,
        "(" << First::toString() << " / " << Second::toString() << ")")

/// Raise power.
SEMT_DEFINE_BINARY_TYPE(EPower,
        pow(First::apply(x),Second::apply(x)),
        Times<EPower<First COMMA Second> COMMA Plus<Times<Second_d COMMA Ln_t<First> > COMMA Divide<Times<Second COMMA First_d> COMMA First > > >,
        "(" << First::toString() << ")^(" << Second::toString() << ")")

/// @}

namespace intern
{
/*!
 * Calculate base ^ expo by shifting the exponent at compile-time.
 * Returns itself instantiated with shifted power and squared base,
 * multiplied by base if last bit of expo is set.
 * You may consider using this for parameter calculations etc.
 * @tparam  epxo    A static, positive integer.
 * @param   base    Is used at run-time.
 * @return  base ^ expo
 */
template<size_t expo> SEMT_INLINE SEMT_PRECISION PowerHelper(SEMT_PRECISION base)
{
    // if the last bit is set, we need to multiply the result with this base,
    // otherwise we just return the next helper applied with the squared base
    return ((expo & 1) ? (base * PowerHelper<(expo >> 1)>(base * base))
                       : (PowerHelper<(expo >> 1)>(base * base)));
}

/// Specialization for base ^ 1, the recursive instantiation stops.
template<> SEMT_INLINE SEMT_PRECISION PowerHelper<1>(SEMT_PRECISION base)
{
    return base;
}

/// Specialization for base ^ 0, needed if simplifications are disabled.
template<> SEMT_INLINE SEMT_PRECISION PowerHelper<0>(SEMT_PRECISION base)
{
    return 1.0;
}
}

/*!
 * Compile-time fixed exponent power type.
 * @ingroup bintypes
 * @sa      intern::PowerHelper
 */
template<typename expr, int expo>
struct Power
{
    typedef SEMT_SIMPLE_TYPE(expr) First;

    static const int LastVar = First::LastVar;

    template<typename var> struct partial
    {
        static const bool dependent = First::template partial<var>::dependent;
        typedef SEMT_SIMPLE_TYPE2(First::template partial<var>::deriv) First_d;
        typedef SEMT_SIMPLE_TYPE(Times<Times<Integer<expo> COMMA Power<First COMMA expo - 1> > COMMA First_d>) deriv;
    };

    static SEMT_INLINE_STATIC SEMT_PRECISION apply(CAR x)
    {
        return ((expo >= 0) ? intern::PowerHelper<expo>(First::apply(x))
                : (1.0 / intern::PowerHelper< -expo>(First::apply(x))));
    }

    SEMT_DEFINE_TOSTRING("(" << First::toString() << ")^(" << expo << ")");
};

} // namespace SEMT

#endif // __SEMT_BINARY_TYPES__H__
