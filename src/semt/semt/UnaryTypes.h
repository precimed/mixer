#ifndef __SEMT_UNARY_TYPES__H__
#define __SEMT_UNARY_TYPES__H__

#include "Forwards.h"

namespace SEMT
{

/*!
 * @addtogroup  untypes
 * @{
 */

/// Directions, please?
template<typename expr>
struct Sgn_t
{
    typedef SEMT_SIMPLE_TYPE(expr) First;
    typedef Sgn_t<First> simple_type;

    static const int FirstVar = First::FirstVar;

    template<typename var> struct partial
    {
        static const bool dependent = First::template partial<var>::dependent;
        // This is zero, but only if First is not zero or this var isn't part of this expression.
        typedef typename Loki::Select<dependent,
                Defined_if<Integer<0>, NotZero, First>,
                Integer<0>
                >::Result deriv;
    };

    static SEMT_INLINE_STATIC SEMT_PRECISION apply(CAR x)
    {
        const SEMT_PRECISION f_value = First::apply(x);
        return (f_value > SEMT_EPS) ? 1.0 : ((f_value < -SEMT_EPS) ? -1.0 : 0.0);
    }

    SEMT_DEFINE_TOSTRING("sgn(" << First::toString() << ')')
};

/*!
 * Absolute value.
 * We use abs(f) = sgn(f) * f and get: abs'(f) = sgn'(f)*f + sgn(f)*f'.
 * Since sgn'(f) is zero unless sgn(f) == 0, we place the checker around sgn(f).
 * @todo Simplify to Zero if Zero argument.
 */
SEMT_DEFINE_UNARY_TYPE(Abs_t,
        std::fabs(First::apply(x)),
        Times<Defined_if<Sgn_t<First> COMMA NotZero > COMMA First_d>,
        "abs(" << First::toString() << ")");

/// Exponential function.
SEMT_DEFINE_UNARY_TYPE(Exp_t,
        exp(First::apply(x)),
        Times<Exp_t<First> COMMA First_d>,
        "exp(" << First::toString() << ")");

/// Natural logarithm.
SEMT_DEFINE_UNARY_TYPE(Ln_t,
        log(First::apply(x)),
        Divide<First_d COMMA First >,
        "ln(" << First::toString() << ")");

/// Sine.
SEMT_DEFINE_UNARY_TYPE(Sin_t,
        sin(First::apply(x)),
        Times<Cos_t<First> COMMA First_d>,
        "sin(" << First::toString() << ")");

/// Inverse sine.
SEMT_DEFINE_UNARY_TYPE(Arcsin_t,
        asin(First::apply(x)),
        Divide<First_d COMMA EPower<Minus<Integer<1> COMMA Power<First COMMA 2 > > COMMA Rational<1 COMMA 2> > >,
        "arcsin(" << First::toString() << ")");

/// Cosine.
SEMT_DEFINE_UNARY_TYPE(Cos_t,
        cos(First::apply(x)),
        Times<Integer<-1> COMMA Times<Sin_t<First> COMMA First_d> >,
        "cos(" << First::toString() << ")");

/// Inverse cosine.
SEMT_DEFINE_UNARY_TYPE(Arccos_t,
        acos(First::apply(x)),
        Divide<Times<Integer<-1> COMMA First_d> COMMA EPower<Minus<Integer<1> COMMA Power<First COMMA 2 > > COMMA Rational<1 COMMA 2> > >,
        "arccos(" << First::toString() << ")");

/// Tangent.
SEMT_DEFINE_UNARY_TYPE(Tan_t,
        tan(First::apply(x)),
        Times<First_d COMMA Plus<Integer<1> COMMA Power<Tan_t<First> COMMA 2> > >,
        "tan(" << First::toString() << ")");

/// Inverse tangent.
SEMT_DEFINE_UNARY_TYPE(Arctan_t,
        atan(First::apply(x)),
        Divide<First_d COMMA Plus<Integer<1> COMMA Power<First COMMA 2 > > >,
        "arctan(" << First::toString() << ")");

/// Hyperbolic sine.
SEMT_DEFINE_UNARY_TYPE(Sinh_t,
        sinh(First::apply(x)),
        Times<Cosh_t<First> COMMA First_d>,
        "sinh(" << First::toString() << ")");

/// Inverse hyperbolic sine.
SEMT_DEFINE_UNARY_TYPE(Arsinh_t,
        asinh(First::apply(x)),
        Divide<First_d COMMA EPower<Plus<Integer<1> COMMA Power<First COMMA 2 > > COMMA Rational<1 COMMA 2> > >,
        "arsinh(" << First::toString() << ")");

/// Hyperbolic cosine.
SEMT_DEFINE_UNARY_TYPE(Cosh_t,
        cosh(First::apply(x)),
        Times<Sinh_t<First> COMMA First_d>,
        "cosh(" << First::toString() << ")");

/// Inverse hyperbolic cosine.
SEMT_DEFINE_UNARY_TYPE(Arcosh_t,
        acosh(First::apply(x)),
        Divide<First_d COMMA EPower<Minus<Power<First COMMA 2 > COMMA Integer<1> > COMMA Rational<1 COMMA 2> > >,
        "arcosh(" << First::toString() << ")");

/// Hyperbolic tangent.
SEMT_DEFINE_UNARY_TYPE(Tanh_t,
        tanh(First::apply(x)),
        Times<First_d COMMA Minus<Integer<1> COMMA Power<Tanh_t<First> COMMA 2> > >,
        "tanh(" << First::toString() << ")");

/// Inverse hyperbolic tangent.
SEMT_DEFINE_UNARY_TYPE(Artanh_t,
        atanh(First::apply(x)),
        Divide<First_d COMMA Minus<Integer<1> COMMA Power<First COMMA 2 > > >,
        "artanh(" << First::toString() << ")");

/// @} // group untypes

} // namespace SEMT

#endif // __SEMT_UNARY_TYPES__H__
