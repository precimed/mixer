#ifndef __SEMT_CONDITIONS__H__
#define __SEMT_CONDITIONS__H__

#include "Forwards.h"

namespace SEMT
{

/*!
 * @addtogroup cond
 * @{
 */

/// Check if less then SEMT_EPS.
SEMT_DEFINE_CONDITION(LessThanZero,
        x < -SEMT_EPS,
        "< 0");

/// Check if greater then SEMT_EPS.
SEMT_DEFINE_CONDITION(GreaterThanZero,
        x > SEMT_EPS,
        "> 0");

/// Check if absolute value is less then SEMT_EPS.
SEMT_DEFINE_CONDITION(AlmostZero,
        std::fabs(x) < SEMT_EPS,
        "== 0");

/// Check if absolute value is greater then SEMT_EPS.
SEMT_DEFINE_CONDITION(NotZero,
        std::fabs(x) > SEMT_EPS,
        "!= 0");

///@}

/*!
 * Throw a runtime error if Condition is not satisfied for cond_arg_ on evaluation.
 * Because the derivative of this type is
 * this type again, all derivatives are protected against the offending evaluation, too.
 * An example is Sgn_t, where this is used to throw a runtime_error if Sgn_t is differentiated at zero.
 * Of course, evaluate to zero is not to be taken literally, see SEMT_EPS.
 * @note    You can disable the whole mechanism by setting SEMT_DISABLE_CHECKERS to 1.
 * @tparam  expr    Expression that needs restriction.
 * @tparam  Condition   A conditional type, supplying \c bool is_true_for(SEMT_PRECISION) and <tt> std::string toString()</tt>.
 * @tparam  cond_arg_   Expression that must fulfill Condition.
 * @ingroup cond
 * @ingroup untypes
 */
template<class expr, class Condition, class cond_arg_>
struct Defined_if
#if SEMT_DISABLE_CHECKERS
: public SEMT_SIMPLE_TYPE(expr)
{
    typedef SEMT_SIMPLE_TYPE(expr) First;
    typedef First simple_type;
#else
{
    typedef SEMT_SIMPLE_TYPE(expr) First;
    typedef SEMT_SIMPLE_TYPE(cond_arg_) cond_arg;
    typedef Defined_if<First, Condition, cond_arg> simple_type;

    static const int LastVar = First::LastVar;

    template<class var> struct partial
    {
        // We may change, if cond_arg changes.
        static const bool dependent = (First::template partial<var>::dependent
                || cond_arg::template partial<var>::dependent);
        typedef typename First::template partial<var>::deriv First_d;

        // If our expression does not depend on this variable, differentiate straight to zero,
        // otherwise just differentiate and preserve checking.
        typedef Defined_if<First_d, Condition, cond_arg> First_d_checked;
        typedef typename Loki::Select<dependent, First_d_checked, Integer<0> >::Result deriv;
    };

    static SEMT_INLINE_STATIC SEMT_PRECISION apply(CAR x)
    {
        if (!Condition::is_true_for(cond_arg::apply(x)))
        {
            std::stringstream msg;
            msg << "Condition " << cond_arg::toString() << " " << Condition::toString()
                    << ((!Loki::IsSameType<First, cond_arg>::value) ?
                            "" : " not satisfied for " + First::toString()) << " at x = " << x;
            throw std::runtime_error(msg.str());
        }
        return First::apply(x);
    }

    // Don't print the expression twice.
    SEMT_DEFINE_TOSTRING("{ " << First::toString() << ", if "
            << ((Loki::IsSameType<First, cond_arg>::value) ? "" : (cond_arg::toString() + " "))
            << Condition::toString() <<" }")
    ;
#endif // SEMT_DISABLE_CHECKERS
};

} // namespace SEMT

#endif // __SEMT_CONDITIONS__H__
