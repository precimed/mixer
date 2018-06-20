#ifndef __SEMT_VARIABLE__H__
#define __SEMT_VARIABLE__H__

#include "Forwards.h"

namespace SEMT
{

/*!
 * Represents the unkown in an expression.
 * @attention   It's virtually a restricion to begin with index 0 and count upwards, because
 *              apply just returns v[index] for a CAR v.
 * @tparam  index   The index of this variable, counting from zero.
 * @ingroup nodes
 */
template<int index> struct Variable
{
    typedef Variable<index> simple_type;

    static const int LastVar = index;

    template<typename var> struct partial
    {
        static const bool dependent = (Loki::IsSameType<simple_type, var>::value) ? true : false;
        typedef typename Loki::Select<Loki::IsSameType<simple_type, var>::value,
                Integer<1>,
                Integer<0>
                >::Result deriv;
    };

    static SEMT_INLINE_STATIC SEMT_PRECISION apply(CAR v)
    {
        return v[index];
    }

    SEMT_DEFINE_TOSTRING('x' << index)
};

} // namespace SEMT

#endif // __SEMT_VARIABLE__H__
