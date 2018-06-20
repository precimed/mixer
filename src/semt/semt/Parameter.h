#ifndef __SEMT_PARAMETER__H__
#define __SEMT_PARAMETER__H__

#include "Forwards.h"

namespace SEMT
{

/*!
 * Use this for inputs for which no derivatives are required.
 * @note They are initialized with 0.0. Use set_parameters to change all parameters at once.
 * @attention   This is definitly not thread safe. At least you have to make sure that
 *              all involved parameters are updated before the an expression is evaluated.
 * @tparam index Counting from zero.
 * @tparam name Divide parameters into groups.
 * @ingroup nodes
 */
template<int index, char name>
class Parameter
{
public:
    typedef Parameter<index> simple_type;

    static const int LastVar = 0;

    template<typename var> struct partial
    {
        static const bool dependent = false;
        typedef Integer<0> deriv;
    };

    /// Assign a new value to this parameter only.
    static void set_value(SEMT_PRECISION new_val)
    {
        value = new_val;
    }

    static SEMT_PRECISION get_value()
    {
        return value;
    }

    static SEMT_PRECISION apply(CAR x0)
    {
        return value;
    }

    SEMT_DEFINE_TOSTRING(name << index)

private:
    static SEMT_PRECISION value;
};

template<int index, char name>
SEMT_PRECISION Parameter<index, name>::value = 0.0;

/*!
 * Update all parameters with a specific identifier.
 * @tparam  name    Character the parameters are identifed by.
 * @tparam  count   Number of parameters to be updated.
 * @attention   Array bounds are not checked for implementation reasons (we don't know
 *              if we're in the first instantiation).
 * @ingroup api
 */
template<int count, char name> struct set_parameters
{
    /*!
     * Assign <tt>t[i]</tt> to <tt>Parameter<i, name></tt> for <tt>0 <= i < count</tt>.
     * @param   t       Array with the new values, must fullfil size() >= count.
     */
    static void to(CAR t)
    {
        Parameter<count - 1, name>::set_value(t[count - 1]);
        set_parameters<count - 1, name>::to(t);
    }
};

/// Stop recursive update of parameters.
template<char name> struct set_parameters<1, name>
{
    static void to(CAR t)
    {
        Parameter<0, name>::set_value(t[0]);
    }
};

} // namespace SEMT

#endif // __SEMT_PARAMETER__H__
