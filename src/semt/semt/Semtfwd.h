/*!
 * SEMT forward header file, users of a function defined via SEMT should include this file only.
 * All external headers used across SEMT are included here.\n
 * @file
 * @ingroup advapi
 */

#ifndef __SEMT_SEMTFWD__H__
#define __SEMT_SEMTFWD__H__

/*!
 * @addtogroup settings
 * @{
 */

/*!
 * The used value type.
 * If it does not support the standard library's math functions, suitable alternatives
 * must be available. In this case, you may wish to redefine SEMT_USE_STDMATH.
 */
#ifndef SEMT_PRECISION
    typedef float SEMT_PRECISION;
#endif

/*!
 * Include standard math header.
 */
#ifndef SEMT_USE_STDMATH
    #define SEMT_USE_STDMATH 1
#endif

/*!
 * Use std::vector as Array.
 * If you change this, you have to redefine Array
 */
#ifndef SEMT_USE_STD_VECTOR
    #define SEMT_USE_STD_VECTOR 1
#endif

/*!
 * Threshold for numeric zero.
 */
#ifndef SEMT_EPS
    #define SEMT_EPS 2.22045e-16
#endif

/*!
 * Throw std::runtime_error on forbidden evaluations (like differentiate abs at zero).
 * @see SEMT::Defined_if
 */
#ifndef SEMT_DISABLE_CHECKERS
    #define SEMT_DISABLE_CHECKERS 0
#endif

/*!
 * Define the Precision of the generated streams.
 */
#ifndef SEMT_STRPREC
    #define SEMT_STRPREC 6
#endif
#define SEMT_STRLEN (SEMT_STRPREC + 9)

/*!
 * Inline policies: methods with no qualifier will be inlined or not.
 */
#ifndef SEMT_INLINE
    #define SEMT_INLINE inline
#endif
/*!
 * Inline policies: methods with static qualifier will be inline or not.
 */
#ifndef SEMT_INLINE_STATIC
    #define SEMT_INLINE_STATIC inline
#endif
/*!
 * Inline policies: methods with virtual qualifier will be inline or not.
 */
#ifndef SEMT_INLINE_VIRTUAL
    #define SEMT_INLINE_VIRTUAL inline
#endif

/*!
 * Define if binary simplifications should be disabled.
 */
#ifndef SEMT_DISABLE_BINARY_SIMPLIFICATIONS
    #define SEMT_DISABLE_BINARY_SIMPLIFICATIONS 0
#endif
/*!
 * Define if unary simplifications should be disabled.
 */
#ifndef SEMT_DISABLE_UNARY_SIMPLIFICATIONS
    #define SEMT_DISABLE_UNARY_SIMPLIFICATIONS 0
#endif
/*!
 * Define if more complex simplifications should be disabled.
 */
#ifndef SEMT_DISABLE_COMPLEX_SIMPLIFICATIONS
    #define SEMT_DISABLE_COMPLEX_SIMPLIFICATIONS 0
#endif

/*!
 * Disable all simplifications at once.
 */
#ifndef SEMT_DISABLE_ALL_SIMPLIFICATIONS
    #define SEMT_DISABLE_ALL_SIMPLIFICATIONS 0
#endif

/*!
 * @} // group settings
 */

/* Hold the promise and disable them all, #undef gives +1 reputation with the compiler.*/
#if SEMT_DISABLE_ALL_SIMPLIFICATIONS
    #undef SEMT_DISABLE_BINARY_SIMPLIFICATIONS
    #undef SEMT_DISABLE_UNARY_SIMPLIFICATIONS
    #undef SEMT_DISABLE_COMPLEX_SIMPLIFICATIONS
    #define SEMT_DISABLE_BINARY_SIMPLIFICATIONS     1
    #define SEMT_DISABLE_UNARY_SIMPLIFICATIONS      1
    #define SEMT_DISABLE_COMPLEX_SIMPLIFICATIONS    1
#endif


#include <cassert>
#include <memory>
#include <stdexcept>

#if SEMT_USE_STD_VECTOR
    #include <vector>
#endif

#if SEMT_USE_STDMATH
    #include <cmath>
#endif

    #include <string>
    #include <sstream>
    #include <iomanip>

    #if SEMT_USE_STD_VECTOR

    /*!
     * Pass a vector of T's to a stream.
     * @tparam T The value type in \c v.
     * @param   os  Output stream.
     * @param   v   Values to print.
     * @return  Reference to \c os.
     */
    template<class T>
    std::ostream& operator<<(std::ostream &os, const std::vector<T>& v)
    {
        size_t l = v.size();
        os << "[ ";
        for (size_t i = 1; i < l; i++)
            os << std::setw(SEMT_STRLEN) << v[i - 1] << ", ";
        os << std::setw(SEMT_STRLEN) << v[l - 1] << " ]";
        return os;
    }

    /*!
     * Split a vector in rows & columns.
     * @tparam rows How many rows to layout.
     * @tparam cols How many columns to layout.
     * @param ex The vector to be formatted as a (rows, cols)-matrix.
     * @return The formatted string.
     * @todo rows could be exchanged for ex.size, would make it a little more flexible.
     */
    template<size_t rows, size_t cols>
    std::string matrix_str(const std::vector<SEMT_PRECISION>& ex)
    {
        size_t values = ex.size();
        if (values != (rows * cols))
            throw std::runtime_error("invalid size in matrix_str.");

        std::ostringstream os;
        os.precision(SEMT_STRPREC);
        for (size_t i = 0; i < rows; ++i)
        {
            os << "| ";
            for (size_t j = 0; j < cols - 1; ++j)
            {
                os << std::setw(SEMT_STRLEN) << ex[i * cols + j] << ", ";
            }
            os << std::setw(SEMT_STRLEN) << ex[(i + 1) * cols - 1] << " |\n";
        }
        return os.str();
    }
    #endif

    #include "loki/NullType.h"
    #include "loki/TypeManip.h"
    #include "loki/Typelist.h"

    namespace SEMT
    {

    /*!
     * Convenient typedefs, the underlying types may change in the future.
     * Must supply SEMT_PRECISION operator[](size_t).
     * @ingroup api
     * @{
     */
    #if SEMT_USE_STD_VECTOR
    typedef std::vector<SEMT_PRECISION> Array;
    #endif
    typedef const Array& CAR;
    /*!
     *   @}
     */

    /*!
     * The abstract base-class for atomic function objects resulting from SEMT operations,
     * used by VectorExpr to store its components.
     * @ingroup dev
     */
    struct AbstractFunctionObject
    {
        virtual ~AbstractFunctionObject()
        {
        }
        virtual SEMT_PRECISION operator()(CAR x) const = 0;
        virtual std::auto_ptr<AbstractFunctionObject> duplicate() const = 0;

        virtual std::string toString() const = 0;
    };

} // namespace SEMT

#endif // __SEMT_SEMTFWD__H__
