#ifndef __SEMT_TRAITS__H__
#define __SEMT_TRAITS__H__

#include "Forwards.h"

namespace SEMT
{

/*!
 * This is from Boost, see boost.org.
 * @{
 */
template<class T, T val>
struct integral_constant
{
    typedef integral_constant<T, val> type;
    typedef T value_type;
    static const T value = val;
};

typedef integral_constant<bool, true> true_type;
typedef integral_constant<bool, false> false_type;

/// @}

/*!
 * This is Boost::enable_if, see boost.org.
 * @{
 */
template<bool B, class T = void>
struct enable_if_c
{
    typedef T type;
};

template<class T>
struct enable_if_c<false, T>
{
};

template<class Cond, class T = void>
struct enable_if : public enable_if_c<Cond::value, T>
{
};

/// @}

/*!
 * This is Boost::disable_if, see boost.org.
 * @{
 */
template<bool B, class T = void>
struct disable_if_c
{
    typedef T type;
};

template<class T>
struct disable_if_c<true, T>
{
};

template<class Cond, class T = void>
struct disable_if : public disable_if_c<Cond::value, T>
{
};

/// @}

/*!
 * @addtogroup nodetraits
 * @{
 */

SEMT_DEFINE_TRAIT(isVar, int i, Variable<i>)
SEMT_DEFINE_TRAIT(isParam, int i, Parameter<i>)
SEMT_DEFINE_TRAIT(isInt, int i, Integer<i>)
SEMT_DEFINE_TRAIT(isZero, , Integer<0>)
SEMT_DEFINE_TRAIT(iSOne, , Integer<1>)
SEMT_DEFINE_TRAIT(isRational, int i COMMA size_t j, Rational<i COMMA j>)
SEMT_DEFINE_TRAIT(isLiteral, class T, Literal<T>)
SEMT_DEFINE_TRAIT(isFunc, SEMT_PRECISION (f)(), Func<f>)

template<class T, class Enable = void>
class isConstExpr : public false_type
{
};
template<class T>
class isConstExpr<T, typename enable_if_c<(
        isInt<T>::value
                || isRational<T>::value
                || isLiteral<T>::value
                || isFunc<T>::value
                || isParam<T>::value
        )>::type>: public true_type
{
};

template<class T, class Enable = void>
struct isLeaf : public false_type
{
};
template<class T>
struct isLeaf<T, typename enable_if_c<(
        isConstExpr<T>::value
                || isVar<T>::value
        )>::type>: public true_type
{
};

/// @}

/*!
 * @addtogroup bintraits
 * @{
 */
SEMT_DEFINE_TRAIT(isPlus, class LHS COMMA class RHS, Plus<LHS COMMA RHS>)
SEMT_DEFINE_TRAIT(isTimes, class LHS COMMA class RHS, Times<LHS COMMA RHS>)
SEMT_DEFINE_TRAIT(isMinus, class LHS COMMA class RHS, Minus<LHS COMMA RHS>)
SEMT_DEFINE_TRAIT(isDivide, class LHS COMMA class RHS, Divide<LHS COMMA RHS>)
SEMT_DEFINE_TRAIT(isEPower, class base COMMA class exp, EPower<base COMMA exp>)
SEMT_DEFINE_TRAIT(isPower, class base COMMA int exp, Power<base COMMA exp>)

template<class T, class Enable = void>
struct isBinOp : public false_type
{
};
template<class T>
struct isBinOp<T, typename enable_if_c<(
        isPlus<T>::value
                || isTimes<T>::value
                || isMinus<T>::value
                || isDivide<T>::value
                || isEPower<T>::value
                || isPower<T>::value
        )>::type>: public true_type
{
};
/// @}

/*!
 * @addtogroup un_traits
 * @{
 */

SEMT_DEFINE_TRAIT(isExp, class T, Exp_t<T>)
SEMT_DEFINE_TRAIT(isLn, class T, Ln_t<T>)

SEMT_DEFINE_TRAIT(isSin, class T, Sin_t<T>)
SEMT_DEFINE_TRAIT(isCos, class T, Cos_t<T>)
SEMT_DEFINE_TRAIT(isTan, class T, Tan_t<T>)
SEMT_DEFINE_TRAIT(isArcsin, class T, Arcsin_t<T>)
SEMT_DEFINE_TRAIT(isArccos, class T, Arccos_t<T>)
SEMT_DEFINE_TRAIT(isArctan, class T, Arctan_t<T>)

SEMT_DEFINE_TRAIT(isSinh, class T, Sinh_t<T>)
SEMT_DEFINE_TRAIT(isCosh, class T, Cosh_t<T>)
SEMT_DEFINE_TRAIT(isTanh, class T, Tanh_t<T>)
SEMT_DEFINE_TRAIT(isArsinh, class T, Arsinh_t<T>)
SEMT_DEFINE_TRAIT(isArcosh, class T, Arcosh_t<T>)
SEMT_DEFINE_TRAIT(isArtanh, class T, Artanh_t<T>)

SEMT_DEFINE_TRAIT(isCond, class T COMMA class C COMMA class A, Defined_if<T COMMA C COMMA A>)
SEMT_DEFINE_TRAIT(isSimpleCond, class arg COMMA class cond, Defined_if<arg COMMA cond COMMA arg>)

template<class T>
struct isComplexCond
{
    static const bool value = isCond<T>::value && !isSimpleCond<T>::value;
};

template<class T, class Enable = void>
struct isUnOp : public false_type
{
};
template<class T>
struct isUnOp<T, typename enable_if_c<(
        isExp<T>::value
                || isLn<T>::value
                || isSin<T>::value
                || isCos<T>::value
                || isTan<T>::value
                || isArcsin<T>::value
                || isArccos<T>::value
                || isArctan<T>::value
                || isSinh<T>::value
                || isCosh<T>::value
                || isTanh<T>::value
                || isArsinh<T>::value
                || isArcosh<T>::value
                || isArtanh<T>::value
                || isCond<T>::value
        )>::type>: public true_type
{
};

/// @}

/*!
 * @addtogroup traits
 * @{
 */

template<class T, class Enable = void>
struct isSemtType : public false_type
{
};

template<class T>
struct isSemtType<T, typename enable_if_c<(
        isLeaf<T>::value
                || isBinOp<T>::value
                || isUnOp<T>::value
        )>::type>: public true_type
{
};

/// @}

}// namespace SEMT

#endif // __SEMT_TRAITS__H__
