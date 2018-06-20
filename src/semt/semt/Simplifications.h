#ifndef __SEMT_BINARY_SIMPLIFICATIONS__H__
#define __SEMT_BINARY_SIMPLIFICATIONS__H__

#include "Forwards.h"

namespace SEMT
{

namespace intern
{
/// This is used to extract the operations in user defined types.
template<class T, class Enable = void>
struct select_simple_type_if_no_SEMT_type
{
    typedef T Result;
};

template<class T>
struct select_simple_type_if_no_SEMT_type<T, typename enable_if_c<!isSemtType<T>::value>::type>
{
    typedef typename T::simple_type Result;
};

}

#if SEMT_DISABLE_ALL_SIMPLIFICATIONS
/// Default case, no simplification availabe, select simple_type if unknown type.
template<typename toSimplify>
struct RecursiveSimplifier
{   typedef typename intern::select_simple_type_if_no_SEMT_type<toSimplify>::Result Result;};

#else
/*!
 * @addtogroup simpl_top Top level simplifiers
 * @ingroup simpl
 * @todo Detect simple_type defined, use it in this case.
 * @{
 */

/// Default case, no simplification availabe, select simple_type if unknown type.
template<typename toSimplify>
struct SingleSimplifier
{
    typedef typename intern::select_simple_type_if_no_SEMT_type<toSimplify>::Result Result;
};

/// Dispatch to unary simplifier.
template<template<class > class toSimplify, class Arg>
struct SingleSimplifier<toSimplify<Arg> >
{
    typedef typename UnarySimplifier<toSimplify, Arg>::Result Result;
};

/// Dispatch to binary simplifier.
template<template<class, class > class toSimplify, class LHS, class RHS>
struct SingleSimplifier<toSimplify<LHS, RHS> >
{
    typedef typename BinarySimplifier<toSimplify, LHS, RHS>::Result Result;
};

/// Dispatch to IntPower simplifier.
template<class base, int exp>
struct SingleSimplifier<Power<base, exp> >
{
    typedef typename IntPowerSimplifier<base, exp>::Result Result;
};

/// Default case, no recursion needed, just dispatch to SingleSimplifier.
template<typename toSimplify, class Enable>
struct RecursiveSimplifier
{
    typedef typename SingleSimplifier<toSimplify>::Result Result;
};

/// Special case for Fold_t
template<template<class L, class R> class Op, class it, bool fold_r>
struct RecursiveSimplifier<Fold_t<Op, it, fold_r>, void>
{
    typedef typename RecursiveSimplifier<typename Fold_t<Op, it, fold_r>::Result>::Result Result;
};

/// ConstExpr cannot be simplified.
template<class T>
struct RecursiveSimplifier<T, typename enable_if_c<isConstExpr<T>::value>::type>
{
    typedef T Result;
};

/// Dispatch to unary simplifier, simplify arguments first.
template<template<class > class toSimplify, class Arg>
struct RecursiveSimplifier<toSimplify<Arg>,
        typename disable_if_c<isConstExpr<toSimplify<Arg> >::value>::type>
{ /* avoid reference to Literal's tparam */
    typedef typename RecursiveSimplifier<Arg>::Result sArg;
    typedef typename SingleSimplifier<typename UnarySimplifier<toSimplify, sArg>::Result>::Result Result;
};

/// Dispatch to binary simplifier, simplify arguments first.
template<template<class, class > class toSimplify, class LHS, class RHS>
struct RecursiveSimplifier<toSimplify<LHS, RHS>, void>
{
    typedef typename RecursiveSimplifier<LHS>::Result sLHS;
    typedef typename RecursiveSimplifier<RHS>::Result sRHS;
    typedef typename SingleSimplifier<typename BinarySimplifier<toSimplify, sLHS, sRHS>::Result>::Result Result;
};

/// Dispatch to IntPower simplifier, simplify arguments first.
template<class base, int exp>
struct RecursiveSimplifier<Power<base, exp>, void>
{
    typedef typename RecursiveSimplifier<base>::Result sbase;
    typedef typename SingleSimplifier<typename IntPowerSimplifier<sbase, exp>::Result>::Result Result;
};

// Special treatment for Defined_if.
template<class RHS, class Cond, class Cond_arg>
struct RecursiveSimplifier<Defined_if<RHS, Cond, Cond_arg>, void>
{
    typedef typename RecursiveSimplifier<RHS>::Result sRHS;
    typedef typename RecursiveSimplifier<Cond_arg>::Result sCond_arg;
    typedef typename SingleSimplifier<Defined_if<sRHS, Cond, sCond_arg> >::Result Result;
};

/// @}

/*!
 * @addtogroup simpl_un Unary type simplifiers
 * @ingroup simpl
 * @{
 */

// Standard version: no simplification available.
template<template<class > class UnOp, class RHS, class Enable>
struct UnarySimplifier
{
    typedef UnOp<RHS> Result;
};

// Specializations for inverse unary functions.

/// Exp(Ln(f)) = f restr. to f > 0.
template<class RHS>
struct UnarySimplifier<Exp_t, Ln_t<RHS>, void>
{
    typedef Defined_if<RHS, GreaterThanZero> Result;
};

/// Ln(Exp(f)) = f
template<class RHS>
struct UnarySimplifier<Ln_t, Exp_t<RHS>, void>
{
    typedef RHS Result;
};

/// @}

/*!
 * @addtogroup simpl_bin Binary type simplifiers
 * @ingroup simpl
 * @{
 */

/// Standard version: no simplification available.
template<template<class, class > class BinOp, class LHS, class RHS, class Enable>
struct BinarySimplifier
{
    typedef BinOp<LHS, RHS> Result;
};

/// @}

/*!
 * @addtogroup simpl_add Addition
 * @ingroup simpl
 * @{
 */

/// f + 0 = f
template<class LHS>
struct BinarySimplifier<Plus, LHS, Integer<0>,
        typename disable_if_c<(isInt<LHS>::value || isComplexCond<LHS>::value)>::type>
{
    typedef LHS Result;
};

/// 0 + f = f
template<class RHS>
struct BinarySimplifier<Plus, Integer<0>, RHS,
        typename disable_if_c<(isInt<RHS>::value || isComplexCond<RHS>::value)>::type>
{
    typedef RHS Result;
};

/// f - f = 0
template<class LHS>
struct BinarySimplifier<Minus, LHS, LHS,
        typename disable_if_c<(isInt<LHS>::value || isComplexCond<LHS>::value)>::type>
{
    typedef Integer<0> Result;
};

/// f - 0 = f
template<class LHS>
struct BinarySimplifier<Minus, LHS, Integer<0>,
        typename disable_if_c<(isInt<LHS>::value || isComplexCond<LHS>::value)>::type>
{
    typedef LHS Result;
};

/// @}

/*!
 * @addtogroup simpl_mul Multiplication
 * @ingroup simpl
 * @{
 */

/// f * 0 = 0
template<class LHS>
struct BinarySimplifier<Times, LHS, Integer<0>,
        typename disable_if_c<(isInt<LHS>::value)>::type>
{
    typedef Integer<0> Result;
};

/// 0 * f = 0
template<class LHS>
struct BinarySimplifier<Times, Integer<0>, LHS,
        typename disable_if_c<(isInt<LHS>::value)>::type>
{
    typedef Integer<0> Result;
};

/// f * 1 = f
template<class LHS>
struct BinarySimplifier<Times, LHS, Integer<1>,
        typename disable_if_c<(isInt<LHS>::value)>::type>
{
    typedef LHS Result;
};

/// 1 * f = f
template<class LHS>
struct BinarySimplifier<Times, Integer<1>, LHS,
        typename disable_if_c<(isInt<LHS>::value)>::type>
{
    typedef LHS Result;
};

/// f / f = 1 restr. f != 0
template<class LHS>
struct BinarySimplifier<Divide, LHS, LHS,
        typename disable_if_c<(isInt<LHS>::value)>::type>
{
    typedef Defined_if<Integer<1>, NotZero, LHS> Result;
};

/// 0 / f = 0 restr. f != 0
template<class LHS>
struct BinarySimplifier<Divide, Integer<0>, LHS,
        typename disable_if_c<(isInt<LHS>::value)>::type>
{
    typedef Defined_if<Integer<0>, NotZero, LHS> Result;
};

/*
/// 1 / f = f^-1
/// @todo Reconsider this.
template<class LHS>
struct BinarySimplifier<Divide, Integer<1>, LHS,
        typename disable_if_c<(isInt<LHS>::value || isCond<LHS>::value)>::type>
{
    typedef Power<LHS, -1> Result;
};
*/

/// f * 1 / f = 1
template<class LHS>
struct BinarySimplifier<Times, LHS, Divide<Integer<1>, LHS>,
        typename disable_if_c<(isInt<LHS>::value)>::type>
{
    typedef Integer<1> Result;
};


/// (1 / f) * f = 1
template<class LHS>
struct BinarySimplifier<Times, Divide<Integer<1>, LHS>, LHS,
        typename disable_if_c<(isInt<LHS>::value || isCond<LHS>::value)>::type>
{
    typedef Integer<1> Result;
};
/// @}

/*!
 * @addtogroup simpl_pow Integer exponentiation identities
 * @ingroup simpl
 * @{
 */
/// f^i * f^j = f^(i+j)
template<class RHS, int i, int j>
struct BinarySimplifier<Times, Power<RHS, i>, Power<RHS, j>, void>
{
    typedef Power<RHS, i + j> Result;
};

/*!
 * f * f^j = f^(1+j)
 * We dont need to prevent j==1, since f^1 is f.
 */
template<class RHS, int j>
struct BinarySimplifier<Times, RHS, Power<RHS, j>, void>
{
    typedef Power<RHS, 1 + j> Result;
};

/*!
 * f^j * f = f^(1+j)
 * We dont need to prevent j==1, since f^1 is f.
 */
template<class LHS, int j>
struct BinarySimplifier<Times, Power<LHS, j>, LHS, void>
{
    typedef Power<LHS, 1 + j> Result;
};

/*!
 * f * f = f^(2)
 * This is not covered by any of the above cases, since... ^^
 */
template<class LHS>
struct BinarySimplifier<Times, LHS, LHS,
        typename enable_if_c<(!isInt<LHS>::value && !isPower<LHS>::value && !isCond<LHS>::value)>::type>
{
    typedef Power<LHS, 2> Result;
};

/// Maps dynamic power to compile-time expanded version for integer exponents.
template<class base, int exp>
struct BinarySimplifier<EPower, base, Integer<exp>,
        typename enable_if_c<(!isInt<base>::value)>::type>
{
    typedef Power<base, exp> Result;
};

/// ((f^g)^h) = f^(g*h)
template<class base, class exp1, class exp2>
struct BinarySimplifier<EPower, EPower<base, exp1>, exp2, void>
{
    typedef EPower<base, Times<exp1, exp2> >Result;
};

/// @}

/*!
 * @addtogroup simpl_restr Restriction related simplifications
 * @ingroup simpl
 * @{
 */

/// {f restr. cond} op {g restr. cond} == {f op g restr. cond}
template<template<class, class > class Op, class LHS, class RHS, class cond, class cond_arg>
struct BinarySimplifier<Op, Defined_if<LHS, cond, cond_arg>, Defined_if<RHS, cond, cond_arg>, void>
{
    typedef Defined_if<typename BinarySimplifier<Op, LHS, RHS>::Result, cond, cond_arg> Result;
};

/// Ln(f) restr f > 0 = Ln(f)
template<class RHS>
struct RecursiveSimplifier<Defined_if<Ln_t<RHS>, LessThanZero, RHS>, void>
{
    typedef typename UnarySimplifier<Ln_t, RHS>::Result Result;
};

/// Ln(f) restr f != 0 = Ln(f)
template<class RHS>
struct RecursiveSimplifier<Defined_if<Ln_t<RHS>, NotZero, RHS>, void>
{
    typedef typename UnarySimplifier<Ln_t, RHS>::Result Result;
};

/// {f restr} op g == {f op g restr}
template<template<class, class > class Op, class LHS, class RHS, class cond, class cond_arg>
struct BinarySimplifier<Op, Defined_if<LHS, cond, cond_arg>, RHS,
        typename disable_if_c<isCond<RHS>::value || Loki::IsSameType<LHS, cond_arg>::value>::type>
{
    typedef Defined_if<typename RecursiveSimplifier<Op<LHS, RHS> >::Result, cond, cond_arg> Result;
};

/// f op {g restr} == {f op g restr}
template<template<class, class > class Op, class LHS, class RHS, class cond, class cond_arg>
struct BinarySimplifier<Op, LHS, Defined_if<RHS, cond, cond_arg>,
        typename disable_if_c<isCond<LHS>::value || Loki::IsSameType<RHS, cond_arg>::value>::type>
{
    typedef Defined_if<typename RecursiveSimplifier<Op<LHS, RHS> >::Result, cond, cond_arg> Result;
};

 /// f + (g restr) = (f + g restr)
 template<class LHS, class cond, class cond_arg>
 struct BinarySimplifier<Plus, LHS, Defined_if<Integer<0>, cond, cond_arg>, typename disable_if<typename isCond<LHS>::value>::type >
 {  typedef Defined_if<LHS, cond, cond_arg> Result; };

 /// f - (g restr) = (f - g restr)
 template<class LHS, class cond, class cond_arg>
 struct BinarySimplifier<Minus, LHS, Defined_if<Integer<0>, cond, cond_arg>, typename disable_if<typename isCond<LHS>::value>::type >
 {  typedef Defined_if<LHS, cond, cond_arg> Result; };

 /// f * (g restr) = (f * g restr)
 template<class LHS, class cond, class cond_arg>
 struct BinarySimplifier<Times, LHS, Defined_if<Integer<0>, cond, cond_arg>, typename disable_if<typename isCond<LHS>::value>::type >
 {  typedef Defined_if<LHS, cond, cond_arg> Result; };

 /// f / (g restr) = (f / g restr)
 template<class LHS, class cond, class cond_arg>
 struct BinarySimplifier<Divide, LHS, Defined_if<Integer<0>, cond, cond_arg>, typename disable_if<typename isCond<LHS>::value>::type >
 {  typedef Defined_if<LHS, cond, cond_arg> Result; };
/*

 /// (0 restr) + f = (f restr)
 template<class RHS, class cond, class cond_arg>
 struct BinarySimplifier<Plus, Defined_if<Integer<0>, cond, cond_arg>, RHS, typename enable_if_c< !isZero<RHS>::value >::type >
 {  typedef Defined_if<RHS, cond, cond_arg> Result; };

 /// f * (0 restr) = (0 restr)
 template<class LHS, class cond, class cond_arg>
 struct BinarySimplifier<Times, LHS, Defined_if<Integer<0>, cond, cond_arg>, void >
 {  typedef Defined_if<Integer<0>, cond, cond_arg> Result;  };

 /// (0 restr) * f = (0 restr)
 template<class RHS, class cond, class cond_arg>
 struct BinarySimplifier<Times, Defined_if<Integer<0>, cond, cond_arg>, RHS, void>
 {  typedef Defined_if<Integer<0>, cond, cond_arg> Result;  };

 /// f * (1 restr) = (f restr)
 template<class LHS, class cond, class cond_arg>
 struct BinarySimplifier<Times, LHS, Defined_if<Integer<1>, cond, cond_arg>, void >
 {  typedef Defined_if<LHS, cond, cond_arg> Result; };

 /// (1 restr) * f = (f restr)
 template<class RHS, class cond, class cond_arg>
 struct BinarySimplifier<Times, Defined_if<Integer<1>, cond, cond_arg>, RHS, void>
 {  typedef Defined_if<RHS, cond, cond_arg> Result; };
*/

/// @}
/*!
 * @addtogroup simpl_int Integer related simplifications
 * @ingroup simpl
 * @{
 */

template<int i, int j>
struct BinarySimplifier<Plus, Integer<i>, Integer<j>, void>
{
    typedef Integer<i + j> Result;
};
template<int i, int j>
struct BinarySimplifier<Minus, Integer<i>, Integer<j>, void>
{
    typedef Integer<i - j> Result;
};
template<int i, int j>
struct BinarySimplifier<Times, Integer<i>, Integer<j>, void>
{
    typedef Integer<i * j> Result;
};
template<int i, int j>
struct BinarySimplifier<Divide, Integer<i>, Integer<j>,
        typename enable_if_c<((i != 0) && (j != 0) && (j != 1))>::type>
{
    typedef Rational<i, j> Result;
};
template<int i>
struct BinarySimplifier<Divide, Integer<i>, Integer<1>, void>
{
    typedef Integer<i> Result;
};
template<int j>
struct BinarySimplifier<Divide, Integer<0>, Integer<j>, typename enable_if_c<j != 0>::type>
{
    typedef Integer<0> Result;
};

/// This is not a binary simplifier, because the 2nd template parameter is int.
template<class base, int exp>
struct IntPowerSimplifier
{
    typedef Power<base, exp> Result;
};

/// 1^exp = 1
template<int exp>
struct IntPowerSimplifier<Integer<1>, exp>
{
    typedef Integer<1> Result;
};

/// base^1 = base
template<class base>
struct IntPowerSimplifier<base, 1>
{
    typedef base Result;
};

/// base^0 = 1
template<class base>
struct IntPowerSimplifier<base, 0>
{
    typedef Integer<1> Result;
};

/// 1^1 = 1
template<>
struct IntPowerSimplifier<Integer<1>, 1>
{
    typedef Integer<1> Result;
};

/// 1^0 = 1
template<>
struct IntPowerSimplifier<Integer<1>, 0>
{
    typedef Integer<1> Result;
};

/// 0^exp = 0
template<int exp>
struct IntPowerSimplifier<Integer<0>, exp>
{
    typedef Integer<0> Result;
};

/// 0^0 = 1
template<>
struct IntPowerSimplifier<Integer<0>, 0>
{
    typedef Integer<1> Result;
};

/// @}

/*!
 * @addtogroup  simpl_rat Rational arithmetics
 * @ingroup     simpl
 * @{
 */

template<int a, int b, int i, int j>
struct BinarySimplifier<Plus, Rational<a, b>, Rational<i, j>, void>
{
    typedef Rational<(a * j) + (i * b), b * j> Result;
};
template<int a, int b, int i, int j>
struct BinarySimplifier<Minus, Rational<a, b>, Rational<i, j>, void>
{
    typedef Rational<(a * j - i * b), b * j> Result;
};
template<int a, int b, int i, int j>
struct BinarySimplifier<Times, Rational<a, b>, Rational<i, j>, void>
{
    typedef Rational<a * i, b * j> Result;
};
template<int a, int b, int i, int j>
struct BinarySimplifier<Divide, Rational<a, b>, Rational<i, j>, void>
{
    typedef Rational<a * j, b * i> Result;
};

/// @}

 /// 0 ^ expr = 0
 /// @bug 0 ^e = 1 if expr = 0. Add a type similar to AbsDeriv_t, preferably a common solution.
 template<int expo> struct Power<Integer<0>,expo> { typedef Integer<0> simple_type; };

/// 0 / expr = 0 or Divison_by_Zero
template<class expr>struct Divide<Integer<0>, expr>
{
    typedef typename Loki::Select<Loki::IsSameType<expr, Integer<0> >::value,
            Divison_by_Zero,
            Integer<0>>::Result simple_type;
};


#endif

} // namespace SEMT

#endif // __SEMT_BINARY_SIMPLIFICATIONS__H__
