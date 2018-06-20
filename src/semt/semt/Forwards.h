/*!
 * This file contains forward declarations and is included in every other header.
 * @file
 * @ingroup dev
 */

#ifndef __SEMT_FORWARDS__H__
#define __SEMT_FORWARDS__H__

#ifndef __SEMT_SEMTFWD__H__
    #include "Semtfwd.h"
#endif
#ifndef __SEMT_MACROS__H__
    #include "Macros.h"
#endif

/*!
 * Namespace for the whole project, everything except operators is defined in here.
 */
namespace SEMT
{

/// Compile time minimum calculation.
template< int l, int r> struct Min
{   static const int value = (l < r) ? l : r;   };

/// Compile time maximum calculation.
template< int l, int r> struct Max
{   static const int value = (l > r) ? l : r;   };

/*!
 * Namespace for internal helpers.
 * @ingroup dev
 */
namespace intern
{}

/*!
 * Detected division by zero at compile-time.
 * @struct Divison_by_Zero
 * @ingroup errtype
 */
struct Divison_by_Zero;

// constant expression & variable forwards
template<int n> struct Integer;
template<int nom, size_t denom> struct Rational;
template<typename T> struct Literal;
template<SEMT_PRECISION (f)()> struct Func;
template<SEMT_PRECISION (f)(CAR)> struct VFunc;
template<int n> struct Variable;
template<int n, char name = 't'> struct Parameter;
template<int count, char name = 't'> struct set_parameters;

// conditions
struct GreaterThanZero;
struct LessThanZero;
struct AlmostZero;
struct NotZero;
template<typename expr, typename Condition, typename cond_arg = expr> struct Defined_if;

// binary type forwards
template<typename LHS, typename RHS> struct Plus;
template<typename LHS, typename RHS> struct Minus;
template<typename LHS, typename RHS> struct Times;
template<typename LHS, typename RHS> struct Divide;
template<typename LHS, typename RHS> struct EPower;
template<typename expr, int expo> struct Power;

// unary type forwards
template<typename ARG> struct Sgn_t;
template<typename ARG> struct Abs_t;
template<typename ARG> struct Exp_t;
template<typename ARG> struct Ln_t;
template<typename ARG> struct Sin_t;
template<typename ARG> struct Arcsin_t;
template<typename ARG> struct Cos_t;
template<typename ARG> struct Arccos_t;
template<typename ARG> struct Tan_t;
template<typename ARG> struct Arctan_t;
template<typename ARG> struct Sinh_t;
template<typename ARG> struct Arsinh_t;
template<typename ARG> struct Cosh_t;
template<typename ARG> struct Arcosh_t;
template<typename ARG> struct Tanh_t;
template<typename ARG> struct Artanh_t;

template<template<int> class it, int i, int end, int skip = 1> struct IntIterator;
template<class typelist> struct TLIterator;
template<template<class L, class R> class Op, class it, bool fold_right = true> struct Fold_t;

#if SEMT_DISABLE_ALL_SIMPLIFICATIONS
template<typename toSimplify> struct RecursiveSimplifier;
#else
template<typename toSimplify> struct SingleSimplifier;
template<typename toSimplify, class Enable = void> struct RecursiveSimplifier;
template<template<class> class UnOp, class RHS, class Enable = void> struct UnarySimplifier;
template<template<class, class> class BinOp, class LHS, class RHS, class Enable = void> struct BinarySimplifier;
template<class base, int exp> struct IntPowerSimplifier;
#endif

template<typename e> struct Expr;

} // namespace SEMT

#endif // __SEMT_FORWARDS__H__
