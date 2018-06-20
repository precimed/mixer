#ifndef __SEMT_SEQUENCE__H__
#define __SEMT_SEQUENCE__H__

#include "Forwards.h"

/*!
 * Serves as placeholder for Iterator-based compile-time loops.
 * Iterator::next should be this type, if Iterator::last ist true.
 * Can be used for template specialization of types defining compile-time loops.
 * @ingroup loops
 * @class Loki::NullType
 */

namespace SEMT
{
/*!
 * @addtogroup loops
 * @{
 */

/*!
 * A simple iterator for a integer-templated list of types.
 * The list this iterators traverses is given as a template template argument, expecting
 * an integer value.\n
 * This is quite different from the usual iterator pattern, since here the iterator is
 * in fact the traversed list itself.
 * @tparam  nodes   An integer parameterized type, musst supply a typedef simple_type.
 * @tparam  i       Starting or current value to instantiate nodes with.
 * @tparam  end     Either upper or lower bound for this iterator to stop, depending on sgn(skip).
 * @tparam  skip    Stepsize, may be either positve or negative. Zero is fine, fool.
 * @attention       Make sure that i + k * skip == end for some k >= 0 if you want nodes<end> to be the
 *                  last instantiated type. nodes<end> is likely not even in the list otherwise.
 */
template<template<int> class nodes, int i, int end, int skip>
struct IntIterator
{
    /// Is this the last type in the list (like list.end() - 1)?
    static const bool last = (skip > 0) ? i >= end : i <= end;

    /// The current item in the list (like *).
    typedef typename nodes<i>::simple_type actual_type;

    /// The next iterator (like ++).
    typedef typename Loki::Select<last, Loki::NullType, IntIterator<nodes, i + skip, end, skip> >::Result next;
};

/*!
 * Iterator for Loki::Typelist.
 */
template<class Head, class Tail>
struct TLIterator<Loki::Typelist<Head, Tail> >
{
    /// Is this the last type in the list (like list.end() - 1)?
    static const bool last = Loki::IsSameType<Tail, Loki::NullType>::value;

    /// The current item in the list (like *).
    typedef typename intern::select_simple_type_if_no_SEMT_type<Head>::Result actual_type;
    //typedef typename Loki::Select<isSemtType<Head>::value, Head, typename Head::simple_type>::Result actual_type;

    /// The next iterator (like ++).
    typedef typename Loki::Select<last, Loki::NullType, TLIterator<Tail> >::Result next;
};

template<> struct TLIterator<Loki::NullType>
{
};

template<class expr>
Loki::Typelist<expr, Loki::NullType> start_tl(const Expr<expr>&)
{
    return Loki::Typelist<expr, Loki::NullType>();
}

/*!
 * Recursively instantiate a binary type over a list of types.
 * @tparam  Op      A type expecting two template parameters, must supply a typedef simple_type,
 *                  or even conform to the expression interface if simplification is disabled.
 * @tparam  it      A type-iterator.
 * @tparam  fold_r  Paranthesize in postorder.
 */
template<template<class L, class R> class Op, class it, bool fold_r>
struct Fold_t
{
    // The recursion stops if it::next is Loki::NullType.
    typedef typename Fold_t<Op, typename it::next, fold_r>::Result inner_loop;

    // If we fold right, the loop goes right, if we go left, the loop goes left.
    typedef typename Loki::Select<fold_r,
            Op<typename it::actual_type, inner_loop>,
            Op<inner_loop, typename it::actual_type>
            >::Result tail;

    // If the iterator is consumed, return the final type,
    // otherwise we instantiate the tail.
    typedef typename Loki::Select<it::last,
            typename it::actual_type,
            tail
    >::Result Result;
};

/// Partial specialization to stop recursive instantiation.
template<template<class L, class R> class Op, bool fold_r>
struct Fold_t<Op, Loki::NullType, fold_r>
{
    typedef Loki::NullType Result;
};

/*!
 * Apply an unary operator to a list of types.
 * @tparam Op   The operator to apply.
 * @tparam it   A type iterator.
 */
template<template<class > class Op, class it>
struct For_each
{
    // The recursion stops if it::next is Loki::NullType.
    typedef typename For_each<Op, typename it::next>::Result inner_loop;

    typedef Op<typename it::actual_type> act_node;

    // If the iterator is consumed, return the final type,
    // otherwise we instantiate the tail.
    typedef Loki::Typelist<act_node, inner_loop> Result;
};

/// Partial specialization to stop recursive instantiation.
template<template<class > class Op>
struct For_each<Op, Loki::NullType>
{
    typedef Loki::NullType Result;
};

/*!
 * Fold a list of nodes with a given operator.
 * @tparam Op       Operation that is performed on the types.
 * @tparam nodes    A type expecting an integral value as a template argument,
 *                  must supply a typedef simple_type.
 * @tparam start    First value to instantiate nodes with.
 * @tparam end      The upper bound of the integral arguments for nodes.
 * @tparam skip     Stepsize, any integer is fine.
 */
SEMT_DEFINE_FOLD(r, true, start COMMA end COMMA skip);

/// Fold left. See foldr.
SEMT_DEFINE_FOLD(l, false, end COMMA start COMMA -skip);

/*!
 * Fold a list of nodes with a given operator.
 * @tparam Op       Operation that is performed on the types.
 * @tparam nodes    A Loki::Typelist with SEMT-compatible types, see \c start_tl().
 */
template<template<class L, class R> class Op, class types>
Expr<typename Fold_t<Op, TLIterator<types>, true>::Result>
fold_tl_r(const types& ts)
{
    return Expr<typename Fold_t<Op, TLIterator<types>, true>::Result>();
}

/// Fold left. See fold_tl_r.
template<template<class L, class R> class Op, class types>
Expr<typename Fold_t<Op, TLIterator<types>, false>::Result>
fold_tl_l(const types& ts)
{
    return Expr<typename Fold_t<Op, TLIterator<types>, true>::Result>();
}

/*!
 * Apply Op to every type in types.
 * @todo example
 */
template<template<class > class Op, class types>
typename For_each<Op, TLIterator<types> >::Result
for_each(const types& ts)
{
    return typename For_each<Op, TLIterator<types> >::Result();
}

/// @} // group loops

}// namesapce SEMT

/*!
 * Continue a list of expressions started \c via start_tl().
 * @ingroup loops
 */
template<class TL, class expr>
typename Loki::TL::Append<TL, expr>::Result operator ,(const TL&, const SEMT::Expr<expr>&)
{
    return typename Loki::TL::Append<TL, expr>::Result();
}

#endif // __SEMT_SEQUENCE__H__
