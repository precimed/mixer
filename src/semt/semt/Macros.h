/*!
 * All macros unrelated to user settings are defined here.
 * You may want to include this when implementing new SEMT-compatible operations.
 * @file
 * @attention   If this file is included a second time, all macros in this file are \#undef'd.
 * @ingroup dev
 */

#ifndef __SEMT_MACROS__H__
#define __SEMT_MACROS__H__

/*!
 * @addtogroup macros
 * @{
 */

/*!
 * Define a static method toString() if printing is enabled.
 * @param   STRING      Is passed to std::ostringstream via operator<< .
 */
#define SEMT_DEFINE_TOSTRING(STRING)                            \
    static SEMT_INLINE_STATIC std::string toString()            \
    {                                                           \
        std::ostringstream ret;                                 \
        ret << STRING;                                          \
        return ret.str();                                       \
    };

#define SEMT_SIMPLE_TYPE(t) t
#define SEMT_SIMPLE_TYPE2(t) typename t
/*
 #ifndef SEMT_DISABLE_ALL_SIMPLIFICATIONS
 /// Forwards to simple_type if SEMT_DISABLE_ALL_SIMPLIFICATIONS is false.
 #define SEMT_SIMPLE_TYPE(t) typename t::simple_type
 /// Like SEMT_SIMPLE_TYPE but typename is put in front anyway.
 #define SEMT_SIMPLE_TYPE2(t) typename t::simple_type
 #else
 #endif
 */

/*!
 * Define a type representing a constant.
 * @param L_NAME    Name of the type.
 * @param T_ARGS    Template parameter(s) of this type.
 * @param T_NAMES   Repeat the template parameter(s) without type specifiers.
 * @param APPLY     Return statement of the apply method.
 * @param STR       This is passed SEMT_DEFINE_TOSTRING.
 */
#define SEMT_DEFINE_CONSTEXPR(L_NAME, T_ARGS, T_NAMES, APPLY, STR)         \
template<T_ARGS> struct L_NAME {                                           \
    typedef L_NAME<T_NAMES> simple_type;                                   \
    static const int LastVar = 0;                                          \
    template<typename var> struct partial                                  \
    {                                                                      \
        static const bool dependent = false;                               \
        typedef Integer<0> deriv;                                          \
    };                                                                     \
    static SEMT_INLINE_STATIC SEMT_PRECISION apply(CAR x)                  \
    {   return APPLY;   }                                                  \
    SEMT_DEFINE_TOSTRING(STR);                                             \
};

/*!
 * Define types that represent unary operations.
 * Types available in the macro parameters: First and First_d.
 * @param TNAME     Name of the type.
 * @param APPLY     Return statement of the apply method, most certainly using First::apply(x),
 *                  but any operation is possible.
 * @param PARTIAL   Type of the partial deriv of this type in terms of TNAME, First, First_d
 *                  and any other known type.
 * @param STRING    Is passed SEMT_DEFINE_TOSTRING.
 */
#define SEMT_DEFINE_UNARY_TYPE(TNAME, APPLY, PARTIAL, STRING)                               \
template<typename RHS>                                                                      \
struct TNAME                                                                                \
{                                                                                           \
    typedef SEMT_SIMPLE_TYPE(RHS) First;                                                    \
    static const int LastVar = First::LastVar;                                              \
    template<typename var> struct partial                                                   \
    {                                                                                       \
        static const bool dependent = First::template partial<var>::dependent;              \
        typedef SEMT_SIMPLE_TYPE2(First::template partial<var>::deriv) First_d;             \
        typedef PARTIAL d_tmp;                                                              \
        typedef typename Loki::Select<dependent, SEMT_SIMPLE_TYPE(d_tmp), Integer<0> >::Result deriv;\
                                                                                            \
    };                                                                                      \
    static SEMT_INLINE_STATIC SEMT_PRECISION apply(CAR x)                                   \
    {                                                                                       \
        return APPLY;                                                                       \
    }                                                                                       \
    SEMT_DEFINE_TOSTRING(STRING);                                                           \
};

/*!
 * Define a templated function, that instantiates a class template with the type of it's argument.
 * @param   op_t A class template expecting one template parameter.
 * @param   op_n The name of the function.
 * @note    All types will be wrapped in Expr.
 */
#define SEMT_DEFINE_UNARYOP(op_t, op_n)                         \
template<typename RHS>                                          \
SEMT_INLINE SEMT::Expr<SEMT::op_t<RHS > >                       \
op_n(const SEMT::Expr<RHS>)                                     \
{                                                               \
    return SEMT::Expr<SEMT::op_t<RHS > >();                     \
};

/*!
 * Define types that represent binary operations.
 * Types available in the macro parameters: First, First_d, Second and Second_d.
 * You can also use LastVar.
 * @param TNAME     Name of the type.
 * @param APPLY     Return statement of the apply method, most certainly using:
 *                  First::apply(x) and Second::apply(x), but any operation is possible.
 * @param PARTIAL   Type of the partial deriv of this type in terms of
 *                  TNAME, First, First_d, Second, Second_d and any other known type.
 * @param STRING    Is passed SEMT_DEFINE_TOSTRING.
 */
#define SEMT_DEFINE_BINARY_TYPE(TNAME, APPLY, PARTIAL, STRING)                              \
template<typename LHS, typename RHS>                                                        \
struct TNAME                                                                                \
{                                                                                           \
    typedef SEMT_SIMPLE_TYPE(LHS) First;                                                    \
    typedef SEMT_SIMPLE_TYPE(RHS) Second;                                                   \
    static const int LastVar = Max<First::LastVar, Second::LastVar>::value;                 \
    template<typename var> struct partial                                                   \
    {                                                                                       \
        static const bool dependent = First::template partial<var>::dependent               \
                                        || Second::template partial<var>::dependent;        \
        typedef SEMT_SIMPLE_TYPE2(First::template partial<var>::deriv) First_d;             \
        typedef SEMT_SIMPLE_TYPE2(Second::template partial<var>::deriv) Second_d;           \
        typedef PARTIAL d_tmp;                                                              \
        typedef typename Loki::Select<dependent, SEMT_SIMPLE_TYPE(d_tmp), Integer<0> >::Result deriv;\
    };                                                                                      \
    static SEMT_INLINE_STATIC SEMT_PRECISION apply(CAR x)                                   \
    {                                                                                       \
        return APPLY;                                                                       \
    }                                                                                       \
    SEMT_DEFINE_TOSTRING(STRING);                                                           \
};

/*!
 * Define a templated function, that instantiates a class template with the types of it's arguments.
 * @param   op_t A class template expecting two template parameters.
 * @param   op_n The name of the defined function.
 * @note    All types will be wrapped in Expr.
 */
#define SEMT_DEFINE_BINARYOP(op_t, op_n)                        \
template<typename LHS, typename RHS>                            \
SEMT_INLINE SEMT::Expr<SEMT::op_t<LHS, RHS > >                  \
        op_n(SEMT::Expr<LHS>, SEMT::Expr<RHS>)                  \
{                                                               \
    return SEMT::Expr<SEMT::op_t<LHS, RHS > >();                \
};

/*!
 * Define a type that represents some condition.
 * @param TNAME     Name of the type.
 * @param TEST      Condition a value x must satisfy.
 * @param STRING    Is passed to SEMT_DEFINE_TOSTRING().
 */
#define SEMT_DEFINE_CONDITION(TNAME, TEST, STRING)              \
struct TNAME                                                    \
{                                                               \
    static SEMT_INLINE bool is_true_for(SEMT_PRECISION x)       \
    { return (TEST); }                                          \
                                                                \
    SEMT_DEFINE_TOSTRING(STRING);                               \
};

/*!
 * Define foldr and foldl.
 * @todo description
 */
#define SEMT_DEFINE_FOLD(token, right, directions)                                              \
template<   template<class L, class R> class Op,                                                \
            template <int i> class nodes,                                                       \
            int start, int end, int skip = 1>                                                   \
Expr<typename Fold_t<Op, IntIterator<nodes, directions>, right >::Result >                      \
    fold##token(    Expr<nodes<start> > = Expr<nodes<start> >(),                                \
                    Expr<Integer<end> > = Expr<Integer<end> >(),                                \
                    Expr<Integer<skip> > = Expr<Integer<skip> >())                              \
{                                                                                               \
    return Expr<typename Fold_t<Op, IntIterator<nodes, directions>, right>::Result >();         \
};

/*!
 * Define a type trait.
 */
#define SEMT_DEFINE_TRAIT(name, t_param, type)          \
        template<class T>                               \
        struct name : public false_type {};             \
        template<t_param>                               \
        struct name<type> : public true_type {};

/*!
 * Avoid cutting of macro parameters.
 */
#define COMMA ,

/// @} // group macros

#else
#undef SEMT_DEFINE_TOSTRING

#undef SEMT_SIMPLE_TYPE
#undef SEMT_SIMPLE_TYPE2

#undef SEMT_DEFINE_CONSTEXPR

#undef SEMT_DEFINE_UNARY_TYPE
#undef SEMT_DEFINE_UNARYOP

#undef SEMT_DEFINE_BINARY_TYPE
#undef SEMT_DEFINE_BINARYOP

#undef SEMT_DEFINE_CONDITION

#undef SEMT_DEFINE_FOLD

#undef COMMA
#endif // __SEMT_MACROS__H__
