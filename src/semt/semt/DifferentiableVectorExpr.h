#ifndef __SEMT_DIFFERENTIABLE_VECTOR_EXPR__H__
#define __SEMT_DIFFERENTIABLE_VECTOR_EXPR__H__

#include "Forwards.h"

namespace SEMT
{
namespace intern
{
/*!
 * @addtogroup dev
 * @{
 */

/*!
 * Bridge between template expressions and homogeneous function objects.
 * Employed by VectorExpr to store different expressions in a homogeneous container.
 * Instantiated by DVEHelper and DifferentiableVectorExpr.\ņ
 * Carries no value, only the virtual table.
 * @tparam expr The expression this object will evaluate.
 */
template<typename expr>
struct ConcreteFunctionObject : public AbstractFunctionObject
{
    /// Virtual cloning.
    virtual SEMT_INLINE_VIRTUAL std::auto_ptr<AbstractFunctionObject> duplicate() const
    {
        return std::auto_ptr<AbstractFunctionObject>(new ConcreteFunctionObject<expr>);
    }

    /// Evaluation is forwarded to expr's Expr::apply(Array).
    virtual SEMT_INLINE_VIRTUAL SEMT_PRECISION operator()(CAR x) const
                                                          {
        return expr::apply(x);
    }

    /// Forward to expr's static Expr::toString().
    virtual SEMT_INLINE_VIRTUAL std::string toString() const
    {
        return expr::toString();
    }
};

/*!
 * Helper to instantiate all derivatives of a single expression.
 * This is done by a two-step recursion and here is the main case:
 * recurse over vars and differentiation depth.
 * @tparam  e       Expression to differentiate.
 * @tparam  vars    Number of variables to be differentiated for.
 * @tparam  deriv_order Maximum differentiation order.
 * @attention   If current_depth > deriv_order or current_var > vars,
 *              we'll end up in an infinite recursion.
 * @see     DifferentiableVectorExpr
 */
template<typename e, size_t vars, size_t deriv_order,
        size_t current_var = 0, size_t current_depth = 1>
struct DVEHelper
{
    /// Differentiate e in respect to Variable<current_var>.
    typedef SEMT_SIMPLE_TYPE2(e::template partial<Variable<current_var> >::deriv) e_d_;
    typedef typename RecursiveSimplifier<e_d_>::Result e_d;
    /*!
     * Append derivatives in appropriate order into data.
     * @param   data    Size must be greater or equal deriv_order.
     */
    static SEMT_INLINE_STATIC void push_back(std::vector<VectorExpr>& data)
    {
        // push back current derivative
        data[current_depth].push_back(new ConcreteFunctionObject<e_d>);
        // start recursion to push back higher derivatives
        DVEHelper<e_d, vars, deriv_order, 0, current_depth + 1 >::push_back(data);
        // start recursion to push back derivatives with respect to next variables
        DVEHelper<e, vars, deriv_order, current_var + 1, current_depth>::push_back(data);
    }
};

/*!
 * Stop recursing the variables.
 * @see DVEHelper
 */
template<typename e, size_t vars, size_t deriv_order, size_t current_depth>
struct DVEHelper<e, vars, deriv_order, vars, current_depth>
{
    static SEMT_INLINE_STATIC void push_back(std::vector<VectorExpr>& data)
    {
    }
    //The Variables are counted from 0 to end-1, while the user will count them from
    //1 to end, so at this point there is nothing left to be done.
};

/*!
 * Stop recursing differentiation depth.
 * @see DVEHelper
 */
template<typename e, size_t vars, size_t deriv_order, size_t current_var>
struct DVEHelper<e, vars, deriv_order, current_var, deriv_order>
{
    typedef SEMT_SIMPLE_TYPE2(e::template partial<Variable<current_var> >::deriv) e_d_;
    typedef typename RecursiveSimplifier<e_d_>::Result e_d;
    static SEMT_INLINE_STATIC void push_back(std::vector<VectorExpr>& data)
    {
        // push back current derivative
        data[deriv_order].push_back(new ConcreteFunctionObject<e_d>);
        // this is the highest differentiation order, so we only need to push back
        // the derivatives with respect to the remaining variables
        DVEHelper<e, vars, deriv_order, current_var + 1, deriv_order>::push_back(data);
    }
};

/*!
 * Stop whole recursion.
 * @see DVEHelper
 */
template<typename e, size_t vars, size_t deriv_order>
struct DVEHelper<e, vars, deriv_order, vars, deriv_order>
{
    static SEMT_INLINE_STATIC void push_back(std::vector<VectorExpr>& data)
    {
    }
};

/*!
 * Special case if differentiation order is zero.
 * @see DVEHelper
 */
template<typename e, size_t vars>
struct DVEHelper<e, vars, 0, 0>
{
    static SEMT_INLINE_STATIC void push_back(std::vector<VectorExpr>& data)
    {
    }
};

/// @}

}// namespace intern

/*!
 * Construct mappings R^n --> R^m and their derivatives step by step.
 * This class is the gateway to the non-template realms of VectorExpr.\n
 * The first derivative is stored in the form: (Jacobian expanded row by row)\n
 * (df0/dx0, df0/dx1, ..., df0/dxn, df1/dx0, df2/dx1. ..., df2/dxn, ... ..., dfm/dx0, ..., dfm/dxn).\n
 * Take this as a map R^n --> R^(n*m) and obtain the form of 2nd derivative by the same rule and so on.
 * @tparam  vars    Number of variables to differentiate for. ( = n)
 * @tparam  deriv_order Number of derivatives to calculate.
 * @note    Be aware of the complexity you demand: template instantiation depth may increase
 *          rapidly for certain expressions.
 * @see     @ref example
 * @ingroup api
 */
template<size_t vars, size_t deriv_order = 1>
class DifferentiableVectorExpr
{
    // _data[i] contains the i-th derivative of this function
    std::vector<VectorExpr> _data;
    DifferentiableVectorExpr& operator=(const DifferentiableVectorExpr& rhs);
    public:
    DifferentiableVectorExpr() :
            _data(deriv_order + 1, VectorExpr(vars))
    {
    }

    /*!
     * Set to a single expression, also used for comma initialization.
     */
    template<typename expr>
    DifferentiableVectorExpr& operator=(const Expr<expr>& ex)
    {
        _data = std::vector<VectorExpr>(deriv_order + 1, VectorExpr(vars));
        push_back(ex);
        return *this;
    }

    /*!
     * Append a component to this vector, as a side effect, the derivatives are also appended.
     * @tparam  e   The type of the component.
     * @param   ex  An instance of this type.
     * @attention   If you put more variables in the expressions then you told me,
     *              I don't care. It's your fault.
     */
    template<typename expr>
    SEMT_INLINE void push_back(const Expr<expr>& ex = Expr<expr>())
    {
        typedef typename RecursiveSimplifier<expr>::Result e_simpl;
        //push back this expression:
        _data[0].push_back(new intern::ConcreteFunctionObject<SEMT_SIMPLE_TYPE(e_simpl)>);
        //start recursion to push back the derivatives:
        intern::DVEHelper<SEMT_SIMPLE_TYPE(e_simpl), vars, deriv_order>::push_back(_data);
    }

    /*!
     * After you're done appending components, get all the derivatives you want.
     * @param   i   The derivatives order, default is the function itself.
     * @return  A reference to a VectorExpr, which hould be valid as long as this object.
     * @attention   Thou shall not ask for more than deriv_order.
     * @note    You are leaving the template realms.
     */
    const VectorExpr& get_derivative(size_t i = 0)
    {
        assert((i <= deriv_order) && "Increase template argument to obtain higher derivatives.");
        return _data[i];
    }

    const VectorExpr& get_function()
    {
        return _data[0];
    }

    /*!
     * Calculate all the expressions for a given vector of unknowns.
     * @param   x       Values for the unknown.
     * @param   res     Vector to store results. Size must be at least size().
     */
    SEMT_INLINE void eval(CAR x, Array& res) const
    {
        _data[0].eval(x, res);
    }

    /*!
     * Calculate all the expressions for a given vector of unknowns.
     * @param   x       Values for the unknown.
     * @return  Vector with the results.
     */
    SEMT_INLINE Array operator()(CAR x) const
    {
        return _data[0](x);
    }

    /*!
     * Syntactic sugar, see the Example for usage.
     */
    template<typename expr>
    DifferentiableVectorExpr& operator ,(const Expr<expr>& e)
    {
        push_back(e);
        return *this;
    }
};

} // namesapce SEMT

#endif // __SEMT_DIFFERENTIABLE_VECTOR_EXPR__H__
