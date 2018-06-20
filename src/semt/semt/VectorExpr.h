#ifndef __SEMT_VECTOREXPR__H__
#define __SEMT_VECTOREXPR__H__

#include "Semtfwd.h"

namespace SEMT
{

/*!
 * Implementation of a vector of expressions, this is the final object representing a map R^n --> R^m.
 * Implements a vector of expression with complete value semantics: you can pass/return it to/from
 * functions either per value or per reference. Construction is usually done by DifferentiableVectorExpr.
 * In fact, this is much like a array of function pointers.\n
 * Evaluation is performed on all stored expressions, but single component evaluation is possible via
 * operator[] and the interface AbstractFunctionObject.\n
 * Use operator<< or the matrix_str template to obtain a string representation.
 *
 * @attention No range checking is performed if NDEBUG is defined.
 * @note        You have left the template area.
 * @todo        Consider reference counting for the data.
 * @ingroup api
 */
class VectorExpr
{
public:

    /// Strong exception safe c-tor.
    VectorExpr(size_t vars) :
            _data(), _fill(0), _vars(vars)
    {
    }

    /// Strong exception safe c-tor.
    VectorExpr(const VectorExpr& rhs);

    /// No-fail destructor.
    ~VectorExpr();

    /// Strong exception safe assignment.
    VectorExpr& operator=(const VectorExpr& other);

    /*!
     * Calculate all the expressions for a given vector of unknowns.
     * @param   x       Values for the unknown.
     * @param   res     Vector to store results. Size must be at least size().
     */
    SEMT_INLINE void eval(CAR x, Array& res) const
    {
        assert(res.size() >= _fill && "eval: pre cond failed");
        assert(x.size() >= _vars && "eval: pre cond failed");
        for (size_t i = 0; i < _fill; ++i)
        {
            res[i] = (*(_data[i]))(x);
        }
    }

    /*!
     * Calculate all the expressions for a given vector of unknowns.
     * @param   x       Values for the unknown.
     * @return  Vector with the results.
     */
    SEMT_INLINE Array operator()(CAR x) const
    {
        Array res(_fill);
        eval(x, res);
        return res;
    }

    /// Get the k-th expression of this vector.
    const AbstractFunctionObject* operator[](size_t k) const;

    /*!
     * Insert another expression in this vector.
     * @param   ex      A pointer to an instance of ConcreteFunctionObject.
     * @details This is usually not called by the user.
     *          DifferentiableVectorExpr uses this while constructing the derivatives.
     * @attention Takes ownership of the pointed-to expression.
     */
    void push_back(const AbstractFunctionObject* ex);

    /// How many expressions are in this vector.
    size_t size() const
    {
        return _fill;
    }

    /// How many variables are assumed in the expressions.
    size_t vars() const
    {
        return _vars;
    }
private:
    std::vector<const AbstractFunctionObject*> _data;
    size_t _fill, _vars;

    void free_data();
    void Swap(VectorExpr& other);
    VectorExpr();
};

} // namespace SEMT

std::ostream& operator<<(std::ostream& os, const SEMT::VectorExpr& ex);

template<size_t rows, size_t cols>
std::string matrix_str(const SEMT::VectorExpr& ex)
{
    const size_t values = ex.size();
    if (values != (rows * cols))
        throw std::runtime_error("invalid size in matrix_str.");

    std::ostringstream os;
    for (size_t i = 0; i < rows; ++i)
    {
        os << "|\t";
        for (size_t j = 0; j < cols; ++j)
        {
            os << ex[i * cols + j]->toString() << ",\t";
        }
        os << "|\n";

    }
    return os.str();
}

#endif // __SEMT_VECTOREXPR__H__
