#include "VectorExpr.h"

namespace SEMT
{
namespace intern
{
/*!
 * @addtogroup dev
 * @{
 */

/// No-fail if the first count pointers in data are valid.
template<class T> void free_data(const std::vector<T*>& data, size_t count)
{
    while (count != 0)
    {
        --count;
        assert((data[count] != 0) && "VectorExpr class invariant broken: null pointer found.");
        delete data[count];
    }
}

/// Exception-neutral swap.
template<class T> void swap(T& t1, T& t2)
{
    const T tmp(t1);
    t1 = t2;
    t2 = tmp;
}

/// @} // group dev

}// namespace intern

VectorExpr::VectorExpr(const VectorExpr& rhs) :
        _data(rhs.size(), 0), _fill(0), _vars(rhs.vars())
{
    const size_t fill = rhs.size();
    try
    {
        // use _fill as counter so we can just call free_data() on error.
        for (_fill = 0; _fill < fill; ++_fill)
        {
            // since rhs is valid, we don't need to check.
            _data[_fill] = (((rhs[_fill])->duplicate()).release());
        }
    } catch (...)
    {
        intern::free_data(_data, _fill);
        _fill = 0;
        _vars = 0;
        throw;
    }
    assert((rhs.size() == (fill)) && "copy ctor: post cond failed");
}

VectorExpr::~VectorExpr()
{
    intern::free_data(_data, _fill);
}

VectorExpr& VectorExpr::operator=(const VectorExpr& other)
{
    VectorExpr tmp(other);
    Swap(tmp);
    return *this;
}

const AbstractFunctionObject* VectorExpr::operator[](size_t i) const
                                                     {
    assert(i < _fill && "out-of-bounds access");
    return _data[i];
}

void VectorExpr::push_back(const AbstractFunctionObject* const ex)
{
    assert(ex && "push_back: pre cond failed");
    _data.push_back(ex);
    ++_fill;
}

void VectorExpr::Swap(VectorExpr& other)
{
    intern::swap(this->_data, other._data);
    intern::swap(this->_fill, other._fill);
    intern::swap(this->_vars, other._vars);
}

} // namespace SEMT

std::ostream& operator<<(std::ostream& os, const SEMT::VectorExpr& ex)
{
    os << "[\t";
    const size_t values = ex.size() - 1;
    for (size_t i = 0; i < values; ++i)
    {
        os << ex[i]->toString() << ",\n\t";
    }
    os << ex[values]->toString() << "\t]";
    return os;
}

