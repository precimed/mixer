#ifndef __SEMT_FUNC_IMPL__H__
#define __SEMT_FUNC_IMPL__H__

#define SEMT_DISABLE_PRINT 0
#include "semt/Semtfwd.h"
#include "semt/VectorExpr.h"

class my_semt_func
{
public:
    my_semt_func();
    const SEMT::VectorExpr& get_func(size_t i) const
    {
        return bunch_of_functions[i];
    }
private:
    std::vector<SEMT::VectorExpr> bunch_of_functions;
};

#endif // __SEMT_FUNC_IMPL__H__
