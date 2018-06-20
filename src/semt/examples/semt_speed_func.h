/*!
 * @file
 * @todo annotations
 */

#ifndef __SEMT_FUNC_IMPL__H__
#define __SEMT_FUNC_IMPL__H__

#define SEMT_DISABLE_PRINT 0
#include "semt/Semtfwd.h"
#include "semt/VectorExpr.h"

class my_semt_vector
{
public:
    my_semt_vector();

    const SEMT::VectorExpr& get_fn(size_t i) const { return functions_[i];}
    void set_parameters(SEMT::CAR new_values) const;
private:
    std::vector<SEMT::VectorExpr> functions_;
};

#endif // __SEMT_FUNC_IMPL__H__
