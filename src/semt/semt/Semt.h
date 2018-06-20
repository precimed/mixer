/*!
 * SEMT main header file, include this to access all features of Semt.
 * @file
 * @ingroup api
 */

#ifndef __SEMT_SEMT__H__
#define __SEMT_SEMT__H__

#include "Semtfwd.h"
#include "Macros.h"
#include "Forwards.h"

#include "Constexpr.h"
#include "Variable.h"
#include "Parameter.h"
#include "Conditions.h"

#include "Traits.h"
#include "Simplifications.h"

#include "Expression.h"

#include "BinaryTypes.h"
#include "BinaryOperators.h"

#include "UnaryTypes.h"
#include "UnaryOperator.h"

#include "Sequence.h"

#include "VectorExpr.h"
#include "DifferentiableVectorExpr.h"

// include again to undef everything
#include "Macros.h"

#endif // __SEMT_SEMT__H__
