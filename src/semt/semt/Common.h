/*!
 * Provides commonly used definitions.
 * @file
 */

#ifndef __SEMT_COMMON__H__
#define __SEMT_COMMON__H__

#include "Shortcuts.h"

/*
 * It doesn't matter if all the SEMT types are instantiated at global or local scope,
 * none of them carry any non-static value.
 */

//@{
//! Variables
DVAR(x0, 0);
DVAR(x1, 1);
DVAR(x2, 2);
DVAR(x3, 3);
DVAR(x4, 4);
DVAR(x5, 5);
DVAR(x6, 6);
DVAR(x7, 7);
DVAR(x8, 8);
DVAR(x9, 9);
DVAR(x10, 10);
//@}

//@{
//! Integers
DINT(mTwo, -2);
DINT(mOne, -1);
DINT(Zero, 0);
DINT(One, 1);
DINT(Two, 2);
DINT(Three, 3);
DINT(Four, 4);
DINT(Five, 5);
DINT(Six, 6);
DINT(Seven, 7);
DINT(Eight, 8);
DINT(Nine, 9);
DINT(Ten, 10);
DINT(Hundred, 100);
DINT(Tousand, 1000);
//@}

//@{
//! Rational numbers.
DRAT(ThreeQuarters, 3, 4);
DRAT(Quarter, 1, 4);
//@}

/// @}
#endif // __SEMT_COMMON__H__
