/*!
 * Include this file to access very evil macros.
 * Actually, I don't expect them to do any harm.
 * Just saves time to write those tedious Expr<BlaType<...> >.
 * @file
 * @ingroup api
 * @par Usage
 *   @dontinclude semt_examples.cpp
 *   @skip int macros
 *   @until }
 */

#ifndef __SEMT_SHORTCUTS__H__
#define __SEMT_SHORTCUTS__H__

namespace SEMT
{

/*!
 * Define a instance of Expr<Variable<INDEX> > with NAME.
 */
#define DVAR(NAME, INDEX) SEMT::Expr<SEMT::Variable<INDEX> > NAME;

/*!
 * Shortcut to Variables.
 * @param INDEX The index of the instantiated variable.
 */
#define VAR(INDEX) SEMT::Expr<SEMT::Variable<INDEX> >()

/*!
 * Define a instance of Expr<Integer<NUMBER> > with NAME.
 */
#define DINT(NAME, value) SEMT::Expr<SEMT::Integer<value> > NAME;

/*!
 * Shortcut to Integers.
 * @param value The value of the instantiated Integer.
 */
#define INT(value) SEMT::Expr<SEMT::Integer<value> >()

/*!
 * Define a instance of Expr<Rational<nom, denom> > with NAME.
 */
#define DRAT(NAME, nom, denom) SEMT::Expr<SEMT::Rational<nom, denom> > NAME;

/*!
 * Shortcut to rational numbers.
 * @param nom The nominator of the instantiated rational.
 * @param denom The denominator of the instantiated rational.
 */
#define RAT(nom, denom) SEMT::Expr<SEMT::Rational<nom, denom> >()

/*!
 * Define a instance of Expr<Parameter<NUMBER, 't'> > with NAME.
 */
#define DPARAM(NAME, INDEX) SEMT::Expr<SEMT::Parameter<INDEX> > NAME;

/*!
 * Shortcut to Parameters named t.
 * @param INDEX The index of the instantiated Parameter.
 */
#define PARAM(INDEX) SEMT::Expr<SEMT::Parameter<INDEX> >()

} // namespace SEMT

#endif // __SEMT_SHORTCUTS__H__
