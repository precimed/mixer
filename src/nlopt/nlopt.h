/* Copyright (c) 2007-2014 Massachusetts Institute of Technology
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
* LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
* OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
* WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef NLOPT_H
#define NLOPT_H

#include <stdarg.h>             /* for va_list */
#include <stddef.h>             /* for ptrdiff_t and size_t */

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

  typedef enum {
    NLOPT_FAILURE = -1,         /* generic failure code */
    NLOPT_INVALID_ARGS = -2,
    NLOPT_OUT_OF_MEMORY = -3,
    NLOPT_ROUNDOFF_LIMITED = -4,
    NLOPT_FORCED_STOP = -5,
    NLOPT_SUCCESS = 1,          /* generic success code */
    NLOPT_STOPVAL_REACHED = 2,
    NLOPT_FTOL_REACHED = 3,
    NLOPT_XTOL_REACHED = 4,
    NLOPT_MAXEVAL_REACHED = 5,
    NLOPT_MAXTIME_REACHED = 6
  } nlopt_result;
#define NLOPT_MINF_MAX_REACHED NLOPT_STOPVAL_REACHED

  /* stopping criteria */
  typedef struct {
    unsigned n;
    double minf_max;
    double ftol_rel;
    double ftol_abs;
    double xtol_rel;
    const double *xtol_abs;
    int *nevals_p, maxeval;
    double maxtime, start;
    int *force_stop;
    char **stop_msg;        /* pointer to msg string to update */
  } nlopt_stopping;

  extern int nlopt_stop_f(const nlopt_stopping * stop, double f, double oldf);
  extern int nlopt_stop_ftol(const nlopt_stopping * stop, double f, double oldf);
  extern int nlopt_stop_x(const nlopt_stopping * stop, const double *x, const double *oldx);
  extern int nlopt_stop_dx(const nlopt_stopping * stop, const double *x, const double *dx);
  extern int nlopt_stop_xs(const nlopt_stopping * stop, const double *xs, const double *oldxs, const double *scale_min, const double *scale_max);
  extern int nlopt_stop_evals(const nlopt_stopping * stop);
  extern int nlopt_stop_time_(double start, double maxtime);
  extern int nlopt_stop_time(const nlopt_stopping * stop);
  extern int nlopt_stop_evalstime(const nlopt_stopping * stop);
  extern int nlopt_stop_forced(const nlopt_stopping * stop);

  typedef double(*nlopt_func) (unsigned n, const double *x,
    double *gradient, /* NULL if not needed */
    void *func_data);

  /* like vsprintf, but reallocs p to whatever size is needed */
  extern char *nlopt_vsprintf(char *p, const char *format, va_list ap);
  extern void nlopt_stop_msg(const nlopt_stopping * s, const char *format, ...)
#ifdef __GNUC__
    __attribute__((format(printf, 2, 3)))
#endif
    ;

  int nlopt_isinf(double x);
  int nlopt_isfinite(double x);
  int nlopt_istiny(double x);
  int nlopt_isnan(double x);

  /* re-entrant qsort */
  extern void nlopt_qsort_r(void *base_, size_t nmemb, size_t size, void *thunk, int(*compar) (void *, const void *, const void *));

  /* seconds timer */
  extern double nlopt_seconds(void);
  extern unsigned long nlopt_time_seed(void);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif

