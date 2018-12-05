// gcc -Wall -I/home/alexeas/.local/include -O2 -c _cmmcost.c
// gcc -L/home/alexeas/.local/lib _cmmcost.o -lgsl -lgslcblas -lm -o _cmmcost
// time ./_cmmcost
#ifndef _CMMCOST
#define _CMMCOST

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>

typedef struct Parameters {
    double * p;
    double * s2;
    size_t n;
    double s02;
    double omega;
} Parameters;

double iff(double x, void * params);

double costdirect(double * z, _Bool * z2use, size_t nz, float * s2,
    unsigned long long int * is2, double * p, double  * sb2, double s02,
    unsigned char * annot);

#endif
