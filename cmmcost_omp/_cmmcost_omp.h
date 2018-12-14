// gcc -I/home/alexeas/.local/include -L/home/alexeas/.local/lib -Wall -O2 -ffast-math -fopenmp -lgsl -lgslcblas -lm -o _cmmcost_omp _cmmcost_omp.c
// time ./_cmmcost_omp
// gcc on Abel:
// gcc -I/cluster/software/VERSIONS/gsl-1.16/include -L/cluster/software/VERSIONS/gsl-1.16/lib -Wall -O2 -ffast-math -fopenmp -lgsl -lgslcblas -lm -o _cmmcost_omp_gcc _cmmcost_omp.c
// time ./_cmmcost_omp_gcc
// icc on Abel (https://software.intel.com/en-us/articles/step-by-step-optimizing-with-intel-c-compiler):
// icc _cmmcost_omp.c -o _cmmcost_omp_icc -L/cluster/software/VERSIONS/gsl-2.2/lib -I/cluster/software/VERSIONS/gsl-2.2/include -std=c99 -O2 -no-prec-div -qopenmp -lgsl -lgslcblas -lm
#ifndef _CMMCOST_OMP
#define _CMMCOST_OMP

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <omp.h>

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

double costsampling(double * z, _Bool * z2use, size_t nz, float * s2,
    unsigned long long int * is2, double * p, double * sb2, double s02,
    unsigned char * annot, size_t n_samples);

void get_s2sample(double * s2sample, int n_s2sample, double * p, double * sb2,
    double s02, float * s2_block, unsigned char * annot_block, int block_size,
    gsl_rng * r);

void cdfsampling(double * zcdf, double * zcdf_qq_annot, double * zgrid,
    _Bool * z2use, size_t nz, size_t n_zgrid, size_t n_qq_annot, float * s2,
    unsigned long long int * is2, double * p, double * sb2, double s02,
    unsigned char * annot, _Bool * qq_template_annot, size_t n_samples);

#endif
