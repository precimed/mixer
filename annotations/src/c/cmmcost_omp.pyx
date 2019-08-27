# cython: language_level=3
cimport cython
from libcpp cimport bool
from libc.stdlib cimport malloc, free
import numpy as np
cimport numpy as np

cdef extern from "_cmmcost_omp.h":
    double costdirect(double * z, bool * z2use, size_t nz, float * s2,
        unsigned long long int * is2, double * p, double  * sb2, double s02,
        unsigned char * annot)

    void costdirect_derivative(double * z, bool * z2use, size_t nz, float * s2,
        unsigned long long int * is2, double * p, double * sb2, double s02,
        unsigned char * annot, size_t n_annot, double * gradient)

    double costsampling(double * z, bool * z2use, size_t nz, float * s2,
        unsigned long long int * is2, double * p, double * sb2, double s02,
        unsigned char * annot, size_t n_samples)

    void cdfsampling(double * zcdf, double * zcdf_qq_annot, double * zgrid,
        bool * z2use, size_t nz, size_t n_zgrid, size_t n_qq_annot, float * s2,
        unsigned long long int * is2, double * p, double * sb2, double s02,
        unsigned char * annot, bool * qq_template_annot, size_t n_samples)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def get_cost(double[:] z, bool[:] z2use, float[:] s2, unsigned long long int[:] is2,
        double[:] p, double[:] sb2, double s02, unsigned char[:] annot, unsigned char n_annot):

    cdef size_t nz = z.shape[0]
    cdef double cost = costdirect(&z[0], &z2use[0], nz, &s2[0], &is2[0], &p[0],
        &sb2[0], s02, &annot[0])

    return cost

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def get_cost_der(double[:] z, bool[:] z2use, float[:] s2, unsigned long long int[:] is2,
        double[:] p, double[:] sb2, double s02, unsigned char[:] annot, unsigned char n_annot):

    cdef size_t nz = z.shape[0]
    cdef size_t na = n_annot
    cdef double[::1] gradient = np.empty(2*n_annot+1, order='C')
    costdirect_derivative(&z[0], &z2use[0], nz, &s2[0], &is2[0], &p[0], &sb2[0],
        s02, &annot[0], na, &gradient[0])

    return np.asarray(gradient)



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def get_costsampling(double[:] z, bool[:] z2use, float[:] s2,
        unsigned long long int[:] is2, double[:] p, double[:] sb2, double s02,
        unsigned char[:] annot, unsigned int n_samples):

    cdef size_t nz = z.shape[0]
    cdef double cost = costsampling(&z[0], &z2use[0], nz, &s2[0], &is2[0], &p[0],
        &sb2[0], s02, &annot[0], n_samples)

    return cost


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def get_cdfsampling( double[::1]                 zgrid,
                     bool[::1]                   z2use,
                     float[::1]                  s2,
                     unsigned long long int[::1] is2,
                     double[::1]                 p,
                     double[::1]                 sb2,
                     double                      s02,
                     unsigned char[::1]          annot,
                                                 qq_template_annot,
                     unsigned int                n_samples ):

    # for arr in [zgrid, z2use, s2, is2, p, sb2, annot, qq_template_annot]:
    if not qq_template_annot.flags['C_CONTIGUOUS']:
        qq_template_annot = np.ascontiguousarray(qq_template_annot)
    cdef bool[:,::1] qq_template_annot_memview = qq_template_annot

    cdef size_t nz = z2use.shape[0]
    cdef size_t n_zgrid = zgrid.shape[0]
    cdef size_t n_qq_annot = qq_template_annot.shape[1]

    cdef double[::1] zcdf = np.empty(n_zgrid, order='C')
    cdef double[:,::1] zcdf_qq_annot = np.empty((n_zgrid, n_qq_annot), order='C')

    cdfsampling(&zcdf[0], &zcdf_qq_annot[0,0], &zgrid[0], &z2use[0], nz,
            n_zgrid, n_qq_annot, &s2[0], &is2[0], &p[0], &sb2[0], s02,
            &annot[0], &qq_template_annot_memview[0,0],n_samples)

    return np.asarray(zcdf), np.asarray(zcdf_qq_annot)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def get_nonzero_min_mean_max(float[:] x):
    cdef float xi = 0, x_max = -1.E9, x_min = 1.E9
    cdef double x_sum = 0
    cdef size_t i, nnz = 0, N = x.shape[0]
    for i in range(0,N):
        xi = x[i]
        if xi != 0:
            nnz += 1
            x_sum += xi
            if xi > x_max: x_max = xi
            if xi < x_min: x_min = xi
    return x_min, x_sum/nnz, x_max
