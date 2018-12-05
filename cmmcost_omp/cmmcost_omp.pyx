# cython: language_level=3
cimport cython
from libcpp cimport bool

cdef extern from "_cmmcost_omp.h":
    double costdirect(double * z, bool * z2use, size_t nz, float * s2,
        unsigned long long int * is2, double * p, double  * sb2, double s02,
        unsigned char * annot)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def get_cost(double[:] z, bool[:] z2use, float[:] s2, unsigned long long int[:] is2,
    double[:] p, double[:] sb2, double s02, unsigned char[:] annot):

    cdef size_t nz = z.shape[0]
    cdef double cost = costdirect(&z[0], &z2use[0], nz, &s2[0], &is2[0], &p[0],
        &sb2[0], s02, &annot[0])

    return cost


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
