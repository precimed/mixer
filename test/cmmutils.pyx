# cython: language_level=3
cimport cython
from cython.parallel import prange


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def fill_s2(double[:]             het,
            unsigned char[:]      annot,
            unsigned long long[:] is2,

            float[:]              r2,
            unsigned int[:]       i21,

            unsigned int[:]       ind_2use,
            double[:]             ssize,
            unsigned long long[:] is2_2use,

            float[:]              s2_2use,
            unsigned char[:]      annot_s2_2use):
    """
    het           [N]  : het for all template SNPs on current chromosome  
    annot         [N]  : annotations for SNPs in template
    is2           [N+1]: cummulative block size for tempalte SNPs, is2[0] = 0

    r2            [L]  : complete r2 vector from ld2npz, L >= N
    i21           [L]  : index of the secod SNP for correlations in r2

    ind_2use      [K]  : indices of SNPs to use, K <= N
    ssize         [K]  : sample size for SNPs to use
    is2_2use      [K+1]: cummulative block size for SNPs to use, is2_2use[0] = 0
    
    s2_2use       [M]  : s2 array corresponding to current chr, M = is2_2use[-1]
    annot_s2_2use [M]  : annot_s2 array corresponding to current chr
    """
    cdef Py_ssize_t i, j, K = ind_2use.shape[0]
    cdef unsigned long long r2_i, s2_2use_i
    cdef double ss
    for i in prange(K, nogil=True, schedule='dynamic'): # this can be converted to prange
        ss = ssize[i]
        s2_2use_i = is2_2use[i]
        for r2_i in range(is2[ind_2use[i]], is2[ind_2use[i]+1]):
            #print(len(s2_2use), len(r2), len(het), len(i21))
            #print(s2_2use_i, r2_i, i21[r2_i]) 
            s2_2use[s2_2use_i] = r2[r2_i]*het[i21[r2_i]]*ss
            annot_s2_2use[s2_2use_i] = annot[i21[r2_i]]
            s2_2use_i = s2_2use_i + 1

