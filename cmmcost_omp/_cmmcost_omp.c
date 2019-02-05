#include "_cmmcost_omp.h"


double iff(double x, void * params) {
    Parameters pars = *(Parameters *) params;
    double neg_half_x2 = -0.5*x*x;
    double res = cos(pars.omega*x)*exp(neg_half_x2*pars.s02)/M_PI;
    for(size_t i=0; i<pars.n; i++) {
        res *= 1 - pars.p[i] + pars.p[i]*exp(neg_half_x2*pars.s2[i]);
    }
    // return GSL_NAN;
    return res;
}


double costdirect(double * z, _Bool * z2use, size_t nz, float * s2,
    unsigned long long int * is2, double * p, double * sb2, double s02,
    unsigned char * annot) {

    size_t nthreads = omp_get_max_threads();
    // printf("Max number of threads = %zu\n", nthreads);

    size_t i, j, k, tid;
    double result[nthreads];
    double error[nthreads];
    Parameters pars[nthreads];
    gsl_function FF[nthreads];

    gsl_integration_workspace ** ww =
        (gsl_integration_workspace **) malloc(nthreads * sizeof(gsl_integration_workspace *));
    for (i=0; i<nthreads; i++) {
        ww[i] = gsl_integration_workspace_alloc(1000);
        pars[i].s02 = s02; // this parameter is always the same
        FF[i].function = &iff;
        FF[i].params = &pars[i];
    }
    
    double cost = 0.;
    size_t nnz = 0;

    gsl_set_error_handler_off();

    #pragma omp parallel for default(shared) private(i,j,k,tid) schedule(static) reduction(+:cost,nnz)
    for (i=0; i<nz; i++) {
        if (z2use[i]) {
            tid = omp_get_thread_num();
            nnz++;
            size_t block_size = is2[i+1] - is2[i];
            //TODO: don't need to  allocate p_block and s2_block. Similar approach as in costsampling can be applied. 
            double * p_block = (double *) malloc(sizeof(double) * block_size);
            double * s2_block = (double *) malloc(sizeof(double) * block_size);

            for (j=0; j<block_size; j++) {
                k = is2[i] + j;
                p_block[j] = p[annot[k]];
                s2_block[j] = sb2[annot[k]] * s2[k];
            }

            pars[tid].p = p_block;
            pars[tid].s2 = s2_block;
            pars[tid].n = block_size;
            pars[tid].omega = z[i];

            int status = gsl_integration_qagiu (&FF[tid], 0, 1.E-7, 1.E-5, 1000,
                ww[tid], &result[tid], &error[tid]);

            if (status) {
                printf ("gsl error (status code: %d): %s\n", status, gsl_strerror (status));
                cost += 25.;
            }
            else if (result[tid] > 0) {
                cost += -log(result[tid]);
            }
            else {
                cost += 25.;
            }

            free(p_block);
            free(s2_block);
        }
    }

    for (i=0; i<nthreads; i++)
        gsl_integration_workspace_free(ww[i]);
    free(ww);

    return cost/(double)nnz;
}


void get_s2sample(double * s2sample, int n_s2sample, double * p, double * sb2,
    double s02, float * s2_block, unsigned char * annot_block, int block_size,
    gsl_rng * r) {

    int i, j;
    double pj, s2eff;

    for (i=0; i<n_s2sample; i++) {
        s2sample[i] = s02;
    }

    for (j=0; j<block_size; j++) {
        s2eff = s2_block[j]*sb2[annot_block[j]];
        pj = p[annot_block[j]];
        for (i=0; i<n_s2sample; i++) {
            if (gsl_rng_uniform(r) < pj)
                s2sample[i] += s2eff;
        }
    }
}



double costsampling(double * z, _Bool * z2use, size_t nz, float * s2,
    unsigned long long int * is2, double * p, double * sb2, double s02,
    unsigned char * annot, size_t n_samples) {

    size_t nthreads = omp_get_max_threads();
    // printf("Max number of threads in costsampling = %zu\n", nthreads);

    size_t i, j, tid;
    double cost = 0.;
    size_t nnz = 0;
    // s2sample_t [nthreads * n_samples] 1D array
    double *s2sample_t;
    // rngs [nthreads] 1D array
    gsl_rng ** rngs;

    s2sample_t = (double *) malloc(nthreads*n_samples*sizeof(double*));
    rngs = (gsl_rng **) malloc(nthreads * sizeof(gsl_rng *));
    for (i=0; i<nthreads; i++) {
        gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
        // 0 seed sets to default seed, which is the same as 1 for gsl_rng_taus, so should start seeding with 1
        gsl_rng_set(r, i+1);
        rngs[i] = r;
    }

    #pragma omp parallel for default(shared) private(i,j,tid) schedule(static) reduction(+:cost,nnz)
    for (i=0; i<nz; i++) {
        if (z2use[i]) {
            tid = omp_get_thread_num();
            nnz++;
            int block_size = is2[i+1] - is2[i];
            double mix_pdf = 0;

            get_s2sample(&s2sample_t[tid*n_samples], n_samples, p, sb2, s02,
                &s2[is2[i]], &annot[is2[i]], block_size, rngs[tid]);

            for (j=0; j<n_samples; j++) {
                // if (s2sample_t[tid*n_samples + j] > 0) - this should not be the case since s20 is always > 0
                mix_pdf += gsl_ran_gaussian_pdf(z[i], sqrt(s2sample_t[tid*n_samples + j]));
            }
            mix_pdf /= (double)n_samples;
            
            if (mix_pdf > 0) {
                cost += -log(mix_pdf);
            }
            else {
                cost += 25.;
            }
        }
    }

    // free memory
    for(i=0; i<nthreads; i++)
        free(rngs[i]);
    free(rngs);
    free(s2sample_t);

    return cost/((double)nnz); 
}



void cdfsampling(double * zcdf, double * zcdf_qq_annot, double * zgrid,
    _Bool * z2use, size_t nz, size_t n_zgrid, size_t n_qq_annot, float * s2,
    unsigned long long int * is2, double * p, double * sb2, double s02,
    unsigned char * annot, _Bool * qq_template_annot, size_t n_samples) {
    /*
    Estimates model CDF on zgrid for all SNPs (zcdf) and subsets of SNPs given
    by qq_template_annot array (zcdf_qq_annot).
    Args:
        zcdf              [n_zgrid]              1D arr
        zcdf_qq_annot     [n_zgrid * n_qq_annot] 1D arr
        zgrid             [n_zgrid]              1D arr
        z2use             [nz]                   1D arr
        s2                [is2[-1]]              1D arr
        is2               [nz+1]                 1D arr
        p                 [n_annot]              1D arr, n_annot = len(unique(annot))
        sb2               [n_annot]              1D arr
        s02
        annot             [is2[-1]]              1D arr
        qq_template_annot [nz * n_qq_annot]      1D arr
        n_samples
    Return:
        Nothing
    Modify:
        zcdf and zcdf_qq_annot arrays
    */
    size_t nthreads = omp_get_max_threads();
    //printf("Max number of threads in get_cdfsampling = %zu\n", nthreads);
    size_t i, j, k, tid, nnz = 0;
    // annot_nnz [n_qq_annot] 1D array
    size_t *annot_nnz;
    // s2sample_t [nthreads * n_samples] 1D array
    // zcdf_t [nthreads * n_zgrid] 1D array
    // zcdf_qq_annot_t [threads * n_zgrid * n_qq_annot] 1D array
    double *zcdf_t, *s2sample_t, *zcdf_qq_annot_t;
    // rngs [nthreads] 1D array
    gsl_rng ** rngs;

    // calloc automatically initializes array with zeros (in contrast to malloc)
    annot_nnz = calloc(n_qq_annot, sizeof(size_t));
    zcdf_t = calloc(nthreads*n_zgrid, sizeof(double));
    zcdf_qq_annot_t = calloc(nthreads*n_zgrid*n_qq_annot, sizeof(double));
    // fill zcdf and zcdf_qq_annot with zeros
    memset(zcdf, 0, n_zgrid*sizeof(double));
    memset(zcdf_qq_annot, 0, n_zgrid*n_qq_annot*sizeof(double));

    s2sample_t = (double *) malloc(nthreads*n_samples*sizeof(double*));
    rngs = (gsl_rng **) malloc(nthreads*sizeof(gsl_rng *));
    for (i=0; i<nthreads; i++) {
        gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
        // 0 seed sets to default seed, which is the same as 1 for gsl_rng_taus, so should start seeding with 1
        gsl_rng_set(r, i+1);
        rngs[i] = r;
    }

    #pragma omp parallel for default(shared) private(i,j,k,tid) schedule(static) reduction(+:nnz)
    for (i=0; i<nz; i++) {
        if (z2use[i]) {
            tid = omp_get_thread_num();
            nnz++;
            double sigmasq;
            int block_size = is2[i+1] - is2[i];
            double * mix_cdf = calloc(n_zgrid, sizeof(double));

            get_s2sample(&s2sample_t[tid*n_samples], n_samples, p, sb2, s02,
                &s2[is2[i]], &annot[is2[i]], block_size, rngs[tid]);


            for (j=0; j<n_samples; j++) {
                sigmasq = sqrt(s2sample_t[tid*n_samples + j]);
                for (k=0; k<n_zgrid; k++)
                    mix_cdf[k] += gsl_cdf_gaussian_P(zgrid[k], sigmasq);
            }

            for (j=0; j<n_zgrid; j++) {
                mix_cdf[j] /= ((double)n_samples);
                zcdf_t[tid*n_zgrid + j] += mix_cdf[j];
                for (k=0; k<n_qq_annot; k++) {
                    if ( qq_template_annot[i*n_qq_annot + k] )
                        zcdf_qq_annot_t[tid*n_zgrid*n_qq_annot + j*n_qq_annot + k] += mix_cdf[j];
                }
            }
            free(mix_cdf);
        }
    }

    // count number of variants in each qq annotation category
    for (i=0; i<nz; i++) {
        if ( z2use[i] ) {
            for (j=0; j<n_qq_annot; j++) {
                if ( qq_template_annot[i*n_qq_annot + j] )
                    annot_nnz[j] += 1;
            }
        }
    }

    // make reduction across nthreads dimention
    for (i=0; i<nthreads; i++) {
        for (j=0; j<n_zgrid; j++) {
            zcdf[j] += zcdf_t[i*n_zgrid + j];
            for (k=0; k<n_qq_annot; k++)
                zcdf_qq_annot[j*n_qq_annot + k] += zcdf_qq_annot_t[i*n_zgrid*n_qq_annot + j*n_qq_annot + k];
        }
    }

    // account for the second tail and number of variants in each category
    for (j=0; j<n_zgrid; j++) {
        zcdf[j] *= 2./((double)nnz);
        for (k=0; k<n_qq_annot; k++)
            zcdf_qq_annot[j*n_qq_annot + k] *= 2./((double)(annot_nnz[k]));
    }

    // free memory
    for(i=0; i<nthreads; i++)
        free(rngs[i]);
    free(zcdf_qq_annot_t);
    free(s2sample_t);
    free(zcdf_t);
    free(rngs);
    free(annot_nnz);
}



int main() {
    size_t i, j;

    // which tests to run
    _Bool test_cost = true;
    _Bool test_cdf = true;

    // parameters for costdirect and costsampling functions
    #define nz    3
    // size_t nz = 3;
    double z[nz] = {1.275, -0.496, 2.983};
    _Bool z2use[nz] = {true, true, true}; // {true, true, false} cost = 1.917154
                                         // {true, true, true} cost = 2.121391
    float s2[9] = {1.593, 0.934, 2.463, 0.719, 1.847, 3.012, 1.927, 0.896, 1.401}; // 3 2 4
    unsigned long long int is2[4] = {0, 3, 5, 9};
    double p[3] = {1., 1., 1.};
    double sb2[3] = {1.17, 2.03, 0.954};
    double s02 = 1.03;
    unsigned char annot[9] = {0,1,0,2,1,1,0,0,2};
    size_t n_samples = 1000;
    double cost;

    /*
    SNP   S2_eff
    0     1.593*1.17 + 0.934*2.03 + 2.463*1.17 + 1.03 = 7.671539999999999
    1     0.719*0.954 + 1.847*2.03 + 1.03 = 5.465336
    2     3.012*2.03 + 1.927*1.17 + 0.896*1.17 + 1.401*0.954 + 1.03 = 11.783824
    */

    // additional parameters for cdfsampling function
    #define n_zgrid   4
    #define n_qq_annot   3
    double zgrid[n_zgrid] = {0.0, -0.1, -0.01, -0.001};
    // qq_template_annot [nz x n_qq_annot] 2D arr
    _Bool qq_template_annot_tmp[nz][n_qq_annot] = { {true,  true,  true},
                                                    {false, true,  true},
                                                    {false,  false, true} };
    _Bool * qq_template_annot = (_Bool *)malloc(nz*n_qq_annot * sizeof(_Bool));
    for (i=0; i<nz; i++) {
        for (j=0; j<n_qq_annot; j++)
            qq_template_annot[i*n_qq_annot + j] = qq_template_annot_tmp[i][j];
    }
    

    // allocate empty arrays for cdfsampling function
    double * zcdf = (double *) malloc(n_zgrid * sizeof(double));
    double * zcdf_qq_annot = (double *) malloc(n_zgrid*n_qq_annot * sizeof(double));

    /*
    z2use = {true, true, true}
    qq_template_annot = { {true,  true,  true}, {false, true,  true}, {false,  false, true} }
    zgrid = {0.0, -0.1, -0.01, -0.001}
    zcdf = {1.0, 0.9712800037085961, 0.9971273420458443, 0.999712733546111}
    zcdf_qq_annot = { 1.0                 1.0                 1.0
                      0.971199207323047   0.9685399893037775  0.9712800037085961
                      0.9971193012706704  0.9968531741508169  0.9971273420458443
                      0.9997119295074843  0.9996853165901003  0.999712733546111 }

    z2use = {true, true, false}
    qq_template_annot = { {true,  true,  true}, {false, true,  true}, {false,  false, true} }
    zgrid = {0.0, -0.1, -0.01, -0.001}
    zcdf = {1.0, 0.9685399893037775, 0.9968531741508169, 0.9996853165901003}
    zcdf_qq_annot = { 1.0                 1.0                 1.0
                      0.971199207323047   0.9685399893037775  0.9685399893037775
                      0.9971193012706704  0.9968531741508169  0.9968531741508169
                      0.9997119295074843  0.9996853165901003  0.9996853165901003 }
    */


    size_t nthreads = omp_get_max_threads();
    printf("Max number of threads = %zu\n", nthreads);

    printf("\n");
    printf("z2use = [%d", z2use[0]);
    for (i=1; i<nz; i++)
        printf("  %d", z2use[i]);
    printf("]\n");
    printf("\n");

    if ( test_cost ) {
        printf("Test cost functions\n");
        cost = costsampling(z, z2use, nz, s2, is2, p, sb2, s02, annot, n_samples);
        printf("sampling cost: %f\n", cost);
        cost =  costdirect(z, z2use, nz, s2, is2, p, sb2, s02, annot);
        printf("direct cost: %f\n", cost);
    }

    printf("\n");

    if ( test_cdf ) {
        printf("Test cdf function\n");
        cdfsampling(zcdf, zcdf_qq_annot, zgrid, z2use, nz, n_zgrid, n_qq_annot, s2,
            is2, p, sb2, s02, annot,  qq_template_annot, n_samples);

        printf("zgrid = [%f", zgrid[0]);
        for (i=1; i<n_zgrid; i++)
            printf("  %f", zgrid[i]);
        printf("]\n");

        printf("zcdf = [%f", zcdf[0]);
        for (i=1; i<n_zgrid; i++)
            printf("  %f", zcdf[i]);
        printf("]\n");

        printf("zcdf_qq_annot = [%f", zcdf_qq_annot[0]);
        for (i=0; i<n_zgrid; i++) {
            for (j=1; j<n_qq_annot; j++)
                printf("  %f", zcdf_qq_annot[i*n_qq_annot + j]);
            if (i < n_zgrid-1)
                printf("\n                 %f", zcdf_qq_annot[(i+1)*n_qq_annot]);
        }
        printf("]\n");
    }


    return 0;
}

