#include "_cmmcost_omp.h"


double ift(double x, void * params) {
    Parameters pars = *(Parameters *) params;
    double neg_half_x2 = -0.5*x*x;
    double res = cos(pars.omega*x)*exp(neg_half_x2*pars.s02)/M_PI; // move /M_PI out
    for(size_t i=0; i<pars.n; i++) {
        res *= 1 + pars.p[i]*(exp(neg_half_x2*pars.s2[i]) - 1);
    }
    // return GSL_NAN;
    return res;
}

double ift_p_sb2_der(double x, void * params) {
    /*
    (f*g*h*q)' = f'*g*h*q + f*g'*h*q + f*g*h'*q + f*g*h*q'
       der          i0         i1         i2         i3
    */
    Parameters_der pars = *(Parameters_der *) params;
    double neg_half_x2 = -0.5*x*x;
    size_t i;
    double der_const = cos(pars.omega*x)*exp(neg_half_x2*pars.s02)/M_PI; // aggregate terms from annot categories other than current derivative (pars.der_ind)
    double der_var = 0; // aggregate terms from current annot category

    /* estimate der_var - part of function which depends on current annot category (pars.der_ind)
                 f  g  h  q
                 _  _  _  _
        forward  f' g  h  q
           |     f  g' h  q      ^
           v     f  g  h' q      |
                 f  g  h  q'  backward
    forward = [1, f, f*g, f*g*h]
    backward = [g*h*q, h*q, q, 1]
    */
    double * forward = (double *) malloc(sizeof(double)*pars.n_in_annot[pars.der_ind]);
    // go forward
    forward[0] = 1;
    for (i=0; i<pars.n_in_annot[pars.der_ind]-1; i++)
        forward[i+1] = forward[i]*(1 + pars.p[pars.der_ind]*(exp(neg_half_x2*pars.s2[pars.der_ind][i]) - 1));
    // go backward and aggreagte
    double tmp_der;
    double backward = 1;
    for (int k=pars.n_in_annot[pars.der_ind]; k--;) {
        if (pars.der_par == P_DER)
            tmp_der = exp(neg_half_x2*pars.s2[pars.der_ind][k]) - 1;
        else if (pars.der_par == SB2_DER)
            tmp_der = pars.p[pars.der_ind]*exp(neg_half_x2*pars.s2[pars.der_ind][k])*(neg_half_x2*pars.s2_no_sb2[k]);
        else {
            printf("Wrong DER code: %d\n", pars.der_par);
            exit(1);
        }
        der_var += forward[k]*backward*tmp_der;
        backward *= backward*(1 + pars.p[pars.der_ind]*(exp(neg_half_x2*pars.s2[pars.der_ind][k]) - 1));
    }
    free(forward);

    // estimate der_const - part of function which doesn't depend on current annot category
    for (i=0; i<pars.n_annot; i++) {
        if (i != pars.der_ind) {
            for (size_t j=0; j<pars.n_in_annot[i]; j++)
                der_const *= 1 + pars.p[i]*(exp(neg_half_x2*pars.s2[i][j]) - 1);
        }
    }

    return der_const*der_var;
}


double ift_s02_der(double x, void * params) {
    Parameters pars = *(Parameters *) params; // here Parameters_der are not required
    double neg_half_x2 = -0.5*x*x;
    double der = cos(pars.omega*x)*exp(neg_half_x2*pars.s02)*neg_half_x2/M_PI; // move /M_PI out
    for(size_t i=0; i<pars.n; i++) {
        der *= 1 + pars.p[i]*(exp(neg_half_x2*pars.s2[i]) - 1);
    }
    return der;
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
        FF[i].function = &ift;
        FF[i].params = &pars[i];
    }
    
    double cost = 0.;
    size_t nnz = 0;

    gsl_set_error_handler_off();
    // static dynamic
    #pragma omp parallel for default(shared) private(i,j,k,tid) schedule(dynamic) reduction(+:cost,nnz)
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


void costdirect_derivative(double * z, _Bool * z2use, size_t nz, float * s2,
    unsigned long long int * is2, double * p, double * sb2, double s02,
    unsigned char * annot, size_t n_annot, double * gradient) {
    /*
    assume we have n anootation categories and f is our likelihood function, then:
    gradient = [f_p1, f_p2, ..., f_pn, f_s21, f_s22, ..., f_s2n, f_s02],
    where f = -log_likelihood,
    f' = [-log(g)]' = -(1/g * g'), where g = likelihood
    */

    size_t i, j, k, tid;
    size_t nthreads = omp_get_max_threads();
    // printf("Max number of threads = %zu\n", nthreads);
    double result[nthreads];
    double error[nthreads];
    Parameters pars[nthreads];
    Parameters_der pars_der[nthreads];
    gsl_function FF[nthreads];
    gsl_function FF_der_p_sb2[nthreads];
    gsl_function FF_der_s02[nthreads];

    gsl_integration_workspace ** ww =
        (gsl_integration_workspace **) malloc(nthreads * sizeof(gsl_integration_workspace *));
    for (i=0; i<nthreads; i++) {
        ww[i] = gsl_integration_workspace_alloc(1000); // if size is too small raises "corrupted size vs. prev_size" error when trying to free ww[i]
        // these parameters are always the same
        pars[i].s02 = s02;
        pars_der[i].p = p;
        pars_der[i].n_annot = n_annot;
        pars_der[i].s02 = s02;

        FF[i].function = &ift; // compare here result and performance of iff and ift
        FF[i].params = &pars[i];

        FF_der_p_sb2[i].function = &ift_p_sb2_der;
        FF_der_p_sb2[i].params = &pars_der[i];

        FF_der_s02[i].function = &ift_s02_der;
        FF_der_s02[i].params = &pars[i]; // use same parameters as in simple ift
    }

    size_t nnz = 0;

    size_t gradient_len = (2*n_annot+1);
    // fill gradient with zeros
    memset(gradient, 0, gradient_len*sizeof(double));
    // array to accumulate gradient in each thread
    double * gradient_t = (double *) calloc(nthreads*gradient_len, sizeof(double));

    gsl_set_error_handler_off();

    #pragma omp parallel for default(shared) private(i,j,k,tid) schedule(dynamic) reduction(+:nnz)
    for (i=0; i<nz; i++) {
        if (z2use[i]) {
            tid = omp_get_thread_num();
            nnz++;
            int status;
            size_t block_size = is2[i+1] - is2[i];
            double * p_block = (double *) malloc(sizeof(double) * block_size);
            double * s2_block = (double *) malloc(sizeof(double) * block_size);
            // n_in_annot: number of variants in each annotation category among all variants in LD with the current
            size_t * n_in_annot = (size_t *) calloc(n_annot, sizeof(size_t));

            for (j=0; j<block_size; j++) {
                k = is2[i] + j;
                p_block[j] = p[annot[k]];
                s2_block[j] = sb2[annot[k]] * s2[k];
                n_in_annot[annot[k]]++;
            }
            pars[tid].p = p_block;
            pars[tid].s2 = s2_block;
            pars[tid].n = block_size;
            pars[tid].omega = z[i];

            // estimate likelihood
            status = gsl_integration_qagiu (&FF[tid], 0, 1.E-7, 1.E-5, 1000,
                ww[tid], &result[tid], &error[tid]);

            double likelihood_z = 0;
            if (status) {
                printf ("gsl error (status code: %d) in likelihood_z: %s\n", status, gsl_strerror (status));
            }
            else if (result[tid] > 0) {
                likelihood_z = result[tid];
            }

            // don't change gradient if likelihood = 0
            if ( likelihood_z == 0) continue;

            // estimate s02 derivative of the likelihood. Parameters are the same as in likelihood estimation
            status = gsl_integration_qagiu (&FF_der_s02[tid], 0, 1.E-7, 1.E-5, 1000,
                ww[tid], &result[tid], &error[tid]);
            double likelihood_z_s02 = 0;
            if (status) {
                printf ("gsl error (status code: %d) in likelihood_z: %s\n", status, gsl_strerror (status));
            }
            else {
                likelihood_z_s02 = result[tid];
            }
            gradient_t[tid*gradient_len+2*n_annot] += -likelihood_z_s02/likelihood_z;

            /*
            Estimate p and sb2 derivatives of likelihood.

            Construct p_in_annot and s2_in_annot containing p and s2 values for each annotation separately
            n_in_annot = [2, 0, 3], 2 - number of variants belonging to the 1st annot, 0 - 2nd annot, 3 - 3rd annot
            s2_in_annot = [ [s2_1, s2_2],      : s2 for varints in the first annot
                            NULL,              : NULL, since current variant is not in LD with any variant from the second annot
                            [s2_1, s2_2, s2_3] : s2 for varints in the third annot
                          ]
            */
            double ** s2_in_annot = (double **) malloc(sizeof(double *)*n_annot);
            double ** s2_no_sb2_annot = (double **) malloc(sizeof(double *)*n_annot);
            for (j=0; j<n_annot; j++) {
                if (n_in_annot[j] > 0) {
                    s2_in_annot[j] = (double *) malloc(sizeof(double)*n_in_annot[j]);
                    s2_no_sb2_annot[j] = (double *) malloc(sizeof(double)*n_in_annot[j]);
                }
                else {
                    s2_in_annot[j] = NULL;
                    s2_no_sb2_annot[j] = NULL;
                }
                
            }
            size_t * ind_in_annot = calloc(n_annot, sizeof(size_t)); // init all elelments with 0
            for (j=0; j<block_size; j++) {
                k = is2[i] + j;
                s2_in_annot[annot[k]][ind_in_annot[annot[k]]] = sb2[annot[k]] * s2[k];
                s2_no_sb2_annot[annot[k]][ind_in_annot[annot[k]]] = s2[k];
                ind_in_annot[annot[k]]++;
            }

            pars_der[tid].s2 = s2_in_annot;
            pars_der[tid].n_in_annot = n_in_annot;
            pars_der[tid].omega = z[i];

            for (j=0; j<n_annot; j++) {
                if (n_in_annot[j] > 0) { // if n_in_annot[j] == 0, current LD block doesn't depend on j-th annotation
                    pars_der[tid].s2_no_sb2 = s2_no_sb2_annot[j];
                    pars_der[tid].der_ind = j;

                    // estimate p derivative
                    pars_der[tid].der_par = P_DER;
                    status = gsl_integration_qagiu (&FF_der_p_sb2[tid], 0, 1.E-7, 1.E-5, 1000,
                        ww[tid], &result[tid], &error[tid]);
                    double likelihood_z_p = 0;
                    if (status) {
                        printf ("gsl error (status code: %d) in likelihood_z_p: %s\n", status, gsl_strerror (status));
                    }
                    else {
                        likelihood_z_p = result[tid];
                    }

                    gradient_t[tid*gradient_len+j] += -likelihood_z_p/likelihood_z;

                    // estimate s2 derivative
                    pars_der[tid].der_par = SB2_DER;
                    status = gsl_integration_qagiu (&FF_der_p_sb2[tid], 0, 1.E-7, 1.E-5, 1000,
                        ww[tid], &result[tid], &error[tid]);
                    double likelihood_z_sb2 = 0;
                    if (status) {
                        printf ("gsl error (status code: %d) in likelihood_z_sb2: %s\n", status, gsl_strerror (status));
                    }
                    else {
                        likelihood_z_sb2 = result[tid];
                    }

                    gradient_t[tid*gradient_len+j+n_annot] += -likelihood_z_sb2/likelihood_z;
                }
            }

            free(ind_in_annot);
            for (j=0; j<n_annot; j++) {
                free(s2_in_annot[j]);
                free(s2_no_sb2_annot[j]);
            }
            free(s2_no_sb2_annot);
            free(s2_in_annot);
            free(n_in_annot);
            free(p_block);
            free(s2_block);
         }
    }

    // free memory
    for (i=0; i<nthreads; i++) {
        gsl_integration_workspace_free(ww[i]);
    }
    free(ww);

    // make reduction across nthreads dimention
    for (i=0; i<nthreads; i++) {
        for (j=0; j<gradient_len; j++)
            gradient[j] += gradient_t[i*gradient_len + j];
    }

    for (i=0; i<2*n_annot+1; i++)
        gradient[i] /= (double)nnz;

    free(gradient_t);
}



int main() {
    size_t i, j;

    // which tests to run
    _Bool test_cost = true;
    _Bool test_cdf = true;
    _Bool test_der = true;

    // parameters for costdirect and costsampling functions
    #define nz    3
    #define number_of_annot 3
    // size_t nz = 3;
    double z[nz] = {1.275, -0.496, 2.983};
    _Bool z2use[nz] = {true, true, true}; // {true, true, false} cost = 1.917154
                                         // {true, true, true} cost = 2.121391
    float s2[9] = {1.593, 0.934, 2.463, 0.719, 1.847, 3.012, 1.927, 0.896, 1.401}; // 3 2 4
    unsigned long long int is2[4] = {0, 3, 5, 9};
    double p[number_of_annot] = {1., 1., 1.};
    double sb2[number_of_annot] = {1.17, 2.03, 0.954};
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


    if (test_der) {
        /*
        z2use = {true, true, true}, should result in (checked analytically for sb2 and s02 derivatives, but not p):
        gradient = [0.108350  0.257204  0.025955  0.079223  0.080213  0.025791  0.049708]
        */
        printf("Running gradient test\n");
        double * gradient = (double *) calloc(2*number_of_annot+1, sizeof(double));
        costdirect_derivative(z, z2use, nz, s2, is2, p, sb2, s02, annot, number_of_annot, gradient);
        printf("gradient = [%f", gradient[0]);
        for (i=1; i<2*number_of_annot+1; i++)
            printf("  %f", gradient[i]);
        printf("]\n");
        free(gradient);
    }


    return 0;
}

