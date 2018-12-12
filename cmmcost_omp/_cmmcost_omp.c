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

    for (i=0; i<nthreads; i++) {
        gsl_integration_workspace_free(ww[i]);
    }
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
    unsigned char * annot, int n_samples) {

    size_t nthreads = omp_get_max_threads();
    // printf("Max number of threads in costsampling = %zu\n", nthreads);

    size_t i, j, tid;
    double cost = 0.;
    size_t nnz = 0;

    // s2samples [nthreads x n_samples] 2D array
    double **s2samples = (double **)malloc(nthreads * sizeof(double *)); 
    gsl_rng ** rngs = (gsl_rng **) malloc(nthreads * sizeof(gsl_rng *));
    for (i=0; i<nthreads; i++) {

        s2samples[i] = (double *)malloc(n_samples * sizeof(double)); 

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

            get_s2sample(s2samples[tid], n_samples, p, sb2, s02, &s2[is2[i]],
                &annot[is2[i]], block_size, rngs[tid]);

            for (j=0; j<n_samples; j++) {
                // if (s2samples[tid][j] > 0) - this should not be the case since s20 is always > 0
                mix_pdf += gsl_ran_gaussian_pdf(z[i], sqrt(s2samples[tid][j]));
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
    for(i=0; i<nthreads; i++) {
        free(s2samples[i]);
        free(rngs[i]);
    }
    free(rngs);
    free(s2samples);

    return cost/(double)nnz; 
}


int main() {
    
    double z[3] = {1.275, -0.496, 2.983};
    _Bool z2use[3] = {true, true, false}; // {true, true, false} cost = 1.917154
                                         // {true, true, true} cost = 2.121391
    size_t nz = 3;
    float s2[9] = {1.593, 0.934, 2.463, 0.719, 1.847, 3.012, 1.927, 0.896, 1.401}; // 3 2 4
    unsigned long long int is2[4] = {0, 3, 5, 9};
    double p[3] = {1., 1., 1.};
    double sb2[3] = {1.17, 2.03, 0.954};
    double s02 = 1.03;
    unsigned char annot[9] = {0,1,0,2,1,1,0,0,2};
    int n_samples = 1000;
    double cost;

    size_t nthreads = omp_get_max_threads();
    printf("Max number of threads = %zu\n", nthreads);

    cost = costsampling(z, z2use, nz, s2, is2, p, sb2, s02, annot, n_samples);
    printf("sampling cost: %f\n", cost);

    cost =  costdirect(z, z2use, nz, s2, is2, p, sb2, s02, annot);
    printf("direct cost: %f\n", cost);

    return 0;
}

