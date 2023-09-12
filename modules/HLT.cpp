#define HLT_C
#include "HLT.hpp"
#include "global.hpp"
#include "tower.hpp"
#include <gsl/gsl_integration.h>

#include <stdio.h>
// #include <mpir.h>


HLT_type::HLT_type(int tmax, int L0, double E0_, int njack, HLT_b type_b, double alpha) {
    Tmax = tmax;
    T = L0;
    Njack = njack;
    E0 = E0_;
    R = (double*)malloc(sizeof(double) * Tmax);
    A = malloc_2<double>(Tmax, Tmax);
    type = type_b;
    if (type == HLT_EXP_b) {
        for (size_t t = 0; t < Tmax; t++) {
            R[t] = 1.0 / (t + 1.0) + 1.0 / (T - t - 1.0);
            for (size_t r = 0; r < Tmax; r++) {
                A[t][r] = exp(-(r + t + 2 - alpha) * E0) / (r + t + 2 - alpha);
                A[t][r] += exp(-(T - r + t - alpha) * E0) / (T - r + t - alpha);
                A[t][r] += exp(-(T + r - t - alpha) * E0) / (T + r - t - alpha);
                A[t][r] += exp(-(2 * T - r - t - 2 - alpha) * E0) / (2 * T - r - t - 2 - alpha);
            }
        }
    }
}

HLT_type::~HLT_type() {
    free(R);
    for (size_t i = 0; i < Tmax; i++) {
        free(A[i]);
    }
    free(A);
}



double integrand_f(double x, void* params) {
    wrapper_smearing* p = (wrapper_smearing*)params;
    int t = p->t;
    int T = p->HLT->T;
    return p->function(x, p->params) * (exp(-(t + 1) * x) + exp(-(T - (t + 1)) * x));
}


void HLT_type::compute_f_EXP_b(wrapper_smearing  Delta, double epsrel) {

    f = (double*)malloc(sizeof(double) * Tmax);
    int Maxiter = 1e+6;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(Maxiter);
    double result, error;
    gsl_function F;
    Delta.HLT = this;
    for (int t = 0;t < Tmax;t++) {
        F.function = &integrand_f;
        Delta.t = t;
        F.params = &Delta;

        // gsl_integration_qags(&F, 0, 1, 0, epsrel, Maxiter, w, &result, &error);
        // int_E0^\infty  b* Delta
        gsl_integration_qagiu(&F, E0, 0, epsrel, Maxiter, w, &result, &error);
        f[t] = result;
    }
    gsl_integration_workspace_free(w);

}

double** HLT_type::HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
    FILE* outfile,   const char* description, wrapper_smearing  Delta, FILE* file_jack) {

    double** r = malloc_2<double>(Tmax, Njack);
    for (int t = 0;t < Tmax;t++)
        for (int j = 0;j < Njack;j++)
            r[t][j] = conf_jack[j][id][t][0];


    double** cov = myres->comp_cov(Tmax, r);
    if (type == HLT_EXP_b) {
        compute_f_EXP_b(Delta);
    }

    HLT_out res(Delta.lambdas);
    
    for (int il =0;il<Delta.lambdas.size(); il++){


    }
    // mpf_t **W=(mpf_t**) malloc(sizeof(mpf_t*)* Tmax);
    double **p;
    return p;
};
