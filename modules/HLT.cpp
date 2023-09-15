#define HLT_C
#include "HLT.hpp"
#include "global.hpp"
#include "tower.hpp"
#include <gsl/gsl_integration.h>

#include <stdio.h>
// #include <mpir.h>

#include "arb_calc.h"


HLT_type_d::HLT_type_d(int tmax, int L0, double E0_, int njack, HLT_b type_b, double alpha) {
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

HLT_type_d::~HLT_type_d() {
    free(R);
    for (size_t i = 0; i < Tmax; i++) {
        free(A[i]);
    }
    free(A);
    if (f_allocated) free(f);
}



double integrand_f(double x, void* params) {
    wrapper_smearing_d* p = (wrapper_smearing_d*)params;
    int t = p->t;
    int T = p->HLT->T;
    return p->function(x, p->params) * (exp(-(t + 1) * x) + exp(-(T - (t + 1)) * x));
}


void HLT_type_d::compute_f_EXP_b(wrapper_smearing_d& Delta, double epsrel) {

    if (f_allocated == true) {
        printf("HLT: recomputing f\n");
        free(f);
    }
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
    f_allocated = true;
}

//////////////////////////////////////////////////////////////////
// HLT  smearing functions
//////////////////////////////////////////////////////////////////
double theta_s1_d(double x, double* p) {
    return 1.0 / (1 + exp(-x / p[0]));
};

double gaussian_for_HLT_d(double x, double* p) {
    return (exp(-(x - p[0]) * (x - p[0]) / (2 * p[1] * p[1]))) / (p[1] * sqrt(2 * M_PIl));
};

//////////////////////////////////////////////////////////////////


double** HLT_type_d::HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
    FILE* outfile, const char* description, wrapper_smearing_d  Delta, FILE* file_jack) {

    double** r = malloc_2<double>(Tmax, Njack);
    for (int t = 0;t < Tmax;t++)
        for (int j = 0;j < Njack;j++)
            r[t][j] = conf_jack[j][id][t][0];


    double** cov = myres->comp_cov(Tmax, r);
    if (type == HLT_EXP_b) {
        compute_f_EXP_b(Delta);
    }
    else { printf("error HLT: type not supported\n");exit(1); }

    HLT_out res(Delta.lambdas);

    for (int il = 0;il < Delta.lambdas.size(); il++) {


    }
    // mpf_t **W=(mpf_t**) malloc(sizeof(mpf_t*)* Tmax);
    double** p;
    return p;
};

////////////////////////////////////////////////////////////////////////////////////////////////
// using ARB
////////////////////////////////////////////////////////////////////////////////////////////////



HLT_type::HLT_type(int tmax, int L0, double E0_, HLT_b type_b, int prec_, double alpha) {
    Tmax = tmax;
    T = L0;
    prec = prec_;
    arb_init(E0);
    arb_set_d(E0, E0_);

    arb_mat_init(R, Tmax, 1);
    arb_mat_init(A, Tmax, Tmax);

    type = type_b;

    arb_t at, aT, arb_alpha;
    arb_init(at);
    arb_init(aT);
    arb_init(arb_alpha);
    arb_set_d(arb_alpha, alpha);

    if (type == HLT_EXP_b) {
        for (size_t t = 0; t < Tmax; t++) {

            arb_set_ui(at, t + 1);
            arb_inv(at, at, prec);

            arb_set_ui(aT, T - t - 1);
            arb_inv(aT, aT, prec);

            // R[t]= 1.0 / (t + 1.0) + 1.0 / (T - t - 1.0);
            arb_add(arb_mat_entry(R, t, 0), at, aT, prec);
            for (size_t r = 0; r < Tmax; r++) {

                // A[t][r] = exp(-(r + t + 2 - alpha) * E0) / (r + t + 2 - alpha);
                arb_set_ui(at, r + t + 2);
                arb_sub(at, at, arb_alpha, prec);

                arb_mul(aT, at, E0, prec);
                arb_neg(aT, aT);
                arb_exp(aT, aT, prec);
                arb_div(arb_mat_entry(A, t, r), aT, at, prec);
                //     A[t][r] += exp(-(T - r + t - alpha) * E0) / (T - r + t - alpha);
                arb_set_ui(at, T - r + t);
                arb_sub(at, at, arb_alpha, prec);

                arb_mul(aT, at, E0, prec);
                arb_neg(aT, aT);
                arb_exp(aT, aT, prec);
                arb_div(aT, aT, at, prec);
                arb_add(arb_mat_entry(A, t, r), arb_mat_entry(A, t, r), aT, prec);
                //     A[t][r] += exp(-(T + r - t - alpha) * E0) / (T + r - t - alpha);
                arb_set_ui(at, T + r - t);
                arb_sub(at, at, arb_alpha, prec);

                arb_mul(aT, at, E0, prec);
                arb_neg(aT, aT);
                arb_exp(aT, aT, prec);
                arb_div(aT, aT, at, prec);
                arb_add(arb_mat_entry(A, t, r), arb_mat_entry(A, t, r), aT, prec);
                //     A[t][r] += exp(-(2 * T - r - t - 2 - alpha) * E0) / (2 * T - r - t - 2 - alpha);
                arb_set_ui(at, 2 * T - r - t - 2);
                arb_sub(at, at, arb_alpha, prec);

                arb_mul(aT, at, E0, prec);
                arb_neg(aT, aT);
                arb_exp(aT, aT, prec);
                arb_div(aT, aT, at, prec);
                arb_add(arb_mat_entry(A, t, r), arb_mat_entry(A, t, r), aT, prec);

            }
        }
    }
}

HLT_type::~HLT_type() {
    // for (size_t i = 0; i < Tmax; i++) {
    //     arb_clear(R[i]);
    //     if (f_allocated)arb_clear(f[i]);
    //     for (size_t r = 0; r < Tmax; r++)         arb_clear(A[i][r]);
    //     free(A[i]);
    // }
    // free(R);
    // free(A);
    arb_mat_clear(A);
    arb_mat_clear(R);
    if (f_allocated) arb_mat_clear(f);
}


int gaussian_for_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    // acb_sin(res, z, prec);
    arb_t* p = (arb_t*)param;
    arb_t sqrt2pi;
    arb_init(sqrt2pi);
    arb_set_d(sqrt2pi, sqrt(2 * M_PIl));


    acb_sub_arb(res, z, p[0], prec);
    acb_mul(res, res, res, prec);
    acb_div_ui(res, res, 2, prec);
    acb_div_arb(res, res, p[1], prec);
    acb_div_arb(res, res, p[1], prec);
    acb_neg(res, res);
    acb_exp(res, res, prec);
    acb_div_arb(res, res, p[1], prec);
    acb_div_arb(res, res, sqrt2pi, prec);

    // return (exp(-(x - p[0]) * (x - p[0]) / (2 * p[1] * p[1]))) / (p[1] * sqrt(2 * M_PIl));
    return 0;
}

int theta_s1_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    //  return 1.0 / (1 + exp(-x / p[0]));
    arb_t* p = (arb_t*)param;
    acb_div_arb(res, z, p[0], prec);
    acb_neg(res, res);
    acb_exp(res, res, prec);
    acb_add_ui(res, res, 1, prec);
    acb_inv(res, res, prec);

    return 0;
}



int no_smearing_for_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */
    acb_set_ui(res, 1);
    return 0;
}


int integrand_f(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    wrapper_smearing* p = (wrapper_smearing*)param;

    p->function(res, z, p->params, order, prec);

    acb_t b; acb_init(b);
    acb_set_ui(b, p->t + 1);
    acb_mul(b, b, z, prec);
    acb_neg(b, b);
    acb_exp(b, b, prec);

    acb_t bT; acb_init(bT);
    acb_set_ui(bT, p->HLT->T - p->t - 1);
    acb_mul(bT, bT, z, prec);
    acb_neg(bT, bT);
    acb_exp(bT, bT, prec);

    acb_add(b, b, bT, prec);
    acb_mul(res, res, b, prec);

    acb_clear(b);
    acb_clear(bT);
    return 0;
    // return p->function(x, p->params) * (exp(-(t + 1) * x) + exp(-(T - (t + 1)) * x));
}

void HLT_type::compute_f_EXP_b(wrapper_smearing& Delta, double epsrel) {

    if (f_allocated == true) {
        printf("HLT: recomputing f\n");
        arb_mat_clear(f);
    }
    arb_mat_init(f, Tmax, 1);
    int Maxiter = 1e+6;
    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

    options->deg_limit = 100;
    options->eval_limit = 100000;
    options->depth_limit = 10000;
    options->verbose = 1;
    // options->use_heap = 1;

    mag_t tol;
    mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -prec);
    slong  goal = prec;
    // Delta.HLT = this;
    acb_t a, b, s;
    acb_init(a);
    acb_init(b);
    acb_init(s);

    for (int t = 0;t < Tmax;t++) {
        // F.function = &integrand_f;
        // Delta.t = t;
        // F.params = &Delta;

        // // gsl_integration_qags(&F, 0, 1, 0, epsrel, Maxiter, w, &result, &error);
        // // int_E0^\infty  b* Delta
        // gsl_integration_qagiu(&F, E0, 0, epsrel, Maxiter, w, &result, &error);
        // f[t] = result;
        Delta.t = t;

        void* param = (void*)&Delta;
        wrapper_smearing* p = (wrapper_smearing*)param;

        acb_set_arb(a, E0);
        acb_set_d(b, 10000);
        acb_calc_integrate(s, integrand_f, param, a, b, goal, tol, options, prec);
        acb_get_real(arb_mat_entry(f, t, 0), s);
    }
    f_allocated = true;
}



double** HLT_type::HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
    FILE* outfile, const char* description, wrapper_smearing  Delta, FILE* file_jack, fit_type_HLT fit_info) {

    int Njack = fit_info.Njack;
    int id = fit_info.corr_id[0];
    double** r = malloc_2<double>(Tmax, Njack);
    for (int t = 0;t < Tmax;t++)
        for (int j = 0;j < Njack;j++)
            r[t][j] = conf_jack[j][id][t][0];


    double** cov = myres->comp_cov(Tmax, r);

    myres->comp_cov_arb(W, Tmax, r,prec);

    for (int t = 0;t < Tmax;t++) {
        for (int r = 0;r < Tmax;r++) {
        // arb_set_d(arb_mat_entry(W,t,r), cov[t][r]);
        printf("cov %d %d = %.12g\n",t,r,cov[t][r]);
        arb_printn(arb_mat_entry(W,t,r), prec / 3.33, 0); flint_printf("\n");
        }
    }

    if (f_allocated == false) {
        printf("HLT: recomputing f\n");
        compute_f_EXP_b(Delta);
    }

    HLT_out res(fit_info.lambdas);

    for (int il = 0;il < fit_info.lambdas.size(); il++) {


    }
    // mpf_t **W=(mpf_t**) malloc(sizeof(mpf_t*)* Tmax);
    double** p;
    return p;
};

