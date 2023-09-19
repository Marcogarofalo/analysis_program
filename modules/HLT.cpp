#define HLT_C
#include "HLT.hpp"
#include "global.hpp"
#include "tower.hpp"
#include <gsl/gsl_integration.h>

#include <stdio.h>
// #include <mpir.h>
#include "myarb.hpp"
#include "arb_calc.h"
#include "mutils.hpp"


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


double gaussian_for_HLT_d(double x, double* p) {
    return (exp(-(x - p[0]) * (x - p[0]) / (2 * p[1] * p[1]))) / (p[1] * sqrt(2 * M_PIl));
};
double theta_for_HLT_d(double x, double* p) {
    return 1.0 / (1.0 + exp(-(p[0] - x) / p[1]));
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

HLT_type::HLT_type(HLT_type_input info_) {
    info = info_;
    error(info.tmax < 0, 1, "HLT_type", "tmax not set");
    error(info.T < 0, 1, "HLT_type", "T not set");
    error(info.prec < 0, 1, "HLT_type", "prec (precision) not set");
    error(info.E0 < 0, 1, "HLT_type", "E0 not set");
    error(info.type_b == HLT_INVALID_b, 1, "HLT_type", "type_b not set");
    arb_init(E0_arb);
    arb_set_d(E0_arb, info.E0);

    if (info.normalize_kernel)  arb_mat_init(R, info.tmax, 1);
    arb_mat_init(A, info.tmax, info.tmax);


    arb_t at, aT, arb_alpha;
    arb_init(at);
    arb_init(aT);
    arb_init(arb_alpha);
    arb_set_d(arb_alpha, info.alpha);

    int prec = info.prec;
    for (size_t t = 0; t < info.tmax; t++) {
        if (info.normalize_kernel) {
            if (info.type_b == HLT_EXP_b) {
                arb_set_ui(at, t + 1);
                arb_inv(arb_mat_entry(R, t, 0), at, prec);
            }
            else if (info.type_b == HLT_EXP_bT) {
                arb_set_ui(at, t + 1);
                arb_inv(arb_mat_entry(R, t, 0), at, prec);
                arb_set_ui(aT, info.T - t - 1);
                arb_inv(aT, aT, prec);
                arb_add(arb_mat_entry(R, t, 0), at, aT, prec);
            }
            else { printf("HLT b_T type not supported\n"); }

        }
        for (size_t r = 0; r < info.tmax; r++) {
            if (info.type_b == HLT_EXP_b) {
                // A[t][r] = exp(-(r + t + 2 - alpha) * E0) / (r + t + 2 - alpha);
                arb_set_ui(at, r + t + 2);
                arb_sub(at, at, arb_alpha, prec);
                arb_mul(aT, at, E0_arb, prec);
                arb_neg(aT, aT);
                arb_exp(aT, aT, prec);
                arb_div(arb_mat_entry(A, t, r), aT, at, prec);
            }
            else if (info.type_b == HLT_EXP_bT) {
                int& T = info.T;
                // A[t][r] = exp(-(r + t + 2 - alpha) * E0) / (r + t + 2 - alpha);
                arb_set_ui(at, r + t + 2);
                arb_sub(at, at, arb_alpha, prec);
                arb_mul(aT, at, E0_arb, prec);
                arb_neg(aT, aT);
                arb_exp(aT, aT, prec);
                arb_div(arb_mat_entry(A, t, r), aT, at, prec);
                //     A[t][r] += exp(-(T - r + t - alpha) * E0) / (T - r + t - alpha);
                arb_set_ui(at, T - r + t);
                arb_sub(at, at, arb_alpha, prec);
                arb_mul(aT, at, E0_arb, prec);
                arb_neg(aT, aT);
                arb_exp(aT, aT, prec);
                arb_div(aT, aT, at, prec);
                arb_add(arb_mat_entry(A, t, r), arb_mat_entry(A, t, r), aT, prec);
                //     A[t][r] += exp(-(T + r - t - alpha) * E0) / (T + r - t - alpha);
                arb_set_ui(at, T + r - t);
                arb_sub(at, at, arb_alpha, prec);
                arb_mul(aT, at, E0_arb, prec);
                arb_neg(aT, aT);
                arb_exp(aT, aT, prec);
                arb_div(aT, aT, at, prec);
                arb_add(arb_mat_entry(A, t, r), arb_mat_entry(A, t, r), aT, prec);
                //     A[t][r] += exp(-(2 * T - r - t - 2 - alpha) * E0) / (2 * T - r - t - 2 - alpha);
                arb_set_ui(at, 2 * T - r - t - 2);
                arb_sub(at, at, arb_alpha, prec);
                arb_mul(aT, at, E0_arb, prec);
                arb_neg(aT, aT);
                arb_exp(aT, aT, prec);
                arb_div(aT, aT, at, prec);
                arb_add(arb_mat_entry(A, t, r), arb_mat_entry(A, t, r), aT, prec);
            }
        }
    }

    arb_mat_init(W, info.tmax, info.tmax);
    arb_clear(at);
    arb_clear(aT);
}

HLT_type::~HLT_type() {
    arb_mat_clear(A);
    if (info.normalize_kernel)  arb_mat_clear(R);
    arb_mat_clear(W);
    if (f_allocated) arb_mat_clear(f);
}

void HLT_type::compute_b(acb_t b, int  t, const acb_t E0) {
    int& prec = info.prec;
    if (info.type_b == HLT_EXP_b) {
        acb_set_ui(b, t + 1);
        acb_mul(b, b, E0, prec);
        acb_neg(b, b);
        acb_exp(b, b, prec);
    }
    if (info.type_b == HLT_EXP_bT) {
        acb_set_ui(b, t + 1);
        acb_mul(b, b, E0, prec);
        acb_neg(b, b);
        acb_exp(b, b, prec);

        acb_t bT; acb_init(bT);
        acb_set_ui(bT, info.T - t - 1);
        acb_mul(bT, bT, E0, prec);
        acb_neg(bT, bT);
        acb_exp(bT, bT, prec);
        acb_add(b, b, bT, prec);
        acb_clear(bT);
    }

}


int gaussian_for_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    // acb_sin(res, z, prec);
    arb_t* p = (arb_t*)param;
    // arb_t sqrt2pi;
    // arb_init(sqrt2pi);
    // arb_set_d(sqrt2pi, sqrt(2 * M_PIl));


    acb_sub_arb(res, z, p[0], prec);
    acb_mul(res, res, res, prec);
    acb_div_ui(res, res, 2, prec);
    acb_div_arb(res, res, p[1], prec);
    acb_div_arb(res, res, p[1], prec);
    acb_neg(res, res);
    acb_exp(res, res, prec);
    // acb_div_arb(res, res, p[1], prec);
    // acb_div_arb(res, res, sqrt2pi, prec);

    // return (exp(-(x - p[0]) * (x - p[0]) / (2 * p[1] * p[1]))) / (p[1] * sqrt(2 * M_PIl));
    return 0;
}

int theta_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    //  return 1.0 / (1 + exp(-(x-p[0]) / p[1]));
    arb_t* p = (arb_t*)param;
    acb_sub_arb(res, z, p[0], prec);
    acb_neg(res, res);
    acb_div_arb(res, res, p[1], prec);
    acb_neg(res, res);
    acb_exp(res, res, prec);
    acb_add_ui(res, res, 1, prec);
    acb_inv(res, res, prec);

    return 0;
}



int no_smearing_for_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */
    return 0;
}


int integrand_f(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    wrapper_smearing* p = (wrapper_smearing*)param;
    acb_set_ui(res, 1);
    p->function(res, z, p->params, order, prec);
    if (p->HLT->info.normalize_kernel) acb_mul_arb(res, res, p->Norm, prec);

    acb_t b; acb_init(b);
    p->HLT->compute_b(b, p->t, z);

    acb_mul(res, res, b, prec);

    // * e^(alpha E)
    acb_set_d(b, p->HLT->info.alpha);
    acb_mul(b, b, z, prec);
    acb_exp(b, b, prec);
    acb_mul(res, res, b, prec);

    acb_clear(b);

    // we could use the c++ struct but it is slightly slower
    // myacb E(z, prec);
    // E = exp(-E * (p->t + 1)) + exp(-E * (p->HLT->T - p->t - 1));
    // acb_mul(res, res, E.a, prec);
    return 0;
}

wrapper_smearing::wrapper_smearing(int  (*f)(acb_ptr, const acb_t, void*, slong, slong), std::vector<double> p, HLT_type* HLT_) {
    function = f;
    Np = p.size();
    // params = (arb_t*)malloc(sizeof(arb_t) * Np);
    params = _arb_vec_init(Np);
    for (int i = 0;i < Np;i++) {
        // arb_init(params[i]);
        arb_set_d((params + i), p[i]);
    }
    HLT = HLT_;
    arb_init(Norm);
    if (HLT->info.normalize_kernel)
        normilise_smearing();
    else
        arb_set_ui(Norm, 1);
};

void wrapper_smearing::normilise_smearing() {
    int Maxiter = 1e+6;
    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

    options->deg_limit = 100;
    options->eval_limit = 100000;
    options->depth_limit = 10000;
    options->verbose = 0;
    mag_t tol;
    mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -HLT->info.prec);
    slong  goal = HLT->info.prec;
    acb_t a, b, s;
    acb_init(a);
    acb_init(b);
    acb_init(s);

    acb_set_ui(a, 0);
    acb_set_d(b, 10000);
    acb_calc_integrate(s, function, (void*)params, a, b, goal, tol, options, HLT->info.prec);
    acb_get_real(Norm, s);
    printf("Normalization for smearing function "); arb_printn(Norm, HLT->info.prec / 3.33, 0); flint_printf("\n");

    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
};


void HLT_type::compute_f_EXP_b(wrapper_smearing& Delta) {

    if (f_allocated == true) {
        printf("HLT: recomputing f\n");
        arb_mat_clear(f);
    }
    arb_mat_init(f, info.tmax, 1);
    int Maxiter = 1e+6;
    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

    options->deg_limit = 100;
    options->eval_limit = 100000;
    options->depth_limit = 10000;
    options->verbose = 0;
    // options->use_heap = 1;

    mag_t tol;
    mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -info.prec);
    slong  goal = info.prec;
    // Delta.HLT = this;
    acb_t a, b, s;
    acb_init(a);
    acb_init(b);
    acb_init(s);

    for (int t = 0;t < info.tmax;t++) {

        void* param = (void*)&Delta;
        Delta.t = t;
        acb_set_arb(a, E0_arb);
        acb_set_d(b, 10000);
        acb_calc_integrate(s, integrand_f, param, a, b, goal, tol, options, info.prec);
        acb_get_real(arb_mat_entry(f, t, 0), s);
    }
    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
    f_allocated = true;
}


void HLT_type::check_reconstruction(wrapper_smearing& Delta, arb_mat_t g, std::array<double, 3> range) {
    acb_t res;
    acb_init(res);
    arb_t res_re;
    arb_init(res_re);
    arb_t res_HLT;
    arb_init(res_HLT);
    acb_t E;
    acb_init(E);
    arb_t E_re;
    arb_init(E_re);
    arb_t b; arb_init(b);
    arb_t bT; arb_init(bT);
    double dh = (range[1] - range[0]) / range[2];
    printf("check HLT  g: \n");
    int& prec = info.prec;
    for (int t = 0;t < info.tmax;t++) {
        printf("t= %d   g= ", t);
        arb_printn(arb_mat_entry(g, t, 0), prec / 3.33, 0); flint_printf("\n");
    }
    for (int i = 0;i < range[2] + 1;i++) {
        arb_set_d(E_re, range[0] + i * dh);
        acb_set_arb(E, E_re);
        Delta.function(res, E, (void*)Delta.params, 0, prec);
        if (info.normalize_kernel) acb_mul_arb(res, res, Delta.Norm, prec);

        acb_get_real(res_re, res);
        double d = std::stod(arb_get_str(res_re, 15, 10));
        // printf("smearing function: ");
        // arb_printn(res_re, prec / 3.33, 0); flint_printf("\n");
        arb_set_ui(res_HLT, 0);
        for (int t = 0;t < info.tmax;t++) {
            arb_set_ui(b, t + 1);
            arb_mul(b, b, E_re, prec);
            arb_neg(b, b);
            arb_exp(b, b, prec);

            arb_set_ui(bT, Delta.HLT->info.T - t - 1);
            arb_mul(bT, bT, E_re, prec);
            arb_neg(bT, bT);
            arb_exp(bT, bT, prec);

            arb_add(b, b, bT, prec);
            arb_addmul(res_HLT, arb_mat_entry(g, t, 0), b, prec);

            // myarb mE(E_re, prec);
            // mE = exp(-mE * (t + 1)) + exp(-mE * (T - t - 1));
            // arb_addmul(res_HLT, arb_mat_entry(g, t, 0), mE.a, prec);

        }
        // printf("recostructed function: ");
        // arb_printn(res_HLT, prec / 3.33, 0); flint_printf("\n");
        double d1 = std::stod(arb_get_str(res_HLT, 15, 10));
        arb_sub(b, res_re, res_HLT, prec);
        double diff = std::stod(arb_get_str(b, 15, 10));
        printf("%.12g  %.12g %.12g  %.12g\n", range[0] + i * dh, d, d1, diff);
    }
    acb_clear(res);
    arb_clear(res_re);
    acb_clear(E);
    arb_clear(E_re);
    arb_clear(res_HLT);
    arb_clear(b);
    arb_clear(bT);

}

double** HLT_type::HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
    FILE* outfile, const char* description, wrapper_smearing& Delta, FILE* file_jack, fit_type_HLT fit_info) {

    int& Tmax = info.tmax;
    int& prec = info.prec;
    int Njack = fit_info.Njack;
    int id = fit_info.corr_id[0];
    double** r = malloc_2<double>(info.tmax, Njack);
    for (int t = 0;t < Tmax;t++)
        for (int j = 0;j < Njack;j++)
            r[t][j] = conf_jack[j][id][t][0];

    // double** cov = myres->comp_cov(Tmax, r);
    myres->comp_cov_arb(W, Tmax, r, prec);

    if (f_allocated == false) {
        printf("HLT: recomputing f\n");
        compute_f_EXP_b(Delta);
    }

    HLT_out res(fit_info.lambdas);

    arb_t lam;
    arb_init(lam);
    arb_mat_t Wl;
    arb_mat_init(Wl, Tmax, Tmax);
    arb_mat_t Wf;
    arb_mat_init(Wf, Tmax, 1);
    arb_mat_t g;
    arb_mat_init(g, Tmax, 1);
    arb_mat_t RT;
    arb_mat_t RTWR;
    if (info.normalize_kernel) {
        arb_mat_init(RT, 1, Tmax);
        arb_mat_init(RTWR, 1, 1);
        for (int t = 0;t < Tmax;t++) {
            arb_set(arb_mat_entry(RT, 0, t), arb_mat_entry(R, t, 0));
        }
    }
    for (int il = 0;il < fit_info.lambdas.size(); il++) {
        // cov/C0^2
        arb_set_d(lam, r[0][Njack - 1]);
        arb_mul(lam, lam, lam, prec);
        arb_mat_scalar_div_arb(Wl, W, lam, prec);
        // lam*cov/C0^2
        arb_set_d(lam, fit_info.lambdas[il]);
        arb_mat_scalar_mul_arb(Wl, Wl, lam, prec);
        //W=(1-lam)A+ lam *Cov/c0^2
        arb_sub_ui(lam, lam, 1, prec);
        arb_neg(lam, lam);
        arb_mat_scalar_addmul_arb(Wl, A, lam, prec);
        // g=W^-1 f
        arb_mat_inv(Wl, Wl, prec);
        arb_mat_mul(Wf, Wl, f, prec);
        arb_mat_scalar_mul_arb(Wf, Wf, lam, prec);// *(1-lambda) is not inside f
        arb_mat_set(g, Wf);
        if (info.normalize_kernel) {
            // lam=1- R^T W^-1 f
            arb_mat_mul(RTWR, RT, Wf, prec);
            arb_sub_ui(lam, arb_mat_entry(RTWR, 0, 0), 1, prec);
            arb_neg(lam, lam);
            // RTWR= R^T W^-1 R
            arb_mat_mul(Wf, Wl, R, prec);
            arb_mat_mul(RTWR, RT, Wf, prec);
            // g+=  W^-1 R (1- R^T W^-1 f ) /(R^T W^-1 R)
            arb_div(lam, lam, arb_mat_entry(RTWR, 0, 0), prec);
            arb_mat_scalar_mul_arb(Wf, Wf, lam, prec);
            arb_mat_add(g, g, Wf, prec);
        }
        check_reconstruction(Delta, g, { info.E0, fit_info.maxE_check_reconstuct ,fit_info.stepsE_check_reconstuct });
    }
    arb_clear(lam);
    arb_mat_clear(Wl);
    arb_mat_clear(Wf);
    arb_mat_clear(g);
    if (info.normalize_kernel) {
        arb_mat_clear(RT);
        arb_mat_clear(RTWR);
    }
    // mpf_t **W=(mpf_t**) malloc(sizeof(mpf_t*)* Tmax);
    double** p;
    return p;
};

