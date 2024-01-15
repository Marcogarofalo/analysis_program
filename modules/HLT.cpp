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
#include "correlators_analysis.hpp"


// HLT_type_d::HLT_type_d(int tmax, int L0, double E0_, int njack, HLT_b type_b, double alpha) {
//     Tmax = tmax;
//     T = L0;
//     Njack = njack;
//     E0 = E0_;
//     R = (double*)malloc(sizeof(double) * Tmax);
//     A = malloc_2<double>(Tmax, Tmax);
//     type = type_b;
//     if (type == HLT_EXP_b) {
//         for (size_t t = 0; t < Tmax; t++) {
//             R[t] = 1.0 / (t + 1.0) + 1.0 / (T - t - 1.0);
//             for (size_t r = 0; r < Tmax; r++) {
//                 A[t][r] = exp(-(r + t + 2 - alpha) * E0) / (r + t + 2 - alpha);
//                 A[t][r] += exp(-(T - r + t - alpha) * E0) / (T - r + t - alpha);
//                 A[t][r] += exp(-(T + r - t - alpha) * E0) / (T + r - t - alpha);
//                 A[t][r] += exp(-(2 * T - r - t - 2 - alpha) * E0) / (2 * T - r - t - 2 - alpha);
//             }
//         }
//     }
// }

// HLT_type_d::~HLT_type_d() {
//     free(R);
//     for (size_t i = 0; i < Tmax; i++) {
//         free(A[i]);
//     }
//     free(A);
//     if (f_allocated) free(f);
// }



// double integrand_f(double x, void* params) {
//     wrapper_smearing_d* p = (wrapper_smearing_d*)params;
//     int t = p->t;
//     int T = p->HLT->T;
//     return p->function(x, p->params) * (exp(-(t + 1) * x) + exp(-(T - (t + 1)) * x));
// }


// void HLT_type_d::compute_f_EXP_b(wrapper_smearing_d& Delta, double epsrel) {

//     if (f_allocated == true) {
//         printf("HLT: recomputing f\n");
//         free(f);
//     }
//     f = (double*)malloc(sizeof(double) * Tmax);
//     int Maxiter = 1e+6;
//     gsl_integration_workspace* w = gsl_integration_workspace_alloc(Maxiter);
//     double result, error;
//     gsl_function F;
//     Delta.HLT = this;
//     for (int t = 0;t < Tmax;t++) {
//         F.function = &integrand_f;
//         Delta.t = t;
//         F.params = &Delta;

//         // gsl_integration_qags(&F, 0, 1, 0, epsrel, Maxiter, w, &result, &error);
//         // int_E0^\infty  b* Delta
//         gsl_integration_qagiu(&F, E0, 0, epsrel, Maxiter, w, &result, &error);
//         f[t] = result;
//     }
//     gsl_integration_workspace_free(w);
//     f_allocated = true;
// }

// //////////////////////////////////////////////////////////////////
// // HLT  smearing functions
// //////////////////////////////////////////////////////////////////


// double gaussian_for_HLT_d(double x, double* p) {
//     return (exp(-(x - p[0]) * (x - p[0]) / (2 * p[1] * p[1]))) / (p[1] * sqrt(2 * M_PIl));
// };
// double theta_for_HLT_d(double x, double* p) {
//     return 1.0 / (1.0 + exp(-(p[0] - x) / p[1]));
// };

// //////////////////////////////////////////////////////////////////


// double** HLT_type_d::HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
//     FILE* outfile, const char* description, wrapper_smearing_d  Delta, FILE* file_jack) {

//     double** r = malloc_2<double>(Tmax, Njack);
//     for (int t = 0;t < Tmax;t++)
//         for (int j = 0;j < Njack;j++)
//             r[t][j] = conf_jack[j][id][t][0];


//     double** cov = myres->comp_cov(Tmax, r);
//     if (type == HLT_EXP_b) {
//         compute_f_EXP_b(Delta);
//     }
//     else { printf("error HLT: type not supported\n");exit(1); }

//     HLT_out res(Delta.lambdas);

//     for (int il = 0;il < Delta.lambdas.size(); il++) {


//     }
//     // mpf_t **W=(mpf_t**) malloc(sizeof(mpf_t*)* Tmax);
//     double** p;
//     return p;
// };

////////////////////////////////////////////////////////////////////////////////////////////////
// using ARB
////////////////////////////////////////////////////////////////////////////////////////////////

// exp(-(r+t+alpha)*E0)/(r+t)
void compute_A_HLT_EXP_b(arb_t A, slong t, slong r, arb_t E0, arb_t alpha, int prec) {
    arb_t at, aT;
    arb_init(at);
    arb_init(aT);

    arb_set_ui(at, r + t);
    arb_add(at, at, alpha, prec);
    arb_mul(aT, at, E0, prec);
    arb_neg(aT, aT);
    arb_exp(aT, aT, prec);
    arb_div(A, aT, at, prec);

    arb_clear(at);arb_clear(aT);
}

// exp(-(r+t+alpha)*E0)/((r+t)*r*t)
void compute_A_HLT_EXP_boverT(arb_t A, slong t, slong r, arb_t E0, arb_t alpha, int prec) {
    arb_t at, aT;
    arb_init(at);
    arb_init(aT);

    arb_set_ui(at, r + t);
    arb_add(at, at, alpha, prec);
    arb_mul(aT, at, E0, prec);
    arb_neg(aT, aT);
    arb_exp(aT, aT, prec);
    arb_div(A, aT, at, prec);

    arb_div_ui(A, A, r, prec);
    arb_div_ui(A, A, t, prec);
    arb_clear(at);arb_clear(aT);
}
// exp(-(r+t+alpha)*E0)/((r+t)*r^2*t^2)
void compute_A_HLT_EXP_boverT2(arb_t A, slong t, slong r, arb_t E0, arb_t alpha, int prec) {
    arb_t at, aT;
    arb_init(at);
    arb_init(aT);

    arb_set_ui(at, r + t);
    arb_add(at, at, alpha, prec);
    arb_mul(aT, at, E0, prec);
    arb_neg(aT, aT);
    arb_exp(aT, aT, prec);
    arb_div(A, aT, at, prec);

    arb_div_ui(A, A, r*r, prec);
    arb_div_ui(A, A, t*t, prec);
    arb_clear(at);arb_clear(aT);
}

void compute_A_HLT_EXP_bT(arb_t A, slong t, slong r, arb_t E0, arb_t alpha, int T, int prec) {
    arb_t at, aT;
    arb_init(at);
    arb_init(aT);

    // A[t][r] = exp(-(r + t + 2 + alpha) * E0) / (r + t + 2 + alpha);
    arb_set_ui(at, r + t);
    arb_add(at, at, alpha, prec);
    // arb_sub(at, at, arb_alpha, prec);
    arb_mul(aT, at, E0, prec);
    arb_neg(aT, aT);
    arb_exp(aT, aT, prec);
    arb_div(A, aT, at, prec);
    //     A[t][r] += exp(-(T - r + t + alpha) * E0) / (T - r + t + alpha);
    arb_set_ui(at, T - r + t);
    arb_add(at, at, alpha, prec);
    // arb_sub(at, at, arb_alpha, prec);
    arb_mul(aT, at, E0, prec);
    arb_neg(aT, aT);
    arb_exp(aT, aT, prec);
    arb_div(aT, aT, at, prec);
    arb_add(A, A, aT, prec);
    //     A[t][r] += exp(-(T + r - t + alpha) * E0) / (T + r - t + alpha);
    arb_set_ui(at, T + r - t);
    arb_add(at, at, alpha, prec);
    // arb_sub(at, at, arb_alpha, prec);
    arb_mul(aT, at, E0, prec);
    arb_neg(aT, aT);
    arb_exp(aT, aT, prec);
    arb_div(aT, aT, at, prec);
    arb_add(A, A, aT, prec);
    //     A[t][r] += exp(-(2 * T - r - t - 2 + alpha) * E0) / (2 * T - r - t - 2 + alpha);
    arb_set_ui(at, 2 * T - r - t);
    arb_add(at, at, alpha, prec);
    // arb_sub(at, at, arb_alpha, prec);
    arb_mul(aT, at, E0, prec);
    arb_neg(aT, aT);
    arb_exp(aT, aT, prec);
    arb_div(aT, aT, at, prec);
    arb_add(A, A, aT, prec);

    arb_clear(at);arb_clear(aT);
}

void compute_A_HLT_deriv_EXP_b(arb_t A, slong t, slong r, arb_t E0, arb_t alpha, int tmin, int prec) {
    arb_t at, aT;
    arb_init(at);
    arb_init(aT);
    int n = t + r - 2 * tmin;
    int sign = n % 2 == 0 ? 1 : -1; // (-1)^(n-k)
    int binomial = 1;
    // doing the binomial expansion
    for (int k = 0;k <= n;k++) {

        arb_set_ui(at, r + t - k);
        arb_add(at, at, alpha, prec);
        arb_mul(aT, at, E0, prec);
        arb_neg(aT, aT);
        arb_exp(aT, aT, prec);
        arb_div(aT, aT, at, prec);

        arb_mul_si(aT, aT, sign, prec);
        arb_mul_ui(aT, aT, binomial, prec);


        arb_add(A, A, aT, prec);
        // update sign and binomial (which is always an integer)
        sign *= -1;
        binomial *= n - k;
        binomial /= k + 1;

    }

    arb_clear(at);arb_clear(aT);
}

HLT_type::HLT_type(HLT_type_input info_) {
    info = info_;
    error(info.tmax < 0, 1, "HLT_type", "tmax not set");
    error(info.T < 0, 1, "HLT_type", "T not set");
    error(info.prec < 0, 1, "HLT_type", "prec (precision) not set");
    error(info.E0 < 0, 1, "HLT_type", "E0 not set");
    error(info.type_b == HLT_INVALID_b, 1, "HLT_type", "type_b not set");

    arb_init(E0_arb);
    arb_set_d(E0_arb, info.E0);

    // init function to compute A
    if (info.A_strategy == HLT_A_fast)
        compute_A = std::bind(&HLT_type::compute_A_fast, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    else if (info.A_strategy == HLT_A_INTEGRAL) {
        // compute_A = std::bind(&HLT_type::compute_A_integral, this, std::placeholders::_1, std::placeholders::_2);
        printf("HLT_type: HLT_A_INTEGRAL not enymore supported\n");    exit(1);
    }
    else {
        printf("HLT_type: A_strategy not set\n A_strategy=%d\n", info.A_strategy);    exit(1);
    }
    // compute f strategy
    if (info.type_b == HLT_EXP_b)
        compute_b_re = std::bind(&HLT_type::compute_bt_re, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    else if (info.type_b == HLT_EXP_bT)
        compute_b_re = std::bind(&HLT_type::compute_btT_re, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    else if (info.type_b == HLT_deriv_EXP_b)
        compute_b_re = std::bind(&HLT_type::compute_deriv_bt_re, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    else if (info.type_b == HLT_EXP_bovert)
        compute_b_re = std::bind(&HLT_type::compute_btovert_re, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    else if (info.type_b == HLT_EXP_bovert2)
        compute_b_re = std::bind(&HLT_type::compute_btovert2_re, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    else {
        printf("HLT_type: type_b not set\n type_b=%d\n", info.type_b);    exit(1);
    }
    //


    if (info.normalize_kernel)  arb_mat_init(R, info.tmax - info.tmin, 1);
    arb_mat_init(A, info.tmax - info.tmin, info.tmax - info.tmin);
    arb_mat_init(A_ref, info.tmax - info.tmin, info.tmax - info.tmin);


    arb_t at, aT, zero_arb;
    arb_init(at);
    arb_init(aT);
    arb_init(zero_arb);
    arb_init(alpha_arb);
    arb_set_d(alpha_arb, info.alpha);
    arb_set_ui(zero_arb,0);

    double alpha = info.alpha;

    int prec = info.prec;

    // arb_t sinhE0;
    // arb_init(sinhE0);
    // arb_div_ui(sinhE0, E0_arb, 2, prec);
    // arb_sinh(sinhE0, sinhE0, prec);
    // arb_mul_ui(sinhE0, sinhE0, 2, prec);
    // arb_t binomial;
    // arb_init(binomial);

    for (size_t t = info.tmin; t < info.tmax; t++) {
        if (info.normalize_kernel) {
            if (info.type_b == HLT_EXP_b) {
                arb_set_ui(at, t);
                arb_inv(arb_mat_entry(R, t - info.tmin, 0), at, prec);
            }
            else if (info.type_b == HLT_EXP_bT) {
                arb_set_ui(at, t);
                arb_inv(arb_mat_entry(R, t - info.tmin, 0), at, prec);
                arb_set_ui(aT, info.T - t);
                arb_inv(aT, aT, prec);
                arb_add(arb_mat_entry(R, t - info.tmin, 0), at, aT, prec);
            }
            else { printf("HLT b_T type not supported\n"); }

        }
        for (size_t r = info.tmin; r < info.tmax; r++) {
            arb_ptr Atr = arb_mat_entry(A, t - info.tmin, r - info.tmin);
            arb_ptr Atr_ref = arb_mat_entry(A_ref, t - info.tmin, r - info.tmin);
            if (info.type_b == HLT_EXP_b) {

                compute_A_HLT_EXP_b(Atr, t, r, E0_arb, alpha_arb, prec);
                compute_A_HLT_EXP_b(Atr_ref, t, r, E0_arb, zero_arb, prec);
            }
            else if (info.type_b == HLT_EXP_bT) {
                int& T = info.T;

                compute_A_HLT_EXP_bT(Atr, t, r, E0_arb, alpha_arb, T, prec);
                compute_A_HLT_EXP_bT(Atr_ref, t, r, E0_arb, zero_arb, T, prec);
            }
            else if (info.type_b == HLT_deriv_EXP_b) {

                compute_A_HLT_deriv_EXP_b(Atr, t, r, E0_arb, alpha_arb, info.tmin, prec);
                compute_A_HLT_deriv_EXP_b(Atr_ref, t, r, E0_arb, zero_arb, info.tmin, prec);

            }
            else if (info.type_b == HLT_EXP_bovert) {
                compute_A_HLT_EXP_boverT(Atr, t, r, E0_arb, alpha_arb, prec);
                compute_A_HLT_EXP_boverT(Atr_ref, t, r, E0_arb, zero_arb, prec);
            }
            else if (info.type_b == HLT_EXP_bovert2) {
                compute_A_HLT_EXP_boverT2(Atr, t, r, E0_arb, alpha_arb, prec);
                compute_A_HLT_EXP_boverT2(Atr_ref, t, r, E0_arb, zero_arb, prec);
            }
            else { printf("HLT b_T type not supported\n"); }
        }
    }

    if (info.normalize_kernel) {
        arb_mat_init(RT, 1, info.tmax - info.tmin);
        arb_mat_init(RTWR, 1, 1);
        for (int t = 0;t < info.tmax - info.tmin;t++) {
            arb_set(arb_mat_entry(RT, 0, t), arb_mat_entry(R, t, 0));
        }
    }

    arb_mat_init(W, info.tmax - info.tmin, info.tmax - info.tmin);
    arb_init(lam);
    arb_mat_init(Wl, info.tmax - info.tmin, info.tmax - info.tmin);
    arb_mat_init(Wf, info.tmax - info.tmin, 1);

    arb_clear(at);
    arb_clear(aT);
    arb_clear(zero_arb);
}

HLT_type::~HLT_type() {
    arb_mat_clear(A);
    arb_clear(lam);
    arb_mat_clear(Wl);
    arb_mat_clear(Wf);

    arb_mat_clear(W);

    if (f_allocated) {
        arb_mat_clear(f);
        arb_clear(K2);
        arb_clear(K2_ref);
    }

    if (info.normalize_kernel) {
        arb_mat_clear(R);
        arb_mat_clear(RT);
        arb_mat_clear(RTWR);
    }

    if (f_allocated == true) {
        arb_mat_clear(f);
        arb_mat_clear(f_ref);
        arb_clear(K2);
        arb_clear(K2_ref);
    }
}


void HLT_type::compute_b(acb_t b, int  t, const acb_t E0) {
    int& prec = info.prec;
    if (info.type_b == HLT_EXP_b) {
        acb_set_si(b, -t);
        acb_mul(b, b, E0, prec);
        acb_exp(b, b, prec);
    }
    if (info.type_b == HLT_EXP_bT) {
        acb_set_ui(b, t);
        acb_mul(b, b, E0, prec);
        acb_neg(b, b);
        acb_exp(b, b, prec);

        acb_t bT; acb_init(bT);
        acb_set_ui(bT, info.T - t);
        acb_mul(bT, bT, E0, prec);
        acb_neg(bT, bT);
        acb_exp(bT, bT, prec);
        acb_add(b, b, bT, prec);
        acb_clear(bT);
    }
}

// void HLT_type::compute_b_re(arb_t b, int  t, const arb_t E0) {
//     int& prec = info.prec;
//     if (info.type_b == HLT_EXP_b) {
//         arb_set_si(b, -t);
//         arb_mul(b, b, E0, prec);
//         // arb_neg(b, b);
//         arb_exp(b, b, prec);
//     }
//     if (info.type_b == HLT_EXP_bT) {
//         arb_set_ui(b, t);
//         arb_mul(b, b, E0, prec);
//         arb_neg(b, b);
//         arb_exp(b, b, prec);

//         arb_t bT; arb_init(bT);
//         arb_set_ui(bT, info.T - t);
//         arb_mul(bT, bT, E0, prec);
//         arb_neg(bT, bT);
//         arb_exp(bT, bT, prec);
//         arb_add(b, b, bT, prec);
//         arb_clear(bT);
//     }

// }


void HLT_type::compute_deriv_bt_re(arb_t b, arb_t  apt, int t, const arb_t E0) {
    int& prec = info.prec;
    // (2(alpha+t)-t+tmin)/2
    arb_mul_ui(b, apt, 2, prec);
    arb_add_si(b, b, info.tmin - t, prec);
    arb_div_ui(b, b, 2, prec);
    arb_neg(b, b);
    arb_mul(b, b, E0, prec);
    arb_exp(b, b, prec);
    arb_t tmp;
    arb_init(tmp);
    arb_div_ui(tmp, E0, 2, prec);
    arb_sinh(tmp, tmp, prec);
    arb_mul_ui(tmp, tmp, 2, prec);
    for (int i = info.tmin;i < t;i++) {
        arb_mul(b, b, tmp, prec);
    }
    arb_clear(tmp);
    // int neg = (t - info.tmin) % 2;
}


void HLT_type::compute_bt_re(arb_t b, arb_t  apt, int t, const arb_t E0) {
    int& prec = info.prec;
    arb_neg(b, apt);
    arb_mul(b, b, E0, prec);
    arb_exp(b, b, prec);
}

void HLT_type::compute_btovert_re(arb_t b, arb_t  apt, int t, const arb_t E0) {
    int& prec = info.prec;
    arb_neg(b, apt);
    arb_mul(b, b, E0, prec);
    arb_exp(b, b, prec);
    arb_div_ui(b,b,t,prec);
}

void HLT_type::compute_btovert2_re(arb_t b, arb_t  apt, int t, const arb_t E0) {
    int& prec = info.prec;
    arb_neg(b, apt);
    arb_mul(b, b, E0, prec);
    arb_exp(b, b, prec);
    arb_div_ui(b,b,t*t,prec);
}


void HLT_type::compute_btT_re(arb_t b, arb_t  apt, int t, const arb_t E0) {
    int& prec = info.prec;

    arb_t bT; arb_init(bT);
    arb_sub_ui(bT, apt, info.T, prec);
    arb_neg(bT, bT);
    // arb_set_ui(bT, info.T - t);
    arb_mul(bT, bT, E0, prec);
    arb_neg(bT, bT);
    arb_exp(bT, bT, prec);

    arb_neg(b, apt);
    arb_mul(b, b, E0, prec);
    arb_neg(b, b);
    arb_exp(b, b, prec);

    arb_add(b, b, bT, prec);

    arb_clear(bT);

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
    // acb_neg(res, res);
    acb_div_arb(res, res, p[1], prec);
    // acb_neg(res, res);
    acb_exp(res, res, prec);
    acb_add_ui(res, res, 1, prec);
    acb_inv(res, res, prec);

    return 0;
}

int c_theta_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    // if (order > 1)
    //     flint_abort();  /* Would be needed for Taylor method. */

    //  return 1.0 / (1 + exp((x-p[0]) / p[1]));
    arb_t* p = (arb_t*)param;

    arb_sub(acb_realref(res), acb_realref(z), p[0], prec);
    arb_div(acb_realref(res), acb_realref(res), p[1], prec);
    arb_exp(acb_realref(res), acb_realref(res), prec);
    arb_add_ui(acb_realref(res), acb_realref(res), 1, prec);
    // arb_inv(acb_realref(res), acb_realref(res), prec);
    // arb_mul(acb_realref(res), acb_realref(res), p[2], prec);
    arb_div(acb_realref(res), p[2], acb_realref(res), prec);

    return 0;
}

int deriv_c_theta_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    arb_t tmp;    arb_init(tmp);
    arb_t* p = (arb_t*)param;
    arb_ptr r=acb_realref(res);
    arb_ptr  Es=p[0];
    arb_ptr  sigma=p[1];


    arb_div(tmp,acb_realref(z),sigma,prec);
    arb_exp(r, tmp, prec);

    arb_div(tmp,Es,sigma,prec);
    arb_exp(tmp, tmp, prec);

    arb_add(r,r,tmp,prec);
    arb_log(r,r,prec);
    arb_mul(r,r,sigma,prec);

    arb_sub(r,r,acb_realref(z),prec);
    arb_mul(r, p[2], r, prec);

    arb_clear(tmp);
    return 0;
}



int c1_theta_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    //  return 1.0 / (1 + exp(-(x-p[0]) / p[1]));
    arb_t* p = (arb_t*)param;

    arb_sub(acb_realref(res), acb_realref(z), p[0], prec);
    arb_div(acb_realref(res), acb_realref(res), p[1], prec);
    arb_exp(acb_realref(res), acb_realref(res), prec);
    arb_add_ui(acb_realref(res), acb_realref(res), 1, prec);
    // arb_inv(acb_realref(res), acb_realref(res), prec);
    // arb_mul(acb_realref(res), acb_realref(res), p[2], prec);
    arb_t tmp;
    arb_init(tmp);
    arb_set(tmp, p[3]);
    arb_mul(tmp, acb_realref(z), tmp, prec);
    arb_add(tmp, tmp, p[2], prec);
    arb_div(acb_realref(res), tmp, acb_realref(res), prec);

    arb_clear(tmp);
    return 0;
}

int c2_theta_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    arb_t* p = (arb_t*)param;
    arb_set_ui(acb_imagref(res), 0);

    arb_t tmp;    arb_init(tmp);
    // P[2]+P[3]*z
    arb_set(tmp, p[3]);
    arb_mul(tmp, acb_realref(z), tmp, prec);
    arb_add(acb_realref(res), tmp, p[2], prec);
    // P[4]*z^2
    arb_set(tmp, p[4]);
    arb_mul(tmp, acb_realref(z), tmp, prec);
    arb_mul(tmp, acb_realref(z), tmp, prec);
    arb_add(acb_realref(res), tmp, acb_realref(res), prec);

    //  (1 + exp((x-p[0]) / p[1]));
    arb_sub(tmp, acb_realref(z), p[0], prec);
    arb_div(tmp, tmp, p[1], prec);
    arb_exp(tmp, tmp, prec);
    arb_add_ui(tmp, tmp, 1, prec);

    arb_div(acb_realref(res), acb_realref(res), tmp, prec);

    arb_clear(tmp);
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

    arb_set_ui(acb_imagref(res), 0);
    p->function(res, z, p->params, order, prec);
    if (p->HLT->info.normalize_kernel) acb_mul_arb(res, res, p->Norm, prec);

    arb_add_ui(p->b_tmp, p->HLT->alpha_arb, p->t, prec);
    p->HLT->compute_b_re(p->b_tmp, p->b_tmp, p->t, acb_realref(z));


    arb_mul(acb_realref(res), acb_realref(res), p->b_tmp, prec);

    // * e^(-(alpha+t) E)
    // arb_add_ui(p->b_tmp, p->HLT->alpha_arb, p->t, prec);
    // arb_mul(p->b_tmp, p->b_tmp, acb_realref(z), prec);
    // arb_exp(p->b_tmp, p->b_tmp, prec);
    // arb_mul(acb_realref(res), acb_realref(res), p->b_tmp, prec);

    return 0;
}

int integrand_f0(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    wrapper_smearing* p = (wrapper_smearing*)param;

    arb_set_ui(acb_imagref(res), 0);
    p->function(res, z, p->params, order, prec);
    if (p->HLT->info.normalize_kernel) acb_mul_arb(res, res, p->Norm, prec);

    // arb_add_ui(p->b_tmp, p->HLT->alpha_arb, p->t, prec);
    arb_set_ui(p->b_tmp, p->t);
    p->HLT->compute_b_re(p->b_tmp, p->b_tmp, p->t, acb_realref(z));

    arb_mul(acb_realref(res), acb_realref(res), p->b_tmp, prec);

    return 0;
}


int integrand_f_pinf(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    // if (order > 1)
    //     flint_abort();  /* Would be needed for Taylor method. */
    // w=e^(E0-z)
    // \int_E0^\infty dz =\int_0^1 dw/w
    // z(w)=E0-log w

    wrapper_smearing* p = (wrapper_smearing*)param;
    arb_t* par = (arb_t*)p->params;

    printf("function  t=%d in  ", p->t); arb_printn(acb_realref(z), prec / 3.33, 0); flint_printf("\n");
    arb_set_ui(acb_imagref(res), 0);

    // p->function(res, zw, p->params, order, prec);
    // p->function_mapped_pinf(res, zw, p->params, order, prec);
    // e^(E0-E*)/sigma
    arb_sub(acb_realref(res), p->HLT->E0_arb, par[0], prec);
    arb_div(acb_realref(res), acb_realref(res), par[1], prec);
    arb_exp(acb_realref(res), acb_realref(res), prec);
    printf("function  e^(E0-E*)/sigma =  "); arb_printn(acb_realref(res), prec / 3.33, 0); flint_printf("\n");
    // w^(1/simga)
    arb_inv(p->b_tmp, par[1], prec);
    arb_pow(p->b_tmp, acb_realref(z), p->b_tmp, prec);
    // arb_log(p->b_tmp, acb_realref(z), prec);
    // arb_div(p->b_tmp,p->b_tmp,par[1], prec);
    // arb_exp(p->b_tmp, p->b_tmp, prec);
    printf("function  ( w^(1/simga) =  "); arb_printn(p->b_tmp, prec / 3.33, 0); flint_printf("\n");
    // res =p[2]/( w^(1/simga)+e^(E0-E*)/sigma)
    arb_add(acb_realref(res), acb_realref(res), p->b_tmp, prec);
    printf("function  ( w^(1/simga)+e^(E0-E*)/sigma) =  "); arb_printn(acb_realref(res), prec / 3.33, 0); flint_printf("\n");
    arb_div(acb_realref(res), par[2], acb_realref(res), prec);
    printf("function  p[2]/( w^(1/simga)+e^(E0-E*)/sigma) =  "); arb_printn(acb_realref(res), prec / 3.33, 0); flint_printf("\n");

    // w^(t-1+alpha-1/sigma)
    arb_inv(p->b_tmp, par[1], prec);
    arb_add_ui(p->b_tmp, p->b_tmp, p->t - 1, prec);
    arb_add(p->b_tmp, p->b_tmp, p->HLT->alpha_arb, prec);
    arb_pow(p->b_tmp, acb_realref(z), p->b_tmp, prec);
    arb_mul(acb_realref(res), acb_realref(res), p->b_tmp, prec);

    // * e^(-(alpha+t) E0)
    arb_set_d(p->b_tmp, -p->HLT->info.alpha - p->t);
    arb_mul(p->b_tmp, p->b_tmp, p->HLT->E0_arb, prec);
    arb_exp(p->b_tmp, p->b_tmp, prec);
    arb_mul(acb_realref(res), acb_realref(res), p->b_tmp, prec);
    printf("function out "); arb_printn(acb_realref(res), prec / 3.33, 0); flint_printf("\n");

    return 0;
}


// int integrand_f_pinf(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
//     // if (order > 1)
//     //     flint_abort();  /* Would be needed for Taylor method. */
//     // w=e^(E0-z)
//     // \int_E0^\infty dz =\int_0^1 dw/w
//     // z(w)=E0-log w

//     wrapper_smearing* p = (wrapper_smearing*)param;
//     arb_t* par = (arb_t*)p->params;

//     printf("function  t=%d in  ", p->t); arb_printn(acb_realref(z), prec / 3.33, 0); flint_printf("\n");
//     arb_set_ui(acb_imagref(res), 0);

//     // p->function(res, zw, p->params, order, prec);
//     // p->function_mapped_pinf(res, zw, p->params, order, prec);
//     // e^(E0-E*)/sigma
//     arb_sub(acb_realref(res), p->HLT->E0_arb, par[0], prec);
//     arb_div(acb_realref(res), acb_realref(res), par[1], prec);
//     arb_exp(acb_realref(res), acb_realref(res), prec);

//     //  z^{-(alpha+t-1+1/sigma)} e^{(E0-E*)/sigma}
//     arb_inv(p->b_tmp, par[1], prec);
//     arb_add_ui(p->b_tmp, p->b_tmp, p->t - 1, prec);
//     arb_add(p->b_tmp, p->b_tmp, p->HLT->alpha_arb, prec);
//     arb_neg(p->b_tmp, p->b_tmp);
//     arb_pow(p->b_tmp, acb_realref(z), p->b_tmp, prec);
//     arb_mul(acb_realref(res), acb_realref(res), p->b_tmp, prec);

//     // z^{-(alpha+t-1)} +z^{-(alpha+t-1+1/sigma)} e^{(E0-E*)/sigma}
//     arb_add_ui(p->b_tmp, p->HLT->alpha_arb, p->t - 1, prec);
//     arb_neg(p->b_tmp, p->b_tmp);
//     arb_pow(p->b_tmp, acb_realref(z), p->b_tmp, prec);
//     arb_add(acb_realref(res), acb_realref(res), p->b_tmp, prec);


//     arb_div(acb_realref(res), par[2], acb_realref(res), prec);

//     arb_add_ui(p->b_tmp, p->HLT->alpha_arb, p->t, prec);
//     arb_mul(p->b_tmp, p->b_tmp, p->HLT->E0_arb, prec);
//     arb_neg(p->b_tmp, p->b_tmp);
//     arb_exp(p->b_tmp, p->b_tmp, prec);
//     arb_mul(acb_realref(res), acb_realref(res), p->b_tmp, prec);

//     printf("function out "); arb_printn(acb_realref(res), prec / 3.33, 0); flint_printf("\n");

//     return 0;
// }


int integrand_K2(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    wrapper_smearing* p = (wrapper_smearing*)param;
    acb_set_ui(res, 1);
    p->function(res, z, p->params, order, prec);
    if (p->HLT->info.normalize_kernel) acb_mul_arb(res, res, p->Norm, prec);
    acb_mul(res, res, res, prec);


    // * e^(-alpha E)
    acb_t b; acb_init(b);
    acb_set_d(b, -p->HLT->info.alpha);
    acb_mul(b, b, z, prec);
    acb_exp(b, b, prec);
    acb_mul(res, res, b, prec);

    acb_clear(b);

    return 0;
}

int integrand_K20(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    wrapper_smearing* p = (wrapper_smearing*)param;
    acb_set_ui(res, 1);
    p->function(res, z, p->params, order, prec);
    if (p->HLT->info.normalize_kernel) acb_mul_arb(res, res, p->Norm, prec);
    acb_mul(res, res, res, prec);

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
    arb_init(b_tmp);
};

void wrapper_smearing::normilise_smearing() {
    int Maxiter = 1e+6;
    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

    options->deg_limit = HLT->info.integration_deg_limit;
    options->eval_limit = HLT->info.integration_eval_limit;
    options->depth_limit = HLT->info.integration_depth_limit;
    options->verbose = HLT->info.integration_verbose;

    mag_t tol;
    mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -HLT->info.prec);
    slong  goal = HLT->info.prec;
    acb_t a, b, s;
    acb_init(a);
    acb_init(b);
    acb_init(s);

    acb_set_ui(a, 0);
    acb_set_d(b, HLT->info.integration_maxE);
    int arb_calc_result = acb_calc_integrate(s, function, (void*)params, a, b, goal, tol, options, HLT->info.prec);
    error(arb_calc_result == ARB_CALC_NO_CONVERGENCE, 1, "normilise_smearing", "acb_calc_integrate returned ARB_CALC_NO_CONVERGENCE");

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
        arb_mat_clear(f_ref);
        arb_clear(K2);
        arb_clear(K2_ref);
    }
    arb_mat_init(f, info.tmax - info.tmin, 1);
    arb_mat_init(f_ref, info.tmax - info.tmin, 1);
    arb_init(K2);
    arb_init(K2_ref);
    int Maxiter = 1e+6;
    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

    options->deg_limit = info.integration_deg_limit;
    options->eval_limit = info.integration_eval_limit;
    options->depth_limit = info.integration_depth_limit;
    options->verbose = info.integration_verbose;
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
    acb_set_arb(a, E0_arb);
    acb_set_d(b, info.integration_maxE);
    // acb_set_ui(a, 0);
    // acb_set_ui(b, 1);
    void* param = (void*)&Delta;

    // std::vector<double>  v(info.tmax - info.tmin);
    for (int t = info.tmax - 1;t >= info.tmin;t--) {
        // double time = timestamp();
        Delta.t = t;
        // int arb_calc_result = acb_calc_integrate(s, integrand_f_pinf, param, a, b, goal, tol, options, info.prec);
        int arb_calc_result = acb_calc_integrate(s, integrand_f, param, a, b, goal, tol, options, info.prec);
        error(arb_calc_result == ARB_CALC_NO_CONVERGENCE, 1, "compute_f_EXP_b", "acb_calc_integrate f(t=%d)  returned ARB_CALC_NO_CONVERGENCE", t);
        acb_get_real(arb_mat_entry(f, t - info.tmin, 0), s);
        // f_ref
        arb_calc_result = acb_calc_integrate(s, integrand_f0, param, a, b, goal, tol, options, info.prec);
        error(arb_calc_result == ARB_CALC_NO_CONVERGENCE, 1, "compute_f0_EXP_b", "acb_calc_integrate f(t=%d)  returned ARB_CALC_NO_CONVERGENCE", t);
        acb_get_real(arb_mat_entry(f_ref, t - info.tmin, 0), s);

    }

    int arb_calc_result = acb_calc_integrate(s, integrand_K2, param, a, b, goal, tol, options, info.prec);
    error(arb_calc_result == ARB_CALC_NO_CONVERGENCE, 1, "compute_f_EXP_b", "acb_calc_integrate K2 returned ARB_CALC_NO_CONVERGENCE");
    acb_get_real(K2, s);
    // K2_ref
    arb_calc_result = acb_calc_integrate(s, integrand_K20, param, a, b, goal, tol, options, info.prec);
    error(arb_calc_result == ARB_CALC_NO_CONVERGENCE, 1, "compute_f_EXP_b", "acb_calc_integrate K2 returned ARB_CALC_NO_CONVERGENCE");
    acb_get_real(K2_ref, s);


    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
    f_allocated = true;
}



int integrand_A(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
    // if (order > 1)
    //     flint_abort();  /* Would be needed for Taylor method. */

    wrapper_smearing* p = (wrapper_smearing*)param;
    arb_set_ui(acb_imagref(res), 0);
    p->function(res, z, p->params, order, prec);
    if (p->HLT->info.normalize_kernel) acb_mul_arb(res, res, p->Norm, prec);

    for (int t = p->HLT->info.tmin;t < p->HLT->info.tmax;t++) {
        arb_set_ui(p->b_tmp, t);
        p->HLT->compute_b_re(p->b_tmp, p->b_tmp, t, acb_realref(z));
        arb_submul(acb_realref(res), arb_mat_entry(p->HLT->g, t - p->HLT->info.tmin, 0), p->b_tmp, prec);
    }
    arb_mul(acb_realref(res), acb_realref(res), acb_realref(res), prec);

    // * e^(-alpha E)
    arb_set_d(p->b_tmp, -p->HLT->info.alpha);
    arb_mul(p->b_tmp, p->b_tmp, acb_realref(z), prec);
    arb_exp(p->b_tmp, p->b_tmp, prec);
    arb_mul(acb_realref(res), acb_realref(res), p->b_tmp, prec);

    return 0;
}



void HLT_type::compute_A_integral(arb_t Ag, wrapper_smearing& Delta) {
    int Maxiter = 1e+6;
    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

    options->deg_limit = info.integration_deg_limit;
    options->eval_limit = info.integration_eval_limit;
    options->depth_limit = info.integration_depth_limit;
    options->verbose = info.integration_verbose;

    mag_t tol;
    mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -info.prec);
    slong  goal = info.prec;
    acb_t a, b, s;
    acb_init(a);
    acb_init(b);
    acb_init(s);

    acb_set_arb(a, E0_arb);
    acb_set_d(b, info.integration_maxE);

    int arb_calc_result = acb_calc_integrate(s, integrand_A, (void*)&Delta, a, b, goal, tol, options, info.prec);
    error(arb_calc_result == ARB_CALC_NO_CONVERGENCE, 1, "compute_A", "acb_calc_integrate returned ARB_CALC_NO_CONVERGENCE");
    acb_get_real(Ag, s);


    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
}


void HLT_type::compute_A_fast(arb_t Ag, arb_t Ag_ref, wrapper_smearing& Delta) {
    arb_mat_t gAg, gt, v;
    int& prec = info.prec;
    arb_mat_init(gAg, 1, 1);
    arb_ptr gAg00=arb_mat_entry(gAg, 0, 0);
    arb_mat_init(v, info.tmax - info.tmin, 1);
    arb_mat_init(gt, 1, info.tmax - info.tmin);
    arb_mat_transpose(gt, g);
    arb_set_ui(Ag, 0);

    arb_mat_mul(v, A, g, prec);
    arb_mat_mul(gAg, gt, v, prec);
    arb_add(Ag, Ag, gAg00, prec);
    // printf("gAg=");arb_printn(arb_mat_entry(gAg, 0, 0), prec / 3.33, 0); flint_printf("\n");

    arb_mat_mul(gAg, gt, f, prec);
    arb_mul_si(gAg00, gAg00, -2, prec);
    arb_add(Ag, Ag, gAg00, prec);
    // printf("gf=");arb_printn(arb_mat_entry(gAg, 0, 0), prec / 3.33, 0); flint_printf("\n");

    arb_add(Ag, Ag, K2, prec);
    ///////////////////////////
    // A_ref
    //////////////////////////
    arb_set_ui(Ag_ref, 0);

    arb_mat_mul(v, A_ref, g, prec);
    arb_mat_mul(gAg, gt, v, prec);
    arb_add(Ag_ref, Ag_ref, gAg00, prec);

    arb_mat_mul(gAg, gt, f_ref, prec);
    arb_mul_si(gAg00, gAg00, -2, prec);
    arb_add(Ag_ref, Ag_ref, gAg00, prec);

    arb_add(Ag_ref, Ag_ref, K2_ref, prec);

    // clear
    arb_mat_clear(gAg);
    arb_mat_clear(gt);
    arb_mat_clear(v);
}

void HLT_type::check_reconstruction(wrapper_smearing& Delta, const char* description,
    int il, fit_type_HLT fit_info, std::array<double, 3> range) {
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
    double dh = (range[1] - range[0]) / range[2];
    int& prec = info.prec;
    // printf("check HLT  g: \n");
    // for (int t = 0;t < info.tmax;t++) {
    //     printf("t= %d   g= ", t);
    //     arb_printn(arb_mat_entry(g, t, 0), prec / 3.33, 0); flint_printf("\n");
    // }
    FILE* outfile = fit_info.outfile_kernel;
    fprintf(outfile, " \n\n");
    fprintf(outfile, "# HLT of correlator id=%d\n", fit_info.corr_id[0]);
    for (int i = 0;i < range[2] + 1;i++) {
        arb_set_d(E_re, range[0] + i * dh);
        acb_set_arb(E, E_re);
        Delta.function(res, E, (void*)Delta.params, 0, prec);
        if (info.normalize_kernel) acb_mul_arb(res, res, Delta.Norm, prec);

        acb_get_real(res_re, res);
        double d = arbtod(res_re);
        arb_set_ui(res_HLT, 0);
        for (int t = info.tmin;t < info.tmax;t++) {
            arb_set_ui(b, t);
            compute_b_re(b, b, t, E_re);
            arb_addmul(res_HLT, arb_mat_entry(g, t - info.tmin, 0), b, prec);

        }


        double d1 = arbtod(res_HLT);
        // arb_sub(b, res_re, res_HLT, prec);
        fprintf(outfile, "%-20.12g  %-25.12g  0    %-25.12g  0  %-20.12g  %-25.12g  0   %-20.12g  %-25.12g  0  \n", range[0] + i * dh, d, d1,
            range[0] + i * dh, d1, range[0] + i * dh, d1);
    }
    char name[NAMESIZE];
    mysprintf(name, NAMESIZE, "%s_lam%.8g", description, lambdas[il]);
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", name, info.tmax, info.tmax, lambdas[il], 0.0);
    fprintf(outfile, "%.12g  %.12g \n", Ag[il] / A0, A0);
    acb_clear(res);
    arb_clear(res_re);
    acb_clear(E);
    arb_clear(E_re);
    arb_clear(res_HLT);
    arb_clear(b);

}

// void HLT_type::compute_B() {

// }

void HLT_type::compute_A_and_B(wrapper_smearing& Delta, int  il) {
    //////////////////////////////////// A
    arb_t A_arg, A_ref;
    arb_init(A_arg);arb_init(A_ref);
    compute_A(A_arg, A_ref, Delta);
    Ag[il] = arbtod(A_arg);
    Ag_ref[il] = arbtod(A_ref);
    //////////////////////////////////// B
    int& Tmax = info.tmax;
    int& prec = info.prec;
    arb_mat_t Wg, gBg, gt;
    arb_mat_init(Wg, Tmax - info.tmin, 1);
    arb_mat_init(gBg, 1, 1);
    arb_mat_init(gt, 1, Tmax - info.tmin);
    arb_mat_transpose(gt, g);
    arb_mat_mul(Wg, W, g, prec);
    arb_mat_mul(gBg, gt, Wg, prec);
    Bg[il] = arbtod(arb_mat_entry(gBg, 0, 0));
    //////////////////////////////////// clear
    arb_mat_clear(gt);
    arb_mat_clear(Wg);
    arb_mat_clear(gBg);
    arb_clear(A_arg);
    arb_clear(A_ref);
}

// r[t]=c(t+1)
void  HLT_type::compute_tilderho(double* rho, double** r, fit_type_HLT& fit_info) {
    const int Njack = fit_info.Njack;
    // myarb a(info.prec), ar(info.prec);
    arb_t a, ar;
    arb_init(a);
    arb_init(ar);
    for (int j = 0;j < Njack;j++) {
        arb_set_ui(a, 0);
        for (int t = 0;t < info.tmax - info.tmin;t++) {
            // r[t]= c(t+1)
            arb_set_d(ar, r[t][j]);
            arb_addmul(a, arb_mat_entry(g, t, 0), ar, info.prec);
        }
        // if (j == Njack - 1) { printf("a="); arb_printn(a, info.prec / 3.33, 0); flint_printf("\n"); }
        rho[j] = arbtod(a);

    }
    arb_clear(ar);
    arb_clear(a);


}

void HLT_type::compute_g(double lambda) {
    int& Tmax = info.tmax;
    int& prec = info.prec;
    // lam*cov/C0^2
    arb_set_d(lam, lambda);
    arb_mat_scalar_mul_arb(Wl, W, lam, prec);
    // //W=(1-lam)A+ lam *Cov/c0^2
    // // arb_sub_ui(lam, lam, 1, prec);
    // // arb_neg(lam, lam);
    // W=A/A0+ lam *Cov/c0^2
    arb_inv(lam, A0_arb, prec);
    arb_mat_scalar_addmul_arb(Wl, A, lam, prec);
    int isolve = arb_mat_solve(Wf, Wl, f, prec);
    error(isolve == 0, 1, "HLT_of_corr was not able to invert the matrix W\n TIP: increase precision or lambda_min",
        "prec=%d lambda=%g", prec, lambda);
    arb_mat_scalar_mul_arb(g, Wf, lam, prec);// *(1-lambda) is not inside f
    // arb_mat_set(g, Wf);
    if (info.normalize_kernel) {
        // lam=1- R^T W^-1 f
        arb_mat_mul(RTWR, RT, g, prec);
        arb_sub_ui(lam, arb_mat_entry(RTWR, 0, 0), 1, prec);
        arb_neg(lam, lam);
        // RTWR= R^T W^-1 R
        // arb_mat_mul(Wf, Wl, R, prec);
        arb_mat_solve(Wf, Wl, R, prec);
        arb_mat_mul(RTWR, RT, Wf, prec);
        // g+=  W^-1 R (1- R^T W^-1 f ) /(R^T W^-1 R)
        arb_div(lam, lam, arb_mat_entry(RTWR, 0, 0), prec);
        arb_mat_scalar_mul_arb(Wf, Wf, lam, prec);
        arb_mat_add(g, g, Wf, prec);
    }


}


fit_result HLT_type::HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
    const char* description, wrapper_smearing& Delta, FILE* file_jack, fit_type_HLT fit_info) {

    error(fit_info.outfile_kernel == NULL, 1, "HLT_of_corr", "outfile kernel not set");
    error(fit_info.outfile_AoverB == NULL, 1, "HLT_of_corr", "outfile AoverB not set");
    error(fit_info.corr_id.size() != 1, 1, "HLT_of_corr", "corr_id must be of size 1, instead %d", fit_info.corr_id.size());
    error(info.tmax <= 0, 1, "HLT_of_corr", "tmax not set: tmax= %d", info.tmax);
    error(info.tmin < 0, 1, "HLT_of_corr", "tmin not set: tmin= %d", info.tmin);
    error(info.tmin >= info.tmax, 1, "HLT_of_corr", "tminmust be smaller than tmax: tmin= %d  tmax=%d", info.tmin, info.tmax);
    int& Tmax = info.tmax;
    int& Tmin = info.tmin;
    int& prec = info.prec;
    int Njack = fit_info.Njack;
    int id = fit_info.corr_id[0];
    Ag.resize(fit_info.nlambda_max);
    Ag_ref.resize(fit_info.nlambda_max);
    Bg.resize(fit_info.nlambda_max);
    lambdas.resize(fit_info.nlambda_max);
    double** r = malloc_2<double>(Tmax - Tmin, Njack);
    double** rho = malloc_2<double>(fit_info.nlambda_max, Njack);
    for (int t = Tmin;t < Tmax;t++)
        for (int j = 0;j < Njack;j++)
            r[t - Tmin][j] = conf_jack[j][id][t][0];

    // for (int d = 0;d < Tmax - Tmin;d++) {
    //     printf(" r= %.12g\n", r[d][Njack - 1]);
    // }
    if (info.type_b == HLT_deriv_EXP_b) {
        // replace c(1), c(2), c(3),...
        // with    c(1), -c^1(1+1/2),  c^2(2), -c^3(2+1/2)
        // minus sign compensate the derivative of the exponential
        for (int d = 1;d < Tmax - Tmin;d++) {
            for (int t = Tmax - Tmin - 1;t >= d;t--) {
                for (int j = 0;j < Njack;j++)
                    r[t][j] = r[t - 1][j] - r[t][j];
            }
        }
    }
    // for (int d = 0;d < Tmax - Tmin;d++) {
    //     printf(" r after= %.12g\n", r[d][Njack - 1]);
    // }

    myres->comp_cov_arb(W, Tmax - Tmin, r, prec);
    //// test covariance
    // double** cov = myres->comp_cov(Tmax - Tmin, r);
    // for (int t = 0;t < Tmax - Tmin;t++) {
    //     for (int r = 0;r < Tmax - Tmin;r++) {
    //         double tmp = arbtod(arb_mat_entry(W, t, r));
    //         double diff = fabs(tmp - cov[t][r]);
    //         printf("%-12.5g",tmp);
    //         // if (diff > 1e-6) {
    //         //     printf("%-4d  %-4d  %-20.12g %-20.12g diff= %-20.12g  ", t, r, tmp, cov[t][r], diff);
    //         //     printf(" error in the covariance\n");
    //         //     exit(1);
    //         // }
    //     }
    //     printf("\n");
    // }
    if (fit_info.diag_cov == true) {
        // printf("diagonal covariance\n");
        for (int t = 0;t < Tmax - Tmin;t++) {
            for (int r = 0;r < Tmax - Tmin;r++) {
                if (t == r) continue;
                else arb_set_ui(arb_mat_entry(W, t, r), 0);
            }
        }
    }
    // read cov
    // FILE* file = open_file("cov_NT.txt", "r");
    // for (int t = 0;t < Tmax - Tmin;t++) {
    //     for (int r = 0;r < Tmax - Tmin;r++) {
    //         double d=arbtod(arb_mat_entry(W, t, r));
    //         printf("%d  %d   %.12g\n",t,r,d);
    //         // int tt,rr;
    //         // char nn[NAMESIZE];
    //         // double cov_read;
    //         // fscanf(file,"%s  %d  %d   %lf\n",nn, &tt, &rr, &cov_read);
    //         // arb_set_d(arb_mat_entry(W, t, r), cov_read);
    //     }
    // }

    if (f_allocated == false) {
        printf("HLT: recomputing f\n");
        compute_f_EXP_b(Delta);
    }

    HLT_out res(fit_info.lambdas);

    arb_mat_init(g, Tmax - Tmin, 1);
    arb_init(A0_arb); arb_init(A0_arb_ref);
    compute_A(A0_arb, A0_arb_ref, Delta);
    A0 = arbtod(A0_arb);
    A0_ref = arbtod(A0_arb_ref);
    printf("A0=%.12g\n", A0);


    double Bnorm = conf_jack[Njack - 1][id][Tmin][0] * conf_jack[Njack - 1][id][Tmin][0];

    // cov/C0^2
    // arb_set_d(lam, conf_jack[Njack - 1][id][Tmin][0]);
    // arb_mul(lam, lam, lam, prec);
    arb_set_d(lam, Bnorm);
    arb_mat_scalar_div_arb(W, W, lam, prec);

    lambdas[0] = fit_info.lambda_start;
    int  same = 0;
    double* diff = (double*)calloc(Njack, sizeof(double));

    // first computation
    compute_g(lambdas[same]);
    compute_A_and_B(Delta, same);
    check_reconstruction(Delta, description, same,
        fit_info, { info.E0, fit_info.maxE_check_reconstuct ,fit_info.stepsE_check_reconstuct });
    compute_tilderho(rho[same], r, fit_info);

    printf("%-20.12g %-20.12g  %-20.12g   %-20.12g  %-20.12g    %d\n", lambdas[same], rho[same][Njack - 1],
        myres->comp_error(rho[same]), diff[Njack - 1], myres->comp_error(diff), same);

    // update

    int same_end = 0;
    int same_max = 1;
    int same_start = 0;
    same++;
    int iter = 1;
    while (iter < fit_info.nlambda_max) {
        lambdas[iter] = lambdas[iter - 1] * fit_info.reduce_lambda;
        compute_g(lambdas[iter]);
        compute_A_and_B(Delta, iter);
        check_reconstruction(Delta, description, iter,
            fit_info, { info.E0, fit_info.maxE_check_reconstuct ,fit_info.stepsE_check_reconstuct });
        compute_tilderho(rho[iter], r, fit_info);

        for (int j = 0;j < Njack;j++) {
            diff[j] = fabs(rho[iter - same][j] - rho[iter][j]);
        }
        // fprintf(fit_info.outfile_AoverB, "%-20.12g %-20.12g  %-20.12g  %-20.12g  %-20.12g  %-20.12g  %-20.12g  %s\n",
        //     lambdas[iter], Ag[iter], Bg[iter], A0, Ag[iter] / A0 + lambdas[iter] * Bg[iter],
        //     rho[iter][Njack - 1], myres->comp_error(rho[iter]), description);
        printf("%-20.12g %-20.12g  %-20.12g   %-20.12g  %-20.12g    %d %d\n", lambdas[iter], rho[iter][Njack - 1],
            myres->comp_error(rho[iter]), diff[Njack - 1], myres->comp_error(diff), iter, same);

        if (diff[Njack - 1] < myres->comp_error(rho[iter])) {// if we found that the result is compatible with the previuos
            same++;
        }
        else if (same > 1) {//we need to try shifting the plateaux of 1 reset same
            for (int ib = 1; ib < same;ib++) {
                for (int j = 0;j < Njack;j++) {
                    diff[j] = fabs(rho[iter - same + ib][j] - rho[iter][j]);
                }
                if (diff[Njack - 1] < myres->comp_error(rho[iter])) {
                    // printf("shift the pateaux forward\n");
                    break;
                }
                else { same--; }
            }

        }
        if (same >= same_max) {
            same_max = same;
            same_start = iter - same + 1;
        }
        iter++;
    }
    same_end = same_start - 1 + same_max;
    // same_start++;// there is an offset
    printf("same [%d,%d]  n= %d\n", same_start, same_end, same_max);

    // struct fit_type fit_info;
    struct fit_type const_fit_info;
    const_fit_info.Nvar = 1;
    const_fit_info.Npar = 1;
    const_fit_info.N = 1;
    const_fit_info.Njack = Njack;
    const_fit_info.function = constant_fit;
    const_fit_info.n_ext_P = 0;
    const_fit_info.linear_fit = true;
    const_fit_info.T = fit_info.nlambda_max;
    struct kinematic kinematic_2pt;
    const_fit_info.codeplateaux = true;
    const_fit_info.tmin = same_start;
    const_fit_info.tmax = same_end;
    double** mt = malloc_2<double>(fit_info.nlambda_max, 2);
    for (int i = 0;i < fit_info.nlambda_max;i++) {
        mt[i][0] = rho[i][Njack - 1];
        mt[i][1] = myres->comp_error(rho[i]);
    }
    int store_l0 = file_head.l0;
    file_head.l0 = fit_info.nlambda_max * 2;

    char name[NAMESIZE] = "HLT";
    double* chi2;
    struct fit_result fit_out = try_fit(option, const_fit_info.tmin, const_fit_info.tmax, 1, mt, rho, Njack, &chi2, const_fit_info);
    file_head.l0 = store_l0;

    //////////////////////////////////// print Aover B
    fprintf(fit_info.outfile_AoverB, "\n\n%-20s %-20s  %-20s  %-20s   %-20s  %-20s   %-20s  %-20s  %-20s  %-20s  %-20s  %-20s\n",
        "#lambda", "A", "B/Bnorm", "A0", "W", "rho", "drho", "label", "fit", "dfit", "Aref", "A0ref");
    for (int i = 0;i < fit_info.nlambda_max;i++) {
        fprintf(fit_info.outfile_AoverB, "%-20.12g %-20.12g  %-20.12g  %-20.12g  %-20.12g  %-20.12g  %-20.12g   %s  %-20.12g  %-20.12g  %-20.12g  %-20.12g\n",
            lambdas[i], Ag[i], Bg[i], A0, Ag[i] / A0 + lambdas[i] * Bg[i],
            rho[i][Njack - 1], myres->comp_error(rho[i]), description, fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]),
            Ag_ref[i], A0_ref);
    }
    fprintf(fit_info.outfile_AoverB, "\n\n #%s fit in [%.12g,%.12g] chi2=%.5g  %.5g\n", description, lambdas[const_fit_info.tmax], lambdas[const_fit_info.tmin], fit_out.chi2[Njack - 1], 0.0);
    fprintf(fit_info.outfile_AoverB, "%.12g  %.12g \n", fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]));

    //////////////////////////////////// print
    //////
    fprintf(fit_info.outfile, "\n\n%-10s %-20s  %-20s  %-20s  %-20s\n", "#lambda", "rho", "drho", "fit", "dfit");
    for (int i = 0;i < fit_info.nlambda_max;i++) {
        fprintf(fit_info.outfile, "%-20.12g %-20.12g  %-20.12g   %-20.12g  %-20.12g\n", lambdas[i], rho[i][Njack - 1],
            mt[i][1], fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]));
    }


    fprintf(fit_info.outfile, "\n\n #%s fit in [%.12g,%.12g] chi2=%.5g  %.5g\n", description, lambdas[const_fit_info.tmax], lambdas[const_fit_info.tmin], fit_out.chi2[Njack - 1], 0.0);
    fprintf(fit_info.outfile, "%.12g  %.12g \n", fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]));

    printf("%s= %.12g  %.12g \n", description, fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]));
    // fprintf(fit_info.outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 0, info.tmax, 0.0, 0.0);
    // fprintf(fit_info.outfile, "0  0 \n");

    arb_mat_clear(g);
    error(same_max < fit_info.nsame && fit_info.strict_checks, 1, "HLT_of_corr", "we could not determine lambda\n try increasing fit_info.nlambda_max",
        "or decresing fit_info.nsame");

    // we do not write the result
    // fwrite(fit_out.P[0], sizeof(double), Njack, file_jack);
    // corr_counter++;
    arb_clear(A0_arb);
    arb_clear(A0_arb_ref);
    free_2(Tmax - Tmin, r);
    return fit_out;
};

