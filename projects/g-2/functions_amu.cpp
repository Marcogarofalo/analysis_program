#include "functions_amu.hpp"
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <algorithm>
#include <string>

#include "correlators_analysis.hpp"
#include "fit_all.hpp"
#include "global.hpp"
#include "linear_fit.hpp"
#include "mutils.hpp"
#include "non_linear_fit.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "tower.hpp"

extern "C" {
    // #include "../external/rzeta/src/dzeta_function.h"
    // #include "../external/rzeta/src/qzeta_function.h"
#include "dzeta_function.h"
#include "qzeta_function.h"

}

double integrand_K(double x, void* params) {
    double z = *(double*)params;
    double arg = (z / 2) * x / (sqrt(1 - x));
    double j0 = sin(arg) / arg;
    double f = (1 - x) * (1 - j0 * j0);
    return f;
}

double kernel_K(double z, double epsrel) {
    int Maxiter = 1e+6;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(Maxiter);
    double result, error;
    gsl_function F;
    F.function = &integrand_K;
    F.params = &z;

    gsl_integration_qags(&F, 0, 1, 0, epsrel, Maxiter, w, &result, &error);
    // printf("K(%g)=%g\n", z, result);
    gsl_integration_workspace_free(w);

    return result;
}

double gm2_step_function(double t_d, double t1_d) {
    return 1.0 / (1.0 + exp(-2 * (t_d - t1_d)));
}

double integrand_K_W(double x, void* params) {
    double z = *(double*)params;
    double arg = (z / 2) * x / (sqrt(1 - x));
    double j0 = sin(arg) / arg;
    double f = (1 - x) * (1 - j0 * j0);

    double  t = z / muon_mass_MeV;
    constexpr double d = 0.15;
    constexpr double t0_d = 0.4 / d;
    constexpr double t1_d = 1.0 / d;
    double theta = gm2_step_function((t * 197.326963) / 0.15, t0_d) - gm2_step_function((t * 197.326963) / 0.15, t1_d);
    f *= theta;

    return f;
}

double kernel_K_W(double z, double epsrel) {
    int Maxiter = 1e+8;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(Maxiter);
    double result, error;
    gsl_function F;
    F.function = &integrand_K_W;
    F.params = &z;

    gsl_integration_qags(&F, 0, 1, 0, epsrel, Maxiter, w, &result, &error);
    // printf("K(%g)=%g\n", z, result);
    gsl_integration_workspace_free(w);

    return result;
}


double integrate_simpson38(int lower, int upper, double* f) {
    double integration = 0;
    integration = f[lower] + f[upper];
    for (int i = 1; i <= (upper - lower); i++) {
        int k = lower + i;

        if (i % 3 == 0) {
            integration += 2 * f[k];
        }
        else {
            integration += 3 * f[k];
        }

    }

    integration *= 3.0 / 8.0;
    return integration;
}

double integrate_reinman(int lower, int upper, double* f) {
    double integration = 0;
    for (int i = 0; i <= (upper - lower); i++) {
        int k = lower + i;
        integration += f[k];
    }

    return integration;
}



// double compute_amu_sd(double* VV, double Z, double a, double q2, double (*int_scheme)(int, int, double*)) {
//     constexpr double d = 0.15;
//     constexpr double t1_d = 0.4 / d;
//     int T = file_head.l0;

//     double* fi = (double*)malloc(sizeof(double) * T / 2);
//     fi[0] = 0;
//     for (int t_a = 1; t_a < T / 2; t_a++) {
//         double t = t_a * a; // time in fm.
//         double z = muon_mass_MeV * (t / 197.326963);
//         // printf("time = %d, z=%g\n", t_a, z);
//         double K = z * z * kernel_K(z);
//         double theta = gm2_step_function(t / 0.15, t1_d);
//         double VV_sub = Z * Z * VV[t_a] - (1.0 / (2.0 * M_PI * M_PI * pow(t_a, 5)));
//         fi[t_a] = K * VV_sub * (1 - theta);
//         // printf("t=%d  K=%g   VV=%g  1-theta=%g  amu=%g\n",t_a,K, VV[t_a],1-theta, amu);
//     }

//     double amu = int_scheme(0, T / 2 - 1, fi);

//     free(fi);
//     amu *= 4 * alpha_em * alpha_em *
//         q2 / (muon_mass_MeV * muon_mass_MeV * (a / 197.326963) * (a / 197.326963));

//     return amu;
// }


double* compute_amu_full(double**** in, int id, int Njack, double* Z, double* a, double q2,
    double (*int_scheme)(int, int, double*), FILE* outfile, const char* description,
    const char* resampling, int isub) {
    constexpr double d = 0.15;
    constexpr double t1_d = 0.4 / d;
    int T = file_head.l0;
    double** fi = double_malloc_2(T / 2, Njack);
    double* amu = (double*)malloc(sizeof(double) * Njack);
    double** amut = malloc_2<double>(T / 2, Njack);
    double* ft = (double*)malloc(sizeof(double) * T / 2);

    double** Kt = double_malloc_2(T / 2, Njack);
    // double** thetat = double_malloc_2(T / 2, Njack);
    double** corr_sub = double_malloc_2(T / 2, Njack);


    for (int j = 0;j < Njack;j++) {

        ft[0] = 0;
        for (int t_a = 1; t_a < T / 2; t_a++) {
            double t = t_a * a[j]; // time in fm.
            double z = muon_mass_MeV * (t / 197.326963);
            double K = z * z * kernel_K(z);
            // double theta = gm2_step_function(t / 0.15, t1_d);
            double VV_sub;
            // if (isub==-3)
            //     VV_sub = Z[j] * Z[j] * in[j][id][t_a][0] + in[j][isub][t_a][0]*(1 - theta);// it does not work, we need an extra param isub need to be the id of the tree-correlator
            if (isub == -2)
                VV_sub = Z[j] * Z[j] * in[j][id][t_a][0];
            else if (isub == -1)
                VV_sub = Z[j] * Z[j] * in[j][id][t_a][0] - (1.0 / (2.0 * M_PI * M_PI * pow(t_a, 5)));// perturbative
            else
                VV_sub = Z[j] * Z[j] * in[j][id][t_a][0] + in[j][isub][t_a][0];// free-theory
            ft[t_a] = K * VV_sub;

            fi[t_a][j] = ft[t_a];
            Kt[t_a][j] = K;
            // thetat[t_a][j] = (1 - theta);
            corr_sub[t_a][j] = VV_sub;
            // printf("t=%d  K=%g   VV=%g  1-theta=%g  amu=%g\n",t_a,K, VV[t_a],1-theta, amu);
        }
        // }
        // int tmax=0;
        // for ( ; tmax < T / 2; tmax++) {
        //     double err=myres->comp_error(fi[tmax]);
        //     if (fi[tmax][Njack-1]<err) break;
        // }
        // printf("tmax=%d\n", tmax);
        // for (int j = 0;j < Njack;j++) {
        // amu[j] = int_scheme(0, T / 2 - 1, ft);
        for (int t_a = 1;t_a < T / 2;t_a++) {
            amut[t_a][j] = int_scheme(0, t_a, ft);
            amut[t_a][j] *= 4 * alpha_em * alpha_em *
                q2 / (muon_mass_MeV * muon_mass_MeV * (a[j] / 197.326963) * (a[j] / 197.326963));
        }
    }
    myres->copy(amu, amut[T / 2 - 1]);
    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < T / 2; t++) {
        fprintf(outfile, "%d   %.15g   %.15g\t", t, amut[t][Njack - 1], error_jackboot(resampling, Njack, amut[t]));
        fprintf(outfile, "%d   %.15g   %.15g\t", t, amu[Njack - 1], error_jackboot(resampling, Njack, amu));
        fprintf(outfile, "%.15g   %.15g\t", fi[t][Njack - 1], error_jackboot(resampling, Njack, fi[t]));
        fprintf(outfile, "%.15g   %.15g\t", corr_sub[t][Njack - 1], error_jackboot(resampling, Njack, corr_sub[t]));
        fprintf(outfile, "%.15g   %.15g\n", Kt[t][Njack - 1], error_jackboot(resampling, Njack, Kt[t]));
    }
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 1, T / 2 - 1, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", amu[Njack - 1], error_jackboot(resampling, Njack, amu));

    free(ft);
    free_2(T / 2, fi);
    free_2(T / 2, Kt);
    // free_2(T / 2, thetat);
    free_2(T / 2, corr_sub);
    free_2(T / 2, amut);
    return amu;
}



double* compute_amu_sd(double**** in, int id, int Njack, double* Z, double* a, double q2,
    double (*int_scheme)(int, int, double*), FILE* outfile, const char* description,
    const char* resampling, int isub, int tmin) {

    constexpr double d = 0.15;
    constexpr double t1_d = 0.4 / d;
    int T = file_head.l0;
    double** fi = double_malloc_2(T / 2, Njack);
    double* amu = (double*)malloc(sizeof(double) * Njack);
    double* ft = (double*)malloc(sizeof(double) * T / 2);

    double** Kt = double_malloc_2(T / 2, Njack);
    double** thetat = double_malloc_2(T / 2, Njack);
    double** corr_sub = double_malloc_2(T / 2, Njack);


    for (int j = 0;j < Njack;j++) {

        ft[0] = 0;
        for (int t_a = 1; t_a < T / 2; t_a++) {
            double t = t_a * a[j]; // time in fm.
            double z = muon_mass_MeV * (t / 197.326963);
            double K = z * z * kernel_K(z);
            double theta = gm2_step_function(t / 0.15, t1_d);
            double VV_sub;
            if (isub == -2)
                VV_sub = Z[j] * Z[j] * in[j][id][t_a][0];
            else if (isub == -1)
                VV_sub = Z[j] * Z[j] * in[j][id][t_a][0] - (1.0 / (2.0 * M_PI * M_PI * pow(t_a, 5)));// perturbative
            else
                VV_sub = Z[j] * Z[j] * in[j][id][t_a][0] + in[j][isub][t_a][0];// free-theory
            ft[t_a] = K * VV_sub * (1 - theta);

            fi[t_a][j] = ft[t_a];
            Kt[t_a][j] = K;
            thetat[t_a][j] = (1 - theta);
            corr_sub[t_a][j] = VV_sub;
            // printf("t=%d  K=%g   VV=%g  1-theta=%g  amu=%g\n",t_a,K, VV[t_a],1-theta, amu);
        }

        amu[j] = int_scheme(tmin, T / 2 - 1, ft);

        amu[j] *= 4 * alpha_em * alpha_em *
            q2 / (muon_mass_MeV * muon_mass_MeV * (a[j] / 197.326963) * (a[j] / 197.326963));
    }

    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < T / 2; t++) {
        fprintf(outfile, "%d   %.15g   %.15g\t", t, fi[t][Njack - 1], error_jackboot(resampling, Njack, fi[t]));
        fprintf(outfile, "%.15g   %.15g\t", Kt[t][Njack - 1], error_jackboot(resampling, Njack, Kt[t]));
        fprintf(outfile, "%.15g   %.15g\t", thetat[t][Njack - 1], error_jackboot(resampling, Njack, thetat[t]));
        fprintf(outfile, "%.15g   %.15g\n", corr_sub[t][Njack - 1], error_jackboot(resampling, Njack, corr_sub[t]));
    }
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 0, T / 2, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", amu[Njack - 1], error_jackboot(resampling, Njack, amu));

    free(ft);
    free_2(T / 2, fi);
    free_2(T / 2, Kt);
    free_2(T / 2, thetat);
    free_2(T / 2, corr_sub);
    return amu;
}



double* compute_amu_W(double**** in, int id, int Njack, double* Z, double* a, double q2, double (*int_scheme)(int, int, double*), FILE* outfile, const char* description, const char* resampling) {
    constexpr double d = 0.15;
    constexpr double t0_d = 0.4 / d;
    constexpr double t1_d = 1.0 / d;
    int T = file_head.l0;
    double** fi = double_malloc_2(T / 2, Njack);
    double* amu = (double*)malloc(sizeof(double) * Njack);
    double* ft = (double*)malloc(sizeof(double) * T / 2);

    double** Kt = double_malloc_2(T / 2, Njack);
    double** thetat = double_malloc_2(T / 2, Njack);
    double** corr_sub = double_malloc_2(T / 2, Njack);


    for (int j = 0;j < Njack;j++) {

        ft[0] = 0;
        for (int t_a = 1; t_a < T / 2; t_a++) {
            double t = t_a * a[j]; // time in fm.
            double z = muon_mass_MeV * (t / 197.326963);
            double K = z * z * kernel_K(z);
            double theta = gm2_step_function(t / 0.15, t0_d) - gm2_step_function(t / 0.15, t1_d);
            double VV_sub = Z[j] * Z[j] * in[j][id][t_a][0];//- (1.0 / (2.0 * M_PI * M_PI * pow(t_a, 5)));
            ft[t_a] = K * VV_sub * theta;

            fi[t_a][j] = ft[t_a];
            Kt[t_a][j] = K;
            thetat[t_a][j] = theta;
            corr_sub[t_a][j] = VV_sub;
            // printf("t=%d  K=%g   VV=%g  1-theta=%g  amu=%g\n",t_a,K, VV[t_a],1-theta, amu);
        }

        amu[j] = int_scheme(0, T / 2 - 1, ft);

        amu[j] *= 4 * alpha_em * alpha_em *
            q2 / (muon_mass_MeV * muon_mass_MeV * (a[j] / 197.326963) * (a[j] / 197.326963));
    }

    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < T / 2; t++) {
        fprintf(outfile, "%d   %.15g   %.15g\t", t, fi[t][Njack - 1], error_jackboot(resampling, Njack, fi[t]));
        fprintf(outfile, "%.15g   %.15g\t", Kt[t][Njack - 1], error_jackboot(resampling, Njack, Kt[t]));
        fprintf(outfile, "%.15g   %.15g\t", thetat[t][Njack - 1], error_jackboot(resampling, Njack, thetat[t]));
        fprintf(outfile, "%.15g   %.15g\n", corr_sub[t][Njack - 1], error_jackboot(resampling, Njack, corr_sub[t]));
    }
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 0, T / 2, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %.15g\n", amu[Njack - 1], error_jackboot(resampling, Njack, amu));

    free(ft);
    free_2(T / 2, fi);
    free_2(T / 2, Kt);
    free_2(T / 2, thetat);
    free_2(T / 2, corr_sub);
    return amu;
}

double* compute_amu_LD(double**** in, int id, int Njack, double* Z, double* a, double q2, double (*int_scheme)(int, int, double*), FILE* outfile, const char* description, const char* resampling) {
    constexpr double d = 0.15;
    constexpr double t0_d = 0.4 / d;
    constexpr double t1_d = 1.0 / d;
    int T = file_head.l0;
    double** fi = double_malloc_2(T / 2, Njack);
    double* amu = (double*)malloc(sizeof(double) * Njack);
    double* ft = (double*)malloc(sizeof(double) * T / 2);

    double** Kt = double_malloc_2(T / 2, Njack);
    double** thetat = double_malloc_2(T / 2, Njack);
    double** corr_sub = double_malloc_2(T / 2, Njack);


    for (int j = 0;j < Njack;j++) {

        ft[0] = 0;
        for (int t_a = 1; t_a < T / 2; t_a++) {
            double t = t_a * a[j]; // time in fm.
            double z = muon_mass_MeV * (t / 197.326963);
            double K = z * z * kernel_K(z);
            double theta = gm2_step_function(t / 0.15, t1_d);
            double VV_sub = Z[j] * Z[j] * in[j][id][t_a][0];//- (1.0 / (2.0 * M_PI * M_PI * pow(t_a, 5)));
            ft[t_a] = K * VV_sub * theta;

            fi[t_a][j] = ft[t_a];
            Kt[t_a][j] = K;
            thetat[t_a][j] = theta;
            corr_sub[t_a][j] = VV_sub;
            // printf("t=%d  K=%g   VV=%g  1-theta=%g  amu=%g\n",t_a,K, VV[t_a],1-theta, amu);
        }

        amu[j] = int_scheme(0, T / 2 - 1, ft);

        amu[j] *= 4 * alpha_em * alpha_em *
            q2 / (muon_mass_MeV * muon_mass_MeV * (a[j] / 197.326963) * (a[j] / 197.326963));
    }

    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < T / 2; t++) {
        fprintf(outfile, "%d   %.15g   %.15g\t", t, fi[t][Njack - 1], error_jackboot(resampling, Njack, fi[t]));
        fprintf(outfile, "%.15g   %.15g\t", Kt[t][Njack - 1], error_jackboot(resampling, Njack, Kt[t]));
        fprintf(outfile, "%.15g   %.15g\t", thetat[t][Njack - 1], error_jackboot(resampling, Njack, thetat[t]));
        fprintf(outfile, "%.15g   %.15g\n", corr_sub[t][Njack - 1], error_jackboot(resampling, Njack, corr_sub[t]));
    }
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 0, T / 2, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %.15g\n", amu[Njack - 1], error_jackboot(resampling, Njack, amu));

    free(ft);
    free_2(T / 2, fi);
    free_2(T / 2, Kt);
    free_2(T / 2, thetat);
    free_2(T / 2, corr_sub);
    return amu;
}

double* compute_amu_bounding(double**** in, int id, int Njack, double* Z, double* a, double q2, double (*int_scheme)(int, int, double*), FILE* outfile,
    const char* description, const char* resampling, int isub, int bound, double* Mpi, double* Meff, int L) {
    constexpr double d = 0.15;
    constexpr double t1_d = 0.4 / d;
    int T = file_head.l0;
    double** fi = double_malloc_2(T / 2, Njack);
    double** amu = malloc_2<double>(T / 2 + 10, Njack);
    double** amu_above = malloc_2<double>(T / 2, Njack);
    double** amu_below = malloc_2<double>(T / 2, Njack);
    double* ft = (double*)malloc(sizeof(double) * T / 2);
    double* above = (double*)calloc(T / 2, sizeof(double));
    double* below = (double*)calloc(T / 2, sizeof(double));
    double** mass = malloc_2<double>(T / 2, Njack);

    double** Kt = double_malloc_2(T / 2, Njack);
    // double** thetat = double_malloc_2(T / 2, Njack);
    double** corr_sub = double_malloc_2(T / 2, Njack);
    double* E2 = (double*)malloc(sizeof(double) * Njack);

    int t_end = T / 2 - 1;
    for (int j = 0;j < Njack;j++) {
        // compute the effective mass
        for (int t_c = 1; t_c < T / 2; t_c++) {
            // mass[t_c][j] = M_eff_T_ct_ctp1(t_c, T, in[j][id][t_c][0], in[j][id][(t_c + 1) % T][0]);
            mass[t_c][j] = log(in[j][id][t_c][0] / in[j][id][(t_c + 1) % T][0]);
        }
    }

    for (int t_c = 2; t_c < T / 2; t_c++) {
        // when the effective mass loose signal
        if (mass[t_c][Njack - 1] > mass[t_c - 1][Njack - 1] ||
            myres->comp_error(mass[t_c]) > fabs(mass[t_c][Njack - 1]) * 0.010) {
            t_end = t_c - 1;
            break;
        }
    }

    printf("t_end=%d\n", t_end);

    for (int j = 0;j < Njack;j++) {
        E2[j] = 2 * sqrt(Mpi[j] * Mpi[j] + (2 * M_PI / L) * (2 * M_PI / L));
        ft[0] = 0;
        for (int t_a = 1; t_a < T / 2; t_a++) {
            double t = t_a * a[j]; // time in fm.
            double z = muon_mass_MeV * (t / 197.326963);
            double K = z * z * kernel_K(z);
            // double theta = gm2_step_function(t / 0.15, t1_d);
            double VV_sub;
            if (isub == -2)
                VV_sub = Z[j] * Z[j] * in[j][id][t_a][0];
            else if (isub == -1)
                VV_sub = Z[j] * Z[j] * in[j][id][t_a][0] - (1.0 / (2.0 * M_PI * M_PI * pow(t_a, 5)));// perturbative
            else
                VV_sub = Z[j] * Z[j] * in[j][id][t_a][0] + in[j][isub][t_a][0];// free-theory
            ft[t_a] = VV_sub;

            fi[t_a][j] = ft[t_a];
            Kt[t_a][j] = K;
            // thetat[t_a][j] = (1 - theta);
            corr_sub[t_a][j] = VV_sub;
            // printf("t=%d  K=%g   VV=%g  1-theta=%g  amu=%g\n",t_a,K, VV[t_a],1-theta, amu);

        }
        double coef = 4 * alpha_em * alpha_em *
            q2 / (muon_mass_MeV * muon_mass_MeV * (a[j] / 197.326963) * (a[j] / 197.326963));
        above[0] = 0;
        below[0] = 0;
        for (int t_c = 1; t_c < T / 2; t_c++) {
            above[t_c] = Kt[t_c][j] * ft[t_c];
            below[t_c] = Kt[t_c][j] * ft[t_c];

            double min_mass;
            // we do not propagate the error on the mass
            if (t_c < t_end) min_mass = mass[t_c][j];
            else min_mass = mass[t_end][j];

            for (int tp = t_c + 1; tp < T / 2; tp++) {
                above[tp] = Kt[tp][j] * ft[t_c] * exp(-E2[j] * (tp - t_c));
                if (bound == 1) {
                    above[tp] = Kt[tp][j] * ft[t_c] * exp(-Meff[j] * (tp - t_c));
                    below[tp] = Kt[tp][j] * ft[t_c] * (exp(-min_mass * (tp - t_c))  /* +  exp(-min_mass * (T-tp - t_c))  */);
                }
                else if (bound == 2) {
                    below[tp] = 0;
                }
                else if (bound == 3) {
                    below[tp] = Kt[tp][j] * ft[t_c] * (exp(-min_mass * (tp - t_c))  /* +  exp(-min_mass * (T-tp - t_c)) */);
                }
                else {
                    printf("%s: bound method not valid bound = %d\n", __func__, bound);
                    exit(1);
                }
            }
            amu_above[t_c][j] = int_scheme(0, T / 2 - 1, above);
            amu_below[t_c][j] = int_scheme(0, T / 2 - 1, below);

            amu_above[t_c][j] *= coef;
            amu_below[t_c][j] *= coef;
            amu[t_c][j] = (amu_above[t_c][j] + amu_below[t_c][j]) / 2.0;
        }
        amu[0][j] = 0;
        // for (int t_c = 1; t_c < T / 2; t_c++) {
        //     amu_above[t_c][j] = Kt[t_c][j] * ft[t_c];
        //     amu_below[t_c][j] = Kt[t_c][j] * ft[t_c];
        // }
        // amu[j] = int_scheme(0, T / 2 - 1, ft);
    }

    /////////// find tmin 
    int start_fit = 1;
    for (int t_c = 1; t_c < T / 2; t_c++) {
        if (amu_above[t_c][Njack - 1] - amu_below[t_c][Njack - 1] < myres->comp_error(amu_above[t_c])) {
            start_fit = t_c;
            break;
        }
    }
    int end_fit = start_fit + 0.40926 / a[Njack - 1];
    error(end_fit > T / 2, 1, "compute_amu_bounding", "tmax fit larger than T/2");
    /////// plateau fit
    struct fit_type const_fit_info;
    const_fit_info.Nvar = 1;
    const_fit_info.Npar = 1;
    const_fit_info.N = 1;
    const_fit_info.Njack = Njack;
    const_fit_info.function = constant_fit;
    const_fit_info.n_ext_P = 0;
    const_fit_info.linear_fit = true;
    const_fit_info.T = T;
    struct kinematic kinematic_2pt;
    const_fit_info.codeplateaux = true;
    const_fit_info.tmin = start_fit;
    const_fit_info.tmax = end_fit;
    double** mt = malloc_2<double>(T / 2, 2);
    for (int t = 0;t < T / 2; t++) {
        mt[t][0] = amu[t][Njack - 1];
        mt[t][1] = myres->comp_error(amu[t]);
    }
    int store_l0 = file_head.l0;
    file_head.l0 = T;

    char name[NAMESIZE] = "HLT";
    double* chi2;
    char** option;
    struct fit_result fit_out = try_fit(option, const_fit_info.tmin, const_fit_info.tmax, 1, mt, amu, Njack, &chi2, const_fit_info);// option is not used 
    file_head.l0 = store_l0;

    //////////////////////////////////////////////////////////////
    // printing
    //////////////////////////////////////////////////////////////
    // above
    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < T / 2; t++) {
        fprintf(outfile, "%d   %.15g   %.15g \t", t, amu_above[t][Njack - 1], myres->comp_error(amu_above[t]));
        fprintf(outfile, "%d   %.15g   %.15g\n", t, fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]));
    }
    fprintf(outfile, "\n\n #%s_above fit in [%d,%d] chi2=%.5g  %.5g\n", description, start_fit, end_fit, fit_out.chi2[Njack - 1], 0.0);
    fprintf(outfile, "   %.15g   %.15g   %d\n", fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]), t_end);

    // average

    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < T / 2; t++) {
        fprintf(outfile, "%d   %.15g   %.15g \t", t, amu[t][Njack - 1], myres->comp_error(amu[t]));
        fprintf(outfile, "%d   %.15g   %.15g\n", t, fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]));
    }
    fprintf(outfile, "\n\n #%s_ave fit in [%d,%d] chi2=%.5g  %.5g\n", description, start_fit, end_fit, fit_out.chi2[Njack - 1], 0.0);
    fprintf(outfile, "   %.15g   %.15g   %d\n", fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]), t_end);

    // below

    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < T / 2; t++) {
        fprintf(outfile, "%d   %.15g   %.15g \t", t, amu_below[t][Njack - 1], myres->comp_error(amu_below[t]));
        fprintf(outfile, "%d   %.15g   %.15g\n", t, fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]));
    }
    fprintf(outfile, "\n\n #%s_below fit in [%d,%d] chi2=%.5g  %.5g\n", description, start_fit, end_fit, fit_out.chi2[Njack - 1], 0.0);
    fprintf(outfile, "   %.15g   %.15g   %d\n", fit_out.P[0][Njack - 1], myres->comp_error(fit_out.P[0]), t_end);


    file_head.l0 = store_l0;
    double* P = (double*)malloc(sizeof(double) * Njack);
    for (size_t j = 0; j < Njack; j++) {
        P[j] = fit_out.P[0][j];
    }

    fit_out.clear();
    free_2(T / 2, mt);
    free_2(T / 2, amu);
    free_2(T / 2, amu_above);
    free_2(T / 2, amu_below);
    free_2(T / 2, mass);

    free(ft);
    free(above);
    free(below);
    free_2(T / 2, fi);
    free_2(T / 2, Kt);
    // free_2(T / 2, thetat);
    free_2(T / 2, corr_sub);

    return P;
}


/********************************************************************************
 * Z
 ********************************************************************************/
template<int idn, int idd>
double ZAl_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double M_PS = fit_info.ext_P[0][j];
    double M_PS_OS = fit_info.ext_P[1][j];
    double GPS = fit_info.ext_P[2][j];
    double GPS_OS = fit_info.ext_P[3][j];


    double num = in[j][idn][t][0];
    double den = (in[j][idd][t + 1][0] - in[j][idd][t - 1][0]) / 2.0;
    double RA = (2 * fit_info.mu * num / (den));

    return RA * M_PS_OS * sinh(M_PS_OS) * GPS / (M_PS * sinh(M_PS) * GPS_OS);

}
template double ZVl_lhs<4, 3>(int, double****, int, fit_type);
template double ZVl_lhs<10, 9>(int, double****, int, fit_type);
template double ZVl_lhs<16, 15>(int, double****, int, fit_type);

template<int idn, int idd>
double ZVl_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;

    double num = in[j][idn][t][0];
    double den = (in[j][idd][t + 1][0] - in[j][idd][t - 1][0]) / 2.0;
    return (2 * fit_info.mu * num / (den));

}
template double ZAl_lhs<1, 0>(int, double****, int, fit_type);
template double ZAl_lhs<7, 6>(int, double****, int, fit_type);
template double ZAl_lhs<13, 12>(int, double****, int, fit_type);

template<int id>
double GPS_OS_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double mu = file_head.musea;
    double mass = fit_info.ext_P[0][j];

    double GPS_OS = sqrt(in[j][id][t][0] * 2 * mass / (exp(-mass * t) + exp(-mass * (T - t))));

    return GPS_OS;
}
template double GPS_OS_lhs<4>(int, double****, int, fit_type);
template double GPS_OS_lhs<10>(int, double****, int, fit_type);
template double GPS_OS_lhs<16>(int, double****, int, fit_type);

template<int id>
double GPS_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double mu = file_head.musea;
    double mass = fit_info.ext_P[0][j];

    double GPS_OS = sqrt(in[j][id][t][0] * 2 * mass / (exp(-mass * t) + exp(-mass * (T - t))));

    return GPS_OS;
}
template double GPS_lhs<1>(int, double****, int, fit_type);
template double GPS_lhs<7>(int, double****, int, fit_type);
template double GPS_lhs<13>(int, double****, int, fit_type);

template<int id>
double lhs_ct(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;

    double a = fit_info.ext_P[0][j];
    double Z = fit_info.ext_P[1][j];

    double GPS_OS = Z * in[j][id][t][0] / ((a / 197.326963) * (a / 197.326963) * (a / 197.326963));

    return GPS_OS;
}
template double lhs_ct<2>(int, double****, int, fit_type);
template double lhs_ct<5>(int, double****, int, fit_type);
template double lhs_ct<8>(int, double****, int, fit_type);
template double lhs_ct<14>(int, double****, int, fit_type);

double* interpol_Z(int Nmus, int Njack, double** Meta, double** Z, double* aMetas_exp,
    FILE* outfile, const char* description, const char* resampling) {
    //  Z(s1) = ( 1  Meta(s1) )  (P[0])
    //  Z(s2) = ( 1  Meta(s2) )  (P[1])

    double** matrix = double_malloc_2(Nmus, Nmus);
    double* Zj = (double*)malloc(sizeof(double) * Nmus);
    double* Zint = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++) {
        // printf("######################\njack %d\n",j);
        for (int i = 0;i < Nmus;i++) {
            Zj[i] = Z[i][j];
            // printf("%20.12g  =", Zj[i]);
            for (int k = 0;k < Nmus;k++) {
                matrix[i][k] = pow(Meta[i][j], k);
                // printf("%20.12g ", matrix[i][k]);
            }
            // printf("\n");
        }
        double* P = LU_decomposition_solver(Nmus, matrix, Zj);
        Zint[j] = 0;
        // printf("P:");
        for (int k = 0;k < Nmus;k++) {
            // printf("%20.12g  ", P[k]);
            Zint[j] += pow(aMetas_exp[j], k) * P[k];
        }
        // printf("\n   ms_MK=%20.12g  %20.12g \n",aMetas_exp[j], error_jackboot(resampling, Njack, aMetas_exp));
        // printf("\n   value=%20.12g  \n",Zint[j]);
        if (Nmus == 2) {
            Zint[j] = ((Z[0][j] - Z[1][j]) * aMetas_exp[j] - Z[0][j] * Meta[1][j] + Z[1][j] * Meta[0][j]) / (Meta[0][j] - Meta[1][j]);
            // printf("\n   value exact=%20.12g  \n",Zint[j]);
        }
        // printf("%20.12g  %20.12g\n",Z[0][j], Z[1][j] );
        free(P);
    }
    free(Zj);
    free_2(Nmus, matrix);
    fprintf(outfile, " \n\n# aMetas_exp  Zint  err\n");

    fprintf(outfile, "%.15g   %.15g   %.15g\t", aMetas_exp[Njack - 1], Zint[Njack - 1], error_jackboot(resampling, Njack, Zint));

    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 0, 0, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %.15g\n", Zint[Njack - 1], error_jackboot(resampling, Njack, Zint));
    printf("%s (%.15g) =  %.15g   %.15g\n", description, aMetas_exp[Njack - 1], Zint[Njack - 1], error_jackboot(resampling, Njack, Zint));

    return(Zint);
}
/////////////////////////////////////////////////////////////////////////
// compute_DV
/////////////////////////////////////////////////////////////////////////

double compute_cotd(int Nvar, double* x) {
    double omega = x[0];
    double mass = x[1];
    double L = x[2];
    double Mrho = x[3];
    double grhopipi = x[4];
    double a = x[5];

    double k = sqrt(omega * omega / 4 - mass * mass);

    double ho = (grhopipi * grhopipi * k * k * k * 2) * log((omega + 2 * k) / (2 * mass)) / (6.0 * M_PI * M_PI * omega);
    double kM = sqrt(Mrho * Mrho / 4. - mass * mass);
    double hM = (grhopipi * grhopipi * kM * kM * kM * 2) * log((Mrho + 2 * kM) / (2 * mass)) / (6.0 * M_PI * M_PI * Mrho);
    // double hM = (grhopipi * grhopipi * kM * kM * kM * 2) * log((Mrho + 2 * kM) / (2 * mass)) / (6.0 * M_PI * M_PI * Mrho);
    double h1M = (grhopipi * grhopipi * kM * kM) * (1 + (1 + 2 * mass * mass / (Mrho * Mrho)) * (Mrho / kM) * log((Mrho + 2 * kM) / (2 * mass))) / (6.0 * M_PI * M_PI * Mrho);
    double gamma = (grhopipi * grhopipi * k * k * k) / (6.0 * M_PI * omega * omega);

    double cotd = Mrho * Mrho - omega * omega - hM - (omega * omega - Mrho * Mrho) * h1M / (2 * Mrho) + ho;
    cotd /= omega * gamma;
    return cotd;
}

double compute_delta(double omega, void* params) {
    double* P = (double*)params;
    int Nvar = P[0] + 1e-6;
    double* xh = (double*)malloc(sizeof(double) * Nvar);
    for (int i = 1;i < Nvar;i++) {
        xh[i] = P[i];
    }
    xh[0] = omega;
    // printf("Nvar=%d  omega=%g   P1=%g L=%g Mrho=%g grpp=%g a =%g\n", Nvar, xh[0], xh[1], xh[2], xh[3], xh[4], xh[5]);
    return  std::atan(1. / compute_cotd(Nvar, xh));
}

double compute_deriv_delta(int Nvar, double* x) {

    gsl_function F;
    double result, abserr;

    F.function = &compute_delta;
    double* P;
    P = (double*)malloc(sizeof(double) * Nvar);
    for (int i = 0;i < Nvar;i++) {
        P[i] = x[i];
    }
    P[0] = Nvar;
    F.params = (void*)P;
    gsl_deriv_central(&F, x[0], 1e-3, &result, &abserr);
    // printf("derivdelta= %g %g\n", result, abserr);
    // result=(F.function(x[0]+1e-8,F.params)-F.function(x[0],F.params))/1e-8;
    double k = sqrt(x[0] * x[0] / 4. - x[1] * x[1]);
    result *= 2 * k / (sqrt(k * k + x[1] * x[1]));

    free(P);
    return result;
}


double w_js(int j, int s, double q) {
    double z[2] = { 0, 0 };
    int dvec[3] = { 0, 0, 0 };
    double gamma = 1;
    double A = 1;
    dzeta_function(z, q * q, j, s, dvec, gamma, A, 1.e-3, 1.e+6, 5);
    // std::complex<double>  zc(z[0], z[1]);
    double r = z[0] / (pow(q, j + 1) * sqrt(2. * j + 1.) * pow(pi_greco, 3. / 2.));
    return r;
}

double Z00(double q) {
    double z[2] = { 0, 0 };
    int dvec[3] = { 0, 0, 0 };
    double gamma = 1;
    double A = 1;
    dzeta_function(z, q * q, 0, 0, dvec, gamma, A, 1.e-3, 1.e-6, 12);
    // qzeta_function(z, q * q, 0, 0, dvec, gamma, A, 1.e-3, 1.e-6, 12);
    return z[0];
}

double omega_QC2(int n, int Nvar, double* x, int Npar, double* P) {
    double omega = x[0];
    double mass = x[1];
    double L = x[2];
    double Mrho = x[3];
    double grhopipi = x[4];
    double a = x[5];
    double k = sqrt(omega * omega / 4. - mass * mass);
    double q = k * (L * (a / 197.326963) / (2. * pi_greco));

    // double r = w_js(0, 0, q) - w_js(2, 0, q) - (3. / sqrt(6)) * (w_js(2, -2, q) + w_js(2, 2, q));
    // double r = w_js(0, 0, q) / (pow(pi_greco, 3. / 2.) * q);
    double r = Z00(q) / (pow(pi_greco, 3. / 2.) * q);

    double cotd = compute_cotd(Nvar, x);
    return cotd - r;

}

double compute_phi(double q, void* params) {
    (void)(params); /* avoid unused parameter warning */
    // double M = w_js(0, 0, q) - w_js(2, 0, q) - (3. / sqrt(6)) * (w_js(2, -2, q) + w_js(2, 2, q));
    return std::atan(-pow(pi_greco, 3. / 2.) * q / Z00(q));
}


double compute_deriv_phi(double q) {

    gsl_function F;
    double result, abserr;

    F.function = &compute_phi;
    F.params = 0;
    gsl_deriv_central(&F, q, 1e-3, &result, &abserr);
    // result=(F.function(q+1e-8,F.params)-F.function(q,F.params))/1e-8;
    // printf("derivphi= %g  %g\n", result, abserr);
    return result;
}

double matrix_element_nuA2(int Nvar, double* x) {
    double omega = x[0];
    double mass = x[1];
    double L = x[2];
    double Mrho = x[3];
    double grhopipi = x[4];
    double a = x[5];
    double k = sqrt(omega * omega / 4. - mass * mass);
    double q = k * (L * (a / 197.326963) / (2. * pi_greco));


    double FGS2 = compute_FGS2(Nvar, x);

    double r = 2 * pow(k, 5) * FGS2 / (3. * M_PI * omega * omega);

    double deriv_d = compute_deriv_delta(Nvar, x);
    double deriv_phi = compute_deriv_phi(q);
    // printf(" r= %g / ( %g  %g  %g  %g) \n", r, k, deriv_d, q, deriv_phi);
    r /= (k * deriv_d + q * deriv_phi);

    return r;
}


double integrand_V_infL(double omega, void* params) {


    double* x = (double*)params;
    double mass = x[1];
    double L = x[2];
    double Mrho = x[3];
    double grhopipi = x[4];
    double a = x[5];

    double t = x[6];
    x[0] = omega;
    double FGS2 = compute_FGS2(7, x);
    // if (once == 0) {
    //     x[0] = 500;
    //     printf(" t= %g\n", t);
    //     once = -1;
    //     printf("FGS2=%.15g\n", compute_FGS2(7, x));
    //     once = 1;
    // }

    return omega * omega * pow(1 - pow(2 * mass / omega, 2), 3. / 2.) * FGS2 * exp(-omega * t);
    // return  exp(-omega * t);

}
double* compute_V_infL(int Nvar, double* x, int T) {
    double omega = x[0];
    double mass = x[1];
    double L = x[2];
    double Mrho = x[3];
    double grhopipi = x[4];
    double a = x[5];

    double* P = (double*)malloc(sizeof(double) * (Nvar + 1));
    for (int i = 0;i < Nvar;i++) {
        P[i] = x[i];
    }
    double* V = (double*)malloc(sizeof(double) * T);
    V[0] = 0;
    int Maxiter = 1e+8;
    double epsrel = 1e-8;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(Maxiter);

    for (int t = 1; t < T;t++) {
        P[Nvar] = t * (a / 197.326963);

        double result, error;
        gsl_function F;
        F.function = &integrand_V_infL;
        F.params = (void*)P;
        // printf("here : %g\n",P[Nvar]);
                //  gsl_integration_qags(&F, 2*mass, 10000, 1e+5, epsrel, Maxiter, w, &result, &error);
        gsl_integration_qagiu(&F, 2 * mass, 0, epsrel, Maxiter, w, &result, &error);
        // printf("K(%g)=%g\n", z, result);

        V[t] = result / (48. * M_PI * M_PI);
    }
    gsl_integration_workspace_free(w);
    free(P);
    return V;
}
double compute_V_infL_t(int Nvar, double* x, double t) {
    double omega = x[0];
    double mass = x[1];
    double L = x[2];
    double Mrho = x[3];
    double grhopipi = x[4];
    double a = x[5];

    double* P = (double*)malloc(sizeof(double) * (Nvar + 1));
    for (int i = 0;i < Nvar;i++) {
        P[i] = x[i];
    }
    double V = 0;

    int Maxiter = 1e+6;
    double epsrel = 1e-6;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(Maxiter);

    P[Nvar] = t;

    double result, error;
    gsl_function F;
    F.function = &integrand_V_infL;
    F.params = (void*)P;
    // printf("here : %g\n",P[Nvar]);
            //  gsl_integration_qags(&F, 2*mass, 10000, 1e+5, epsrel, Maxiter, w, &result, &error);
    gsl_integration_qagiu(&F, 2 * mass, 0, epsrel, Maxiter, w, &result, &error);
    // printf("K(%g)=%g\n", z, result);

    V = result / (48. * M_PI * M_PI);

    gsl_integration_workspace_free(w);
    free(P);
    return V;
}


double* compute_VL(int n, double** omega, double** nuA2, double* a, int T, int j) {
    double* VL = (double*)calloc(T, sizeof(double));
    for (int t = 0; t < T;t++) {
        double ta = t * (a[j] / 197.326963);
        for (int in = 0;in < n;in++) {
            VL[t] += nuA2[in][j] * exp(-omega[in][j] * ta);

        }

    }
    return VL;

}

double compute_VL_t(int n, double* omega, double* nuA2, double t) {
    double VL = 0;
    for (int in = 0;in < n;in++) {
        VL += nuA2[in] * exp(-omega[in] * t);
    }
    return VL;
}


double** compute_DVt(int L, int Njack, double* Mpi, double* Mrho, double* a, double* grhopipi, FILE* outfile, const char* description, const char* resampling) {

    int T = file_head.l0;
    constexpr int Nvar = 6;
    constexpr int Nomegas = 25;
    double** omega = double_malloc_2(Nomegas, Njack);
    double** nuA2 = double_malloc_2(Nomegas, Njack);
    double** DVt = double_malloc_2(T, Njack);
    double** VinfL = (double**)malloc(sizeof(double*) * Njack);
    double** VL = (double**)malloc(sizeof(double*) * Njack);


    int p2[Nomegas * Nomegas * Nomegas];
    int count = 0;
    for (int i = 0;i < Nomegas;i++) {
        for (int j = 0;j < Nomegas;j++) {
            for (int k = 0;k < Nomegas;k++) {
                p2[count] = i * i + j * j + k * k;
                count++;
            }
        }
    }
    std::sort(p2, p2 + Nomegas * Nomegas * Nomegas);
    int p2s[Nomegas + 1];
    count = 0;
    p2s[0] = 0;
    int old = -1;
    while (count < Nomegas + 1) {
        for (int i = 0;i < Nomegas * Nomegas * Nomegas;i++) {
            if (p2[i] > old) {
                p2s[count] = p2[i];
                count++;
                old = p2[i];
            }
        }
    }
    // for (int i = 0;i < Nomegas+1;i++)
    //     printf("sorted: %d\n", p2s[i]);


    fprintf(outfile, " \n\n# omega  cotd  M11\n");
    for (int j = Njack - 1;j < Njack;j++) {
        double x[Nvar];
        x[0] = 1;
        x[1] = Mpi[j];
        x[2] = L;
        x[3] = Mrho[j];
        x[4] = grhopipi[j];
        x[5] = a[j];
        // first we plot cotd*k/m
        if (j == Njack - 1) {
            for (int i = 0;i < 1500;i++) {
                x[0] = x[1] * (2. + 30.5 * i / 1500.0);
                double k = sqrt(x[0] * x[0] / 4 - x[1] * x[1]);
                double cotd = compute_cotd(Nvar, x);
                double q = k * (L * (x[5] / 197.326963) / (2. * pi_greco));
                // double r = w_js(0, 0, q) - w_js(2, 0, q) - (3. / sqrt(6)) * (w_js(2, -2, q) + w_js(2, 2, q));
                double r = Z00(q) / (pow(pi_greco, 3. / 2.) * q);
                fprintf(outfile, "%g    %g  %g\n", x[0], cotd, r);
            }

            fprintf(outfile, "\n\n #QC2 fit in [%d,%d] chi2=%.5g  %.5g\n", 0, 0, 0.0, 0.0);

        }
        double p = 2 * M_PI / ((double)L * (x[5] / 197.326963));

        for (int i = 0; i < Nomegas;i++) {
            double xmin = 2 * sqrt(x[1] * x[1] + p2s[i] * p * p) + 1e-8;
            double xmax = 2 * sqrt(x[1] * x[1] + p2s[i + 1] * p * p) - 1e-8;
            // printf("%d   %d\n", p2s[i], p2s[i + 1]);
            double oj = rtsafe(omega_QC2, 0/*n*/, 6, x, 0/* Npar */, nullptr/* Pkcot */, 0/*ivar*/, 0. /*input*/, xmin, xmax, 1e-5, 100, 1e-4);
            int seed = rand();
            omega[i] = fake_sampling(resampling, oj, oj / 1e+6, Njack, seed);
            // if (j==Njack-1) printf("omega[%d]=%f  x=%f  %f\n",i, oj, xmin,xmax);
            x[0] = omega[i][j];
            double A = matrix_element_nuA2(Nvar, x);
            nuA2[i] = fake_sampling(resampling, A, A / 1e+6, Njack, seed + 1);
            count++;
        }
    }

    for (int j = 0;j < Njack;j++) {
        double x[Nvar];
        x[0] = 1;
        x[1] = Mpi[j];
        x[2] = L;
        x[3] = Mrho[j];
        x[4] = grhopipi[j];
        x[5] = a[j];

        VinfL[j] = compute_V_infL(Nvar, x, T);
        VL[j] = compute_VL(Nomegas, omega, nuA2, a, T, j);
        for (int t = 0;t < T;t++) {
            DVt[t][j] = (VinfL[j][t] - VL[j][t]) * pow(a[j] / 197.326963, 3);
        }
    }


    for (int i = 0; i < Nomegas;i++) {
        printf("omega nuA2[%d]=%g  %g\n", i, omega[i][Njack - 1], nuA2[i][Njack - 1]);
        // printf("nuA2[%d]=%g  %g\n", i, nuA2[i][Njack - 1], error_jackboot(resampling, Njack, nuA2[i]));
        fprintf(outfile, "%g  %g\n", omega[i][Njack - 1], error_jackboot(resampling, Njack, omega[i]));
    }

    fprintf(outfile, " \n\n# t  inf  L\n");
    for (int t = 1;t < T;t++) {
        fprintf(outfile, "%d   %g   %g\n", t, VinfL[Njack - 1][t], VL[Njack - 1][t]);
    }
    fprintf(outfile, "\n\n #DVt fit in [%d,%d] chi2=%.5g  %.5g\n", 0, 0, 0.0, 0.0);
    fprintf(outfile, "%.5g  %.5g\n", 0.0, 0.0);

    // fprintf(outfile, " \n\n");
    // fprintf(outfile, "#\n");
    // for (int t = 1; t < T / 2; t++) {
    //      fprintf(outfile, "%d   %.15g   %.15g\t", t, fi[t][Njack - 1], error_jackboot(resampling, Njack, fi[t]));
    // }
    // fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 0, T / 2, 0.0, 0.0);
    // fprintf(outfile, "   %.15g   %15.g\n", amu[Njack - 1], error_jackboot(resampling, Njack, amu));

    // free_2(5, omega);
    // free_2(5, nuA2);
    return DVt;
}


double integrand_DV(double t, void* params) {
    // printf("t=%g\n", t);
    if (t < 1e-6) {
        return 0;
    }

    constexpr int Nvar = 6;
    double* x = (double*)params;
    int Nomegas = (int)(x[Nvar] + 0.001);
    double* omega = (double*)malloc(sizeof(double) * Nomegas);
    double* nuA2 = (double*)malloc(sizeof(double) * Nomegas);
    for (int i = 0; i < Nomegas;i++) {
        omega[i] = x[Nvar + 1 + i];
        nuA2[i] = x[Nvar + 1 + Nomegas + i];

    }

    double VinfL = compute_V_infL_t(Nvar, x, t);
    // printf("VinfL=%g\n", VinfL);
    double VL = compute_VL_t(Nomegas, omega, nuA2, t);
    // printf("VL=%g\n", VL);
    double z = muon_mass_MeV * t;
    double K = z * z * kernel_K(z);
    // printf("K=%g\n", K);

    constexpr double d = 0.15;
    constexpr double t0_d = 0.4 / d;
    constexpr double t1_d = 1.0 / d;
    double theta = gm2_step_function((t * 197.326963) / 0.15, t0_d) - gm2_step_function((t * 197.326963) / 0.15, t1_d);

    double V_sub = (VinfL - VL) * K * theta;
    double q2 = 1;

    V_sub *= 4 * alpha_em * alpha_em * q2 / ((muon_mass_MeV) * (muon_mass_MeV));

    free(omega);
    free(nuA2);
    // printf("%g  %g\n", t, V_sub);
    return V_sub;
}
double integrand_DV_full(double t, void* params) {
    // printf("t=%g\n", t);
    if (t < 1e-6) {
        return 0;
    }

    constexpr int Nvar = 6;
    double* x = (double*)params;
    int Nomegas = (int)(x[Nvar] + 0.001);
    double* omega = (double*)malloc(sizeof(double) * Nomegas);
    double* nuA2 = (double*)malloc(sizeof(double) * Nomegas);
    for (int i = 0; i < Nomegas;i++) {
        omega[i] = x[Nvar + 1 + i];
        nuA2[i] = x[Nvar + 1 + Nomegas + i];

    }

    double VinfL = compute_V_infL_t(Nvar, x, t);
    // printf("VinfL=%g\n", VinfL);
    double VL = compute_VL_t(Nomegas, omega, nuA2, t);
    // printf("VL=%g\n", VL);
    double z = muon_mass_MeV * t;
    double K = z * z * kernel_K(z);
    // printf("K=%g\n", K);

    constexpr double d = 0.15;
    constexpr double t0_d = 0.4 / d;
    constexpr double t1_d = 1.0 / d;
    double theta = 1.0; //gm2_step_function((t * 197.326963) / 0.15, t0_d) - gm2_step_function((t * 197.326963) / 0.15, t1_d);

    double V_sub = (VinfL - VL) * K * theta;
    double q2 = 1;

    V_sub *= 4 * alpha_em * alpha_em * q2 / ((muon_mass_MeV) * (muon_mass_MeV));

    free(omega);
    free(nuA2);
    // printf("%g  %g\n", t, V_sub);
    return V_sub;
}

double integrand_DV_M(double t, void* params) {
    // printf("t=%g\n", t);
    if (t < 1e-6) {
        return 0;
    }

    constexpr int Nvar = 6;
    double* x = (double*)params;
    int Nomegas = (int)(x[Nvar] + 0.001);
    double* omega = (double*)malloc(sizeof(double) * Nomegas);
    double* nuA2 = (double*)malloc(sizeof(double) * Nomegas);
    for (int i = 0; i < Nomegas;i++) {
        omega[i] = x[Nvar + 1 + i];
        nuA2[i] = x[Nvar + 1 + Nomegas + i];

    }

    // double VinfL = compute_V_infL_t(Nvar, x, t);
    // printf("VinfL=%g\n", VinfL);
    double VL = compute_VL_t(Nomegas, omega, nuA2, t);
    // printf("VL=%g\n", VL);
    double z = muon_mass_MeV * t;
    double K = z * z * kernel_K(z);
    // printf("K=%g\n", K);

    constexpr double d = 0.15;
    constexpr double t0_d = 0.4 / d;
    constexpr double t1_d = 1.0 / d;
    double theta = 1.0; //gm2_step_function((t * 197.326963) / 0.15, t0_d) - gm2_step_function((t * 197.326963) / 0.15, t1_d);

    double V_sub = (VL)*K * theta;
    double q2 = 1;

    V_sub *= 4 * alpha_em * alpha_em * q2 / ((muon_mass_MeV) * (muon_mass_MeV));

    free(omega);
    free(nuA2);
    // printf("%g  %g\n", t, V_sub);
    return V_sub;
}

double* compute_DVt_and_integrate(int L, int Njack, double* Mpi, double* Mrho, double* a, double* grhopipi, FILE* outfile, const char* description,
    const char* resampling, int Nomegas, int window) {

    int T = file_head.l0;
    constexpr int Nvar = 6;
    double** omega = double_malloc_2(Nomegas, Njack);
    double** nuA2 = double_malloc_2(Nomegas, Njack);
    double** DVt = double_malloc_2(T, Njack);
    double** VinfL = (double**)malloc(sizeof(double*) * Njack);
    double** VL = (double**)malloc(sizeof(double*) * Njack);


    int p2[Nomegas * Nomegas * Nomegas];
    int count = 0;
    for (int i = 0;i < Nomegas;i++) {
        for (int j = 0;j < Nomegas;j++) {
            for (int k = 0;k < Nomegas;k++) {
                p2[count] = i * i + j * j + k * k;
                count++;
            }
        }
    }
    std::sort(p2, p2 + Nomegas * Nomegas * Nomegas);
    int p2s[Nomegas + 1];
    count = 0;
    p2s[0] = 0;
    int old = -1;
    while (count < Nomegas + 1) {
        for (int i = 0;i < Nomegas * Nomegas * Nomegas;i++) {
            if (p2[i] > old) {
                p2s[count] = p2[i];
                count++;
                old = p2[i];
            }
        }
    }


    fprintf(outfile, " \n\n# omega  cotd  M11\n");
    for (int j = Njack - 1;j < Njack;j++) {
        double x[Nvar];
        x[0] = 1;
        x[1] = Mpi[j];
        x[2] = L;
        x[3] = Mrho[j];
        x[4] = grhopipi[j];
        x[5] = a[j];
        // first we plot cotd*k/m
        if (j == Njack - 1) {
            for (int i = 0;i < 100 /* 1500 */;i++) {
                x[0] = x[1] * (2. + 30.5 * i / 1500.0);
                double k = sqrt(x[0] * x[0] / 4 - x[1] * x[1]);
                double cotd = compute_cotd(Nvar, x);
                double q = k * (L * (x[5] / 197.326963) / (2. * pi_greco));
                // double r = w_js(0, 0, q) - w_js(2, 0, q) - (3. / sqrt(6)) * (w_js(2, -2, q) + w_js(2, 2, q));
                double r = Z00(q) / (pow(pi_greco, 3. / 2.) * q);
                fprintf(outfile, "%g    %g  %g\n", x[0], cotd, r);
            }

            fprintf(outfile, "\n\n #QC2 fit in [%d,%d] chi2=%.5g  %.5g\n", 0, 0, 0.0, 0.0);

        }
        double p = 2 * M_PI / ((double)L * (x[5] / 197.326963));

        for (int i = 0; i < Nomegas;i++) {
            double xmin = 2 * sqrt(x[1] * x[1] + p2s[i] * p * p) + 1e-8;
            double xmax = 2 * sqrt(x[1] * x[1] + p2s[i + 1] * p * p) - 1e-8;
            // printf("%d   %d\n", p2s[i], p2s[i + 1]);
            double oj = rtsafe(omega_QC2, 0/*n*/, 6, x, 0/* Npar */, nullptr/* Pkcot */, 0/*ivar*/, 0. /*input*/, xmin, xmax, 1e-5, 100, 1e-4);
            int seed = rand();
            omega[i] = fake_sampling(resampling, oj, oj / 1e+6, Njack, seed);
            // if (j==Njack-1) printf("omega[%d]=%f  x=%f  %f\n",i, oj, xmin,xmax);
            x[0] = omega[i][j];
            double A = matrix_element_nuA2(Nvar, x);
            nuA2[i] = fake_sampling(resampling, A, A / 1e+6, Njack, seed + 1);
            count++;
            printf("omega nuA2[%d]=%g  %g\n", i, omega[i][Njack - 1], nuA2[i][Njack - 1]);
            fprintf(outfile, "%g  %g\n", omega[i][Njack - 1], error_jackboot(resampling, Njack, omega[i]));
        }
    }

    for (int j = 0;j < Njack;j++) {
        double x[Nvar];
        x[0] = 1;
        x[1] = Mpi[j];
        x[2] = L;
        x[3] = Mrho[j];
        x[4] = grhopipi[j];
        x[5] = a[j];

        VinfL[j] = compute_V_infL(Nvar, x, T);
        VL[j] = compute_VL(Nomegas, omega, nuA2, a, T, j);
        for (int t = 0;t < T;t++) {
            DVt[t][j] = (VinfL[j][t] - VL[j][t]) * pow(a[j] / 197.326963, 3);
        }
    }



    fprintf(outfile, " \n\n# t  inf  L\n");
    for (int t = 1;t < T;t++) {
        fprintf(outfile, "%d   %g   %g\n", t, VinfL[Njack - 1][t], VL[Njack - 1][t]);
    }
    fprintf(outfile, "\n\n #%s_DVt fit in [%d,%d] chi2=%.5g  %.5g\n", description, 0, 0, 0.0, 0.0);
    fprintf(outfile, "%.5g  %.5g\n", 0.0, 0.0);
    printf("--------------integrate GS------------------\n");
    ///// setup integration
    int Maxiter = 1e+6;
    double epsrel = 1e-6;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(Maxiter);
    double* GS;
    double GS_mean;
    for (int j = Njack - 1;j < Njack;j++) {
        // P[Nvar] = t * (a / 197.326963);
        double* P = (double*)malloc(sizeof(double) * (Nvar + 1 + Nomegas * 2));
        P[0] = 1;
        P[1] = Mpi[j];
        P[2] = L;
        P[3] = Mrho[j];
        P[4] = grhopipi[j];
        P[5] = a[j];
        P[6] = Nomegas;
        for (int i = 0; i < Nomegas;i++) {
            P[Nvar + 1 + i] = omega[i][j];
            P[Nvar + 1 + Nomegas + i] = nuA2[i][j];
        }
        double result, error;
        gsl_function F;
        if (window == 0)
            F.function = &integrand_DV;
        else if (window == 1)
            F.function = &integrand_DV_full;
        else if (window == 2) {
            F.function = &integrand_DV_M;
        }
        F.params = (void*)P;

        // gsl_integration_qagiu(&F, 0, 0, epsrel, Maxiter, w, &result, &error);
        gsl_integration_qags(&F, 0, 0.2, 0, epsrel, Maxiter, w, &result, &error);
        printf("GS=%g +- %g\n", result, error);
        GS_mean = result;
        free(P);
    }
    GS = fake_sampling(resampling, GS_mean, GS_mean / 1e+3, Njack, rand());
    gsl_integration_workspace_free(w);


    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < T / 2; t++) {
        fprintf(outfile, "%d   1   1\t", t);
    }
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 0, T / 2, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", GS[Njack - 1], error_jackboot(resampling, Njack, GS));


    return GS;
}



///////////////

double rhs_amu_log_a4(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double Mpi = x[2];
    double Mpiphys = x[3];

    int ilog = x[4];
    int ia4 = x[5];
    int iMpi = x[6];
    int idc = 3;
    double w02 = 0.17383 * 0.17383;


    if (n == 0)      r = P[0] + a2 * P[1];
    else if (n == 1) r = P[0] + a2 * P[2];

    if (ilog == 0) {
        // nothing    
    }
    else if (ilog == 1) {
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(a2 / w02), 3));
    }
    else if (ilog == 2) {
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(a2 / w02), 3));
    }
    else if (ilog == 3) {
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(a2 / w02), 3));
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(a2 / w02), 3));
    }
    else if (ilog == 4) {
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(a2 / w02), 2));
    }
    else if (ilog == 5) {
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(a2 / w02), 2));
    }
    else if (ilog == 6) {
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(a2 / w02), 2));
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(a2 / w02), 2));
    }
    else if (ilog == 7) {
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(a2 / w02), 1));
    }
    else if (ilog == 8) {
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(a2 / w02), 1));
    }
    else if (ilog == 9) {
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(a2 / w02), 1));
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(a2 / w02), 1));
    }


    if (ia4 == 0) {

    }
    if (ia4 == 1) {
        if (n == 0)      r += +a2 * a2 * P[idc];
        idc++;
    }
    if (ia4 == 2) {
        if (n == 1)      r += +a2 * a2 * P[idc];
        idc++;
    }
    if (ia4 == 3) {
        if (n == 0)           r += +a2 * a2 * P[idc];
        else if (n == 1)      r += +a2 * a2 * P[idc + 1];
        idc += 2;
    }

    switch (iMpi) {
    case 0:
        break;
    case 1:
        if (n == 0)  r += P[idc] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
        idc++;
        break;
    case 2:
        if (n == 1)  r += P[idc] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
        idc++;
        break;
    case 3:
        if (n == 0)  r += P[idc] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
        if (n == 1)  r += P[idc] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
        idc++;
        break;

    default:
        break;
    }

    return r;
}


double linear_fit_mu_correction(int n, int Nvar, double* x, int Npar, double* P) {
    double a2 = x[0];
    double r;
    if (n == 0) r = P[0];
    if (n == 1) r = P[1];
    if (Npar > 2) {
        if (n == 0) r += a2 * P[2];
        if (n == 1) r += a2 * P[3];
    }
    return r;
}

double exp_MpiL(int n, int Nvar, double* x, int Npar, double* P) {
    double MpiL = x[0];
    double r;
    if (n == 0) r = P[0] + P[2] * exp(-MpiL);
    if (n == 1) r = P[1] + P[3] * exp(-MpiL);
    return r;
}
double const_A(int n, int Nvar, double* x, int Npar, double* P) {

    return P[n];
}



double rhs_amu_RF(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double Mpi = x[2];
    double Mpiphys = x[3];

    int ilog = x[4];
    int ia4 = x[5];
    int iMpi = x[6];
    double iw = x[7];
    // int plog = x[9];
    int idc = Npar - 1;
    double w02 = 0.17383 * 0.17383 / (iw);


    if (n == 0)      r = P[0] + a2 * P[1];
    else if (n == 1) r = P[0] + a2 * P[2];
    switch (ilog) {
    case 0:
        break;
    case 1:
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(w02 / a2), 3));
        break;
    case 2:
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(w02 / a2), 3));
        break;
    case 3:
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(w02 / a2), 3));
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(w02 / a2), 3));
        break;
    case 4:
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(w02 / a2), 2));
        break;
    case 5:
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(w02 / a2), 2));
        break;
    case 6:
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(w02 / a2), 2));
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(w02 / a2), 2));
        break;
    case 7:
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(w02 / a2), 1));
        break;
    case 8:
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(w02 / a2), 1));
        break;
    case 9:
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(w02 / a2), 1));
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(w02 / a2), 1));
        break;
    case 10:
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(w02 / a2), -0.2));
        break;
    case 11:
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(w02 / a2), -0.2));
        break;
    case 12:
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(w02 / a2), -0.2));
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(w02 / a2), -0.2));
        break;

    case 13:
        if (n == 0)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), 3));
        idc--;
        break;
    case 14:
        if (n == 1)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), 3));
        idc--;
        break;
    case 15:
        if (n == 0)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), 3));
        if (n == 1)      r += a2 * P[idc - 1] * (1. / pow(log(w02 / a2), 3));
        idc--;idc--;
        break;
    case 16:
        if (n == 0)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), 2));
        idc--;
        break;
    case 17:
        if (n == 1)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), 2));
        idc--;
        break;
    case 18:
        if (n == 0)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), 2));
        if (n == 1)      r += a2 * P[idc - 1] * (1. / pow(log(w02 / a2), 2));
        idc--;idc--;
        break;
    case 19:
        if (n == 0)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), 1));
        idc--;
        break;
    case 20:
        if (n == 1)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), 1));
        idc--;
        break;
    case 21:
        if (n == 0)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), 1));
        if (n == 1)      r += a2 * P[idc - 1] * (1. / pow(log(w02 / a2), 1));
        idc--;idc--;
        break;
    case 22:
        if (n == 0)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), -0.2));
        idc--;
        break;
    case 23:
        if (n == 1)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), -0.2));
        idc--;
        break;
    case 24:
        if (n == 0)      r += a2 * P[idc] * (1. / pow(log(w02 / a2), -0.2));
        if (n == 1)      r += a2 * P[idc - 1] * (1. / pow(log(w02 / a2), -0.2));
        idc--;idc--;
        break;
    default:
        break;
    }





    if (ia4 == 0) {
    }
    if (ia4 == 1) {
        if (n == 0)      r += +a2 * a2 * P[idc];
        idc--;
    }
    if (ia4 == 2) {
        if (n == 1)      r += +a2 * a2 * P[idc];
        idc--;
    }
    if (ia4 == 3) {
        if (n == 0)           r += +a2 * a2 * P[idc];
        else if (n == 1)      r += +a2 * a2 * P[idc - 1];
        idc -= 2;
    }
    if (ia4 == 4) {
        if (n == 0)           r += +a2 * a2 * P[idc];
        else if (n == 1)      r += +a2 * a2 * P[idc];
        idc--;
    }

    switch (iMpi) {
    case 0:
        break;
    case 1:
        if (n == 0)  r += P[idc] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
        idc--;
        break;
    case 2:
        if (n == 1)  r += P[idc] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
        idc--;
        break;
    case 3:
        if (n == 0)  r += P[idc] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
        if (n == 1)  r += P[idc] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
        idc--;
        break;

    default:
        break;
    }

    return r;
}

double rhs_amu_separate(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    int il = x[4];
    int ia2 = x[5];
    int iw = x[7];
    double r = 0;
    double w02 = 0.17383 * 0.17383;
    if (iw == 1) w02 *= 9;
    r = P[0];

    if (ia2 == 1) {
        switch (il) {
        case 0:
            r += a2 * P[1];
            break;
        case 1:
            r += a2 * P[1] * (1. / pow(log(w02 / a2), 1));
            break;
        case 2:
            r += a2 * P[1] * (1. / pow(log(w02 / a2), 2));
            break;
        case 3:
            r += a2 * P[1] * (1. / pow(log(w02 / a2), 3));
            break;
        default:
            break;
        }
    }

    return r;
}


double rhs_amu_a4(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    int il = x[4];
    int ia2 = x[5];
    int ial = x[6];
    int iw = x[7];
    double r = 0;
    int iR = x[6];
    double w02 = 0.17383 * 0.17383;
    if (iw == 1) w02 *= 9;
    double OS = 0, TM = 0;

    int in = il / 4;
    int im = il % 4;

    int an = ial / 4;
    int am = ial % 4;


    OS += P[0] + a2 * P[1] * (1. / pow(log(w02 / a2), in));
    TM += P[0] + a2 * P[2] * (1. / pow(log(w02 / a2), im));
    // if (an == 0) OS += a2 * a2 * P[3];
    // if (am == 0) TM += a2 * a2 * P[3];
    // if (ia2 == 0 an > 0) OS += (1. / pow(log(w02 / a2), an)) * a2 * P[3];
    // if (am > 0) TM += (1. / pow(log(w02 / a2), am)) * a2 * P[3];
// printf("ial=%d  ia2=%d  in=%d  im=%d   an =%d   am=%d  iw=%d    \n",ial,ia2,in,im,an,am,iw);
    if (ia2 == 0 && ial == 0) OS += a2 * a2 * P[3];
    if (ia2 == 1 && ial == 0) TM += a2 * a2 * P[3];
    if (ia2 == 0 && ial > 0)  OS += (1. / pow(log(w02 / a2), an)) * a2 * P[3];
    if (ia2 == 1 && ial > 0)  TM += (1. / pow(log(w02 / a2), am)) * a2 * P[3];
    if (ia2 == 2) {
        if (ial == 0) {
            OS += a2 * a2 * P[3];
            TM += a2 * a2 * P[4];
        }
        else {
            OS += (1. / pow(log(w02 / a2), an)) * a2 * P[3];
            TM += (1. / pow(log(w02 / a2), am)) * a2 * P[4];
        }

    }


    if (n == 0) r = OS;
    if (n == 1) r = TM;



    return r;
}

double rhs_amu_a4_charm(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    int il = x[4];
    int ia2 = x[5];
    int ial = x[6];
    int iw = x[7];
    double r = 0;
    int iR = x[6];
    double w02 = 0.17383 * 0.17383;
    if (iw == 1) w02 *= 9;
    double OS = 0, TM = 0;

    int in = il / 4;
    int im = il % 4;

    int an = ial / 4;
    int am = ial % 4;


    OS += P[0] + a2 * P[1] * (1. / pow(log(w02 / a2), in));
    TM += P[0] + a2 * P[2] * (1. / pow(log(w02 / a2), im));
    // if (an == 0) OS += a2 * a2 * P[3];
    // if (am == 0) TM += a2 * a2 * P[3];
    // if (ia2 == 0 an > 0) OS += (1. / pow(log(w02 / a2), an)) * a2 * P[3];
    // if (am > 0) TM += (1. / pow(log(w02 / a2), am)) * a2 * P[3];

    if (ia2 == 0 && ial == 0) OS += a2 * a2 * P[3];
    if (ia2 == 1 && ial == 0) TM += a2 * a2 * P[3];
    if (ia2 == 0 && ial > 0) OS += (1. / pow(log(w02 / a2), an)) * a2 * P[3];
    if (ia2 == 1 && ial > 0) TM += (1. / pow(log(w02 / a2), am)) * a2 * P[3];
    if (ia2 == 2) {
        if (ial == 0) {
            OS += a2 * a2 * P[3];
            TM += a2 * a2 * P[4];
        }
        else {
            OS += (1. / pow(log(w02 / a2), an)) * a2 * P[3];
            TM += (1. / pow(log(w02 / a2), am)) * a2 * P[4];
        }

    }


    if (n == 0) r = OS;
    if (n == 1) r = TM;



    return r;
}



double rhs_amu_diff_ratio_charm(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    int il = x[4];
    int ia2 = x[5];
    int iw = x[7];
    double r = 0;
    int iR = x[6];
    double w02 = 0.17383 * 0.17383;
    if (iw == 1) w02 *= 9;
    double diff = 0, ratio = 1;

    int in = il / 4;
    int im = il % 4;

    if (ia2 == 0) {
        diff += a2 * P[0] * (1. / pow(log(w02 / a2), in));
        ratio += a2 * P[1] * (1. / pow(log(w02 / a2), im));
    }
    else if (ia2 == 1) {
        diff += a2 * P[0] * (1. / pow(log(w02 / a2), in)) + a2 * a2 * P[2];
        ratio += a2 * P[1] * (1. / pow(log(w02 / a2), im)) + a2 * a2 * P[3];
    }
    else if (ia2 == 2) {
        if (in == 0)diff += +a2 * P[0];
        else diff += +a2 * P[0] + a2 * P[2] * (1. / pow(log(w02 / a2), in));
        if (im == 0)ratio += +a2 * P[1];
        else {
            if (in == 0) ratio += +a2 * P[1] + a2 * P[2] * (1. / pow(log(w02 / a2), im));
            else ratio += +a2 * P[1] + a2 * P[3] * (1. / pow(log(w02 / a2), im));
        }
    }


    if (n == 0) r = diff;
    if (n == 1) r = ratio;

    if (iR == 1) {
        if (n == 0) r = -diff;
    }

    return r;
}

double rhs_amu_diff_ratio(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    int il = x[4];
    int ia2 = x[5];
    int iw = x[7];
    double r = 0;
    int iR = x[6];
    double w02 = 0.17383 * 0.17383;
    if (iw == 1) w02 *= 9;
    double diff = 0, ratio = 1;

    // int in = il / 4;
    // int im = il % 4;

    if (ia2 == 0) {
        diff += a2 * P[0] * (1. / pow(log(w02 / a2), il));
        ratio += a2 * P[1] * (1. / pow(log(w02 / a2), il));
    }
    else if (ia2 == 1) {
        diff += a2 * P[0] * (1. / pow(log(w02 / a2), il)) + a2 * a2 * P[2];
        ratio += a2 * P[1] * (1. / pow(log(w02 / a2), il)) + a2 * a2 * P[3];
    }
    else if (ia2 == 2) {
        diff += +a2 * P[0] + a2 * P[2] * (1. / pow(log(w02 / a2), il));
        ratio += +a2 * P[1] + a2 * P[3] * (1. / pow(log(w02 / a2), il));
    }


    if (n == 0) r = diff;
    if (n == 1) r = ratio;

    if (iR == 1) {
        if (n == 0) r = -diff;
    }

    return r;
}

double rhs_amu_cut(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    int il = x[4];
    int ia2 = x[5];
    int iw = x[7];
    double r = 0;
    double w02 = 0.17383 * 0.17383;
    if (iw == 1) w02 *= 9;

    r = P[0];
    if (ia2 == 1) {
        switch (il) {
        case 0:
            r += a2 * P[1 + n];
            break;
        case 1:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 0));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 1));
            break;
        case 2:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 0));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 2));
            break;
        case 3:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 0));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 3));
            break;
        case 4:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 1));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 0));
            break;
        case 5:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 1));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 1));
            break;
        case 6:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 1));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 2));
            break;
        case 7:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 1));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 3));
            break;
        case 8:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 2));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 0));
            break;
        case 9:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 2));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 1));
            break;
        case 10:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 2));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 2));
            break;
        case 11:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 2));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 3));
            break;
        case 12:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 3));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 0));
            break;
        case 13:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 3));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 1));
            break;
        case 14:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 3));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 2));
            break;
        case 15:
            if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 3));
            if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), 3));
            break;
        default:
            break;
        }
    }



    return r;
}


double rhs_amu_cut_charm(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    int il = x[4];
    int ia2 = x[5];
    int iw = x[7];
    int ial = x[6];
    double r = 0;
    double w02 = 0.17383 * 0.17383;
    if (iw == 1) w02 *= 9;

    int in = il / 4;
    int im = il % 4;

    int an = ial / 4;
    int am = ial % 4;

    r = P[0];
    if (ia2 == 1) {
        if (n == 0) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), in));
        if (n == 1) r += a2 * P[1 + n] * (1. / pow(log(w02 / a2), im));
    }
    if (ia2 == 2) {

        r += a2 * P[1 + n];
        if (n == 0) {
            if (in == 0) r += a2 * P[3 + n] * a2;
            else   r += a2 * P[3 + n] * (1. / pow(log(w02 / a2), in));
        }
        if (n == 1) {
            if (im == 0) r += a2 * P[3 + n] * a2;
            else r += a2 * P[3 + n] * (1. / pow(log(w02 / a2), im));
        }
    }



    return r;
}

double rhs_amu_pade(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    int il = x[4];
    int ia2 = x[5];
    int iw = x[7];
    int who_pade = x[6];
    double r = 0;
    double w02 = 0.17383 * 0.17383;
    if (iw == 1) w02 *= 9;

    double logw = (1. / pow(log(w02 / a2), il));

    if (who_pade == 0) {
        if (n == 0) r = (P[0] + P[2] * a2 * logw) / (1 + P[3] * a2);
        else if (n == 1) r = P[0] + P[1] * a2 * logw;
    }
    else if (who_pade == 1) {
        if (n == 0) r = P[0] + P[1] * a2 * logw;
        else if (n == 1) r = (P[0] + P[2] * a2 * logw) / (1 + P[3] * a2);
    }



    return r;
}


double rhs_amu_FVE_RF(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    double r = rhs_amu_RF(n, Nvar, x, Npar, P);
    if (n == 0)      r += -(10.0 / 9.0) * GS * (P[3] * a2);
    else if (n == 1) r += -(10.0 / 9.0) * GS * (P[4] * a2);
    return r;
}

double rhs_amu_diff_RF(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    double Mpi = x[2];
    double Mpiphys = x[3];

    int ilog = x[4];
    int ia4 = x[5];
    int iMpi = x[6];
    double iw = x[7];
    // int plog = x[9];
    int idc = Npar - 1;
    double w02 = 0.17383 * 0.17383 / (iw);

    double r = a2 * P[0] - (10.0 / 9.0) * GS * (P[1] * a2);
    switch (ilog) {
    case 0:
        break;
    case 1:
        r += a2 * P[0] * (-1 + 1. / pow(log(w02 / a2), 3));
        break;
    case 2:
        r += a2 * P[0] * (-1 + 1. / pow(log(w02 / a2), 2));
        break;
    case 3:
        r += a2 * P[0] * (-1 + 1. / pow(log(w02 / a2), 1));
        break;
    case 4:
        r += a2 * P[0] * (-1 + 1. / pow(log(w02 / a2), -0.2));
        break;
    case 5:
        r += a2 * P[idc] * (1. / pow(log(w02 / a2), 3));
        idc--;
        break;
    case 6:
        r += a2 * P[idc] * (1. / pow(log(w02 / a2), 2));
        idc--;
        break;
    case 7:
        r += a2 * P[idc] * (1. / pow(log(w02 / a2), 1));
        idc--;
        break;
    case 8:
        r += a2 * P[idc] * (1. / pow(log(w02 / a2), -0.2));
        idc--;
        break;
    default:
        break;
    }

    if (ia4 == 0) {
    }
    if (ia4 == 1) {
        r += +a2 * a2 * P[idc];
        idc--;
    }

    switch (iMpi) {
    case 0:
        break;
    case 1:
        r += P[idc] * (Mpi - Mpiphys);// all in MeV
        idc--;
        break;
    default:
        break;
    }


    return r;
}

double rhs_amu_diff(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    double Mpi = x[2];
    double Mpiphys = x[3];

    int ilog = x[4];
    int ia4 = x[5];
    int iMpi = x[6];
    double iw = x[7];
    // int plog = x[9];
    int idc = Npar - 1;
    double w02 = 0.17383 * 0.17383 / (iw);

    double r = a2 * P[0];

    if (ia4 == 0) {
    }
    if (ia4 == 1) {
        r += +a2 * a2 * P[1];
        idc--;
    }




    return r;
}


double rhs_amu_ratio(int n, int Nvar, double* x, int Npar, double* P) {
    double GS = x[1];
    double a2 = x[0];
    double Mpi = x[2];
    double Mpiphys = x[3];

    int ilog = x[4];
    int ia4 = x[5];
    int iMpi = x[6];
    double iw = x[7];
    // int plog = x[9];
    int idc = Npar - 1;
    double w02 = 0.17383 * 0.17383 / (iw);

    double r = 1. + a2 * P[0];

    if (ia4 == 0) {
    }
    if (ia4 == 1) {
        r += +a2 * a2 * P[1];
        idc--;
    }




    return r;
}


double rhs_amu_common(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a = x[0];
    if (n == 0)      r = P[0] + a * P[1];
    else if (n == 1) r = P[0] + a * P[2];
    return r;
}



double rhs_amu_common_a2_FVE(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double GS = x[1];
    double Mpi = x[2];
    double Mpiphys = x[3];
    if (n == 0)      r = P[0] + a2 * P[1] - (10.0 / 9.0) * GS * (P[3] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
    else if (n == 1) r = P[0] + a2 * P[2] - (10.0 / 9.0) * GS * (P[4] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
    return r;
}


double rhs_amu_common_a2_FVE_log_eq(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double GS = x[1];
    double Mpi = x[2];
    double Mpiphys = x[3];
    double w02 = 0.17383 * 0.17383;
    if (n == 0)      r = P[0] + a2 * P[1] / pow(log(a2 / w02), 3) - (10.0 / 9.0) * GS * (P[3] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
    else if (n == 1) r = P[0] + a2 * P[2] - (10.0 / 9.0) * GS * (P[4] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
    return r;
}
double rhs_amu_common_a2_FVE_log_op(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double GS = x[1];
    double Mpi = x[2];
    double Mpiphys = x[3];
    double w02 = 0.17383 * 0.17383;
    if (n == 0)      r = P[0] + a2 * P[1] - (10.0 / 9.0) * GS * (P[3] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
    else if (n == 1) r = P[0] + a2 * P[2] / pow(log(a2 / w02), 3) - (10.0 / 9.0) * GS * (P[4] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
    return r;
}

double rhs_amu_common_a2_FVE_log_eq_op(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double GS = x[1];
    double Mpi = x[2];
    double Mpiphys = x[3];
    double w02 = 0.17383 * 0.17383;
    if (n == 0)      r = P[0] + a2 * P[1] / pow(log(a2 / w02), 3) - (10.0 / 9.0) * GS * (P[3] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
    else if (n == 1) r = P[0] + a2 * P[2] / pow(log(a2 / w02), 3) - (10.0 / 9.0) * GS * (P[4] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));

    return r;
}



double rhs_amu_common_a2_FVE_log_a4(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double GS = x[1];
    double Mpi = x[2];
    double Mpiphys = x[3];
    double w02 = 0.17383 * 0.17383;
    int ilog = x[4];
    int a4 = x[5];
    if (n == 0)      r = P[0] + a2 * P[1] - (10.0 / 9.0) * GS * (P[3] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
    else if (n == 1) r = P[0] + a2 * P[2] - (10.0 / 9.0) * GS * (P[4] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));

    if (ilog == 0) {
        // nothing    
    }
    else if (ilog == 1) {
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(a2 / w02), 3));
    }
    else if (ilog == 2) {
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(a2 / w02), 3));
    }
    else if (ilog == 3) {
        if (n == 0)      r += a2 * P[1] * (-1 + 1. / pow(log(a2 / w02), 3));
        if (n == 1)      r += a2 * P[2] * (-1 + 1. / pow(log(a2 / w02), 3));
    }

    if (a4 == 0) {

    }
    if (a4 == 1) {
        if (n == 0)      r += +a2 * a2 * P[6];
    }
    if (a4 == 2) {
        if (n == 1)      r += +a2 * a2 * P[6];
    }



    return r;
}

double rhs_amu_common_a2_FVE_log_eq_a4_eq(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double GS = x[1];
    double Mpi = x[2];
    double Mpiphys = x[3];
    double w02 = 0.17383 * 0.17383;
    if (n == 0)      r = P[0] + a2 * P[1] / pow(log(a2 / w02), 3) - (10.0 / 9.0) * GS * (P[3] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963)) + a2 * a2 * P[6];
    else if (n == 1) r = P[0] + a2 * P[2] - (10.0 / 9.0) * GS * (P[4] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));

    return r;
}
double rhs_amu_common_a2_FVE_log_op_a4_eq(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double GS = x[1];
    double Mpi = x[2];
    double Mpiphys = x[3];
    double w02 = 0.17383 * 0.17383;
    if (n == 0)      r = P[0] + a2 * P[1] - (10.0 / 9.0) * GS * (P[3] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963)) + a2 * a2 * P[6];
    else if (n == 1) r = P[0] + a2 * P[2] / pow(log(a2 / w02), 3) - (10.0 / 9.0) * GS * (P[4] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));

    return r;
}

double rhs_amu_common_a2_FVE_log_eq_op_a4_eq(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double GS = x[1];
    double Mpi = x[2];
    double Mpiphys = x[3];
    double w02 = 0.17383 * 0.17383;
    if (n == 0)      r = P[0] + a2 * P[1] / pow(log(a2 / w02), 3) - (10.0 / 9.0) * GS * (P[3] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
    else if (n == 1) r = P[0] + a2 * P[2] / pow(log(a2 / w02), 3) - (10.0 / 9.0) * GS * (P[4] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));

    return r;
}


double rhs_amu_common_a2_FVE_a4_eq(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double GS = x[1];
    double Mpi = x[2];
    double Mpiphys = x[3];
    if (n == 0)      r = P[0] + a2 * P[1] - (10.0 / 9.0) * GS * (P[3] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963)) + a2 * a2 * P[6];
    else if (n == 1) r = P[0] + a2 * P[2] - (10.0 / 9.0) * GS * (P[4] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));

    return r;
}

double rhs_amu_common_a2_FVE_a4_op(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double GS = x[1];
    double Mpi = x[2];
    double Mpiphys = x[3];
    if (n == 0)      r = P[0] + a2 * P[1] - (10.0 / 9.0) * GS * (P[3] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963));
    else if (n == 1) r = P[0] + a2 * P[2] - (10.0 / 9.0) * GS * (P[4] * a2) + P[5] * (Mpi - Mpiphys * (sqrt(a2) / 197.326963)) + a2 * a2 * P[6];

    return r;
}


double rhs_amu_common_a4(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    if (n == 0)      r = P[0] + a2 * P[1] + a2 * a2 * P[3];
    else if (n == 1) r = P[0] + a2 * P[2] + a2 * a2 * P[4];

    return r;
}
double rhs_amu_common_a4_n0(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    if (n == 0)      r = P[0] + a2 * P[1] + a2 * a2 * P[3];
    else if (n == 1) r = P[0] + a2 * P[2];

    return r;
}
double rhs_amu_common_a4_n1(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    if (n == 0)      r = P[0] + a2 * P[1];
    else if (n == 1) r = P[0] + a2 * P[2] + a2 * a2 * P[3];

    return r;
}
double rhs_amu_common_log_a4_n0(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double w02 = 0.030276;
    if (n == 0)      r = P[0] + (a2 / pow(log(a2 / w02), 3)) * P[1] + a2 * a2 * P[3];
    else if (n == 1) r = P[0] + (a2 / pow(log(a2 / w02), 3)) * P[2];
    return r;
}
double rhs_amu_common_log_a4_n1(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    double w02 = 0.030276;
    if (n == 0)      r = P[0] + (a2 / pow(log(a2 / w02), 3)) * P[1];
    else if (n == 1) r = P[0] + (a2 / pow(log(a2 / w02), 3)) * P[2] + a2 * a2 * P[3];
    return r;
}

template<int ieq, int iop>
double lhs_amu_common(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    double GS = gjack.en[e].jack[58][j];
    if (n == 0)        r = gjack.en[e].jack[ieq][j];
    else if (n == 1)   r = gjack.en[e].jack[iop][j];
    return r;
}

double lhs_amu_common_GS(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    double GS = gjack.en[e].jack[58][j];
    int ieq = fit_info.corr_id[0];
    int iop = fit_info.corr_id[1];
    if (n == 0)        r = gjack.en[e].jack[ieq][j] + (10.0 / 9.0) * GS;
    else if (n == 1)   r = gjack.en[e].jack[iop][j] + (10.0 / 9.0) * GS;
    return r;
}



double lhs_amu_common_GS_diff(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    double GS = gjack.en[e].jack[58][j];
    int ieq = fit_info.corr_id[0];
    int iop = fit_info.corr_id[1];
    r = gjack.en[e].jack[ieq][j] + (10.0 / 9.0) * GS - gjack.en[e].jack[iop][j] - (10.0 / 9.0) * GS;
    // else if (n == 1)   r = gjack.en[e].jack[iop][j] + (10.0 / 9.0) * GS;
    return r;
}

double lhs_amu_diff(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    double GS = gjack.en[e].jack[58][j];
    int ieq = fit_info.corr_id[0];
    int iop = fit_info.corr_id[1];
    r = gjack.en[e].jack[iop][j] - gjack.en[e].jack[ieq][j];
    // else if (n == 1)   r = gjack.en[e].jack[iop][j] + (10.0 / 9.0) * GS;
    return r;
}

double lhs_amu_ratio(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    double GS = gjack.en[e].jack[58][j];
    int ieq = fit_info.corr_id[0];
    int iop = fit_info.corr_id[1];
    r = gjack.en[e].jack[iop][j] / gjack.en[e].jack[ieq][j];
    // else if (n == 1)   r = gjack.en[e].jack[iop][j] + (10.0 / 9.0) * GS;
    return r;
}


double lhs_amu(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    if (fit_info.corr_id.size() != 2) { printf("error lhs_amu called without ids\n"); exit(1); }
    if (n == 0)        r = gjack.en[e].jack[fit_info.corr_id[0]][j];
    else if (n == 1)   r = gjack.en[e].jack[fit_info.corr_id[1]][j];
    return r;
}
double lhs_Acharm(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;

    return gjack.en[e].jack[fit_info.corr_id[n]][j];
}


double lhs_amu_separate(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    r = gjack.en[e].jack[fit_info.corr_id[0]][j];

    return r;
}


double lhs_amu_diff_ratio(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    if (fit_info.corr_id.size() != 2) { printf("error lhs_amu called without ids\n"); exit(1); }
    if (n == 0)        r = gjack.en[e].jack[fit_info.corr_id[0]][j] - gjack.en[e].jack[fit_info.corr_id[1]][j];
    else if (n == 1)   r = gjack.en[e].jack[fit_info.corr_id[0]][j] / gjack.en[e].jack[fit_info.corr_id[1]][j];
    return r;
}
double lhs_mu_corrections(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    double MpiMeV_phys = fit_info.x[1][e][j];// in principle e--> e+n*fit_info.myen.size(), but htey should be the same
    double a_MeV = sqrt(fit_info.x[0][e][j]) / hbarc;
    if (n == 0)        r = (gjack.en[e].jack[130][j] - gjack.en[e].jack[42][j]) / (MpiMeV_phys - gjack.en[e].jack[1][j] / a_MeV);
    else if (n == 1)   r = (gjack.en[e].jack[132][j] - gjack.en[e].jack[43][j]) / (MpiMeV_phys - gjack.en[e].jack[1][j] / a_MeV);
    return r;
}

double lhs_SD_mu_corrections(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    double MpiMeV_phys = fit_info.x[1][e][j];// in principle e--> e+n*fit_info.myen.size(), but htey should be the same
    double a_MeV = sqrt(fit_info.x[0][e][j]) / hbarc;
    if (n == 0)        r = (gjack.en[e].jack[134][j] - gjack.en[e].jack[25][j]) / (MpiMeV_phys - gjack.en[e].jack[1][j] / a_MeV);
    else if (n == 1)   r = (gjack.en[e].jack[135][j] - gjack.en[e].jack[26][j]) / (MpiMeV_phys - gjack.en[e].jack[1][j] / a_MeV);
    return r;
}



template<int ieq, int iop>
double lhs_amu_common_FVE(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    if (n == 0)        r = gjack.en[e].jack[ieq][j] + (10. / 9.) * gjack.en[e].jack[48][j];
    else if (n == 1)   r = gjack.en[e].jack[iop][j] + (10. / 9.) * gjack.en[e].jack[48][j];
    return r;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// print band
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void print_fit_band_amu_W_l(char** argv, data_all gjack, struct fit_type fit_info,
    struct fit_type fit_info_m0, const char* label, const char* dir_name,
    struct fit_result fit_out, struct fit_result fit_out_m0, int var, int en, double h, std::vector<double> xval,
    double lhs_fun(int, int, int, data_all, struct fit_type)) {

    std::vector<int> myen = fit_info.myen;
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar;
    int Njack = gjack.en[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double** tif_m0 = swap_indices(fit_info_m0.Npar, Njack, fit_out_m0.P);
    double* tmpx = (double*)malloc(sizeof(double) * Nvar);

    double min, max;
    if (fit_info.band_range.size() != 2) {
        min = fit_info.x[var][0][Njack - 1];
        max = fit_info.x[var][0][Njack - 1];
        for (int e = 1; e < myen.size() * N;e++) {
            if (fit_info.x[var][e][Njack - 1] < min)
                min = fit_info.x[var][e][Njack - 1];
            if (fit_info.x[var][e][Njack - 1] > max)
                max = fit_info.x[var][e][Njack - 1];
        }
        if (min == max)  printf("no variation in this direction\n");
        min -= (max - min) / 10.0;
        max += (max - min) / 10.0;

    }
    else {
        min = fit_info.band_range[0];
        max = fit_info.band_range[1];
    }

    max += h / 2.0;
    /// y
    int* ens = (int*)malloc(sizeof(int) * fit_info.N);// we need to init en and en_tot to allocate the other 
    for (int e = 0;e < fit_info.N; e++) { ens[e] = myen.size(); }
    int en_tot = 0;      for (int n = 0;n < N;n++) { en_tot += ens[n]; }// total data to fit
    double*** y = double_malloc_3(Njack, en_tot, 2);// 2 is mean/error
    int count = 0;
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < ens[n];e++) {
            double* tmpj = (double*)malloc(sizeof(double) * Njack);
            for (int j = 0;j < Njack;j++) {
                tmpj[j] = lhs_fun(n, myen[e], j, gjack, fit_info);
                // if (j == Njack - 1 && n == 0) printf("y =%g   GS=%g  FVE=%g  e=%d  \n", tmpj[j],
                //     (10.0 / 9.0) * fit_info.x[1][e][j],
                //     (10.0 / 9.0) * fit_info.x[1][e][j] * (1 - fit_out.P[3][j] * fit_info.x[0][e][j]), e);
                for (int i = 0; i < Nvar; i++) {
                    tmpx[i] = fit_info.x[i][e + count][j];
                }
                double DGS = fit_info.function(n, Nvar, tmpx, Npar, tif[j]);
                tmpx[1] = 0;
                DGS -= fit_info.function(n, Nvar, tmpx, Npar, tif[j]);
                tmpx[1] = fit_info.x[1][e + count][j];

                double DMpi = fit_info.function(n, Nvar, tmpx, Npar, tif[j]);
                tmpx[2] = 0; tmpx[3] = 0;
                DMpi -= fit_info.function(n, Nvar, tmpx, Npar, tif[j]);

                // double compare = (10.0 / 9.0) * fit_info.x[1][e][j] * (fit_out.P[4][j] * fit_info.x[0][e][j])
                //     - fit_out.P[5][j] * (fit_info.x[2][e][j] - fit_info.x[3][e][j] * (sqrt(fit_info.x[0][e][j]) / 197.326963));

                tmpj[j] -= (DGS + DMpi);
                // if(compare==-(DGS+DMpi)){ printf("error  compare=%.12g    , %.12g\n", compare, -(DGS+DMpi) )  ; exit(1); }
                // if (n == 0)  tmpj[j] += (10.0 / 9.0) * fit_info.x[1][e][j] * (fit_out.P[3][j] * fit_info.x[0][e][j])
                //     - fit_out.P[5][j] * (fit_info.x[2][e][j] - fit_info.x[3][e][j] * (sqrt(fit_info.x[0][e][j]) / 197.326963));
                // else if (n == 1)  tmpj[j] += (10.0 / 9.0) * fit_info.x[1][e][j] * (fit_out.P[4][j] * fit_info.x[0][e][j]) - fit_out.P[5][j] * (fit_info.x[2][e][j] - fit_info.x[3][e][j] * (sqrt(fit_info.x[0][e][j]) / 197.326963));
            }
            double err = error_jackboot(argv[1], Njack, tmpj);
            for (int j = 0;j < Njack;j++) {
                y[j][e + count][0] = tmpj[j];
                y[j][e + count][1] = err;
            }


            free(tmpj);
        }
        count += ens[n];
    }
    mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_ysub_%s.txt", argv[3], label, dir_name);
    FILE* fy = open_file(namefile, "w+");

    for (int n = 0;n < N; n++) {

        mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_n%d_%s.txt", argv[3], label, n, dir_name);
        f = open_file(namefile, "w+");
        double* tmpy = (double*)malloc(sizeof(double) * Njack);
        printf("writing: %s\n", namefile);


        count = 0;
        for (int n = 0;n < N;n++) {
            for (int e = 0;e < ens[n];e++) {
                for (int v = 0; v < Nvar;v++) {
                    fprintf(fy, " %g\t ", fit_info.x[v][e + count][Njack - 1]);
                }
                fprintf(fy, " %g   %g  \t ", y[Njack - 1][e + count][0], y[Njack - 1][e + count][1]);
                fprintf(fy, " %d   \n ", n);
            }
            count += ens[n];
            fprintf(fy, "\n\n");
        }
        // band range
        double pos = min;
        while (pos < max) {
            for (int j = 0;j < Njack;j++) {
                for (int i = 0; i < xval.size(); i++) {
                    tmpx[i] = xval[i];
                }
                for (int i = xval.size(); i < Nvar; i++) {
                    tmpx[i] = fit_info.x[i][en + n * myen.size()][j];
                }
                tmpx[var] = pos;
                tmpy[j] = fit_info.function(n, Nvar, tmpx, Npar, tif[j]);
            }
            fprintf(f, "%g  \t %g  %g\n", pos, tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));
            pos += h;
        }
        free(tmpy);
        fclose(f);
    }
    free(tmpx);
    fclose(fy);
    free_2(Njack, tif);
    free_2(Njack, tif_m0);
    free_3(Njack, en_tot, y);
    free(ens);

}


double rhs_2exp(int n, int Nvar, double* x, int Npar, double* P) {
    double r = 0;
    double xpower = 1;
    double t = x[0];
    double t0 = x[1];
    double a = x[3];

    return P[0] * (P[0] * exp(P[1] * (t * a)) + P[2] * exp(P[3] * (t * a))) / P[0] * exp(P[1] * (t0)) + P[2] * exp(P[3] * (t0));
}

double rhs_1exp(int n, int Nvar, double* x, int Npar, double* P) {
    double r = 0;
    double xpower = 1;
    double t = x[0];
    double t0 = x[1];
    double a = x[3];
    return P[0] * exp(P[1] * (t * a - t0));
}

double rhs_poly(int n, int Nvar, double* x, int Npar, double* P) {
    double r = 0;
    double xpower = 1;
    double t = x[0];
    double t0 = x[1];
    double a = x[3];
    for (int i = 0;i < Npar;i++) {
        r += P[i] * xpower;
        xpower *= (t * a - t0);
    }
    if (isnan(r)) {
        printf("nan found:\n");
        printf("t= %g  t0=%g   a=%g  \n", t, t0, a);
        for (int i = 0;i < Npar;i++) printf("P[%d]=%g\t", i, P[i]);
        printf("Nvar= %d  Npar=%d  \n", Nvar, Npar);
    }
    return r;
}


double** corr_plus_dm(int j, double**** in, int t, struct fit_type fit_info) {
    double** r = double_calloc_2(fit_info.N, 2);
    double dmu = fit_info.guess[1] - fit_info.guess[0];
    int id0 = fit_info.corr_id[0];
    int id1 = fit_info.corr_id[1];
    int ibolla = 42;
    int ibcorr = fit_info.corr_id[2];
    r[0][0] = in[j][id1][t][0] - dmu * (in[j][ibcorr][t][0] - in[j][ibolla][t][0] * in[j][id0][t][0]);
    return r;
}

double** corr_plus_dm_correlated(int j, double**** in, int t, struct fit_type fit_info) {
    double** r = double_calloc_2(fit_info.N, 2);
    double dmu = fit_info.guess[1] - fit_info.guess[0];
    //{ 4,47, 40, 53,var/*55*/ };//P5P5, P5P5_corr_dmu,  P5P5_mu+dm ,P5P5_corr_bolla, bolla*P5P5dmu
    int id0 = fit_info.corr_id[0];
    int id_corr_dmu = fit_info.corr_id[1];
    int id_dmu = fit_info.corr_id[2];
    int id_corr_bolla = fit_info.corr_id[3];
    int ibolla = 42;
    int ibcorr = fit_info.corr_id[4];
    r[0][0] = in[j][id0][t][0] - (in[j][id_corr_dmu][t][0] - in[j][id_dmu][t][0])
        - dmu * (in[j][ibcorr][t][0] - in[j][ibolla][t][0] * in[j][id_corr_bolla][t][0]);
    return r;
}

double** mu_sea_correction(int j, double**** in, int t, struct fit_type fit_info) {
    double** r = double_calloc_2(fit_info.N, 2);
    double dmu = fit_info.guess[1] - fit_info.guess[0];
    //{ 4,47, 40, 53,var/*55*/ };//P5P5, P5P5_corr_dmu,  P5P5_mu+dm ,P5P5_corr_bolla, bolla*P5P5dmu
    int id0 = fit_info.corr_id[0];
    int id_corr_dmu = fit_info.corr_id[1];
    int id_dmu = fit_info.corr_id[2];
    int id_corr_bolla = fit_info.corr_id[3];
    int ibolla = 42;
    int ibcorr = fit_info.corr_id[4];
    r[0][0] = -dmu * (in[j][ibcorr][t][0] - in[j][ibolla][t][0] * in[j][id_corr_bolla][t][0]);
    r[0][1] = 0;
    return r;
}


void compute_syst_eq28(data_all in, const char* outpath, const char* filename) {
    int N = in.Nfits;
    int Njack = in.fits[0].Njack;
    char name[NAMESIZE];
    mysprintf(name, NAMESIZE, "%s/%s", outpath, filename);
    FILE* f = open_file(name, "w+");
    printf("writing: %s\n", name);
    std::vector<double> aves(N);
    std::vector<double> errors(N);
    for (int i = 0; i < N; i++) {
        aves[i] = in.fits[i].P[0][Njack - 1];
        errors[i] = error_jackboot(in.resampling.c_str(), Njack, in.fits[i].P[0]);
    }

    double ave = 0, err = 0;
    for (int i = 0; i < N; i++) {
        ave += aves[i];
        err += pow(errors[i], 2);
    }
    ave /= (double)N;

    for (int i = 0; i < N; i++) {
        err += pow(ave - aves[i], 2);
    }
    err = sqrt(err / (double)N);

    for (int i = 0; i < N; i++) {
        fprintf(f, "%s    %g     %g   %g   %g  %g  %d  %d\n", in.fits[i].name, aves[i], errors[i], ave, err, in.fits[i].chi2[Njack - 1], in.fits[i].dof, in.fits[i].Npar);
    }
    printf("systematics  %s: N=%d\n", filename, N);
    printf("mean(eq28)= %g  %g \n", ave, err);
    fclose(f);
}
