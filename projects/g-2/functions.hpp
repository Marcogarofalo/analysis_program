#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include "tower.hpp"
#include "fit_all.hpp"
extern "C" {
#include "../external/rzeta/src/dzeta_function.h"
}
#include <algorithm>

using namespace std::complex_literals;

constexpr double muon_mass_MeV = 105.6583755;
constexpr double alpha_em = 0.0072973525693;
constexpr double muon_mass_fm = muon_mass_MeV * 197.326963;


double integrand_K(double x, void* params) {
    double z = *(double*)params;
    double arg = (z / 2) * x / (sqrt(1 - x));
    double j0 = sin(arg) / arg;
    double f = (1 - x) * (1 - j0 * j0);
    return f;
}

double kernel_K(double z, double epsrel = 1e-7) {
    int Maxiter = 1e+3;
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
double* compute_amu_sd(double**** in, int id, int Njack, double* Z, double* a, double q2, double (*int_scheme)(int, int, double*), FILE* outfile, const char* description, const char* resampling) {
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
            double VV_sub = Z[j] * Z[j] * in[j][id][t_a][0] - (1.0 / (2.0 * M_PI * M_PI * pow(t_a, 5)));
            ft[t_a] = K * VV_sub * (1 - theta);

            fi[t_a][j] = ft[t_a];
            Kt[t_a][j] = K;
            thetat[t_a][j] = (1 - theta);
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
    fprintf(outfile, "   %.15g   %15.g\n", amu[Njack - 1], error_jackboot(resampling, Njack, amu));

    free(ft);
    free_2(T / 2, fi);
    free_2(T / 2, Kt);
    free_2(T / 2, thetat);
    free_2(T / 2, corr_sub);
    return amu;
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
template<int idn, int idd>
double ZVl_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;

    double num = in[j][idn][t][0];
    double den = (in[j][idd][t + 1][0] - in[j][idd][t - 1][0]) / 2.0;
    return (2 * fit_info.mu * num / (den));

}

template<int id>
double GPS_OS_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double mu = file_head.musea;
    double mass = fit_info.ext_P[0][j];

    double GPS_OS = sqrt(in[j][id][t][0] * 2 * mass / (exp(-mass * t) + exp(-mass * (T - t))));

    return GPS_OS;
}

template<int id>
double GPS_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double mu = file_head.musea;
    double mass = fit_info.ext_P[0][j];

    double GPS_OS = sqrt(in[j][id][t][0] * 2 * mass / (exp(-mass * t) + exp(-mass * (T - t))));

    return GPS_OS;
}


double* interpol_Z(int Nmus, int Njack, double** Meta, double** Z, double* aMetas_exp,
    FILE* outfile, const char* description, const char* resampling) {
    //  Z(s1) = ( 1  Meta(s1) )  (P[0])
    //  Z(s1) = ( 1  Meta(s2) )  (P[1])

    double** matrix = double_malloc_2(Nmus, Nmus);
    double* Zj = (double*)malloc(sizeof(double) * Nmus);
    double* Zint = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++) {
        for (int i = 0;i < Nmus;i++) {
            Zj[i] = Z[i][j];
            for (int k = 0;k < Nmus;k++) {
                matrix[i][k] = pow(Meta[i][j], k);
            }
        }
        double* P = LU_decomposition_solver(Nmus, matrix, Zj);
        Zint[j] = 0;
        for (int k = 0;k < Nmus;k++) {
            Zint[j] += pow(aMetas_exp[j], k) * P[k];
        }
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
    double hM = (grhopipi * grhopipi * k * k * k * 2) * log((Mrho + 2 * k) / (2 * mass)) / (6.0 * M_PI * M_PI * Mrho);
    double h1 = (grhopipi * grhopipi * k * k) * (1 + (1 + 2 * mass * mass / (Mrho * Mrho) * (Mrho / k) * log((Mrho + 2 * k) / (2 * mass)))) / (6.0 * M_PI * M_PI * Mrho);
    double gamma = (grhopipi * grhopipi * k * k * k) / (6.0 * M_PI * omega * omega);

    double cotd = Mrho * Mrho - omega * omega - hM - (omega * omega - Mrho * Mrho) * h1 / (2 * Mrho) + ho;
    cotd /= omega * gamma;
    return cotd;
}

double compute_delta(double x, void* params) {
    double* P = (double*)params;
    int Nvar = P[0] + 1e-6;
    double* xh = (double*)malloc(sizeof(double) * Nvar);
    for (int i = 0;i < Nvar;i++) {
        xh[i] = P[i];
    }
    xh[0] = x;
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
    gsl_deriv_central(&F, x[0], 1e-8, &result, &abserr);
    double k = sqrt(x[0] * x[0] / 4. - x[1] * x[1]);
    result *= 2 * k / (sqrt(k * k + x[1] * x[1]));

    free(P);
    return result;
}


double w_js(int j, int s, double q) {
    double z[2] = { 0,0 };
    int dvec[3] = { 0, 0, 0 };
    double gamma = 1;
    double A = 1;
    dzeta_function(z, q * q, j, s, dvec, gamma, A, 1.e-3, 1.e+6, 5);
    // std::complex<double>  zc(z[0], z[1]);
    double r = z[0] / (pow(q, j + 1) * sqrt(2. * j + 1.) * pow(pi_greco, 3. / 2.));
    return r;
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

    // printf("%g %g %g\n", omega,k,qsq);exit(1);
    // dzeta_function(z, qsq, 0, 0, dvec, gamma, A, 1.e-3, 1.e+6, 5);
    // std::complex<double>  zc(z[0], z[1]);
    // double r = real(zc * 2. * pi_greco / (pow(pi_greco, 3. / 2.) * L * (a / 197.326963) * gamma * mass));
    double r = w_js(0, 0, q) - w_js(2, 0, q) - (3. / sqrt(6)) * (w_js(2, -2, q) + w_js(2, 2, q));

    double cotd = compute_cotd(Nvar, x);
    return cotd - r;

}

double compute_phi(double q, void* params) {
    (void)(params); /* avoid unused parameter warning */
    double M = w_js(0, 0, q) - w_js(2, 0, q) - (3. / sqrt(6)) * (w_js(2, -2, q) + w_js(2, 2, q));
    return -std::atan(1.0 / M);
}


double compute_deriv_phi(double q) {

    gsl_function F;
    double result, abserr;

    F.function = &compute_phi;
    F.params = 0;
    gsl_deriv_central(&F, q, 1e-8, &result, &abserr);

    return result;
}
inline double compute_FGS2(int Nvar, double* x) {
    double omega = x[0];
    double mass = x[1];
    double L = x[2];
    double Mrho = x[3];
    double grhopipi = x[4];
    double a = x[5];
    double k = sqrt(omega * omega / 4. - mass * mass);
    double q = k * (L * (a / 197.326963) / (2. * pi_greco));

    double ho = (grhopipi * grhopipi * k * k * k * 2) * log((omega + 2 * k) / (2 * mass)) / (6.0 * M_PI * M_PI * omega);
    double hM = (grhopipi * grhopipi * k * k * k * 2) * log((Mrho + 2 * k) / (2 * mass)) / (6.0 * M_PI * M_PI * Mrho);
    double h1 = (grhopipi * grhopipi * k * k) * (1 + (1 + 2 * mass * mass / (Mrho * Mrho) * (Mrho / k) * log((Mrho + 2 * k) / (2 * mass)))) / (6.0 * M_PI * M_PI * Mrho);
    double gamma = (grhopipi * grhopipi * k * k * k) / (6.0 * M_PI * omega * omega);

    double Apipi0 = hM - Mrho * h1 / 2. + pow(grhopipi * mass, 2) / (6 * M_PI * M_PI);
    std::complex<double> Apipi = hM + (omega * omega - Mrho * Mrho) * h1 / (2 * Mrho) - ho + 1i * omega * gamma;
    std::complex<double> FGS = (Mrho * Mrho - Apipi0) / (Mrho * Mrho - omega * omega - Apipi);
    double FGS2 = std::abs(FGS);
    FGS2 *= FGS2;
    return FGS2;
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
    r /= k * deriv_d + q * deriv_phi;

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
    int Maxiter = 1e+3;
    int epsrel = 1e-5;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(Maxiter);

    for (int t = 1; t < T;t++) {
        P[Nvar] = t * (a / 197.326963);

        double result, error;
        gsl_function F;
        F.function = &integrand_V_infL;
        F.params = (void*)P;
        // printf("here : %g\n",P[Nvar]);
                //  gsl_integration_qags(&F, 2*mass, 10000, 1e+5, epsrel, Maxiter, w, &result, &error);
        gsl_integration_qagiu(&F, 2 * mass, 1e-4, epsrel, Maxiter, w, &result, &error);
        // printf("K(%g)=%g\n", z, result);

        V[t] = result / (48. * M_PI * M_PI);
    }
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

double** compute_DVt(int L, int Njack, double* Mpi, double* Mrho, double* a, double* grhopipi, FILE* outfile, const char* description, const char* resampling) {

    int T = file_head.l0;
    constexpr int Nvar = 6;
    constexpr int Nomegas = 5;
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
    int p2s[Nomegas+1];
    count = 0;
    p2s[0]=0;
    while (count < Nomegas+1) {
        double old = p2[0];
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
    for (int j = 0;j < Njack;j++) {
        double x[Nvar];
        x[0] = 1;
        x[1] = Mpi[j];
        x[2] = L;
        x[3] = Mrho[j];
        x[4] = grhopipi[j];
        x[5] = a[j];
        // first we plot cotd*k/m
        if (j == Njack-1) {
            for (int i = 0;i < 1500;i++) {
                x[0] = x[1] * (2. + 30.5 * i / 1500.0);
                double k = sqrt(x[0] * x[0] / 4 - x[1] * x[1]);
                double cotd = compute_cotd(Nvar, x);
                double q = k * (L * (x[5] / 197.326963) / (2. * pi_greco));
                double r = w_js(0, 0, q) - w_js(2, 0, q) - (3. / sqrt(6)) * (w_js(2, -2, q) + w_js(2, 2, q));
                fprintf(outfile, "%g    %g  %g\n", x[0], cotd, r);
            }

            fprintf(outfile, "\n\n #QC2 fit in [%d,%d] chi2=%.5g  %.5g\n", 0, 0, 0.0, 0.0);

        }
        double p = 2 * M_PI / ((double)L * (x[5] / 197.326963));

        for (int i = 0; i < Nomegas;i++) {
            double xmin = 2 * sqrt(x[1] * x[1] + p2s[i] * p * p) + 1e-8;
            double xmax = 2 * sqrt(x[1] * x[1] + p2s[i + 1] * p * p) - 1e-8;
            // printf("%d   %d\n", p2s[i], p2s[i + 1]);
            omega[i][j] = rtsafe(omega_QC2, 0/*n*/, 6, x, 0/* Npar */, nullptr/* Pkcot */, 0/*ivar*/, 0. /*input*/, xmin, xmax, 1e-5, 100, 1e-4);
            x[0] = omega[i][j];
            nuA2[i][j] = matrix_element_nuA2(Nvar, x);
            count++;
        }

        VinfL[j] = compute_V_infL(Nvar, x, T);
        VL[j] = compute_VL(Nomegas, omega, nuA2, a, T, j);
        for (int t = 0;t < T;t++) {
            DVt[t][j] = VinfL[j][t] - VL[j][t];
        }
    }


    for (int i = 0; i < Nomegas;i++) {
        printf("omega[%d]=%g  %g\n", i, omega[i][Njack - 1], error_jackboot(resampling, Njack, omega[i]));
        printf("nuA2[%d]=%g  %g\n", i, nuA2[i][Njack - 1], error_jackboot(resampling, Njack, nuA2[i]));
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



///////////////
double rhs_amu_common(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a = x[0];
    if (n == 0)      r = P[0] + a * P[1];
    else if (n == 1) r = P[0] + a * P[2];

    return r;
}
double rhs_amu_common_a4(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    if (n == 0)      r = P[0] + a2 * P[1] + a2 * a2 * P[3];
    else if (n == 1) r = P[0] + a2 * P[2] + a2 * a2 * P[4];

    return r;
}

template<int ieq, int iop>
double lhs_amu_common(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    if (n == 0)        r = gjack.en[e].jack[ieq][j];
    else if (n == 1)   r = gjack.en[e].jack[iop][j];
    return r;
}