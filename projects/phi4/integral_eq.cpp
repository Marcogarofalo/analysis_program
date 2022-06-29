#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <iostream>
#include "integral_eq_QC3.hpp"
#include "mutils.hpp"
#include "non_linear_fit.hpp"



double rhs_laurent_pole(int n, int Nvar, double* x, int Npar, double* P) {
    error(Npar % 2 != 0, 1, "rhs_laurent_pole:", "Npar=%d but it must be multiple of two since the parameters are complex", Npar);
    std::complex<double> z(x[0], x[1]);

    std::complex<double> p(P[0], P[1]);
    std::complex<double> am1(P[2], P[3]);

    std::complex<double> r = am1 / (z * z - p);

    if (Npar >= 6) {
        std::complex<double> a0(P[4], P[5]);
        r += a0;
    }

    if (n == 0)      return real(r);
    else if (n == 1) return imag(r);
    else { printf("%s\n", __func__);  exit(1); }
}

double compute_kiso(double E3_m, double* P) {
    return P[0] / (E3_m * E3_m - P[1]) + P[2];
}

double compute_kcot(int Nvar, double* x, int Npar, double* P) {
    double r;
    r = 1. / P[0];

    return r;
}

int main(int argc, char** argv) {


    int Npar = 1;
    double* P = (double*)malloc(sizeof(double) * 1);
    P[0] = -0.149458;
    double* PKiso = (double*)malloc(sizeof(double) * 3);
    PKiso[0] = 204.692;
    PKiso[1] = 9.12076;
    PKiso[2] = -2491.55;

    int Njack = 15;
    int seed = 1;
    int tot_parK = 4;
    double* mean = (double*)malloc(sizeof(double) * tot_parK);
    mean[0] = 204.692;
    mean[1] = 9.12076;
    mean[2] = -2491.55;
    mean[3] = -0.149458;
    double** cov = double_calloc_2(tot_parK, tot_parK);
    cov[0][0] = pow(23, 2);
    cov[1][1] = pow(0.0015, 2);
    cov[2][2] = pow(6.7e+2, 2);
    cov[3][3] = pow(0.002, 2);
    // cov[0][1] = -1.77; cov[0][2] = 9.49;     cov[0][3] = 3.95;
    // ;                  cov[1][2] = 0.000297; cov[1][3] = -8.63e-5;
    // ; ;                                      cov[2][3] = -254;

    for (int i = 0; i < tot_parK;i++) {
        for (int j = i + 1; j < tot_parK;j++) {
            cov[i][j] = cov[i][j] * sqrt(cov[i][i] * cov[j][j]);
            cov[j][i] = cov[i][j];
        }
    }

    for (int i = 0; i < tot_parK;i++) {
        cov[i][i] += 1e-7;
    }

    double** tmp = fake_sampling_covariance("jack", mean, Njack, tot_parK, cov, seed);


    int N = 600;
    // int N = 1000;
    double d = 0.005;
    double eps = 0.0005;
    printf("E3    M3_re   M3_im      Kdf_re  kdf_im  Finf_re  Finf_im\n");


    double Emin = 3.0189;
    double Emax = 3.0191;
    int NE = 6;
    double dE = (Emax - Emin) / ((double)NE);

    std::vector<double> E3(NE);
    // std::vector<std::complex<double>> M3(NE);
    double*** M3 = double_malloc_3(NE, 2, Njack);
    printf("resampling parameters:\n");
    for (int i = 0; i < tot_parK;i++)
        printf("%g   %g\n", tmp[i][Njack - 1], error_jackboot("jack", Njack, tmp[i]));

    printf("computing M3:\n");
    for (int i = 0;i < NE;i++) {
        for (int j = 0;j < Njack;j++) {
            E3[i] = Emin + i * dE;
            P[0] = tmp[3][j];
            PKiso[0] = tmp[0][j];
            PKiso[1] = tmp[1][j];
            PKiso[2] = tmp[2][j];

            std::complex<double> m3 = compute_M3_sym(E3[i], N, Npar, P, compute_kcot, PKiso, compute_kiso, eps);
            M3[i][0][j] = m3.real(); M3[i][1][j] = m3.imag();

            std::complex<double> Kdf = compute_kiso(E3[i], PKiso);
            Eigen::MatrixXcd D = compute_D(E3[i], N, Npar, P, compute_kcot, eps);
            std::complex<double> Finf = comput_Finf(E3[i], D, N, Npar, P, compute_kcot, eps);
            // printf("jack =%-4d%-18.8g%-14g%-18g%-14g%-18g||%-25g%-25g%-25g%-25g\n",
                // j, E3[i], real(m3), imag(m3), real(Kdf), imag(Kdf), PKiso[0], PKiso[1], PKiso[2], P[3]);
            printf("%-18.8g%-14g%-18g%-14g%-18g%-18.12g%-20.12g\n", E3[i], real(m3), imag(m3), real(Kdf), imag(Kdf), real(Finf), imag(Finf));
        }
        printf("%-18.8g%-14g%-18g%-14g%-18g\n", E3[i], M3[i][0][Njack - 1], error_jackboot("jack", Njack, M3[i][0]),
            M3[i][1][Njack - 1], error_jackboot("jack", Njack, M3[i][1]));
    }

    exit(1);

    {
        fit_type fit_info;
        fit_info.Njack = 1;
        fit_info.N = 2;
        fit_info.myen = std::vector<int>(NE);
        fit_info.Nvar = 2;
        fit_info.Npar = 4;
        fit_info.function = rhs_laurent_pole;
        ////// allocation
        int* en = (int*)malloc(sizeof(int) * fit_info.N);// we need to init en and en_tot to allocate the other 
        for (int e = 0;e < fit_info.N; e++) { en[e] = fit_info.myen.size(); }
        int en_tot = 0;      for (int n = 0;n < fit_info.N;n++) { en_tot += en[n]; }// total data to fit

        double*** y = double_malloc_3(fit_info.Njack, en_tot, 2);// 2 is mean/error
        double*** x = double_malloc_3(fit_info.Njack, en_tot, fit_info.Nvar);

        //init x
        int count = 0;
        for (int j = 0;j < fit_info.Njack;j++) {
            int count = 0;
            for (int n = 0;n < fit_info.N;n++) {
                for (int e = 0;e < en[n];e++) {
                    x[j][count][0] = E3[e];
                    x[j][count][1] = 0;
                    count++;
                }
            }
        }
        ////////////////////////////////////////// y
        count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < en[n];e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    double yd;
                    y[j][e + count][0] = M3[e][n][j];// y must be zero

                }
                double dy = error_jackboot("jack", Njack, M3[e][n]);
                for (int j = 0;j < fit_info.Njack;j++)
                    y[j][e + count][1] = dy;
            }
            count += en[n];
        }
        //////  init end
        int j = fit_info.Njack - 1.;
        fit_info.guess = { 9.1145,   2.08327e-06,   254.085,   -10.0645 };
        double* guess = (double*)malloc(sizeof(double) * fit_info.Npar);
        for (int i = 0;i < fit_info.Npar; i++) {
            guess[i] = fit_info.guess[i];
        }


        fit_info.noderiv = false;
        fit_info.Prange = { 1e-3,1e-3,100,100 };
        fit_info.verbosity = 0;
        fit_info.acc = 1e-9;
        fit_info.h = 1e-4;
        fit_info.lambda = 1e-4;
        fit_info.precision_sum = 2;
        fit_info.repeat_start = 1;
        // fit_info.second_deriv=true;
        // guess= guess_for_non_linear_fit_Nf(fit_info.N, en, x[j], y[j], fit_info.Nvar, fit_info.Npar, fit_info.function, guess, fit_info);
        double* P = non_linear_fit_Nf(fit_info.N, en, x[j], y[j], fit_info.Nvar, fit_info.Npar, fit_info.function, guess, fit_info);
        double chi2 = compute_chi_non_linear_Nf(fit_info.N, en, x[j], y[j], P, fit_info.Nvar, fit_info.Npar, fit_info.function);
        printf("chi2/dof=%g   dof=%d\n",
            chi2 / (en_tot - fit_info.Npar), (en_tot - fit_info.Npar));
        printf("p=%g   %g   am1=%g   %g   \n", P[0], P[1], P[2], P[3]);

        free(P);
        fit_info.restore_default();
    }

}