#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include <unistd.h>
#include <sys/time.h>
#include <fcntl.h>
#include "global.hpp"


#include "linear_fit.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"

double f(int n, int Nvar, double* x, int Npar, double* P) {

    return  P[0] * exp(x[0] * P[1]);
}

void test_res(non_linear_fit_result fit, int NE, int Npar) {
    printf("chi2=%-18g     chi2/dof=%g\n", fit.chi2, fit.chi2 / (NE - Npar));
    printf("min=%g   %g\n", fit.P[0], fit.P[1]);
    double chi = sqrt(fit.chi2 / (NE - Npar));
    double gchidof = 1.32723;
    double gP0 = 1.75309e-05, dgP1 = 6.99e-06 * gchidof;
    double gP1 = 0.0414701, dgP0 = 0.001104 * gchidof;
    if (fabs(fit.P[0] - gP0) > 1e-4) { printf("the minimum should be ( %g, %g )\n", gP0, gP1); exit(1); }
    if (fabs(fit.P[1] - gP1) > 1e-4) { printf("the minimum should be ( %g, %g )\n", gP0, gP1); exit(1); }
    if (fabs(chi - gchidof) > 1e-4) { printf("chi2/dof should be  %g \n", gchidof * gchidof); exit(1); }
    printf("the fit is correct!\n\n");
}


int main() {
    int NE = 13;
    double** y = double_malloc_2(NE, 2);
    double** x = double_malloc_2(NE, 2);
    double** sigmax = double_malloc_2(NE, 2);

    x[0][0] = 368.5;     y[0][0] = 73.332;        sigmax[0][0] = 0.66;    y[0][1] = 1.5;
    x[1][0] = 364.2;     y[1][0] = 62.668;        sigmax[1][0] = 0.66;    y[1][1] = 1.0;
    x[2][0] = 359.2;     y[2][0] = 52.004;        sigmax[2][0] = 0.66;    y[2][1] = 0.8;
    x[3][0] = 354.5;     y[3][0] = 44.006;        sigmax[3][0] = 0.66;    y[3][1] = 0.7;
    x[4][0] = 348.7;     y[4][0] = 34.675;        sigmax[4][0] = 0.66;    y[4][1] = 1.2;
    x[5][0] = 343.7;     y[5][0] = 28.010;        sigmax[5][0] = 0.66;    y[5][1] = 1.6;
    x[6][0] = 338.7;     y[6][0] = 22.678;        sigmax[6][0] = 0.66;    y[6][1] = 1.2;
    x[7][0] = 334.2;     y[7][0] = 17.346;        sigmax[7][0] = 0.66;    y[7][1] = 1.5;
    x[8][0] = 329.0;     y[8][0] = 14.680;        sigmax[8][0] = 0.66;    y[8][1] = 1.6;
    x[9][0] = 324.0;     y[9][0] = 10.681;        sigmax[9][0] = 0.66;    y[9][1] = 1.2;
    x[10][0] = 319.1;    y[10][0] = 8.015;        sigmax[10][0] = 0.66;   y[10][1] = 0.8;
    x[11][0] = 314.6;    y[11][0] = 6.682;        sigmax[11][0] = 0.66;   y[11][1] = 1.0;
    x[12][0] = 308.7;    y[12][0] = 5.349;        sigmax[12][0] = 0.66;   y[12][1] = 1.5;


    ///////////////////////
    {
        fit_type fit_info;
        fit_info.Njack = 1;
        fit_info.N = 1;
        fit_info.myen = std::vector<int>(NE);
        for (int i = 0;i < fit_info.myen.size();i++) fit_info.myen[i] = i;
        auto myen = fit_info.myen;
        fit_info.Nvar = 1; // not used
        fit_info.Npar = 2; // what we are minimizing
        ////// allocation
        int* en = (int*)malloc(sizeof(int) * fit_info.N);// we need to init en and en_tot to allocate the other 
        for (int e = 0;e < fit_info.N; e++) { en[e] = myen.size(); }
        int en_tot = 0;      for (int n = 0;n < fit_info.N;n++) { en_tot += en[n]; }// total data to fit

        double* guess = (double*)malloc(sizeof(double) * fit_info.Npar);
        guess[0] = 1; guess[1] = 0.001;
        fit_info.verbosity = -1;
        fit_info.acc = 1e-6;
        fit_info.h = { 0.001, 0.001 };
        fit_info.lambda = 1e-4;

        fit_info.function = f;
        printf("testing the LM fit\n");
        fit_info.verbosity = 0;
        non_linear_fit_result fit = non_linear_fit_Nf(fit_info.N, en, x, y, fit_info.Nvar, fit_info.Npar, fit_info.function, guess, fit_info);
        test_res(fit, NE, fit_info.Npar);
        free(fit.P);

        printf("testing the NM fit\n");
        fit_info.NM = true;
        fit_info.verbosity = 0;
        guess[0] = 1.7e-5; guess[1] = 0.041;
        fit_info.h = { 1e-6, 0.0001 };

        fit = non_linear_fit_Nf(fit_info.N, en, x, y, fit_info.Nvar, fit_info.Npar, fit_info.function, guess, fit_info);
        test_res(fit, NE, fit_info.Npar);
        free(fit.P);


        // fit_info.noderiv = true;
        // fit_info.acc = 1e-16;
        // fit_info.Prange = { 1e-5, 0.01 };
        // guess[0] = 1e-5; guess[1] = 0.04;
        // fit_info.verbosity = 0;
        // fit = non_linear_fit_Nf(fit_info.N, en, x, y, fit_info.Nvar, fit_info.Npar, fit_info.function, guess, fit_info);
        // test_res(fit, NE, fit_info.Npar);

// TOO SLOW TO CONVERGE
        // fit_info.second_deriv = true;
        // fit_info.verbosity = 2;
        // fit_info.maxiter=1e+9;
        // fit_info.acc = 1e-12;

        // guess[0] = 1.e-05; guess[1] = 0.04;
        // fit = non_linear_fit_Nf(fit_info.N, en, x, y, fit_info.Nvar, fit_info.Npar, fit_info.function, guess, fit_info);
        // test_res(fit, NE, fit_info.Npar);
        free(guess);free(en);

    }
    free_2(NE, x);
    free_2(NE, y);
    free_2(NE, sigmax);
}