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
#include "minimizer.hpp"

using namespace std::complex_literals;


double rhs_amu_common_a4(int n, int Nvar, double* x, int Npar, double* P) {

    std::complex<double> r(P[0], P[1]);
    std::complex<double> ir(0.5, 0.3);

    r = ((r - ir) * std::conj(r - ir)) + 1. + ((r - ir) * std::conj(r - ir)) * ((r - ir) * std::conj(r - ir));
    return real(r);
}

int main() {


    {
        fit_type fit_info;
        fit_info.Njack = 1;
        fit_info.N = 1;
        fit_info.myen = { 1 };
        auto myen = fit_info.myen;
        fit_info.Nvar = 1; // not used
        fit_info.Npar = 2; // what we are minimizing
        fit_info.function = rhs_amu_common_a4;
        ////// allocation
        int* en = (int*)malloc(sizeof(int) * fit_info.N);// we need to init en and en_tot to allocate the other 
        for (int e = 0;e < fit_info.N; e++) { en[e] = myen.size(); }
        int en_tot = 0;      for (int n = 0;n < fit_info.N;n++) { en_tot += en[n]; }// total data to fit

        double*** y = double_malloc_3(fit_info.Njack, en_tot, 2);// 2 is mean/error
        double*** x = double_malloc_3(fit_info.Njack, en_tot, fit_info.Nvar);

        //init x
        int count = 0;
        for (int j = 0;j < fit_info.Njack;j++) {
            int count = 0;
            for (int n = 0;n < fit_info.N;n++) {
                for (int e = 0;e < en[n];e++) {
                    for (int i = 0; i < fit_info.Nvar; i++) {
                        x[j][count][i] = 1;
                    }
                    count++;
                }
            }
        }
        ////////////////////////////////////////// y
        count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < en[n];e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    y[j][e + count][0] = 0;// y must be zero
                    y[j][e + count][1] = 1.0;
                }
            }
            count += en[n];
        }
        //////  init end
        int j = fit_info.Njack - 1.;
        double* guess = (double*)malloc(sizeof(double) * fit_info.Npar);
        guess[0] = 1;
        guess[1] = 1;
        // fit_info.noderiv=true;
        fit_info.verbosity = -1;
        fit_info.acc = 1e-6;
        fit_info.h = { 0.001, 0.001 };
        fit_info.lambda = 1e-4;
        fit_info.second_deriv = true;
        double* P = non_linear_fit_Nf(fit_info.N, en, x[j], y[j], fit_info.Nvar, fit_info.Npar, fit_info.function, guess, fit_info);
        printf("min=%g   %g\n", P[0], P[1]);
        if (fabs(P[0] - 0.5) > 1e-4) { printf("the minimum should be ( 0.5, 0.3 )\n"); exit(1); }
        if (fabs(P[1] - 0.3) > 1e-4) { printf("the minimum should be ( 0.5, 0.3 )\n"); exit(1); }

        free(P);
        fit_info.noderiv = true;
        P = non_linear_fit_Nf(fit_info.N, en, x[j], y[j], fit_info.Nvar, fit_info.Npar, fit_info.function, guess, fit_info);
        printf("min=%g   %g\n", P[0], P[1]);
        if (fabs(P[0] - 0.5) > 1e-4) { printf("the minimum should be ( 0.5, 0.3 )\n"); exit(1); }
        if (fabs(P[1] - 0.3) > 1e-4) { printf("the minimum should be ( 0.5, 0.3 )\n"); exit(1); }
    }
    ////////////////////////
    //   using the build in function
    {
        fit_type   fit_info;
        fit_info.N = 1;
        fit_info.second_deriv = true;
        fit_info.mean_only = false;
        fit_info.Njack = 1;
        fit_info.Nvar = 0; // it is important that the value is correct since we need to pass all the x
        fit_info.Npar = 2; // what we are minimizing
        fit_info.function = rhs_amu_common_a4;
        fit_result min = minimize_functions_Nf(fit_info);
        printf("min=%g   %g\n", min.P[0][0], min.P[1][0]);
        if (fabs(min.P[0][0] - 0.5) > 1e-4) { printf("the minimum should be ( 0.5, 0.3 )\n"); exit(1); }
        if (fabs(min.P[1][0] - 0.3) > 1e-4) { printf("the minimum should be ( 0.5, 0.3 )\n"); exit(1); }
    }



}