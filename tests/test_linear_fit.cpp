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

//a*1+b*x
double* linear(int M, double x) {
    double* r;

    r = (double*)malloc(sizeof(double) * M);
    r[0] = 1.;// function that multiply a
    r[1] = x;// function that multiply b

    return r;

}



//a*1+b*x+c*x^2
double* parabolic(int M, double x) {
    double* r;

    r = (double*)malloc(sizeof(double) * M);
    r[0] = 1.;// function that multiply a
    r[1] = x;// function that multiply b
    r[2] = x * x;// function that multiply c
    return r;

}


//a*1+b*x+c/x
double* pole(int M, double x) {
    double* r;

    r = (double*)malloc(sizeof(double) * M);
    r[0] = 1.;// function that multiply a
    r[1] = x;// function that multiply b
    r[2] = 1. / x;// function that multiply c
    return r;

}

//a*1+b*y+c/x
double* pol_xy(int M, double* x) {
    double* r;

    r = (double*)malloc(sizeof(double) * M);
    r[0] = 1.;// function that multiply a
    r[1] = x[0];// function that multiply b
    r[2] = x[1];// function that multiply c
    return r;

}


int main() {
    double** y, * x;
    double** fit, chi2;
    int N = 4, i;

    x = (double*)malloc(sizeof(double) * N);
    y = (double**)malloc(sizeof(double*) * N);
    for (i = 0;i < N;i++)
        y[i] = (double*)malloc(sizeof(double) * 2);

    //x       y_mean       y_error
    x[0] = 1;    y[0][0] = 0.9;   y[0][1] = 0.07;
    x[1] = 1.5;  y[1][0] = 1.6;   y[1][1] = 0.08;
    x[2] = 2;    y[2][0] = 2.1;   y[2][1] = 0.09;
    x[3] = 2.5;  y[3][0] = 2.45;  y[3][1] = 0.1;

    printf("fitting data as a+bx\n");
    fit = linear_fit(N, x, y, 2, linear);
    chi2 = compute_chisqr(N, x, y, 2, fit, linear);
    printf("chi2=%f\n", chi2);
    //                                  a       a_err      b         b_err
    printf("a=%g+-%g \t  b=%g+-%g\n", fit[0][0], fit[0][1], fit[1][0], fit[1][1]);

    free(fit);
    /////////////////////////
    printf("fitting data as a+bx+c*x^2\n");
    fit = linear_fit(N, x, y, 3, parabolic);
    chi2 = compute_chisqr(N, x, y, 3, fit, parabolic);
    printf("chi2=%f\n", chi2);
    //                                            a       a_err      b         b_err
    printf("a=%g+-%g \t  b=%g+-%g   c=%g+-%g\n", fit[0][0], fit[0][1], fit[1][0], fit[1][1], fit[2][0], fit[2][1]);

    free(fit);
    /////////////////////////
    printf("fitting data as a+bx+c/x\n");
    fit = linear_fit(N, x, y, 3, pole);
    chi2 = compute_chisqr(N, x, y, 3, fit, pole);
    printf("chi2=%f\n", chi2);
    //                                            a       a_err      b         b_err
    printf("a=%g+-%g \t  b=%g+-%g   c=%g+-%g\n", fit[0][0], fit[0][1], fit[1][0], fit[1][1], fit[2][0], fit[2][1]);


    ///////////////////
    for (i = 0;i < N;i++)
        free(y[i]);
    free(y);
    printf("2 variables\n");

    double** xy;
    xy = (double**)malloc(sizeof(double*) * 5);
    y = (double**)malloc(sizeof(double*) * 5);
    for (i = 0;i < 5;i++) {
        xy[i] = (double*)malloc(sizeof(double) * 2);
        y[i] = (double*)malloc(sizeof(double) * 2);
    }
    xy[0][0] = 1; xy[0][1] = 1.3; y[0][0] = 6.2;  y[0][1] = 0.1;
    xy[1][0] = 2; xy[1][1] = 1.5; y[1][0] = 8.8;  y[1][1] = 0.1;
    xy[2][0] = 3; xy[2][1] = 1.2; y[2][0] = 9.9;  y[2][1] = 0.1;
    xy[3][0] = 4; xy[3][1] = 1.1; y[3][0] = 11.6; y[3][1] = 0.1;
    xy[4][0] = 5; xy[4][1] = 1.7; y[4][0] = 15.4; y[4][1] = 0.1;
    y[2][0] = 9.2;
    fit = linear_fit_Nvar(5, xy, y, 3, pol_xy);
    chi2 = compute_chisqr_Nvar(5, xy, y, 3, fit, pol_xy);
    printf("chi2=%f\n", chi2);
    printf("a=%g+-%g \t  b=%g+-%g   c=%g+-%g\n", fit[0][0], fit[0][1], fit[1][0], fit[1][1], fit[2][0], fit[2][1]);
}



/*

   degrees of freedom    (FIT_NDF)                        : 2
rms of residuals      (FIT_STDFIT) = sqrt(WSSR/ndf)    : 1.47882
variance of residuals (reduced chisquare) = WSSR/ndf   : 2.18691
p-value of the Chisq distribution (FIT_P)              : 0.112263

Final set of parameters            Asymptotic Standard Error
=======================            ==========================
a               = -0.0969224       +/- 0.1889       (194.9%)
b               = 1.06323          +/- 0.1115       (10.49%)

   */




   /*
   *
   degrees of freedom    (FIT_NDF)                        : 1
   rms of residuals      (FIT_STDFIT) = sqrt(WSSR/ndf)    : 0.131126
   variance of residuals (reduced chisquare) = WSSR/ndf   : 0.0171939
   p-value of the Chisq distribution (FIT_P)              : 0.895676

   Final set of parameters            Asymptotic Standard Error
   =======================            ==========================
   a               = -1.012           +/- 0.05988      (5.917%)
   b               = 2.26719          +/- 0.07628      (3.364%)
   c               = -0.353508        +/- 0.02221      (6.282%)

   */


   /*
       * degrees of freedom    (FIT_NDF)                        : 1
   rms of residuals      (FIT_STDFIT) = sqrt(WSSR/ndf)    : 0.250358
   variance of residuals (reduced chisquare) = WSSR/ndf   : 0.0626791
   p-value of the Chisq distribution (FIT_P)              : 0.802311

   Final set of parameters            Asymptotic Standard Error
   =======================            ==========================
   a               = 1.87203          +/- 0.2396       (12.8%)
   b               = 0.464971         +/- 0.07457      (16.04%)
   c               = -1.43875         +/- 0.1735       (12.06%)

       */
