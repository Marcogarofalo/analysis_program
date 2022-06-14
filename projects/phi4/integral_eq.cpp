#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <iostream>
#include "integral_eq_QC3.hpp"



std::complex<double> compute_kiso(double E3_m, double* P) {
    return P[0] / (E3_m * E3_m - P[1]) + P[2];
}

double compute_kcot(int Nvar, double* x, int Npar, double* P) {
    double r;
    r = 1. / P[0];

    return r;
}

int main(int argc, char** argv) {
    double E3_m = 3.0150;

    int Npar = 1;
    double* P = (double*)malloc(sizeof(double) * 1);
    P[0] = -0.149458;
    double* PKiso = (double*)malloc(sizeof(double) * 3);
    PKiso[0] = 204.692;
    PKiso[1] = 9.12076;
    PKiso[2] = -2491.55;

    // test paper
    // double E3_m = 3.3;
    // int Npar = 1;
    // double* P = (double*)malloc(sizeof(double) * 1);
    // P[0] = -0.296;
    // double* PKiso = (double*)malloc(sizeof(double) * 3);
    // PKiso[0] = 0;
    // PKiso[1] = 1e+6;
    // PKiso[2] = 0;

    int N = 100;
    double d = 0.005;
    double eps = 1e-4;
    printf("E3    M3       Kdf\n");
    for (int i = 0;i < 100;i++) {
        double E3 = E3_m +i*0.0001;
        std::complex<double> M3 = compute_M3_sym(E3, N, Npar, P, compute_kcot, PKiso, compute_kiso, d, eps);
        std::complex<double> Kdf= compute_kiso(E3 , PKiso);
        printf("%-16g%-12g%-16g%-12g%-12g\n",E3, real(M3),imag(M3),real(Kdf),imag(Kdf));
        
    }

}