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

    int N =  600;
    double d = 0.005;
    double eps =  0.0005;
    printf("E3    M3_re   M3_im      Kdf_re  kdf_im  Finf_re  Finf_im\n");
    for (int i = 0;i < 1;i++) {
        double E3 = 3.02 + i * 0.000001;
    // for (int i = 0;i < 2;i++) {
    //     double E3 = E3_m + i * 0.005;

        std::complex<double> M3 = compute_M3_sym(E3, N, Npar, P, compute_kcot, PKiso, compute_kiso,  eps);
        std::complex<double> Kdf = compute_kiso(E3, PKiso);
        Eigen::MatrixXcd D = compute_D(E3, N, Npar, P, compute_kcot,  eps);

        std::complex<double> Finf = comput_Finf(E3, D, N, Npar, P, compute_kcot,  eps);

        printf("%-18.8g%-14g%-18g%-14g%-18g%-18.12g%-20.12g\n", E3, real(M3), imag(M3), real(Kdf), imag(Kdf), real(Finf), imag(Finf));

    }

}