#include <math.h>
#include <complex.h>
#include "tower.hpp"
#include <Eigen/Dense>  
#include <iostream>

using namespace std::complex_literals;

double kcot_generic(int Nvar, double* x, int Npar, double* P) {
    double r;
    r = 1. / P[0];

    return r;
}

std::complex<double> M2_from_kcot(int Nvar, double* x, int Npar, double* P, double kcot(int, double*, int, double*)) {
    double k_m = x[0];
    double E3_m = x[1];
    double k_mcotd = kcot(Nvar, x, Npar, P);
    // double k_m = sqrt(E_m * E_m / 4. - 1);
    double omega_k = sqrt(k_m * k_m + 1);
    double Es_mk2 = pow(E3_m - omega_k, 2) - k_m * k_m;

    std::complex<double> r;
    if (Es_mk2 >= 4)
        r *= -1i * sqrt(Es_mk2 / 4. - 1);
    else
        r *= sqrt(fabs(Es_mk2 / 4. - 1));
    r += k_mcotd;
    r /= ( 16. * M_PI * sqrt(Es_mk2));
    return 1. / r;
}

inline double J(double x) {
    if (x <= 0) return 0;
    else if (x > 0 && x < 1) return exp(-(1. / x) * exp(-(1. / (1. - x))));
    else if (x >= 1) return 1;
    else {
        printf("error in J(x): x<0");
        exit(-1);
    }
}

double H(double k, double p, double E3_m) {

    double omega_k = sqrt(k * k + 1);
    double Es_mk2 = pow(E3_m - omega_k, 2) - k * k;
    double omega_p = sqrt(p * p + 1);
    double Es_mp2 = pow(E3_m - omega_p, 2) - p * p;
    // if (k == 0 || p == 0)
    //     printf("k,p=%g %g  omega_k=%g omega_p=%g  Es_mk2=%g  Es_mp2=%g\n", k, p, omega_k, omega_p, Es_mk2, Es_mp2);
    // alpha hard coded to -1
    return J(Es_mk2 / 4.) * J(Es_mp2 / 4.);

}

std::complex<double> Gs_pke(double k, double p, double E3_m, double d, double eps) {

    double pk = p * k;
    double omega_k = sqrt(k * k + 1);
    double omega_p = sqrt(p * p + 1);
    std::complex<double> r;

    if (p < d) {
        r = H(k, p, E3_m) / (1 + E3_m * E3_m + 1i * eps + 2 * omega_k - 2 * E3_m * (1 + omega_k));
    }
    else if (k < d) {
        r = H(k, p, E3_m) / (1 + E3_m * E3_m + 1i * eps + 2 * omega_p - 2 * E3_m * (1 + omega_p));
    }
    else {
        r = log((2 * pk - pow(E3_m - omega_k - omega_p, 2) + p * p + k * k + 1 - 1i * eps) /
            (-2 * pk - pow(E3_m - omega_k - omega_p, 2) + p * p + k * k + 1 - 1i * eps));
        // std::cout << "log=" << r << std::endl;
        r *= -H(k, p, E3_m) / (4 * pk);
        // std::cout << "H=" << -H(k, p, E3_m) / (4 * pk) << std::endl;
        // std::cout << "r=" << r << std::endl;

    }
    return r;
}

double Pk(double k, double d) {
    double omega_k = sqrt(k * k + 1);
    return k * k * d / (M_PI * M_PI * 2 * 2 * omega_k);
}


Eigen::MatrixXcd compute_D(double E3_m, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {

    Eigen::MatrixXcd oneplusMGP = Eigen::MatrixXcd::Identity(N, N);
    Eigen::MatrixXcd MGM = Eigen::MatrixXcd::Zero(N, N);

    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[2];
    double kmax = sqrt(pow((E3_m * E3_m + 1) / (2 * E3_m), 2) - 1);
    double d = kmax / N;
    for (int i = 0; i < N;i++) {
        double p = d * i;
        x[0] = p;
        x[1] = E3_m;
        std::complex<double> m2 = M2_from_kcot(Nvar, x, Npar, P, kcot);
        // std::cout << "m2=" << m2 << std::endl;
        for (int j = 0; j < N;j++) {
            double k = d * j;
            std::complex<double> G = Gs_pke(p, k, E3_m, d, eps);
            // std::cout << "G=" << G << std::endl;

            oneplusMGP(i, j) += m2 * G * Pk(k, d);
            x[0] = k;
            MGM(i, j) += m2 * G * M2_from_kcot(Nvar, x, Npar, P, kcot);
            // std::cout << "m2(p)=" << M2_from_kcot(Nvar, x, Npar, P, kcot) << std::endl;
        }
    }
    Eigen::MatrixXcd oneplusMGP1 = oneplusMGP.inverse();
    Eigen::MatrixXcd D = -1 * oneplusMGP1 * MGM;
    // for (int i = 0; i < N;i++) {
    //     double k = d * i;
    //     x[0] = k;
    //     printf("m2(%d)=%g    %g    \n", i, M2_from_kcot(Nvar, x, Npar, P, kcot).real(), M2_from_kcot(Nvar, x, Npar, P, kcot).imag());
    //     printf("Gs_pke(%d,1)=%g    %g   Gs_pke(1,%d)=%g    %g  \n", i, Gs_pke(k, d, E3_m, d, eps).real(), Gs_pke(k, d, E3_m, d, eps).imag(), i, Gs_pke(d, k, E3_m, d, eps).real(), Gs_pke(d, k, E3_m, d, eps).imag());

    //     printf("MGM(%d,1)=%g    %g   MGM(1,%d)=%g    %g  \n", i, MGM(i, 1).real(), MGM(i, 1).imag(), i, MGM(1, i).real(), MGM(1, i).imag());
    //     printf("oneplusMGP(%d,1)=%g    %g   oneplusMGP(1,%d)=%g    %g  \n", i, oneplusMGP(i, 1).real(), oneplusMGP(i, 1).imag(), i, oneplusMGP(1, i).real(), oneplusMGP(1, i).imag());
    //     printf("oneplusMGP1(%d,1)=%g    %g   oneplusMGP1(1,%d)=%g    %g  \n", i, oneplusMGP1(i, 1).real(), oneplusMGP1(i, 1).imag(), i, oneplusMGP1(1, i).real(), oneplusMGP1(1, i).imag());
    //     printf("D(%d,1)=%g    %g   D(1,%d)=%g    %g  \n", i, D(i, 1).real(), D(i, 1).imag(), i, D(1, i).real(), D(1, i).imag());
    // }
    return D;
}

std::complex<double> rho(double k_m, double E3_m) {
    double omega_k = sqrt(k_m * k_m + 1);
    double Es_mk2 = pow(E3_m - omega_k, 2) - k_m * k_m;
    std::complex<double> r = J(Es_mk2 / 4.) / (16. * M_PI * sqrt(Es_mk2));
    if (Es_mk2 >= 4)
        r *= -1i * sqrt(Es_mk2 / 4. - 1);
    else
        r *= sqrt(fabs(Es_mk2 / 4. - 1));
    return r;
}

Eigen::VectorXcd  cal_L(double E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {
    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[2];
    Eigen::VectorXcd L(N);
    double kmax = sqrt(pow((E3_m * E3_m + 1) / (2 * E3_m), 2) - 1);
    double d = kmax / N;
    for (int i = 0;i < N;i++) {
        double k = i * d;
        L(i) = std::complex<double>(1. / 3., 0);
        x[0] = k;
        x[1] = E3_m;
        L(i) -= M2_from_kcot(Nvar, x, Npar, P, kcot) * rho(k, E3_m);
        // std::cout << "M2="<<M2_from_kcot(Nvar, x, Npar, P, kcot) << "   rho=" << rho(k, E3_m)<< std::endl;
        for (int j = 0;j < N;j++) {
            double p = j * d;
            // std::cout << "D =" << D(i, j) << "  pk" << Pk(p, d) << "  rho" << rho(p, E3_m) << std::endl;
            L(i) -= D(i, j) * Pk(p, d) * rho(p, E3_m);
        }
    }
    return L;
}


Eigen::VectorXcd  cal_L1(double E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {
    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[2];
    Eigen::VectorXcd L(3);

    double k[3];
    k[0] = 0;
    k[1] = sqrt(pow(E3_m - 1, 2) - 4);
    k[2] = sqrt(pow(E3_m - 1, 2) - 4);// this with a minus
    double kmax = sqrt(pow((E3_m * E3_m + 1) / (2 * E3_m), 2) - 1);
    double d = kmax / N;
    for (int i = 0;i < 3;i++) {

        L(i) = std::complex<double>(1. / 3., 0);
        x[0] = k[i];
        x[1] = E3_m;
        L(i) -= M2_from_kcot(Nvar, x, Npar, P, kcot) * rho(k[i], E3_m);

        int i_min = k[i] / d;
        int i_max = k[i] / d + 1;
        if (i_max >= N) { printf("error N*d too small\n"); exit(-1); }
        // std::cout << "M2="<<M2_from_kcot(Nvar, x, Npar, P, kcot) << "   rho=" << rho(k, E3_m)<< std::endl;
        for (int j = 0;j < N;j++) {
            double p = j * d;

            std::complex<double> m = (D(i_max, j) - D(i_min, j)) / d;
            std::complex<double> q = D(i_min, j) - m * ((double)i_min * d);
            std::complex<double> Dave = q + m * k[i];
            // printf("k=%g j=%d  d=%g  imin=%d imax=%d  D(min)=%g+I%g  D(max)=%g+I%g  D=%g+I%g\n", k[i], j, d, i_min, i_max,
            //     real(D(i_min, j)), D(i_min, j).imag(), D(i_max, j).real(), D(i_max, j).imag(), Dave.real(), Dave.imag());
            // std::cout<< "D ="<< D(i, j)<<"  pk" << Pk(p, d)<<"  rho"  << rho(p, E3_m)  << std::endl;
            L(i) -= Dave * Pk(p, d) * rho(p, E3_m);
        }
    }
    return L;
}


std::complex<double> comput_Finf(double E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {
    std::complex<double> F = std::complex<double>(0, 0);
    Eigen::VectorXcd L = cal_L(E3_m, D, N, Npar, P, kcot, eps);
    double kmax = sqrt(pow((E3_m * E3_m + 1) / (2 * E3_m), 2) - 1);
    double d = kmax / N;
    for (int i = 0;i < N;i++) {
        double k = i * d;
        // std::cout << i << "   " << Pk(k, d) << "   " << L(i) << "   " << rho(k, E3_m) << std::endl;
        F += Pk(k, d) * L(i) * rho(k, E3_m);
    }
    return F;
}

std::complex<double> compute_M3_sym(double E3_m, int N, int Npar, double* P,
    double kcot(int, double*, int, double*), double* PKiso, double compute_kiso(double, double*), double eps) {

    std::complex<double> kiso = compute_kiso(E3_m, PKiso);
    Eigen::MatrixXcd D = compute_D(E3_m, N, Npar, P, kcot, eps);
    std::complex<double> Finf = comput_Finf(E3_m, D, N, Npar, P, kcot, eps);
    Eigen::VectorXcd  L = cal_L1(E3_m, D, N, Npar, P, kcot, eps);

    std::complex<double> M3 = std::complex<double>(0, 0);

    for (int i = 0;i < 3; i++) {
        for (int j = 0;j < 3; j++) {
            // std::cout << j <<"   " << L(j)<<"   " << (1. / (1. / kiso + Finf))<<   std::endl;
            M3 += L(i) * (1. / (1. / kiso + Finf)) * L(j);
            // M3 += (1. / (1. / kiso + Finf)) ;
        }
    }
    // printf("M3=%g %g     %g  %g     %g   %g\n",M3.real(),M3.imag(), kiso.real(), kiso.imag(),Finf.real(),Finf.imag());
    return M3;

}

