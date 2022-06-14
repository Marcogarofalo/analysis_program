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
    double k_mcotd = kcot(Nvar, x, Npar, P);
    // double k_m = sqrt(E_m * E_m / 4. - 1);
    double E_m = 2 * sqrt(k_m * k_m + 1);

    std::complex<double> r = 16 * M_PI * E_m;
    r /= k_mcotd - 1i * k_m;
    return r;
}

inline double J(double x) {
    if (x <= 0) return 0;
    if (x > 0 && x <= 1) return exp(-(1. / x) * exp(-(1. / (1. - x))));
    if (x > 1) return 1;
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
        std::complex<double> r = log((2 * pk - pow(E3_m - omega_k - omega_p, 2) + p * p + k * k + 1 - 1i * eps) /
            (-2 * pk - pow(E3_m - omega_k - omega_p, 2) + p * p + k * k + 1 - 1i * eps));
        // std::cout << "log=" << r << std::endl;
        r *= -H(k, p, E3_m) / (4 * pk);
        // std::cout << "H=" << -H(k, p, E3_m) / (4 * pk) << std::endl;
    }
    return r;
}

double Pk(double k, double d) {
    double omega_k = sqrt(k * k + 1);
    return k * k * d / (M_2_PI * M_2_PI * omega_k);
}


Eigen::MatrixXcd compute_D(double E3_m, int N, int Npar, double* P, double kcot(int, double*, int, double*), double d, double eps) {

    Eigen::MatrixXcd oneplusMGP = Eigen::MatrixXcd::Identity(N, N);
    Eigen::MatrixXcd MGM = Eigen::MatrixXcd::Zero(N, N);

    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[1];
    for (int i = 0; i < N;i++) {
        double p = d * i;
        x[0] = p;
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
    Eigen::MatrixXcd D = -1 * oneplusMGP1 * oneplusMGP;

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

Eigen::VectorXcd  cal_L(double E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double d, double eps) {
    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[1];
    Eigen::VectorXcd L(N);
    for (int i = 0;i < N;i++) {
        double k = i * d;
        L(i) = std::complex<double>(1. / 3., 0);
        x[0] = k;
        L(i) -= M2_from_kcot(Nvar, x, Npar, P, kcot) * rho(k, E3_m);
        // std::cout << "M2="<<M2_from_kcot(Nvar, x, Npar, P, kcot) << "   rho=" << rho(k, E3_m)<< std::endl;
        for (int j = 0;j < N;j++) {
            double p = j * d;
            // std::cout<< "D ="<< D(i, j)<<"  pk" << Pk(p, d)<<"  rho"  << rho(p, E3_m)  << std::endl;
            L(i) += D(i, j) * Pk(p, d) * rho(p, E3_m);
        }
    }
    return L;
}


Eigen::VectorXcd  cal_L1(double E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double d, double eps) {
    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[1];
    Eigen::VectorXcd L(3);
    int i_min = 0.1 / d;
    int i_max = 0.1 / d + 1;
    if (i_max>=N) {printf("error N too small\n"); exit(-1);} 
    for (int i = 0;i < 3;i++) {
        double k = fabs(0.1 * (i - 1));//1,0,1
        L(i) = std::complex<double>(1. / 3., 0);
        x[0] = k;
        L(i) -= M2_from_kcot(Nvar, x, Npar, P, kcot) * rho(k, E3_m);
        // std::cout << "M2="<<M2_from_kcot(Nvar, x, Npar, P, kcot) << "   rho=" << rho(k, E3_m)<< std::endl;
        for (int j = 0;j < N;j++) {
            double p = j * d;

            std::complex<double> Dave = (D(i_min, j) + D(i_max, j)) / 2.;
            // std::cout<< "D ="<< D(i, j)<<"  pk" << Pk(p, d)<<"  rho"  << rho(p, E3_m)  << std::endl;
            L(i) += Dave * Pk(p, d) * rho(p, E3_m);
        }
    }
    return L;
}


std::complex<double> Finf(double E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double d, double eps) {
    std::complex<double> F = std::complex<double>(0, 0);
    Eigen::VectorXcd L = cal_L(E3_m, D, N, Npar, P, kcot, d, eps);
    for (int i = 0;i < N;i++) {
        double k = i * d;
        F += Pk(k, d) * L(i) * rho(k, E3_m);
    }
    return F;
}

std::complex<double> compute_M3_sym(double E3_m, int N, int Npar, double* P,
    double kcot(int, double*, int, double*), double* PKiso, std::complex<double> compute_kiso(double, double*), double d, double eps) {

    std::complex<double> kiso = compute_kiso(E3_m, PKiso);
    Eigen::MatrixXcd D = compute_D(E3_m, N, Npar, P, kcot, d, eps);
    std::complex<double> Finf = (E3_m, D, N, Npar, P, kcot, d, eps);
    Eigen::VectorXcd  L = cal_L1(E3_m, D, N, Npar, P, kcot, d, eps);

    std::complex<double> M3 = std::complex<double>(0, 0);

    for (int i = 0;i < 3; i++) {
        for (int j = 0;j < 3; j++) {
            // std::cout << j <<"   " << L(j)<<"   " << (1. / (1. / kiso + Finf))<<   std::endl;
            M3 += L(i) * (1. / (1. / kiso + Finf)) * L(j);
        }
    }
    return M3;

}

