#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <new>
#include <utility>
#include <Eigen/Dense>  
#include "integral_eq_QC3.hpp"

using namespace std::complex_literals;

inline double compute_kmax(double E3_m) {
    double xmin = 0.02;
    return sqrt(pow(((E3_m * E3_m + 1) / (4.) - xmin) * 2 / E3_m, 2) - 1);
}
inline double compute_d(double kmax, int N) {
    return kmax / ((double)N - 1.);
}


double kcot_generic(int Nvar, double* x, int Npar, double* P) {
    double r;
    r = 1. / P[0];

    return r;
}

std::complex<double> M2_from_kcot(int Nvar, double* x, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {
    double k_m = x[0];
    double E3_m = x[1];
    double k_mcotd = kcot(Nvar, x, Npar, P);
    // double k_m = sqrt(E_m * E_m / 4. - 1);
    double omega_k = sqrt(k_m * k_m + 1);
    std::complex<double> Es_mk2 = pow(E3_m - omega_k, 2) - k_m * k_m + 1i * eps;

    // if (Es_mk2 >= 4)
    //     r *= -1i * sqrt((Es_mk2 + 1i * eps) / 4. - 1.);
    // else
    //     r *= sqrt(fabs((Es_mk2 + 1i * eps) / 4. - 1.));
    std::complex<double> q = sqrt(Es_mk2 / 4. - 1.);

    return 16. * M_PI * sqrt(Es_mk2) / (k_mcotd - 1i * q);
}

std::complex<double> M2_from_kcot_complex(int Nvar, double* x, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {
    double k_m = x[0];

    std::complex<double> E3_m_c(x[1], x[2]);

    double k_mcotd = kcot(Nvar, x, Npar, P);
    // double k_m = sqrt(E_m * E_m / 4. - 1);
    double omega_k = sqrt(k_m * k_m + 1);
    std::complex<double> Es_mk2 = pow(E3_m_c - omega_k, 2) - k_m * k_m + 1i * eps;

    double iq = (Es_mk2 / 4. - 1.).imag();
    std::complex<double> q = sqrt(Es_mk2 / 4. - 1.);
    if (iq < 0)  q = -q;
    // std::cout << "Es_mk2 / 4. - 1.  ="<<Es_mk2 / 4. - 1.<<"\n";
    // std::cout << "q  ="<<q<<"\n";
    // q =-conj(q);

    // std::cout << "-conj(q)  ="<<q<<"\n";
    double iE = (Es_mk2).imag();
    std::complex<double> sE = sqrt(Es_mk2);
    // if (iE<0)  sE= -sE;
    // sE=-conj(sE);
    return 16. * M_PI * sE / (k_mcotd - 1i * q);

}

inline double J(double x) {
    double xmin = 0.02 - 1e-12;// fernando set them to these values
    double xmax = 0.97;

    if (x < xmin) return 0;
    else if (x >= xmin && x < xmax) return exp(-(1. / x) * exp(-(1. / (1. - x))));
    else if (x >= xmax) return 1;
    else {
        printf("error in J(x): x<0");
        exit(-1);
    }
}

inline std::complex<double> J(std::complex<double> x) {
    double xmin = 0.02 - 1e-12;// fernando set them to these values
    double xmax = 0.97;


    if (real(x) < xmin) return 0;
    else if (real(x) >= xmin && real(x) < xmax) return exp(-(1. / x) * exp(-(1. / (1. - x))));
    else if (real(x) >= xmax) return 1;
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
    // printf("call H (%.12g, %.12g)=%.12g  %.12g",Es_mk2,Es_mp2,J(Es_mk2 / 4.) , J(Es_mp2 / 4.));
    return J(Es_mk2 / 4.) * J(Es_mp2 / 4.);

}
std::complex<double>  H(double k, double p, std::complex<double> E3_m) {

    double omega_k = sqrt(k * k + 1);
    std::complex<double> Es_mk2 = pow(E3_m - omega_k, 2) - k * k;
    double omega_p = sqrt(p * p + 1);
    std::complex<double> Es_mp2 = pow(E3_m - omega_p, 2) - p * p;
    return J(Es_mk2 / 4.) * J(Es_mp2 / 4.);

}

std::complex<double> Gs_pke(double k, double p, double E3_m, double d, double eps) {

    double omega_k = sqrt(k * k + 1);
    double omega_p = sqrt(p * p + 1);
    std::complex<double> r;

    double th = d / 2.;
    if (p < th) {
        r = H(k, p, E3_m) / (1 + E3_m * E3_m + 1i * eps + 2 * omega_k - 2 * E3_m * (1 + omega_k));
    }
    else if (k < th) {
        r = H(k, p, E3_m) / (1 + E3_m * E3_m + 1i * eps + 2 * omega_p - 2 * E3_m * (1 + omega_p));
    }
    else {
        // r = log((2 * pk - pow(E3_m - omega_k - omega_p, 2) + p * p + k * k + 1 - 1i * eps) /
        //     (-2 * pk - pow(E3_m - omega_k - omega_p, 2) + p * p + k * k + 1 - 1i * eps));
        // r *= -H(k, p, E3_m) / (4 * pk);
        double alpha = pow(E3_m - omega_k - omega_p, 2) - p * p - k * k - 1.;
        r = -1 * H(k, p, E3_m) / (4 * p * k) * log((alpha - 2 * p * k + 1i * eps) / (alpha + 2 * p * k + 1i * eps));
        // printf("alpha=%g   H=%g\n",alpha,H(k, p, E3_m));
        // std::cout << "log=" << r << std::endl;

    // std::cout << "H=" << -H(k, p, E3_m) / (4 * pk) << std::endl;
    // std::cout << "r=" << r << std::endl;

    }
    return r;
}
std::complex<double> Gs_pke(double k, double p, std::complex<double> E3_m, double d, double eps) {

    double omega_k = sqrt(k * k + 1);
    double omega_p = sqrt(p * p + 1);
    std::complex<double> r;

    double th = d / 2.;
    if (p < th) {
        r = H(k, p, E3_m) / (1.0 + E3_m * E3_m + 1i * eps + 2. * omega_k - 2. * E3_m * (1. + omega_k));
    }
    else if (k < th) {
        r = H(k, p, E3_m) / (1.0 + E3_m * E3_m + 1i * eps + 2. * omega_p - 2. * E3_m * (1. + omega_p));
    }
    else {
        std::complex<double> alpha = pow(E3_m - omega_k - omega_p, 2) - p * p - k * k - 1.;
        r = -1. * H(k, p, E3_m) / (4 * p * k) * log((alpha - 2 * p * k + 1i * eps) / (alpha + 2 * p * k + 1i * eps));
        // std::cout<< "   a=" << (alpha - 2 * p * k + 1i * eps) / (alpha + 2 * p * k + 1i * eps)<< std::endl;
        // std::cout<< " arg=" << std::arg((alpha - 2 * p * k + 1i * eps) / (alpha + 2 * p * k + 1i * eps))<< std::endl;
        // std::cout<< " abs=" << std::abs((alpha - 2 * p * k + 1i * eps) / (alpha + 2 * p * k + 1i * eps))<< std::endl;
        // std::cout<< "loga=" << log((alpha - 2 * p * k + 1i * eps) / (alpha + 2 * p * k + 1i * eps))<< std::endl;
    }
    return r;
}


inline double Pk(double k, double d) {
    double omega_k = sqrt(k * k + 1);
    return k * k * d / (M_PI * M_PI * 2 * 2 * omega_k);
}


Eigen::MatrixXcd compute_D(double E3_m, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {

    Eigen::MatrixXcd oneplusMGP = Eigen::MatrixXcd::Identity(N, N);
    Eigen::MatrixXcd MGM = Eigen::MatrixXcd::Zero(N, N);

    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[2];
    double kmax = compute_kmax(E3_m);
    double d = compute_d(kmax, N);
    for (int i = 0; i < N;i++) {
        double p = d * i;
        x[0] = p;
        x[1] = E3_m;
        std::complex<double> m2 = M2_from_kcot(Nvar, x, Npar, P, kcot, eps);
        // std::cout << "m2=" << m2 << std::endl;
        for (int j = 0; j < N;j++) {
            double k = d * j;
            std::complex<double> G = Gs_pke(p, k, E3_m, d, eps);
            // std::cout << "G=" << G << std::endl;

            oneplusMGP(i, j) += m2 * G * Pk(k, d);
            x[0] = k;
            MGM(i, j) += m2 * G * M2_from_kcot(Nvar, x, Npar, P, kcot, eps);
            // std::cout << "m2(p)=" << M2_from_kcot(Nvar, x, Npar, P, kcot) << std::endl;
        }

        // for (int j = 0; j < N;j++) {
        //     double k = d * j;
        //     std::complex<double> G = Gs_pke(p, k, E3_m, d, eps);
        //     printf("1pMGP(%d,%d)=(%.12g,  %.12g)  m2=(%.12g,  %.12g)  G=(%.12g,  %.12g) pk=%g d=%g k=%g\n", i, j, oneplusMGP(i, j).real(), oneplusMGP(i, j).imag(), m2.real(), m2.imag(),G.real(), G.imag(),Pk(k, d) , d, k);
        // }
    }
    // exit(666);
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


Eigen::MatrixXcd compute_D(std::complex<double> E3_m, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {

    Eigen::MatrixXcd oneplusMGP = Eigen::MatrixXcd::Identity(N, N);
    Eigen::MatrixXcd MGM = Eigen::MatrixXcd::Zero(N, N);

    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[3];
    double kmax = compute_kmax(real(E3_m));
    double d = compute_d(kmax, N);
    for (int i = 0; i < N;i++) {
        double p = d * i;
        x[0] = p;
        x[1] = E3_m.real();
        x[2] = E3_m.imag();
        std::complex<double> m2 = M2_from_kcot_complex(Nvar, x, Npar, P, kcot, eps);

        for (int j = 0; j < N;j++) {
            double k = d * j;
            std::complex<double> G = Gs_pke(p, k, E3_m, d, eps);
            oneplusMGP(i, j) += m2 * G * Pk(k, d);
            x[0] = k;
            MGM(i, j) += m2 * G * M2_from_kcot_complex(Nvar, x, Npar, P, kcot, eps);
        }
    }
    Eigen::MatrixXcd oneplusMGP1 = oneplusMGP.inverse();
    Eigen::MatrixXcd D = -1 * oneplusMGP1 * MGM;
    return D;
}

std::complex<double> rho(double k_m, double E3_m) {
    double omega_k = sqrt(k_m * k_m + 1);
    long double Es_mk2 = pow(E3_m - omega_k, 2) - k_m * k_m;
    std::complex<double> r = J(Es_mk2 / 4.) / (16. * M_PI * sqrt(Es_mk2));
    if (Es_mk2 >= 4)
        r *= -1i * ((double)sqrt(Es_mk2 / 4. - 1));
    else
        r *= sqrt(-(Es_mk2 / 4. - 1));
    // printf("rho= %.12g     %.12Lg   %.12g \n",k_m, Es_mk2 / 4., J(Es_mk2 / 4.) );
    return r;
}


std::complex<double> rho(double k_m, std::complex<double> E3_m) {
    double omega_k = sqrt(k_m * k_m + 1);
    std::complex<double> Es_mk2 = pow(E3_m - omega_k, 2) - k_m * k_m;
    std::complex<double> r = J(Es_mk2 / 4.) / (16. * M_PI * sqrt(Es_mk2));
    // if (Es_mk2.imag() < 0) r *= -1.0;


    double a = (Es_mk2 / 4. - 1.0).imag();
    std::complex<double> rhoc = sqrt(Es_mk2 / 4. - 1.0);
    if (a < 0)  rhoc = -(rhoc);

    // std::complex<double> rho1 = -1i * sqrt(Es_mk2 / 4. - 1.0);
    r *= -1i * rhoc;
    // if (abs(rhoc - rho1) > 1e-4) {
    //     printf("\n\n different  E=(%g,%g) rhoc=(%g,%g)   rho1=(%g,%g) \n\n", Es_mk2.real(), Es_mk2.imag(), rhoc.real(), rhoc.imag(), rho1.real(), rho1.imag());
    // }
    return r;
}

Eigen::VectorXcd  cal_L(double E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {
    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[2];
    typedef Eigen::Matrix< std::complex<  long double >, Eigen::Dynamic, 1 > VectorXcld;
    VectorXcld Ld(N);
    Eigen::VectorXcd L(N);
    double kmax = compute_kmax(E3_m);
    double d = compute_d(kmax, N);
    for (int i = 0;i < N;i++) {
        double k = i * d;
        Ld(i) = std::complex<long double>(1. / 3., 0);
        x[0] = k;
        x[1] = E3_m;
        Ld(i) -= M2_from_kcot(Nvar, x, Npar, P, kcot, 0) * rho(k, E3_m);
        // if (i == 100  ) {
        //     printf("M2=(%.12g,%.12g)   rho=(%.12g,%.12g)\n",  M2_from_kcot(Nvar, x, Npar, P, kcot, 0.0).real(),
        //              M2_from_kcot(Nvar, x, Npar, P, kcot, 0.0).imag(), rho(k, E3_m).real(), rho(k, E3_m).imag());
        // }
        // std::cout << "M2="<<M2_from_kcot(Nvar, x, Npar, P, kcot) << "   rho=" << rho(k, E3_m)<< std::endl;
        for (int j = 0;j < N;j++) {
            double p = j * d;
            // if (i == 100  ) {
            //     printf("D(%d,%d)=(%.12g,%.12g)  pk=%.12g  rho=(%.12g,%.12g) L=(%.12Lf, %.12Lf)\n", i,j,D(i, j).real(),
            //         D(i, j).imag(), Pk(p, d), rho(p, E3_m).real(), rho(p, E3_m).imag(),Ld(i).real(),Ld(i).imag() );
            // }
            Ld(i) -= D(i, j) * Pk(p, d) * rho(p, E3_m);
        }
        L(i) = Ld(i);
        // printf("L(%d)=(%.12g  %.12g)\n",i,L[i].real(),L[i].imag());
    }

    return L;
}
Eigen::VectorXcd  cal_L(std::complex<double> E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {
    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[3];
    typedef Eigen::Matrix< std::complex<  long double >, Eigen::Dynamic, 1 > VectorXcld;
    VectorXcld Ld(N);
    Eigen::VectorXcd L(N);
    double kmax = compute_kmax(real(E3_m));
    double d = compute_d(kmax, N);

    for (int i = 0;i < N;i++) {
        double k = i * d;
        Ld(i) = std::complex<long double>(1. / 3., 0);
        x[0] = k;
        x[1] = E3_m.real();
        x[2] = E3_m.imag();
        Ld(i) -= M2_from_kcot_complex(Nvar, x, Npar, P, kcot, 0) * rho(k, E3_m);

        for (int j = 0;j < N;j++) {
            double p = j * d;
            Ld(i) -= D(i, j) * Pk(p, d) * rho(p, E3_m);
        }
        L(i) = Ld(i);

    }

    return L;
}


Eigen::VectorXcd  cal_L1(double E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {
    int Nvar = 1; // kcotd can have only one variable, k_m
    double x[2];
    Eigen::VectorXcd L(3);

    double k[3];
    k[0] = 0;
    k[1] = sqrt(pow(E3_m - 1, 2) / 4.0 - 1);
    k[2] = sqrt(pow(E3_m - 1, 2) / 4.0 - 1);// this with a minus but the modulus is the same
    double kmax = compute_kmax(E3_m);
    double d = compute_d(kmax, N);
    for (int i = 0;i < 3;i++) {

        L(i) = std::complex<double>(1. / 3., 0);
        x[0] = k[i];
        x[1] = E3_m;
        L(i) -= M2_from_kcot(Nvar, x, Npar, P, kcot, 0) * rho(k[i], E3_m);

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
    double kmax = compute_kmax(E3_m);
    double d = compute_d(kmax, N);
    for (int i = 0;i < N;i++) {
        double k = i * d;

        F += Pk(k, d) * L(i) * rho(k, E3_m);
        // if (i == N - 1) printf("%-18.12g%-18.12g%-18.12g%-18.12g%-18.12g%-18.12g\n", k, rho(k, E3_m).real(), rho(k, E3_m).imag(), L(i).real(), L(i).imag(), Pk(k, d));

    }

    return F;
}


std::complex<double> comput_Finf(std::complex<double> E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps) {
    std::complex<double> F = std::complex<double>(0, 0);
    Eigen::VectorXcd L = cal_L(E3_m, D, N, Npar, P, kcot, eps);
    double kmax = compute_kmax(real(E3_m));
    double d = compute_d(kmax, N);
    for (int i = 0;i < N;i++) {
        double k = i * d;
        F += Pk(k, d) * L(i) * rho(k, E3_m);
        // if (i == N - 1) printf("%-18.12g%-18.12g%-18.12g%-18.12g%-18.12g%-18.12g\n", k, rho(k, E3_m).real(), rho(k, E3_m).imag(), L(i).real(), L(i).imag(), Pk(k, d));
    }
    if (E3_m.imag() < 0) return conj(F);
    return (F);
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

