#ifndef integral_eq_QC3_H
#define integral_eq_QC3_H

#include <Eigen/Dense>  


typedef Eigen::Matrix< std::complex<  long double >, Eigen::Dynamic, Eigen::Dynamic > MatrixXcld;
typedef Eigen::Matrix< std::complex<  long double >, Eigen::Dynamic, 1 > VectorXcld;


std::complex<double> compute_M3_sym(double E3_m, int N, int Npar, double* P,
    double kcot(int, double*, int, double*), double* PKiso, double compute_kiso(double, double*), double eps);
Eigen::MatrixXcd compute_D(double E3_m, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps);
std::complex<double> comput_Finf(double E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P,
    double kcot(int, double*, int, double*), double eps);
Eigen::VectorXcd  cal_L1(double E3_m, Eigen::MatrixXcd D, int N, int Npar, double* P, double kcot(int, double*, int, double*), double eps);

#endif