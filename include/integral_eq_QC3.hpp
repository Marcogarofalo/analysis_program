#ifndef integral_eq_QC3_H
#define integral_eq_QC3_H

#include <Eigen/Dense>  


std::complex<double> compute_M3_sym(double E3_m, int N, int Npar, double* P,
    double kcot(int, double*, int, double*), double* PKiso, std::complex<double> compute_kiso(double, double*), double d, double eps);
#endif