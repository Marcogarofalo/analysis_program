#ifndef extra_fun_phi4_H
#define extra_fun_phi4_H


#include "fit_all.hpp"
#include <cstdlib>
#ifdef PYTHON
#include "QC3_interface.hpp"
#endif
// #include "mass_phi4.hpp"

void print_fit_band_phi4(char** argv, data_all gjack, struct fit_type fit_info,
    struct fit_type fit_info_m0, const char* label, const char* dir_name,
    struct fit_result fit_out, struct fit_result fit_out_m0, int var, int en, double h);
// QC2
double rhs_kcotd_m_new(int n, int Nvar, double* x, int Npar, double* P);
double lhs_kcotd_m_deltaE_g_new(int n, int e, int j, data_all gjack, struct fit_type fit_info);
double compute_k_m_g_new(int n, int e, int j, data_all gjack, struct fit_type fit_info);

// QC3
double lhs_E3orE1_m_complex_new(int n, int e, int j, data_all gjack, struct fit_type fit_info);
#ifdef PYTHON
double rhs_E3_m_QC3_pole_new(int n, int Nvar, double* x, int Npar, double* P);


void print_fit_band_QC3_phi4(char** argv, data_all gjack, struct fit_type fit_info,
    struct fit_type fit_info_E3_poly, const char* label, const char* dir_name,
    struct fit_result fit_out, struct fit_result fit_out_poly, int var, int en, double h) ;
#endif //PYTHON

#endif // !extra_fun_phi4_H