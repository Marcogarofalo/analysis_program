#ifndef common_integral_eq_H
#define common_integral_eq_QC3_H
using namespace std::complex_literals;

data_all setup_data_for_fits(int NE, int Njack, double*** M, double*** F);
void write_M3(int NE, int Njack, std::vector<double>& E3, double*** M3, double*** F, std::string namef, const char* prefix);
void read_M3(int& NE, int& Njack, std::vector<double>& E3, double***& M3, double***& F, std::string namef, const char* prefix);
void compute_M3(int NE, double Emin, double dE, int Njack, std::vector<double>& E3, double***& M3, double***& F, int N, int Npar,
    double** P, double compute_kcot(int, double*, int, double*), double** PKiso, double compute_kiso(double, double*),
    double eps);

double rhs_laurent_pole(int n, int Nvar, double* x, int Npar, double* P);
double rhs_BW(int n, int Nvar, double* x, int Npar, double* P);
double rhs_absBW(int n, int Nvar, double* x, int Npar, double* P);
double rhs_F(int n, int Nvar, double* x, int Npar, double* P);
double lhs_M3(int n, int e, int j, data_all gjack, struct fit_type fit_info);
double lhs_absM3(int n, int e, int j, data_all gjack, struct fit_type fit_info);

double lhs_F(int n, int e, int j, data_all gjack, struct fit_type fit_info);
double compute_kiso(double E3_m, double* P);
double denom_M(int n, int Nvar, double* x, int Npar, double* P);
double denom_M_direct(int n, int Nvar, double* x, int Npar, double* P);
double compute_kcot(int Nvar, double* x, int Npar, double* P);

#endif
