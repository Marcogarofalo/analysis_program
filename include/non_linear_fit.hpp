#ifndef non_linear_fit_H
#define non_linear_fit_H


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <vector>
#include "global.hpp"
#include "fit_all.hpp"

struct h_deriv {
    std::vector<double> h;
    // cast h_deriv as std::vector<double>
    operator std::vector<double>() const { return h; }
    // copy constructor
    h_deriv(const h_deriv& rhs) {
        h = rhs.h;
    }
    // move constructor
    h_deriv(h_deriv&& rhs) {
        h = std::move(rhs.h);
    }
    // copy ass
    h_deriv& operator = (const h_deriv& rhs) {
        h = rhs.h;
        return *this;
    }
    // move ass
    h_deriv& operator = (h_deriv&& rhs) noexcept {
        h = std::move(rhs.h);
        return *this;
    }

    h_deriv(const double& rhs) {
        h = std::vector<double>(1);
        h[0] = rhs;
    }
    h_deriv& operator = (const double& rhs) {
        h = std::vector<double>(1);
        h[0] = rhs;
        return *this;
    }
    h_deriv& operator = (const double&& rhs) {
        h = std::vector<double>(1);
        h[0] = rhs;
        return *this;
    }

    h_deriv& operator = (const std::vector<double>& rhs) {
        h = rhs;
        return *this;
    }
};

struct  fit_type {
    double (*function)(int, int, double*, int, double*);//N, Nvar, x ,Npar,P
    double (*f1)(int, int, double*, int, double*);
    double (*f2)(int, int, double*, int, double*);
    double** ext_P;  //parameter to not fit, they will be add to the Nvar,  Nvar=variables+(prameter to not fit)

    int N, Npar, Nvar, Njack;
    int n;// an index to be passed to (*function)( iN ,int,double*,int,double*)

    int n_ext_P = 0; //number of external parameter that will no be fitted
    int custom = false; // false=0 means default fit , 1 custom fit options
    double lambda = 0.001;
    double acc = 0.001;
    h_deriv h = 1e-5;
    int devorder = 4;  // 2 , 4 , -2 adaptive h=h*parameter
    int repeat_start = 1;
    double chi2_gap_jackboot = 1;
    int guess_per_jack = 3;
    int precision_sum = 0;// 0 float, 1 float kahan,  >1 long double
    bool mean_only = false;
    bool unstable = false; // if true avoid thing that may return error
    bool noderiv = false;
    bool covariancey = false;
    bool second_deriv = false;

    std::vector<int> corr_id = {};
    std::vector<double>  Prange = {};
    std::vector<double>  guess = {};
    std::vector<int> myen = {};
    double*** x;
    bool allocated_x = false;

    // GEVP
    int t0_GEVP = 3;
    int value_or_vector = 0;  // eigenvalues=0 or eigenvalues =1
    bool GEVP_tpt0 = false;
    bool GEVP_swap_t_t0 = false;
    int GEVP_ignore_warning_after_t = 1000;

    // HANKEL
    int HENKEL_size = 1;

    int T;
    int L;
    bool plateaux_scan = false;
    FILE* f_plateaux_scan = NULL;
    std::string name_plateaux_scan;

    std::string resampling;
    int verbosity = 0;
    double mu = 0;
    //
    std::vector<double> band_range{};
    double** cov;
    double** cov1;
    void compute_cov_fit(char** argv, data_all gjack, double lhs_fun(int, int, int, data_all, struct fit_type), struct fit_type fit_info);
    void compute_cov1_fit();
    void restore_default();
    void malloc_x();

};

struct fit_result {
    int Njack;
    int Npar;
    double** P;// [par][jack]
    double* chi2;
    double*** C;
    char name[NAMESIZE];

};
EXTERN fit_type default_fit_info;

struct fit_result malloc_fit(struct  fit_type  fit_info);
void free_fit_result(struct  fit_type  fit_info, struct fit_result  out);

////////////////////////////////////////
struct non_linear_fit_result {
    double* P;// [par][jack]
    double chi2;
};

// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(Nvariables,variables[], Nparameters,parameters[])
//*funk(Nvariables,variables[], Nparameters,parameters[]) must return a vector[Nparameters] that contain the derivatives of fun respect to the parameters
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2
//double  *non_linear_fit(int ensemble,double **x, double **y ,int Nvar, int Npar, double fun(int,double*,int,double*) , double *funk(int,double*,int,double*) );
double* non_linear_fit(int ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, double*, int, double*));

double compute_chi_non_linear(int ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, double*, int, double*));
double compute_chi_non_linear_Nf_cov1(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), double** cov1);

//double  *non_linear_fit_Nf(int N, int *ensemble ,double **x, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess , double lambda=0.001, double acc=0.001, double h=1e-5, std::vector<double>  Prange={},  int devorder=4 , int verbosity=1, int precision_sum=0);

non_linear_fit_result non_linear_fit_Nf(int N, int* ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess, fit_type fit_info = default_fit_info);
double* non_linear_fit_Nf_cov(int N, int* ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess, fit_type fit_info, double** cov1);


double compute_chi_non_linear_Nf(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), fit_type fit_info = default_fit_info);
double compute_chi_non_linear_Nf_kahan(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), fit_type fit_info = default_fit_info);
double** covariance_non_linear_fit_Nf(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), fit_type fit_info = default_fit_info);

//double  *guess_for_non_linear_fit_Nf(int N, int *ensemble ,double **x, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess , int repeat=1, double lambda=0.001, double acc=0.001, double h=1e-5, std::vector<double>  Prange={},  int devorder=4 , int verbosity=1, int precision_sum=0);
double* guess_for_non_linear_fit_Nf(int N, int* ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess, fit_type fit_info = default_fit_info);
double* guess_for_non_linear_fit_Nf_covariance(int N, int* ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess);

double* non_linear_fit_Nf_sigmax(int N, int* ensemble, double** x, double** sigmax, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess);
double* non_linear_fit_Nf_sigmax_iterative(int N, int* ensemble, double** x, double** sigmax, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess);
double* non_linear_fit_Nf_covariance(int N, int* ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess, double** cov, int devorder = 4);
double* non_linear_fit_Nf_sigmax_covariance(int N, int* ensemble, double** x, double** sigmax, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess, double** cov);

double rtsafe(void (*funcd)(double, int, double*, double*, double*), int Npar, double* P, double x1, double x2, double xacc);//non funziona rotine di merda!
double rtbis(double (*func)(double, double, int, double*), double input, int Npar, double* P, double x1, double x2, double xacc);
double rtbis_func_eq_input(double (*func)(int, int, double*, int, double*), int n, int Nvar, double* x, int Npar, double* P, int ivar, double input, double x1, double x2, double xacc, int Pedanticness = 0);
double rtnewt(double (*func)(int, int, double*, int, double*), int n, int Nvar, double* x, int Npar, double* P, int ivar, double input, double xstart, double xmin, double xmax, float xacc, int JMAX = 1000, double h = 1e-5);
double  rtsafe(double (*func)(int, int, double*, int, double*), int n, int Nvar, double* x, int Npar, double* P, int ivar, double input, double x1, double x2, float xacc, int JMAX = 1000, double h = 1e-5);

double* der_O4_fun_Nf_h(int n, int Nvar, double* x, int Npar, double* P, double fun(int, int, double*, int, double*), std::vector<double> h);
double* derN_fun_Nf_var_h(int n, int Nvar, double* x, int Npar, double* P, double fun(int, int, double*, int, double*), double h, int N);

struct  fit_result   malloc_copy_fit_result(struct fit_result fit_out);
void  copy_fit_type_into(struct  fit_type* fit_tmp, struct fit_type fit_info);
struct fit_all   save_fit(struct fit_all fit_chi2_good, struct fit_type fit_info, struct fit_result fit_out);




#endif
