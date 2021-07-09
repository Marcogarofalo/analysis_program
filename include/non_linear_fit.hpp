#ifndef non_linear_fit_H
#define non_linear_fit_H


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <vector>


struct fit_result malloc_fit( struct  fit_type  fit_info);
void free_fit_result( struct  fit_type  fit_info,struct fit_result  out);

// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(Nvariables,variables[], Nparameters,parameters[])
//*funk(Nvariables,variables[], Nparameters,parameters[]) must return a vector[Nparameters] that contain the derivatives of fun respect to the parameters
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2
//double  *non_linear_fit(int ensemble,double **x, double **y ,int Nvar, int Npar, double fun(int,double*,int,double*) , double *funk(int,double*,int,double*) );
double  *non_linear_fit(int ensemble,double **x, double **y ,int Nvar, int Npar, double fun(int,double*,int,double*)  );

double compute_chi_non_linear(int ensemble,double **x, double **y, double *P ,int Nvar, int Npar, double fun(int,double*,int,double*));


double  *non_linear_fit_Nf(int N, int *ensemble ,double **x, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess , double lambda=0.001, double acc=0.001, double h=1e-5, std::vector<double>  Prange={},  int devorder=4 , int verbosity=1, int precision_sum=0);



double compute_chi_non_linear_Nf(int N,int *ensemble,double **x, double **y, double *P ,int Nvar, int Npar,  double fun(int,int,double*,int,double*)) ;
double compute_chi_non_linear_Nf_kahan(int N,int *ensemble,double **x, double **y, double *P ,int Nvar, int Npar,  double fun(int,int,double*,int,double*)) ;
double  **covariance_non_linear_fit_Nf(int N, int *ensemble ,double **x, double **y,double *P ,int Nvar, int Npar,  double fun(int,int,double*,int,double*)  );

double  *guess_for_non_linear_fit_Nf(int N, int *ensemble ,double **x, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess , int repeat=1, double lambda=0.001, double acc=0.001, double h=1e-5, std::vector<double>  Prange={},  int devorder=4 , int verbosity=1, int precision_sum=0);
double  *guess_for_non_linear_fit_Nf_covariance(int N, int *ensemble ,double **x, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess );

double  *non_linear_fit_Nf_sigmax(int N, int *ensemble ,double **x, double **sigmax, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess );
double  *non_linear_fit_Nf_sigmax_iterative(int N, int *ensemble ,double **x, double **sigmax, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess );
double  *non_linear_fit_Nf_covariance(int N, int *ensemble ,double **x, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess, double **cov ,int devorder=4);
double  *non_linear_fit_Nf_sigmax_covariance(int N, int *ensemble ,double **x, double **sigmax, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess, double **cov );

double rtsafe(void (*funcd)(double,int,double*, double *, double *),int Npar,double *P , double x1, double x2,double xacc);//non funziona rotine di merda!
double rtbis(double (*func)(double,double,int,double*),double input,int Npar, double *P, double x1, double x2, double xacc);
double rtbis_func_eq_input(double (*func)(int , int , double*,int,double*),int n, int Nvar, double *x,int Npar, double *P, int ivar,double input, double x1, double x2, double xacc, int Pedanticness =0);
double rtnewt(double (*func)(int , int , double*,int,double*),int n, int Nvar, double *x,int Npar, double *P, int ivar,double input, double xstart, double xmin , double xmax,  float xacc, int JMAX=1000, double h=1e-5);
double  rtsafe(double (*func)(int , int , double*,int,double*), int n, int Nvar, double *x,int Npar, double *P, int ivar,double input,   double x1, double x2 ,  float xacc, int JMAX=1000, double h=1e-5);

double *der_O4_fun_Nf_h(int n, int Nvar, double *x,int Npar,double  *P, double fun(int,int,double*,int,double*), double h);
double *derN_fun_Nf_var_h(int n, int Nvar, double *x,int Npar,double  *P, double fun(int,int,double*,int,double*), double h,int N);

struct fit_all   save_fit(struct fit_all fit_chi2_good,struct fit_type fit_info, struct fit_result fit_out);
#endif
