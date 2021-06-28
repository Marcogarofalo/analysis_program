#ifndef linear_fit_H
#define linear_fit_H


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

double **cholesky_decomposition(double **a, int n);
double *cholesky_solver(int n,double **a,double *b);
double *cholesky_solver_if_possible(int n,double **a,double *b);
double **symmetric_matrix_inverse(int N, double **M  );

//return the vector x: the solution of Mx=b
//M is a matrix 
//it uses the LU decomposition see Numerical Receipes
double *LU_decomposition_solver(int N, double **M, double *b);
//return the inverse matrix of M
double **matrix_inverse(int N, double **M  );

double inter_spline(double at, int Npoints, double *x, double *y);
    

//fit N data y[N][value][error]  with the function sum_{i=0;i<M}  a_i f_i(x)
//*fit_functions is a function such that fit_function(int M,double x)[i]=f_i
//it returns a[i][value,error]
double  **linear_fit(int N, double *x, double **y,  int M, double *fit_function(int,double) );
//it needs **a: the result of **linear_fit 
double  compute_chisqr(int N, double *x, double **y,  int M, double **a, double *fit_function(int,double) );
double *constant_fit_to_try(int M, double in);
double *give_jack_linear_fit(int tmin, int tmax, int sep ,double **corr_ave, double **corr_J,int Njack );
double *try_linear_fit(char **option,int tmin, int tmax, int sep ,double **corr_ave, double **corr_J,int Njack ,double **chi2);
double  **global_linear_fit(int N, double **x, double **y,  int M );
double  compute_chisqr_global_fit(int N, double **x, double **y,  int M, double **a );

double  **linear_fit_Nvar(int N, double **x, double **y,  int M, double *fit_function(int,double*) );
double  compute_chisqr_Nvar(int N, double **x, double **y,  int M, double **a, double *fit_function(int,double*) );

double  *linear_fit_Nf(int Nfunc, int *Npoints, double **x, double **y,  int Nvar,int Npar, double *fit_function(int,int,double*,int)  );
double compute_chi_linear_Nf(int N,int *ensemble,double **x, double **y, double *P ,int Nvar, int Npar,  double *fun(int,int,double*,int)) ;

int is_it_positive(double **a, int n);

#endif
