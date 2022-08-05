#define linear_fit_C


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "resampling.hpp"
#include "linear_fit.hpp"
#include "tower.hpp"

//Given a positive-definite symmetric matrix a[1..n][1..n] , this routine constructs its Cholesky
//decomposition, A = L · L T . On input, only the upper triangle of a need be given; it is not
//modified. The Cholesky factor L is returned in the lower triangle of a , except for its diagonal
//elements which are returned in p[1..n] .
double** cholesky_decomposition(double** a, int n) {
    int i, j, k;
    double sum;
    double** L = (double**)malloc(sizeof(double*) * n);
    for (i = 0;i < n;i++)
        L[i] = (double*)calloc(n, sizeof(double));

    for (i = 0;i < n;i++) {
        for (j = i;j < n;j++) {
            sum = a[i][j];

            for (k = i - 1;k >= 0;k--) sum -= L[i][k] * L[j][k];
            if (i == j) {
                if (sum <= 0.0) {
                    printf("cholesky_decomposition failed: the matrix, with rounding errors, is not positive definite,\n here the matrix\n{");
                    printf("{");
                    for (int ii = 0;ii < n;ii++) {

                        for (int jj = 0;jj < ii;jj++)
                            printf("%.5g, ", a[jj][ii]);
                        for (int jj = ii;jj < n - 1;jj++)
                            printf("%.5g ,", a[ii][jj]);
                        printf("%.5g },\n", a[ii][n - 1]);

                    }
                    printf("}\n");
                    /*printf("returning inverse of diag{matrix}");
                    double **L1;
                    for (int ii=0;ii<n;ii++){
                        for (int jj=0;jj<n;jj++)
                            L[ii][jj]=0;

                        L[ii][ii]=a[ii][ii];
                    }
                    L1=cholesky_decomposition(L,  n);
                    return L1;*/
                    exit(0);
                }
                L[i][i] = sqrt(sum);
            }
            else L[j][i] = sum / L[i][i];
        }
    }


    return L;
}


/*check if a matrix is positive define*/
int is_it_positive(double** a, int n) {
    int i, j, k;
    double sum;
    double** L = (double**)malloc(sizeof(double*) * n);
    for (i = 0;i < n;i++)
        L[i] = (double*)calloc(n, sizeof(double));

    for (i = 0;i < n;i++) {
        for (j = i;j < n;j++) {
            sum = a[i][j];

            for (k = i - 1;k >= 0;k--) sum -= L[i][k] * L[j][k];
            if (i == j) {
                if (sum <= 0.0) {
                    /*
                     printf("{");
                     for (int ii=0;ii<n;ii++) {
                            printf("{");
                         for (int jj=0;jj<ii;jj++)
                                 printf("%.10f, ",a[jj][ii]);
                         for (int jj=ii;jj<n-1;jj++)
                                 printf("%.10f ,",a[ii][jj]);
                         printf("%.10f },\n",a[ii][n-1]);

                     }
                     printf("}\n");
                     */
                    return 1;

                }
                L[i][i] = sqrt(sum);
            }
            else L[j][i] = sum / L[i][i];
        }
    }


    return 0;
}

bool is_it_positive_lex_reim(double** a, int n) {
    int i, j, k;
    double sum;
    double** L = (double**)malloc(sizeof(double*) * n);
    for (i = 0;i < n;i++)
        L[i] = (double*)calloc(n, sizeof(double));

    for (i = 0;i < n;i++) {
        for (j = i;j < n;j++) {
            sum = a[i + j * n][0];

            for (k = i - 1;k >= 0;k--) sum -= L[i][k] * L[j][k];
            if (i == j) {
                if (sum <= 0.0) {

                    return false;

                }
                L[i][i] = sqrt(sum);
            }
            else L[j][i] = sum / L[i][i];
        }
    }


    return true;
}

/*Solves the set of n linear equations A · x = b, where a is a positive-definite symmetric matrix.
a[1..n][1..n] and p[1..n] are input as the output of the routine choldc . Only the lower
subdiagonal portion of a is accessed. b[1..n] is input as the right-hand side vector. The
solution vector is returned in x[1..n] . a , n , and p are not modified and can be left in place
for successive calls with different right-hand sides b . b is not modified unless you identify b and
x in the calling sequence, which is allowed.*/
double* cholesky_solver(int n, double** a, double* b)
//(float **a, int n, float p[], float b[], float x[])
{
    int i, k;
    double sum;
    double** L = cholesky_decomposition(a, n);
    double* x = (double*)malloc(sizeof(double) * n);
    for (i = 0;i < n;i++) {
        for (sum = b[i], k = i - 1;k >= 0;k--)
            sum -= L[i][k] * x[k];
        x[i] = sum / L[i][i];
    }
    for (i = n - 1;i >= 0;i--) {
        for (sum = x[i], k = i + 1;k < n;k++)
            sum -= L[k][i] * x[k];
        x[i] = sum / L[i][i];
    }
    for (i = 0;i < n;i++)
        free(L[i]);
    free(L);
    return x;
}

void make_the_matrix_positive(double** M, int N, double eps) {
    double yn = is_it_positive(M, N);
    while (yn > 0) {
        yn *= 2;
        printf("make_the_matrix_positive: covariance matrix not positive defined, adding cov[0][0]*%g*I \n", yn * eps);
        for (int i = 0;i < N;i++)
            M[i][i] += M[0][0] * eps * yn;
        for (int i = 0;i < N;i++) {
            for (int j = 0;j < N;j++)
                printf("%g\t", M[i][j]);
            printf("\n");
        }
        yn *= is_it_positive(M, N);
        if (yn < 0) {
            printf("make_the_matrix_positive:\n Cannot make the covariance matrix positive defined\n");
            exit(1);
        }
    }
    printf("now the matrix is positive defined.  %g\n", yn);
}


//return the vector x: the solution of Mx=b
//M is a matrix 
//it uses the LU decomposition see Numerical Receipes
double* LU_decomposition_solver(int N, double** M, double* b) {
    double** U, ** L, * y, * x;
    int i, j, k;

    x = (double*)malloc(sizeof(double) * N);
    y = (double*)malloc(sizeof(double) * N);
    L = (double**)malloc(sizeof(double*) * N);
    U = (double**)malloc(sizeof(double*) * N);
    for (i = 0;i < N;i++) {
        U[i] = (double*)calloc(N, sizeof(double));
        L[i] = (double*)calloc(N, sizeof(double));
    }

    for (i = 0;i < N;i++)
        L[i][i] = 1;

    for (j = 0;j < N;j++) {
        for (i = 0;i <= j;i++) {
            U[i][j] = M[i][j];
            for (k = 0;k < i;k++)
                U[i][j] -= L[i][k] * U[k][j];
        }
        for (i = j + 1;i < N;i++) {
            L[i][j] = M[i][j];
            for (k = 0;k < j;k++)
                L[i][j] -= L[i][k] * U[k][j];
            L[i][j] /= U[j][j];
        }
    }

    y[0] = b[0] / L[0][0];
    for (i = 0;i < N;i++) {
        y[i] = b[i];
        for (k = 0;k < i;k++)
            y[i] -= L[i][k] * y[k];
        y[i] /= L[i][i];
    }

    x[N - 1] = y[N - 1] / U[N - 1][N - 1];
    for (i = N - 2;i >= 0;i--) {
        x[i] = y[i];
        for (k = i + 1;k < N;k++)
            x[i] -= U[i][k] * x[k];
        x[i] /= U[i][i];
    }

    free(y);
    for (i = 0;i < N;i++) {
        free(L[i]);
        free(U[i]);
    }
    free(U);
    free(L);

    return x;

}

//return the inverse matrix of M
double** matrix_inverse(int N, double** M) {
    double* b, ** r, * a;
    int i, j;

    b = (double*)calloc(N, sizeof(double));
    r = (double**)malloc(sizeof(double*) * N);
    for (i = 0;i < N;i++) {
        r[i] = (double*)malloc(N * sizeof(double));
    }

    for (i = 0;i < N;i++) {
        b[i] = 1.;
        a = LU_decomposition_solver(N, M, b);
        for (j = 0;j < N;j++)
            r[j][i] = a[j];
        free(a);
        b[i] = 0;
    }

    free(b);

    return r;
}


double* cholesky_solver_if_possible(int n, double** a, double* b)
//(float **a, int n, float p[], float b[], float x[])
{
    int i, k, j;
    double sum;
    double* x = (double*)malloc(sizeof(double) * n);

    double** L = (double**)malloc(sizeof(double*) * n);
    for (i = 0;i < n;i++)
        L[i] = (double*)calloc(n, sizeof(double));

    for (i = 0;i < n;i++) {
        for (j = i;j < n;j++) {
            sum = a[i][j];

            for (k = i - 1;k >= 0;k--) sum -= L[i][k] * L[j][k];
            if (i == j) {
                if (sum <= 0.0) {
                    //      printf("cholesky_solver failed: the matrix, with rounding errors, is not positive definite,\n switching to LU_decomposition_solver\n");
                    for (i = 0;i < n;i++)
                        free(L[i]);
                    free(L); free(x);
                    return LU_decomposition_solver(n, a, b);

                }
                L[i][i] = sqrt(sum);
            }
            else L[j][i] = sum / L[i][i];
        }
    }
    for (i = 0;i < n;i++) {
        for (sum = b[i], k = i - 1;k >= 0;k--)
            sum -= L[i][k] * x[k];
        x[i] = sum / L[i][i];
    }
    for (i = n - 1;i >= 0;i--) {
        for (sum = x[i], k = i + 1;k < n;k++)
            sum -= L[k][i] * x[k];
        x[i] = sum / L[i][i];
    }
    for (i = 0;i < n;i++)
        free(L[i]);
    free(L);
    return x;
}

//return the inverse of a symmetric matrix 
double** symmetric_matrix_inverse(int N, double** M) {
    double* b, ** r, * a;
    int i, j;

    b = (double*)calloc(N, sizeof(double));
    r = (double**)malloc(sizeof(double*) * N);
    for (i = 0;i < N;i++) {
        r[i] = (double*)malloc(N * sizeof(double));
    }

    for (i = 0;i < N;i++) {
        b[i] = 1.;
        a = cholesky_solver_if_possible(N, M, b);
        for (j = 0;j < N;j++)
            r[j][i] = a[j];
        free(a);
        b[i] = 0;
    }

    free(b);

    return r;
}



double inter_spline(double at, int Npoints, double* x, double* y) {
    double** M = double_malloc_2(Npoints, Npoints);
    double xn;
    for (int i = 0; i < Npoints;i++) {
        xn = 1;
        for (int j = 0; j < Npoints;j++) {
            M[i][j] = xn;
            xn *= x[i];
        }
    }
    double* P = LU_decomposition_solver(Npoints, M, y);
    free_2(Npoints, M);
    double r = 0;
    xn = 1;
    for (int i = 0; i < Npoints;i++) {
        r += P[i] * xn;
        xn *= at;
    }
    free(P);
    return r;
}


//fit N data (x[N],y[N][value/error])  with the function sum_{i=0;i<M}  a_i f_i(x)
//*fit_functions is a function such that fit_function(int M,double x)[i]=f_i
//it returns a[i][value,error]
double** linear_fit(int N, double* x, double** y, int M, double* fit_function(int, double)) {
    double** alpha, * X, * beta, ** a, ** C, * sigma;
    int i, j, k;


    beta = (double*)calloc(M, sizeof(double));
    a = (double**)malloc(M * sizeof(double*));
    alpha = (double**)malloc(sizeof(double*) * M);
    for (j = 0;j < M;j++) {
        alpha[j] = (double*)calloc(M, sizeof(double));
        a[j] = (double*)calloc(2, sizeof(double));
    }
    for (i = 0;i < N;i++) {
        X = fit_function(M, x[i]);
        for (j = 0;j < M;j++) {
            beta[j] += y[i][0] * X[j] / (y[i][1] * y[i][1]);
            for (k = 0;k < M;k++) {
                alpha[j][k] += X[j] * X[k] / (y[i][1] * y[i][1]);
            }
        }
        free(X);
    }
    if (M == 1) {
        C = (double**)malloc(sizeof(double*) * 1);
        C[0] = (double*)malloc(sizeof(double) * 1);
        C[0][0] = 1. / alpha[0][0];
    }
    else
        C = matrix_inverse(M, alpha);


    for (j = 0;j < M;j++) {
        for (k = 0;k < M;k++) {
            a[j][0] += C[j][k] * beta[k];
        }
        a[j][1] = sqrt(C[j][j]);
    }

    for (j = 0;j < M;j++) {
        free(alpha[j]);
        free(C[j]);
    }
    free(C);free(alpha);free(beta);
    return a;
}



//fit N data (x[N][Nvariables],y[N][value/error])  with the function sum_{i=0;i<M}  a_i f_i(x)
//*fit_functions is a function such that fit_function(int M,double *x)[i]=f_i
//Nvariables is not required in this function but needs to be consistend in the declaration of fit_function and x
//it returns a[i][value,error]
double** linear_fit_Nvar(int N, double** x, double** y, int M, double* fit_function(int, double*)) {
    double** alpha, * X, * beta, ** a, ** C, * sigma;
    int i, j, k;

    beta = (double*)calloc(M, sizeof(double));
    a = (double**)malloc(M * sizeof(double*));
    alpha = (double**)malloc(sizeof(double*) * M);
    for (j = 0;j < M;j++) {
        alpha[j] = (double*)calloc(M, sizeof(double));
        a[j] = (double*)calloc(2, sizeof(double));
    }
    for (i = 0;i < N;i++) {
        X = fit_function(M, x[i]);
        for (j = 0;j < M;j++) {
            beta[j] += y[i][0] * X[j] / (y[i][1] * y[i][1]);
            for (k = 0;k < M;k++) {
                alpha[j][k] += X[j] * X[k] / (y[i][1] * y[i][1]);
            }
        }
        free(X);
    }
    if (M == 1) {
        C = (double**)malloc(sizeof(double*) * 1);
        C[0] = (double*)malloc(sizeof(double) * 1);
        C[0][0] = 1. / alpha[0][0];
    }
    else
        C = matrix_inverse(M, alpha);


    for (j = 0;j < M;j++) {
        for (k = 0;k < M;k++) {
            a[j][0] += C[j][k] * beta[k];
        }
        a[j][1] = sqrt(C[j][j]);
    }

    for (j = 0;j < M;j++) {
        free(alpha[j]);
        free(C[j]);
    }
    free(C);free(alpha);free(beta);
    return a;
}


//fit N data (x[N][Nvariables],y[N][value/error])  with the function sum_{i=0;i<M}  a_i f_i(x)
//*fit_functions is a function such that fit_function(int Nfunc, int Nvar,double *x,int Npar)[i]=f_i
//Nvariables is not required in this function but you need to be consistend in the declaration of fit_function and x
//it returns a[i][value,error]
double* linear_fit_Nf(int Nfunc, int* Npoints, double** x, double** y, int Nvar, int Npar, double* fit_function(int, int, double*, int)) {
    double** alpha, * X, * beta, ** a, ** C, * sigma;
    int i, j, k, e, count, n;
    double* r;

    beta = (double*)calloc(Npar, sizeof(double));
    r = (double*)malloc(Npar * sizeof(double));
    a = (double**)malloc(Npar * sizeof(double*));
    alpha = (double**)malloc(sizeof(double*) * Npar);
    for (j = 0;j < Npar;j++) {
        alpha[j] = (double*)calloc(Npar, sizeof(double));
        a[j] = (double*)calloc(2, sizeof(double));
    }
    count = 0;
    for (n = 0;n < Nfunc;n++) {
        for (e = 0;e < Npoints[n];e++) {
            i = e + count;
            X = fit_function(n, Nvar, x[i], Npar);
            for (j = 0;j < Npar;j++) {
                beta[j] += y[i][0] * X[j] / (y[i][1] * y[i][1]);
                for (k = 0;k < Npar;k++) {
                    alpha[j][k] += X[j] * X[k] / (y[i][1] * y[i][1]);
                }
            }
            free(X);
        }
        count += Npoints[n];
    }
    if (Npar == 1) {
        C = (double**)malloc(sizeof(double*) * 1);
        C[0] = (double*)malloc(sizeof(double) * 1);
        C[0][0] = 1. / alpha[0][0];
    }
    else
        C = matrix_inverse(Npar, alpha);

    /*  for (j=0;j<Npar;j++){
          for (k=0;k<Npar;k++){

          }
          printf("\n");
      }printf("\n");
  */
    for (j = 0;j < Npar;j++) {
        for (k = 0;k < Npar;k++) {
            a[j][0] += C[j][k] * beta[k];
        }
        a[j][1] = sqrt(C[j][j]);
    }
    for (j = 0;j < Npar;j++) {
        r[j] = a[j][0];
    }

    for (j = 0;j < Npar;j++) {
        free(alpha[j]);    free(a[j]);
        free(C[j]);
    }
    free(C);free(alpha);free(beta); free(a);



    return r;
}

double compute_chi_linear_Nf(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double* fun(int, int, double*, int)) {
    double chi2 = 0, * f, tmp;
    int e, j, n, count, k;

    count = 0;
    for (n = 0;n < N;n++) {
        for (e = 0;e < ensemble[n];e++) {
            f = fun(n, Nvar, x[e + count], Npar);
            tmp = 0;
            for (k = 0;k < Npar;k++) {
                tmp += (f[k] * P[k]);
            }
            chi2 += (y[e + count][0] - tmp) * (y[e + count][0] - tmp) / (y[e + count][1] * y[e + count][1]);

        }
        count += ensemble[n];
    }

    return chi2;
}




//it needs **a: the result of **linear_fit 
double  compute_chisqr(int N, double* x, double** y, int M, double** a, double* fit_function(int, double)) {
    double chi = 0, * X, tmp;
    int i, k;

    for (i = 0;i < N;i++) {
        X = fit_function(M, x[i]);
        tmp = 0;
        for (k = 0;k < M;k++) {
            tmp += (X[k] * a[k][0]);
        }
        chi += (y[i][0] - tmp) * (y[i][0] - tmp) / (y[i][1] * y[i][1]);
        free(X);
    }
    return chi;
}

double  compute_chisqr_Nvar(int N, double** x, double** y, int M, double** a, double* fit_function(int, double*)) {
    double chi = 0, * X, tmp;
    int i, k;

    for (i = 0;i < N;i++) {
        X = fit_function(M, x[i]);
        tmp = 0;
        for (k = 0;k < M;k++) {
            tmp += (X[k] * a[k][0]);
        }
        chi += (y[i][0] - tmp) * (y[i][0] - tmp) / (y[i][1] * y[i][1]);
        free(X);
    }
    return chi;
}


double* constant_fit_to_try(int M, double in) {
    double* r;

    r = (double*)malloc(sizeof(double) * M);
    r[0] = 1.;

    return r;
}

double* give_jack_linear_fit(int tmin, int tmax, int sep, double** corr_ave, double** corr_J, int Njack) {

    double*** y, * x, ** tmp, * fit;
    int i, j;

    y = (double***)malloc(sizeof(double**) * Njack);

    for (j = 0;j < Njack;j++) {
        y[j] = (double**)malloc(sizeof(double*) * (tmax - tmin + 1));
        for (i = tmin;i <= tmax;i++) {
            y[j][i - tmin] = (double*)malloc(sizeof(double) * 2);
        }
    }
    x = (double*)malloc(sizeof(double) * (tmax - tmin + 1));
    fit = (double*)malloc(sizeof(double) * Njack);

    for (i = tmin;i <= tmax;i += sep) {
        for (j = 0;j < Njack;j++) {
            y[j][(i - tmin) / sep][0] = corr_J[i][j];
            y[j][(i - tmin) / sep][1] = corr_ave[i][1];
        }
    }
    for (j = 0;j < Njack;j++) {
        tmp = linear_fit((tmax - tmin) / sep + 1, x, y[j], 1, constant_fit_to_try);
        fit[j] = tmp[0][0];
        free(tmp[0]);free(tmp);
    }

    for (j = 0;j < Njack;j++) {
        for (i = tmin;i <= tmax;i++) {
            free(y[j][i - tmin]);
        }
        free(y[j]);
    }
    free(x);
    free(y);
    return fit;
}

double* try_linear_fit(char** option, int tmin, int tmax, int sep, double** corr_ave, double** corr_J, int Njack, double** chi2) {

    double*** y, * x, * r, ** tmp, * fit;
    int i, j;
    double* chi2j;

    y = (double***)malloc(sizeof(double**) * Njack);
    chi2j = (double*)malloc(sizeof(double) * Njack);

    for (j = 0;j < Njack;j++) {
        y[j] = (double**)malloc(sizeof(double*) * (tmax - tmin + 1));
        for (i = tmin;i <= tmax;i++) {
            y[j][i - tmin] = (double*)malloc(sizeof(double) * 2);
        }
    }
    x = (double*)malloc(sizeof(double) * (tmax - tmin + 1));
    fit = (double*)malloc(sizeof(double) * Njack);

    for (i = tmin;i <= tmax;i += sep) {
        for (j = 0;j < Njack;j++) {
            y[j][(i - tmin) / sep][0] = corr_J[i][j];
            y[j][(i - tmin) / sep][1] = corr_ave[i][1];
        }
    }
    for (j = 0;j < Njack;j++) {
        tmp = linear_fit((tmax - tmin) / sep + 1, x, y[j], 1, constant_fit_to_try);
        chi2j[j] = compute_chisqr((tmax - tmin) / sep + 1, x, y[j], 1, tmp, constant_fit_to_try) / ((tmax - tmin) / sep + 1 - 1);
        fit[j] = tmp[0][0];
        free(tmp[0]);free(tmp);
    }
    if (strcmp(option[4], "jack") == 0) {
        r = mean_and_error_jack_biased(Njack, fit);
        (*chi2) = mean_and_error_jack_biased(Njack, chi2j);
        //r=mean_and_error_jack(Njack, fit);
        //(*chi2)=mean_and_error_jack(Njack, chi2j);
    }
    if (strcmp(option[4], "boot") == 0) {
        r = mean_and_error_boot(Njack, fit);
        (*chi2) = mean_and_error_boot(Njack, chi2j);;
    }

    free(fit);
    for (j = 0;j < Njack;j++) {
        for (i = tmin;i <= tmax;i++) {
            free(y[j][i - tmin]);
        }
        free(y[j]);
    }
    free(x);free(chi2j);
    free(y);
    return r;
}
double** global_linear_fit(int N, double** x, double** y, int M) {
    double** alpha, * X, * beta, ** a, ** C, * sigma;
    int i, j, k;


    beta = (double*)calloc(M, sizeof(double));
    a = (double**)malloc(M * sizeof(double*));
    alpha = (double**)malloc(sizeof(double*) * M);
    for (j = 0;j < M;j++) {
        alpha[j] = (double*)calloc(M, sizeof(double));
        a[j] = (double*)calloc(2, sizeof(double));
    }
    for (i = 0;i < N;i++) {
        X = x[i];
        for (j = 0;j < M;j++) {
            beta[j] += y[i][0] * X[j] / (y[i][1] * y[i][1]);
            for (k = 0;k < M;k++) {
                alpha[j][k] += X[j] * X[k] / (y[i][1] * y[i][1]);
            }
        }
        //  free(X);  
    }
    if (M == 1) {
        C = (double**)malloc(sizeof(double*) * 1);
        C[0] = (double*)malloc(sizeof(double) * 1);
        C[0][0] = 1. / alpha[0][0];
    }
    else
        C = matrix_inverse(M, alpha);


    for (j = 0;j < M;j++) {
        for (k = 0;k < M;k++) {
            a[j][0] += C[j][k] * beta[k];
        }
        a[j][1] = sqrt(C[j][j]);
    }

    for (j = 0;j < M;j++) {
        free(alpha[j]);
        free(C[j]);
    }
    free(C);free(alpha);free(beta);
    return a;
}

double  compute_chisqr_global_fit(int N, double** x, double** y, int M, double** a) {
    double chi = 0, * X, tmp;
    int i, k;

    for (i = 0;i < N;i++) {
        X = x[i];
        tmp = 0;
        for (k = 0;k < M;k++) {
            tmp += (X[k] * a[k][0]);
        }
        chi += (y[i][0] - tmp) * (y[i][0] - tmp) / (y[i][1] * y[i][1]);
        // free(X);
    }
    return chi;
}
