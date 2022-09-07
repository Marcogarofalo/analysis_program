#define non_linear_fit_C


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "mutils.hpp"

#define MAXIT 10000

//#include "jacknife.h"
//#include "bootstrap.h"
#include "linear_fit.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
#include "global.hpp"
#include <random>
#include "fit_all.hpp"
#include "sorting.hpp"
#include "resampling.hpp"


double* fit_type::linear_function(int n, int  Nvar, double* x, int Npar) {
    double* P = (double*)calloc(Npar, sizeof(double));
    double* r = (double*)malloc(Npar * sizeof(double));
    for (int i = 0;i < Npar;i++) {
        P[i] = 1.;
        r[i] = this->function(n, Nvar, x, Npar, P);
        P[i] = 0.;
    }
    free(P);
    return r;
}

void fit_type::compute_cov_fit(char** argv, data_all gjack, double lhs_fun(int, int, int, data_all, struct fit_type), struct fit_type fit_info) {
    ////// allocation
    int* en = (int*)malloc(sizeof(int) * N);// we need to init en and en_tot to allocate the other 
    for (int e = 0;e < N; e++) { en[e] = myen.size(); }
    int en_tot = 0;      for (int n = 0;n < N;n++) { en_tot += en[n]; }// total data to fit

    double*** y = double_malloc_3(Njack, en_tot, 2);// 2 is mean/error
    ////////////////////////////////////////// y
    int count = 0;
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < en[n];e++) {
            double* tmpj = (double*)malloc(sizeof(double) * Njack);
            for (int j = 0;j < Njack;j++) {
                tmpj[j] = lhs_fun(n, myen[e], j, gjack, fit_info);
            }
            double err = error_jackboot(argv[1], Njack, tmpj);
            for (int j = 0;j < Njack;j++) {
                y[j][e + count][0] = tmpj[j];
                y[j][e + count][1] = err;
            }
            free(tmpj);
        }
        count += en[n];
    }
    //////  init end
    double** yy = double_malloc_2(en_tot, Njack);
    for (int i = 0;i < en_tot;i++)
        for (int j = 0;j < Njack;j++)
            yy[i][j] = y[j][i][0];


    cov = covariance(argv[1], en_tot, Njack, yy);
    free_2(en_tot, yy);

    int yn = is_it_positive(cov, en_tot);
    while (yn == 1) {
        printf("covariance matrix not positive defined adding eps*cov[0][0]*I \n");
        for (int i = 0;i < en_tot;i++)
            cov[i][i] += cov[0][0] * 1e-12;
        yn = is_it_positive(cov, en_tot);
        printf("now the matrix is positive defined.  %d\n", yn);
    }
    cov1 = symmetric_matrix_inverse(en_tot, cov);
    free(en);
    free_3(Njack, en_tot, y);
    covariancey = true;
}

void fit_type::compute_cov1_fit() {
    int* en = (int*)malloc(sizeof(int) * N);// we need to init en and en_tot to allocate the other 
    for (int e = 0;e < N; e++) { en[e] = myen.size(); }
    int en_tot = 0;      for (int n = 0;n < N;n++) { en_tot += en[n]; }// total data to fit

    if (covariancey) {
        free_2(en_tot, cov1);
        int yn = is_it_positive(cov, en_tot);
        while (yn == 1) {
            printf("covariance matrix not positive defined adding eps*cov[0][0]*I \n");
            for (int i = 0;i < en_tot;i++)
                cov[i][i] += cov[0][0] * 1e-12;
            yn = is_it_positive(cov, en_tot);
            printf("now the matrix is positive defined.  %d\n", yn);
        }
        cov1 = symmetric_matrix_inverse(en_tot, cov);
    }
    free(en);
}

void fit_type::malloc_x() {
    error(Nvar <= 0, 1, "malloc_x for fit", "fit_info.Nvar must be greather than zero: current Nvar=%g\n", Nvar);
    error(myen.size() * N <= 0, 1, "malloc_x for fit",
        "fit_info.myen.size() and N must be greather than zero: current myen.size()=%d  N=%d\n", myen.size(), N);
    error(Njack <= 0, 1, "malloc_x for fit", "fit_info.Njack must be greather than zero: current Nvar=%g\n", Njack);
    x = double_malloc_3(Nvar, myen.size() * N, Njack);
    allocated_x = true;
}

void fit_type::restore_default() {
    custom = false; // 1 means default fit , 0 custom fit options
    lambda = 0.001;
    acc = 0.001;
    h = 1e-5;
    Prange = std::vector<double>();
    guess = std::vector<double>();
    corr_id = std::vector<int>();
    devorder = 4;
    repeat_start = 1;
    maxiter = 200;
    if (allocated_x) {
        for (int i = 0;i < Nvar;i++) {
            for (int j = 0;j < myen.size() * N;j++) {
                free(x[i][j]);
            }
            free(x[i]);
        }
        free(x);
    }
    allocated_x = false;
    // Nvar = 0;
    N = 0;
    myen = std::vector<int>();

    mean_only = false;
    NM = false;
    alpha = 1;
    gamma = 2;
    rho = 0.5;
    sigma = 0.5;

    unstable = false; // if true avoid thing that may return error
    noderiv = false;

    plateaux_scan = false;
    guess_per_jack = 0;
    chi2_gap_jackboot = 1;

    precision_sum = 0;
    verbosity = 0;
    mu = 0;

    t0_GEVP = 3;
    value_or_vector = 0;
    GEVP_tpt0 = false;
    GEVP_swap_t_t0 = false;
    GEVP_ignore_warning_after_t = 1000;

    HENKEL_size = 1;
    if (n_ext_P > 0) {
        for (int i = 0;i < n_ext_P;i++) {
            ext_P[i] = nullptr;
        }
        free(ext_P);
        n_ext_P = 0;
    }

    band_range = std::vector<double>();
    //      f_plateaux_scan=NULL;
    //      name_plateaux_scan="\0";
    if (covariancey) {
        int* en = (int*)malloc(sizeof(int) * N);// we need to init en and en_tot to allocate the other 
        for (int e = 0;e < N; e++) { en[e] = myen.size(); }
        int en_tot = 0;      for (int n = 0;n < N;n++) { en_tot += en[n]; }// total data to fit
        free_2(en_tot, cov);
        free_2(en_tot, cov1);
        free(en);
    }
    covariancey = false;
    second_deriv = false;

    linear_fit = false;
    codeplateaux = false;
    tmin = -1;
    tmax = -1;
}


struct fit_result malloc_fit(struct  fit_type  fit_info) {
    struct fit_result fit_out;
    fit_out.Njack = fit_info.Njack;
    fit_out.P = (double**)malloc(sizeof(double*) * fit_info.Npar);
    fit_out.chi2 = (double*)malloc(sizeof(double) * fit_info.Njack);
    fit_out.C = (double***)malloc(sizeof(double**) * fit_info.Npar);
    for (int i = 0;i < fit_info.Npar;i++) {
        fit_out.P[i] = (double*)malloc(sizeof(double) * fit_info.Njack);

        fit_out.C[i] = (double**)malloc(sizeof(double*) * fit_info.Npar);
        for (int n = 0;n < fit_info.Npar;n++) {
            fit_out.C[i][n] = (double*)malloc(sizeof(double) * fit_info.Njack);

        }
    }
    return fit_out;
}


void free_fit_result(struct  fit_type  fit_info, struct fit_result  out) {
    for (int i = 0; i < fit_info.Npar;i++) {
        free(out.P[i]);
        for (int n = 0;n < fit_info.Npar;n++)
            free(out.C[i][n]);
        free(out.C[i]);
    }
    free(out.P);free(out.C);free(out.chi2);

}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////


struct  fit_result   malloc_copy_fit_result(struct fit_result fit_out) {
    struct  fit_result fit_tmp;
    int i, j, k;
    fit_tmp.Njack = fit_out.Njack;
    fit_tmp.P = (double**)malloc(sizeof(double*) * fit_out.Npar);
    fit_tmp.C = (double***)malloc(sizeof(double**) * fit_out.Npar);
    fit_tmp.chi2 = (double*)malloc(sizeof(double) * fit_out.Njack);

    for (i = 0;i < fit_out.Npar;i++) {
        fit_tmp.P[i] = (double*)malloc(sizeof(double) * fit_out.Njack);
        fit_tmp.C[i] = (double**)malloc(sizeof(double*) * fit_out.Npar);
        for (j = 0;j < fit_out.Njack;j++) {
            fit_tmp.P[i][j] = fit_out.P[i][j];
        }
        for (k = 0;k < fit_out.Npar;k++) {
            fit_tmp.C[i][k] = (double*)malloc(sizeof(double) * fit_out.Njack);
            for (j = 0;j < fit_out.Njack;j++)
                fit_tmp.C[i][k][j] = fit_out.C[i][k][j];
        }

    }
    for (j = 0;j < fit_out.Njack;j++) {
        fit_tmp.chi2[j] = fit_out.chi2[j];
    }
    fit_tmp.Npar = fit_out.Npar;
    mysprintf(fit_tmp.name, NAMESIZE, "%s", fit_out.name);
    return fit_tmp;
}


struct  fit_result   copy_fit_result(struct fit_result fit_out, struct fit_type fit_info) {
    struct  fit_result fit_tmp;
    int i, j, k;
    fit_tmp.Njack = fit_out.Njack;
    fit_tmp.P = (double**)malloc(sizeof(double*) * fit_info.Npar);
    fit_tmp.C = (double***)malloc(sizeof(double**) * fit_info.Npar);
    fit_tmp.chi2 = (double*)malloc(sizeof(double) * fit_out.Njack);

    for (i = 0;i < fit_info.Npar;i++) {
        fit_tmp.P[i] = (double*)malloc(sizeof(double) * fit_out.Njack);
        fit_tmp.C[i] = (double**)malloc(sizeof(double*) * fit_info.Npar);
        for (j = 0;j < fit_out.Njack;j++) {
            fit_tmp.P[i][j] = fit_out.P[i][j];
        }
        for (k = 0;k < fit_info.Npar;k++) {
            fit_tmp.C[i][k] = (double*)malloc(sizeof(double) * fit_out.Njack);
            for (j = 0;j < fit_out.Njack;j++)
                fit_tmp.C[i][k][j] = fit_out.C[i][k][j];
        }

    }
    for (j = 0;j < fit_out.Njack;j++) {
        fit_tmp.chi2[j] = fit_out.chi2[j];
    }
    return fit_tmp;
}

struct  fit_type   copy_fit_type(struct fit_result fit_out, struct fit_type fit_info) {
    struct  fit_type fit_tmp;
    fit_tmp.N = fit_info.N;
    fit_tmp.Npar = fit_info.Npar;
    fit_tmp.Nvar = fit_info.Nvar;
    fit_tmp.function = fit_info.function;

    return fit_tmp;
}

void  copy_fit_type_into(struct  fit_type* fit_tmp, struct fit_type fit_info) {

    fit_tmp->N = fit_info.N;
    fit_tmp->Npar = fit_info.Npar;
    fit_tmp->Nvar = fit_info.Nvar;
    fit_tmp->function = fit_info.function;

}

struct fit_all   save_fit(struct fit_all fit_chi2_good, struct fit_type fit_info, struct fit_result fit_out) {

    struct fit_all fit_tmp;
    int i, j, k;
    int N = fit_chi2_good.Nfits + 1;
    fit_tmp.Nfits = N;
    fit_tmp.info = (struct fit_type*)malloc(sizeof(struct fit_type) * N);
    fit_tmp.out = (struct fit_result*)malloc(sizeof(struct fit_result) * N);
    for (i = 0;i < N - 1;i++) {
        fit_tmp.out[i] = copy_fit_result(fit_chi2_good.out[i], fit_chi2_good.info[i]);
        fit_tmp.info[i] = copy_fit_type(fit_chi2_good.out[i], fit_chi2_good.info[i]);
        for (j = 0;j < fit_chi2_good.info[i].Npar;j++) {
            free(fit_chi2_good.out[i].P[j]);
            for (k = 0;k < fit_chi2_good.info[i].Npar;k++) {
                free(fit_chi2_good.out[i].C[j][k]);
            }
            free(fit_chi2_good.out[i].C[j]);
        }
        free(fit_chi2_good.out[i].chi2);
        free(fit_chi2_good.out[i].P);
        free(fit_chi2_good.out[i].C);

    }
    fit_tmp.out[N - 1] = copy_fit_result(fit_out, fit_info);
    fit_tmp.info[N - 1] = copy_fit_type(fit_out, fit_info);

    return fit_tmp;

}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference 
//derivative 1 accuracy 4
double* der_fun_h(int Nvar, double* x, int Npar, double* P, double fun(int, double*, int, double*), double h) {

    double* Ph = (double*)malloc(sizeof(double) * Npar);
    double* df = (double*)calloc(Npar, sizeof(double));
    int i;

    for (i = 0;i < Npar;i++)
        Ph[i] = P[i];

    for (i = 0;i < Npar;i++) {
        Ph[i] = P[i] - 2. * h;
        df[i] = fun(Nvar, x, Npar, Ph);

        Ph[i] = P[i] - h;
        df[i] -= 8 * fun(Nvar, x, Npar, Ph);

        Ph[i] = P[i] + h;
        df[i] += 8 * fun(Nvar, x, Npar, Ph);

        Ph[i] = P[i] + 2. * h;
        df[i] -= fun(Nvar, x, Npar, Ph);

        Ph[i] = P[i];//you need to leave the parameter as it was before you move to the next parameter
        df[i] /= (12. * h);
    }


    free(Ph);

    return df;
}

/////////////////////////////////////////////// /////////////////////////////////////////////// 
///////////////////////////////////////////////   chi2      /////////////////////////////////////////////// 
/////////////////////////////////////////////// /////////////////////////////////////////////// 

double compute_chi_non_linear_Nf_cov1_double(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), fit_type fit_info) {
    double chi2 = 0, f, f1;
    int e, n, count, e1, n1, count1;

    int en_tot = 0;
    for (n = 0;n < N;n++)
        for (e = 0;e < ensemble[n];e++)
            en_tot++;

    double* tmp = (double*)malloc(sizeof(double) * en_tot);
    count = 0;
    for (n = 0;n < N;n++) {
        for (e = 0;e < ensemble[n];e++) {
            tmp[count] = fun(n, Nvar, x[count], Npar, P) - y[count][0];// f1(n,e,N,N1,Nvar,x1[count],Npar,Npar1,P,fun)-y1[count][0];
            // printf("predicted=%g    latt=%g   n=%d  e=%d\n",fun(n, Nvar, x[count], Npar, P), y[count][0], n, e);
            count++;
        }
    }

    for (int i = 0;i < en_tot;i++)
        chi2 += tmp[i] * fit_info.cov1[i][i] * tmp[i];

    for (int i = 0;i < en_tot;i++)
        for (int j = i + 1;j < en_tot;j++)
            chi2 += 2. * tmp[i] * fit_info.cov1[i][j] * tmp[j];

    free(tmp);
    return chi2;
}


double compute_chi_non_linear_Nf_cov1_long_double(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), fit_type fit_info) {
    double chi2 = 0, f, f1;
    int e, n, count, e1, n1, count1;

    int en_tot = 0;
    for (n = 0;n < N;n++)
        for (e = 0;e < ensemble[n];e++)
            en_tot++;

    double* tmp = (double*)malloc(sizeof(double) * en_tot);
    count = 0;
    for (n = 0;n < N;n++) {
        for (e = 0;e < ensemble[n];e++) {
            tmp[count] = fun(n, Nvar, x[count], Npar, P) - y[count][0];// f1(n,e,N,N1,Nvar,x1[count],Npar,Npar1,P,fun)-y1[count][0];
            // printf("predicted=%g    latt=%g   n=%d  e=%d\n",fun(n, Nvar, x[count], Npar, P), y[count][0], n, e);
            count++;
        }
    }
    long double ld = 0;
    for (int i = 0;i < en_tot;i++)
        ld += tmp[i] * fit_info.cov1[i][i] * tmp[i];

    for (int i = 0;i < en_tot;i++)
        for (int j = i + 1;j < en_tot;j++)
            ld += 2. * tmp[i] * fit_info.cov1[i][j] * tmp[j];

    chi2 = (double)ld;
    free(tmp);
    return chi2;
}

double compute_chi_non_linear(int ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, double*, int, double*)) {
    double chi2 = 0, f;
    int e, j;

    for (e = 0;e < ensemble;e++) {
        f = fun(Nvar, x[e], Npar, P);
        chi2 += (y[e][0] - f) * (y[e][0] - f) / (y[e][1] * y[e][1]);

    }

    return chi2;

}

// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(Nvariables,variables[], Nparameters,parameters[])
//*funk(Nvariables,variables[], Nparameters,parameters[]) must return a vector[Nparameters] that contain the derivatives of fun respect to the parameters
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
//double  *non_linear_fit(int ensemble,double **x, double **y ,int Nvar, int Npar, double fun(int,double*,int,double*) , double *funk(int,double*,int,double*) ){
double* non_linear_fit(int ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, double*, int, double*)) {

    double** alpha, * X, * beta, ** a, ** C, * sigma;
    int i, j, k, e;
    double f, * fk;
    double chi2, chi2_tmp;
    double* P, * P_tmp, lambda, res;
    double h = 0.00001;
    lambda = 0.001;
    res = 1;

    //   double *fkk;
    P = (double*)malloc(Npar * sizeof(double));
    P_tmp = (double*)malloc(Npar * sizeof(double));


    for (j = 0;j < Npar;j++) {
        //  P[j]=0.1;
        P_tmp[j] = P[j];
    }

    beta = (double*)calloc(Npar, sizeof(double));
    alpha = (double**)malloc(sizeof(double*) * Npar);
    for (j = 0;j < Npar;j++) {
        alpha[j] = (double*)calloc(Npar, sizeof(double));
    }
    chi2 = compute_chi_non_linear(ensemble, x, y, P, Nvar, Npar, fun);


    // printf("chi2=%f   res=%.10f Bw=%f   fw=%f\n",chi2,res,P[0],P[1]);
    while (res > 0.001) {
        chi2_tmp = chi2 + 1;

        while (chi2_tmp > chi2) {

            for (e = 0;e < ensemble;e++) {
                f = fun(Nvar, x[e], Npar, P);
                fk = der_fun_h(Nvar, x[e], Npar, P, fun, h);
                //   fk=funk(Nvar,x[e],Npar,P);
                for (j = 0;j < Npar;j++) {
                    beta[j] += (y[e][0] - f) * fk[j] / (y[e][1] * y[e][1]);
                    //      printf("|analitic-numeric|=  |%g -%g|   = %g\n",fk[j],fkk[j],fabs(fk[j]-fkk[j]));
                    for (k = 0;k < Npar;k++) {
                        alpha[j][k] += fk[j] * fk[k] / (y[e][1] * y[e][1]);
                    }
                }
                free(fk);
            }

            for (j = 0;j < Npar;j++) {
                alpha[j][j] *= (lambda + 1.);
            }

            if (Npar == 1) {
                C = (double**)malloc(sizeof(double*) * 1);
                C[0] = (double*)malloc(sizeof(double) * 1);
                C[0][0] = 1. / alpha[0][0];
            }
            else
                C = matrix_inverse(Npar, alpha);

            for (j = 0;j < Npar;j++)
                P_tmp[j] = P[j];
            for (j = 0;j < Npar;j++) {
                for (k = 0;k < Npar;k++) {
                    P_tmp[j] += C[j][k] * beta[k];
                }
            }

            //printf("lambda=%f\n",lambda);
            chi2_tmp = compute_chi_non_linear(ensemble, x, y, P_tmp, Nvar, Npar, fun);
            if (chi2_tmp != chi2_tmp) chi2_tmp = chi2 + 1;
            //       printf("chi2=%f chi2_tmp=%f  res=%f P0=%f P0_tmp=%f  P0=%f P0_tmp=%f\n",chi2,chi2_tmp,res,P[0],P_tmp[0],P[1],P_tmp[1]);


            if (chi2_tmp > chi2)
                lambda *= 10;
            for (j = 0;j < Npar;j++) {
                free(C[j]);
                beta[j] = 0;
                for (k = 0;k < Npar;k++) {
                    alpha[j][k] = 0;
                }
            }
            free(C);


        }
        res = chi2 - chi2_tmp;
        //error(res<0,2,"non_linear_fit","The Levenberg-Marquardt accepted a configuration when the chi2 increased");
        chi2 = chi2_tmp;
        lambda /= 10;
        for (j = 0;j < Npar;j++) {
            P[j] = P_tmp[j];
        }
        // printf("chi2=%f   res=%.10f Bw=%f   fw=%f  P1=%f  P2=%f\n",chi2,res,P[0],P[1],P[2],P[3]);
        //printf("check residue=%f\n",res);
    }


    for (j = 0;j < Npar;j++) {
        free(alpha[j]);
    }
    free(P_tmp);
    free(alpha);free(beta);
    return P;
}

//https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference 
//derivative 1 accuracy 4
double* der_O4_fun_Nf_h(int n, int Nvar, double* x, int Npar, double* P, double fun(int, int, double*, int, double*), std::vector< double > h) {

    double* Ph = (double*)malloc(sizeof(double) * Npar);
    double* df = (double*)calloc(Npar, sizeof(double));
    int i;

    for (i = 0;i < Npar;i++)
        Ph[i] = P[i];

    for (i = 0;i < Npar;i++) {
        Ph[i] = P[i] - 2. * h[i];
        df[i] = fun(n, Nvar, x, Npar, Ph);

        Ph[i] = P[i] - h[i];
        df[i] -= 8 * fun(n, Nvar, x, Npar, Ph);

        Ph[i] = P[i] + h[i];
        df[i] += 8 * fun(n, Nvar, x, Npar, Ph);

        Ph[i] = P[i] + 2. * h[i];
        df[i] -= fun(n, Nvar, x, Npar, Ph);

        Ph[i] = P[i];//you need to leave the parameter as it was before you move to the next parameter
        df[i] /= (12. * h[i]);
    }


    free(Ph);

    return df;
}



//https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference 
//derivative 1 accuracy 2
double* der_O2_fun_Nf_h(int n, int Nvar, double* x, int Npar, double* P, double fun(int, int, double*, int, double*), std::vector< double > h) {

    double* Ph = (double*)malloc(sizeof(double) * Npar);
    double* df = (double*)calloc(Npar, sizeof(double));
    int i;

    for (i = 0;i < Npar;i++)
        Ph[i] = P[i];

    for (i = 0;i < Npar;i++) {

        Ph[i] = P[i] - h[i];
        df[i] = -fun(n, Nvar, x, Npar, Ph);

        Ph[i] = P[i] + h[i];
        df[i] += fun(n, Nvar, x, Npar, Ph);


        Ph[i] = P[i];//you need to leave the parameter as it was before you move to the next parameter
        df[i] /= (2. * h[i]);
    }


    free(Ph);

    return df;
}


double der1_O2_h(int n, int Nvar, double* x, int Npar, double* P, double fun(int, int, double*, int, double*), int dir, std::vector< double > h) {

    double* Ph = (double*)malloc(sizeof(double) * Npar);
    double df;
    int i;

    for (i = 0;i < Npar;i++)
        Ph[i] = P[i];

    i = dir;
    Ph[i] = P[i] - h[i];
    df = -fun(n, Nvar, x, Npar, Ph);

    Ph[i] = P[i] + h[i];
    df += fun(n, Nvar, x, Npar, Ph);

    Ph[i] = P[i];//you need to leave the parameter as it was before you move to the next parameter
    df /= (2. * h[i]);

    free(Ph);
    return df;
}

//https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference 
//derivative 1 accuracy 2
double** der2_O2_fun_Nf_h(int n, int Nvar, double* x, int Npar, double* P, double fun(int, int, double*, int, double*), std::vector< double > h) {

    double* Ph = (double*)malloc(sizeof(double) * Npar);
    double** df = double_calloc_2(Npar, Npar);
    int i, j;

    for (i = 0;i < Npar;i++)
        Ph[i] = P[i];

    for (i = 0;i < Npar;i++) {
        for (j = 0;j < Npar;j++) {
            Ph[i] = P[i] - h[i];
            df[i][j] = -der1_O2_h(n, Nvar, x, Npar, Ph, fun, j, h);

            Ph[i] = P[i] + h[i];
            df[i][j] += der1_O2_h(n, Nvar, x, Npar, Ph, fun, j, h);


            Ph[i] = P[i];//you need to leave the parameter as it was before you move to the next parameter
            df[i][j] /= (2. * h[i]);
        }
    }


    free(Ph);

    return df;
}




double* der_O2_fun_Nf_h_adaptive(int n, int Nvar, double* x, int Npar, double* P, double fun(int, int, double*, int, double*), std::vector< double > h) {

    double* Ph = (double*)malloc(sizeof(double) * Npar);
    double* df = (double*)calloc(Npar, sizeof(double));
    int i;

    for (i = 0;i < Npar;i++)
        Ph[i] = P[i];


    for (i = 0;i < Npar;i++) {
        double hp = h[i] * fabs(P[i]);
        Ph[i] = P[i] - hp;
        df[i] = -fun(n, Nvar, x, Npar, Ph);

        Ph[i] = P[i] + hp;
        df[i] += fun(n, Nvar, x, Npar, Ph);


        Ph[i] = P[i];//you need to leave the parameter as it was before you move to the next parameter
        df[i] /= (2. * hp);
    }


    free(Ph);

    return df;
}


//https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference 
//derivative 1 accuracy 4
// double *der_fun_Nf_h(int n, int Nvar, double *x,int Npar,double  *P, double fun(int,int,double*,int,double*), double h){
//     
//     double *Ph=(double*) malloc(sizeof(double)*Npar);
//     double *df=(double*) calloc(Npar,sizeof(double));
//     int i;
//     
//     for (i=0;i<Npar;i++)
//         Ph[i]=P[i];
//     
//     for (i=0;i<Npar;i++){
//         Ph[i]=P[i]-2.*h;
//         df[i]=fun(n,Nvar,x,Npar,Ph);
//     
//         Ph[i]=P[i]-h;
//         df[i]-=8*fun(n,Nvar,x,Npar,Ph);
//         
//         Ph[i]=P[i]+h;
//         df[i]+=8*fun(n,Nvar,x,Npar,Ph);
//     
//         Ph[i]=P[i]+2.*h;
//         df[i]-=fun(n,Nvar,x,Npar,Ph);
//     
//         Ph[i]=P[i];//you need to leave the parameter as it was before you move to the next parameter
//         df[i]/=(12.*h);
//     }
//     
//     
//     free(Ph);
//     
//     return df;
// }
//https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference 
//derivative 1 accuracy 4
double* derN_fun_Nf_var_h(int n, int Nvar, double* x, int Npar, double* P, double fun(int, int, double*, int, double*), double h, int N) {

    double* xh = (double*)malloc(sizeof(double) * Nvar);
    double* df = (double*)calloc(Nvar, sizeof(double));
    double* tmp;
    int i;

    for (i = 0;i < Nvar;i++)
        xh[i] = x[i];

    if (N == 0) {
        for (i = 0;i < Nvar;i++)
            df[i] = fun(n, Nvar, xh, Npar, P);
    }
    else if (N == 1) {
        for (i = 0;i < Nvar;i++) {
            xh[i] = x[i] - 2. * h;
            df[i] = fun(n, Nvar, xh, Npar, P);

            xh[i] = x[i] - h;
            df[i] -= 8 * fun(n, Nvar, xh, Npar, P);

            xh[i] = x[i] + h;
            df[i] += 8 * fun(n, Nvar, xh, Npar, P);

            xh[i] = x[i] + 2. * h;
            df[i] -= fun(n, Nvar, xh, Npar, P);

            xh[i] = x[i];//you need to leave the parameter as it was before you move to the next parameter
            df[i] /= (12. * h);
        }
    }
    else {
        for (i = 0;i < Nvar;i++) {
            xh[i] = x[i] - 2. * h;
            tmp = derN_fun_Nf_var_h(n, Nvar, xh, Npar, P, fun, h, N - 1);
            df[i] = tmp[i]; free(tmp);

            xh[i] = x[i] - h;
            tmp = derN_fun_Nf_var_h(n, Nvar, xh, Npar, P, fun, h, N - 1);
            df[i] -= 8 * tmp[i]; free(tmp);

            xh[i] = x[i] + h;
            tmp = derN_fun_Nf_var_h(n, Nvar, xh, Npar, P, fun, h, N - 1);
            df[i] += 8 * tmp[i]; free(tmp);

            xh[i] = x[i] + 2. * h;
            tmp = derN_fun_Nf_var_h(n, Nvar, xh, Npar, P, fun, h, N - 1);
            df[i] -= tmp[i]; free(tmp);

            xh[i] = x[i];//you need to leave the parameter as it was before you move to the next parameter
            df[i] /= (12. * h);
        }
    }

    for (i = 1;i <= N;i++)
        df[i] /= ((double)i);
    free(xh);

    return df;
}



double compute_chi_non_linear_Nf(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), fit_type fit_info) {
    double chi2 = 0, f;
    int e, n, count;

    count = 0;
    for (n = 0;n < N;n++) {
        for (e = 0;e < ensemble[n];e++) {
            f = fun(n, Nvar, x[count], Npar, P) - y[count][0];
            // printf("e=%d  n=%d   f=%g  y=%g err=%g    dchi= %g  \n ",e,n,f+y[count][0],y[count][0],y[count][1],f*f/(y[count][1]* y[count][1]));
            //  for(int pp=0;pp<Npar; pp++) printf("P: %g  \t",P[pp]]); printf("\n");
            f /= y[count][1];
            chi2 += f * f;
            count++;
        }
    }
    return chi2;
}


double compute_chi_non_linear_Nf_long(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), fit_type fit_info) {
    long double chi2 = 0, f;
    int e, n, count;

    count = 0;
    for (n = 0;n < N;n++) {
        for (e = 0;e < ensemble[n];e++) {
            f = fun(n, Nvar, x[count], Npar, P) - y[count][0];
            // printf("predicted=%g    latt=%g   n=%d  e=%d\n",fun(n, Nvar, x[count], Npar, P), y[count][0], n, e);
            f /= y[count][1];
            chi2 += f * f;
            count++;
        }
    }
    double chi2_double = chi2;
    return chi2_double;
}


double compute_chi_non_linear_Nf_kahan(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), fit_type fit_info) {
    double chi2 = 0, f;
    int e, n, count;
    double c = 0.0;

    count = 0;
    for (n = 0;n < N;n++) {
        for (e = 0;e < ensemble[n];e++) {
            f = fun(n, Nvar, x[count], Npar, P) - y[count][0];
            //             printf("e=%d  n=%d   f=%g  y=%g err=%g    dchi= %g  \n ",e,n,f+y[count][0],y[count][0],y[count][1],f*f/(y[count][1]* y[count][1]));
            f /= y[count][1];
            double in = f * f - c;
            double t = chi2 + in;
            c = t - chi2 - in;
            chi2 = t;

            count++;
        }
    }
    return chi2;
}



bool compute_alpha(double* beta, double** alpha, int N, int* ensemble,
    double** x, double** y, int Nvar, int Npar, double* P_tmp,
    double fun(int, int, double*, int, double*),
    double* der_fun_Nf_h(int, int, double*, int, double*, double(int, int, double*, int, double*), std::vector< double >), std::vector< double > h,
    fit_type fit_info
) {
    for (int j = 0;j < Npar;j++) {
        beta[j] = 0;
        for (int k = 0;k < Npar;k++) {
            alpha[j][k] = 0;
        }
    }
    int count = 0;
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < ensemble[n];e++) {//printf("e=%d   n=%d   en[%d]=%d\n",e,n,n,ensemble[n]);
            double f = fun(n, Nvar, x[e + count], Npar, P_tmp);
            double* fk = der_fun_Nf_h(n, Nvar, x[e + count], Npar, P_tmp, fun, h);
            if (fit_info.verbosity > 3) {
                printf("n=%d     e=%d   fun=%g  derivative=\t", n, e, f);
                for (int j = 0;j < Npar;j++) { printf("%g \t", fk[j]); }printf("\n");
            }
            for (int j = 0;j < Npar;j++) {
                beta[j] += (y[e + count][0] - f) * fk[j] / (y[e + count][1] * y[e + count][1]);
                for (int k = j;k < Npar;k++) {
                    alpha[j][k] += fk[j] * fk[k] / (y[e + count][1] * y[e + count][1]);
                }
            }
            free(fk);
        }
        count += ensemble[n];
    }
    if (fit_info.second_deriv == true) {
        count = 0;
        for (int n = 0;n < N;n++) {
            for (int e = 0;e < ensemble[n];e++) {
                double f = fun(n, Nvar, x[e + count], Npar, P_tmp);
                double** fkk = der2_O2_fun_Nf_h(n, Nvar, x[e + count], Npar, P_tmp, fun, h);
                for (int j = 0;j < Npar;j++) {
                    for (int k = j;k < Npar;k++) {
                        alpha[j][k] -= (y[e + count][0] - f) * fkk[j][k] / (y[e + count][1] * y[e + count][1]);
                    }
                }
                free_2(Npar, fkk);

            }
            count += ensemble[n];
        }
    }
    return  true;
}

bool compute_alpha_cov1(double* beta, double** alpha, int N, int* ensemble,
    double** x, double** y, int Nvar, int Npar, double* P_tmp,
    double fun(int, int, double*, int, double*),
    double* der_fun_Nf_h(int, int, double*, int, double*, double(int, int, double*, int, double*), std::vector< double >), std::vector< double > h,
    fit_type fit_info
) {
    for (int j = 0;j < Npar;j++) {
        beta[j] = 0;
        for (int k = 0;k < Npar;k++) {
            alpha[j][k] = 0;
        }
    }
    int en_tot = 0;
    for (int n = 0;n < N;n++)
        for (int e = 0;e < ensemble[n];e++)
            en_tot += 1;
    int count = 0;
    double* f_value = (double*)malloc(sizeof(double) * en_tot);
    double** df_value = (double**)malloc(sizeof(double*) * en_tot);
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < ensemble[n];e++) {
            f_value[count] = fun(n, Nvar, x[e], Npar, P_tmp);
            df_value[count] = der_fun_Nf_h(n, Nvar, x[e], Npar, P_tmp, fun, h);
            count++;
        }
    }
    for (int j = 0;j < Npar;j++) {
        for (int i = 0;i < en_tot;i++)
            for (int ii = 0;ii < en_tot;ii++)
                beta[j] += (y[i][0] - f_value[i]) * fit_info.cov1[i][ii] * df_value[ii][j];


        for (int k = j;k < Npar;k++) {
            for (int i = 0;i < en_tot;i++)
                for (int ii = 0;ii < en_tot;ii++)
                    alpha[j][k] += df_value[i][j] * fit_info.cov1[i][ii] * df_value[ii][k];
        }
    }
    if (fit_info.second_deriv == true) {
        count = 0;
        for (int n = 0;n < N;n++) {
            for (int e = 0;e < ensemble[n];e++) {
                printf("second derivative %i  %i\n", n ,e);
                double** fkk = der2_O2_fun_Nf_h(n, Nvar, x[e + count], Npar, P_tmp, fun, h);
                int count1 = 0;
                for (int n1 = 0;n1 < N;n1++) {
                    for (int e1 = 0;e1 < ensemble[n1];e1++) {
                        for (int j = 0;j < Npar;j++) {
                            for (int k = j;k < Npar;k++) {
                                alpha[j][k] -= (y[e1 + count1][0] - f_value[e1 + count1]) * fit_info.cov1[e1 + count1][e + count] * fkk[j][k];
                            }
                        }
                    }
                    count1 += ensemble[n1];
                }
                free_2(Npar, fkk);
            }
            count += ensemble[n];
        }
    }

    for (int i = 0;i < en_tot;i++)
        free(df_value[i]);
    free(df_value);
    free(f_value);
    return  true;
}

// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(index_function,Nvariables,variables[], Nparameters,parameters[])
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
double** covariance_non_linear_fit_Nf(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), fit_type fit_info) {

    double** alpha, ** C;
    int i, j, k, e;
    double f, * fk;
    int n, count;
    std::vector< double > h;
    if (fit_info.h.h.size() == 1) {
        h = std::vector< double >(Npar);
        for (int i = 0; i < Npar;i++)
            h[i] = fit_info.h.h[0];
    }
    else if (fit_info.h.h.size() == Npar) { h = fit_info.h; }
    else { printf("non_linear_fit_Nf: h must me a std::vector<double> of 1 o Npar lenght"); exit(1); }


    int en_tot = 0;
    for (int n = 0;n < N;n++)
        for (int e = 0;e < ensemble[n];e++)
            en_tot += 1;

    alpha = (double**)malloc(sizeof(double*) * Npar);
    for (j = 0;j < Npar;j++) {
        alpha[j] = (double*)calloc(Npar, sizeof(double));
    }
    double* beta = (double*)calloc(Npar, sizeof(double));

    double* (*der_fun_Nf_h)(int, int, double*, int, double*, double(int, int, double*, int, double*), std::vector< double >);
    bool (*alpha_fun)(double*, double**, int, int*, double**, double**, int, int, double*, double(int, int, double*, int, double*),
        double* (int, int, double*, int, double*, double(int, int, double*, int, double*), std::vector< double >), std::vector< double >,
        fit_type);

    if (fit_info.covariancey)        alpha_fun = compute_alpha_cov1;
    else                             alpha_fun = compute_alpha;



    if (fit_info.devorder == -2)     der_fun_Nf_h = der_O2_fun_Nf_h_adaptive;
    else if (fit_info.devorder == 4) der_fun_Nf_h = der_O4_fun_Nf_h;
    else if (fit_info.devorder == 2) der_fun_Nf_h = der_O2_fun_Nf_h;
    else error(true, 1, "non_linear_fit_Nf", "order of the derivative must be 4 (default) , 2,  -2 to set a step different for each parameter h[i]=P[i]*h    ");


    bool computed_alpha = alpha_fun(beta, alpha, N, ensemble,
        x, y, Nvar, Npar, P, fun, der_fun_Nf_h, h, fit_info);

    for (j = 0;j < Npar;j++) {
        for (k = 0;k < j;k++)
            alpha[j][k] = alpha[k][j];
    }

    if (Npar == 1) {
        C = (double**)malloc(sizeof(double*) * 1);
        C[0] = (double*)malloc(sizeof(double) * 1);
        C[0][0] = 1. / alpha[0][0];
    }
    else {
        make_the_matrix_positive(alpha, Npar);
        C = symmetric_matrix_inverse(Npar, alpha);
    }


    for (j = 0;j < Npar;j++) {
        free(alpha[j]);
    }

    free(alpha);
    free(beta);
    return C;
}



void sort_NM(double* chi2s1, double** points1, double* chi2s, double** points, int Npar, int* order) {
    for (int i = 0;i < Npar + 1;i++)
        order[i] = i;
    quickSort(order, chi2s, 0, Npar);
    for (int i = 0;i < Npar + 1;i++) {
        // printf("%d order=%d   chi2s=%g\n", i, order[i], chi2s[i]);
        chi2s1[i] = chi2s[order[i]];
        for (int j = 0;j < Npar;j++)
            points1[i][j] = points[order[i]][j];
    }
}
void compute_centroid_NM(double* centroid, double** points, int Npar) {
    for (int j = 0;j < Npar;j++) {
        centroid[j] = 0;
        for (int i = 0;i < Npar;i++) {// apart the the last point
            centroid[j] += points[i][j];
        }
        centroid[j] /= (double)(Npar);
    }
}

void compute_reflected_NM(double* reflected, double* centroid, double** points, int Npar, double alpha) {
    for (int j = 0;j < Npar;j++) {
        reflected[j] = centroid[j] + alpha * (centroid[j] - points[Npar][j]);
    }
}

void compute_expanded_NM(double* expanded, double* reflected, double* centroid, int Npar, double gamma) {
    for (int j = 0;j < Npar;j++) {
        expanded[j] = centroid[j] + gamma * (reflected[j] - centroid[j]);
    }
}

void compute_contracted_out_NM(double* co, double* reflected, double* centroid, int Npar, double rho) {
    for (int j = 0;j < Npar;j++) {
        co[j] = centroid[j] + rho * (reflected[j] - centroid[j]);
    }
}

void compute_contracted_in_NM(double* ci, double* last, double* centroid, int Npar, double rho) {
    for (int j = 0;j < Npar;j++) {
        ci[j] = centroid[j] + rho * (last[j] - centroid[j]);
    }
}

void shrink_NM(double** points, int Npar, double sigma) {
    for (int i = 1; i < Npar + 1;i++) {
        for (int j = 0;j < Npar;j++) {
            points[i][j] = points[0][j] + sigma * (points[i][j] - points[0][j]);
        }
    }

}

double compute_sd_NM(double* centroid, double** points, int Npar) {
    double sd = 0;
    compute_centroid_NM(centroid, points, Npar);

    for (int j = 0;j < Npar;j++)
        centroid[j] = 0;
    for (int j = 0;j < Npar;j++) {
        for (int i = 0;i < Npar + 1;i++)
            centroid[j] += points[i][j];
        centroid[j] /= ((double)(Npar + 1));
    }
    for (int j = 0;j < Npar;j++) {
        for (int i = 0;i < Npar + 1;i++)
            sd += (centroid[j] - points[i][j]) * (centroid[j] - points[i][j]);
    }
    sd = sqrt(sd / ((double)Npar));
    return sd;
}



//fit N data (x[N][Nvariables],y[N][value/error])  with the function sum_{i=0;i<M}  a_i f_i(x)
//*fit_functions is a function such that fit_function(int Nfunc, int Nvar,double *x,int Npar)[i]=f_i
//Nvariables is not required in this function but you need to be consistend in the declaration of fit_function and x
//it returns a[i][value,error]
double* linear_fit_Nf(int* Npoints, double** x, double** y, fit_type fit_info) {
    int Nfunc = fit_info.N;
    int Nvar = fit_info.Nvar;
    int Npar = fit_info.Npar;

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
    if (fit_info.covariancey == true) {
        int en_tot = 0;
        for (int n = 0;n < Nfunc;n++)
            for (int e = 0;e < Npoints[n];e++)
                en_tot += 1;

        double** df_value = (double**)malloc(sizeof(double*) * en_tot);
        for (int n = 0;n < Nfunc;n++){
            for (int e = 0;e < Npoints[n];e++){
                df_value[count] = fit_info.linear_function(n, Nvar, x[count], Npar);
                count++;
            }
        }

        for (int j = 0;j < Npar;j++) {
            for (int i = 0;i < en_tot;i++)
                for (int ii = 0;ii < en_tot;ii++)
                    beta[j] += (df_value[i][j]) * fit_info.cov1[i][ii] * y[ii][0];

            for (int k = 0;k < Npar;k++) {
                for (int i = 0;i < en_tot;i++)
                    for (int ii = 0;ii < en_tot;ii++)
                        alpha[j][k] += df_value[i][j] * fit_info.cov1[i][ii] * df_value[ii][k];
            }
        }
        
        free_2(en_tot, df_value);
    }
    else {
        for (n = 0;n < Nfunc;n++) {
            for (e = 0;e < Npoints[n];e++) {
                i = e + count;
                X = fit_info.linear_function(n, Nvar, x[i], Npar);
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
    }
    if (Npar == 1) {
        C = (double**)malloc(sizeof(double*) * 1);
        C[0] = (double*)malloc(sizeof(double) * 1);
        C[0][0] = 1. / alpha[0][0];
    }
    else
        C = matrix_inverse(Npar, alpha);

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



/*********************************************************************************
 * x[ensemble][variable number] ,   y[ensemble][0=mean,1=error],
 * fun(index_function,Nvariables,variables[], Nparameters,parameters[])
 * the function return an array[Nparameter]  with the value of the parameters that minimise the chi2
 * devorder can be negative, in that case each parameter has is own h[i]=Parameter[i] *h
 ***********************************************************************************/
non_linear_fit_result non_linear_fit_Nf(int N, int* ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess, fit_type fit_info) {

    double lambda = fit_info.lambda;
    double acc = fit_info.acc;
    std::vector< double > h;
    if (fit_info.h.h.size() == 1) {
        h = std::vector< double >(Npar);
        for (int i = 0; i < Npar;i++)
            h[i] = fit_info.h.h[0];
    }
    else if (fit_info.h.h.size() == Npar) { h = fit_info.h; }
    else { printf("non_linear_fit_Nf: h must me a std::vector<double> of 1 o Npar lenght"); exit(1); }

    std::vector< double > Prange = fit_info.Prange;
    int devorder = fit_info.devorder;
    int verbosity = fit_info.verbosity;
    int precision_sum = fit_info.precision_sum;

    double* (*der_fun_Nf_h)(int, int, double*, int, double*, double(int, int, double*, int, double*), std::vector< double >);
    double (*chi2_fun)(int, int*, double**, double**, double*, int, int, double(int, int, double*, int, double*), fit_type);

    bool (*alpha_fun)(double*, double**, int, int*,
        double**, double**, int, int, double*,
        double(int, int, double*, int, double*),
        double* (int, int, double*, int, double*, double(int, int, double*, int, double*), std::vector< double >), std::vector< double >,
        fit_type
        );

    if (fit_info.covariancey) {
        alpha_fun = compute_alpha_cov1;
        chi2_fun = compute_chi_non_linear_Nf_cov1_double;
        if (precision_sum > 0)    chi2_fun = compute_chi_non_linear_Nf_cov1_long_double;

    }
    else {
        alpha_fun = compute_alpha;
        chi2_fun = compute_chi_non_linear_Nf;

        if (precision_sum == 1) {
            chi2_fun = compute_chi_non_linear_Nf_kahan;
        }
        else if (precision_sum > 1) {
            chi2_fun = compute_chi_non_linear_Nf_long;
        }
    }


    if (devorder == -2) {
        der_fun_Nf_h = der_O2_fun_Nf_h_adaptive;
    }
    else if (devorder == 4) {
        der_fun_Nf_h = der_O4_fun_Nf_h;
    }
    else if (devorder == 2) {
        der_fun_Nf_h = der_O2_fun_Nf_h;
    }
    else {
        error(true, 1, "non_linear_fit_Nf", "order of the derivative must be 4 (default) , 2,  -2 to set a step different for each parameter h[i]=P[i]*h    ");
    }
    if (fit_info.linear_fit) {
        // error(fit_info.covariancey, 2, "non_linear_fit_Nf ", "linear fit with covarince not implemented");
        non_linear_fit_result output;
        // double *(linear_function)(int , int  , double *, int );
        // linear_function=fit_info.linear_function;
        output.P = linear_fit_Nf(ensemble, x, y, fit_info);
        output.chi2 = chi2_fun(N, ensemble, x, y, output.P, Nvar, Npar, fun, fit_info);
        return output;
    }

    double** alpha, * X, * beta, ** a, ** C, * sigma;
    int i, j, k, e;
    double f, * fk;
    double chi2, chi2_tmp;
    double* P, * P_tmp, * P_tmp1, res;
    int n, count, Niter = 0;
    int nerror = 0;
    res = 1;

    P = (double*)malloc(Npar * sizeof(double));
    P_tmp = (double*)malloc(Npar * sizeof(double));

    for (j = 0;j < Npar;j++) {
        P[j] = guess[j];
        P_tmp[j] = P[j];

    }


    beta = (double*)calloc(Npar, sizeof(double));
    alpha = (double**)malloc(sizeof(double*) * Npar);
    for (j = 0;j < Npar;j++) {
        alpha[j] = (double*)calloc(Npar, sizeof(double));
    }
    double** alpha_l = (double**)malloc(sizeof(double*) * Npar);
    for (j = 0;j < Npar;j++) {
        alpha_l[j] = (double*)malloc(Npar * sizeof(double));
    }

    chi2 = chi2_fun(N, ensemble, x, y, P_tmp, Nvar, Npar, fun, fit_info);//printf("chi2 in fit function=%g\n",chi2/(ensemble[0]*N-Npar));

    if (verbosity > 0) {
        printf("non_linear_fit_Nf:\ninitial chi2=%g\n", chi2);
        if (verbosity > 1) {
            for (j = 0;j < Npar;j++)
                printf("  P[%d]=%g", j, P_tmp[j]);
            printf("\n");
        }
    }
    //     printf("chi2=%f   res=%.10f P0=%f   P1=%f\n",chi2,res,P[0],P[1]);

    if (fit_info.NM) {//https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
        double** points = double_malloc_2(Npar + 1, Npar);
        double** points1 = double_malloc_2(Npar + 1, Npar);
        double* chi2s = (double*)malloc((Npar + 1) * sizeof(double));
        double* chi2s1 = (double*)malloc((Npar + 1) * sizeof(double));
        int* order = (int*)malloc((Npar + 1) * sizeof(int));
        double* centroid = (double*)malloc((Npar) * sizeof(double));
        double* reflected = (double*)malloc((Npar) * sizeof(double));
        double* expanded = (double*)malloc((Npar) * sizeof(double));
        double* co = (double*)malloc((Npar) * sizeof(double));
        double* ci = (double*)malloc((Npar) * sizeof(double));

        double  r_chi2, e_chi2, co_chi2, ci_chi2;
        for (i = 0;i < Npar + 1;i++) {
            for (j = 0;j < Npar;j++)
                points[i][j] = P[j];
        }

        for (i = 0;i < Npar;i++) {
            points[i][i] += h[i];
        }
        for (i = 0;i < Npar + 1;i++) {
            chi2s[i] = chi2_fun(N, ensemble, x, y, points[i], Nvar, Npar, fun, fit_info);
        }
        if (verbosity >= 2) {
            printf("%s initial points of NM\n", __func__);
            for (i = 0;i < Npar + 1;i++) {
                printf("chi2=%g   params=", chi2s[i]);
                for (j = 0;j < Npar;j++)  printf("%g  ", points[i][j]);
                printf("\n");
            }
        }
        double sd = fit_info.acc + 10;
        while (sd > fit_info.acc) {

            sort_NM(chi2s1, points1, chi2s, points, Npar, order);
            compute_centroid_NM(centroid, points1, Npar);
            compute_reflected_NM(reflected, centroid, points1, Npar, fit_info.alpha);
            if (verbosity >= 1) {
                printf("%s (NM): best point chi2=%g   params=", __func__, chi2s1[0]);
                for (j = 0;j < Npar;j++)  printf("%g  ", points1[0][j]); printf("\n");
                if (verbosity >= 2) {
                    printf("sd=%g\n", sd);
                    for (i = 1;i < Npar + 1;i++) {
                        printf("point %d chi2=%g   params=", i, chi2s1[i]);
                        for (j = 0;j < Npar;j++)  printf("%g  ", points1[i][j]); printf("\n");
                    }
                    if (verbosity >= 3) {
                        printf("previous step points:\n");
                        printf("centroid: ");for (j = 0;j < Npar;j++)  printf("%g  ", centroid[j]); printf("\n");
                        printf("reflected: ");for (j = 0;j < Npar;j++)  printf("%g  ", reflected[j]); printf("\n");
                    }
                }
            }
            r_chi2 = chi2_fun(N, ensemble, x, y, reflected, Nvar, Npar, fun, fit_info);
            if (r_chi2<chi2s1[Npar - 1] && r_chi2>chi2s1[0]) {
                // printf("reflection\n");
                // substitute the last with the reflected
                for (j = 0;j < Npar;j++)
                    points1[Npar][j] = reflected[j];
                chi2s1[Npar] = r_chi2;

                // go to step 1
                sort_NM(chi2s, points, chi2s1, points1, Npar, order);
                sd = compute_sd_NM(centroid, points, Npar);
                continue;
            }
            else if (r_chi2 < chi2s1[0]) {
                // printf("expansion\n");
                compute_expanded_NM(expanded, reflected, centroid, Npar, fit_info.gamma);
                e_chi2 = chi2_fun(N, ensemble, x, y, expanded, Nvar, Npar, fun, fit_info);
                if (e_chi2 < r_chi2) {
                    for (j = 0;j < Npar;j++)
                        points1[Npar][j] = expanded[j];
                    chi2s1[Npar] = e_chi2;
                }
                else {
                    for (j = 0;j < Npar;j++)
                        points1[Npar][j] = reflected[j];
                    chi2s1[Npar] = r_chi2;
                }
                // go to step 1
                sort_NM(chi2s, points, chi2s1, points1, Npar, order);
                sd = compute_sd_NM(centroid, points, Npar);
                continue;
            }
            else {//for sure r_chi2>=chi2s1[Npar-1]
                if (r_chi2 < chi2s1[Npar]) {
                    // printf("contraction out\n");
                    compute_contracted_out_NM(co, reflected, centroid, Npar, fit_info.rho);
                    co_chi2 = chi2_fun(N, ensemble, x, y, co, Nvar, Npar, fun, fit_info);
                    if (co_chi2 < r_chi2) {
                        for (j = 0;j < Npar;j++)
                            points1[Npar][j] = co[j];
                        chi2s1[Npar] = co_chi2;
                        sort_NM(chi2s, points, chi2s1, points1, Npar, order);
                        sd = compute_sd_NM(centroid, points, Npar);
                        continue;
                    }
                    else {
                        // printf("shrink\n");
                        shrink_NM(points1, Npar, fit_info.sigma);
                        for (i = 1;i < Npar + 1;i++) {
                            chi2s[i] = chi2_fun(N, ensemble, x, y, points1[i], Nvar, Npar, fun, fit_info);
                        }
                        sort_NM(chi2s, points, chi2s1, points1, Npar, order);
                        sd = compute_sd_NM(centroid, points, Npar);
                        continue;
                    }

                }
                else {
                    // printf("contraction in\n");
                    compute_contracted_in_NM(ci, points[Npar], centroid, Npar, fit_info.rho);
                    ci_chi2 = chi2_fun(N, ensemble, x, y, ci, Nvar, Npar, fun, fit_info);
                    if (ci_chi2 < chi2s1[Npar]) {
                        for (j = 0;j < Npar;j++)
                            points1[Npar][j] = ci[j];
                        chi2s1[Npar] = ci_chi2;
                        sort_NM(chi2s, points, chi2s1, points1, Npar, order);
                        sd = compute_sd_NM(centroid, points, Npar);
                        continue;
                    }
                    else {
                        // printf("shrink\n");
                        shrink_NM(points1, Npar, fit_info.sigma);
                        for (i = 1;i < Npar + 1;i++) {
                            chi2s[i] = chi2_fun(N, ensemble, x, y, points1[i], Nvar, Npar, fun, fit_info);
                        }
                        sort_NM(chi2s, points, chi2s1, points1, Npar, order);
                        sd = compute_sd_NM(centroid, points, Npar);
                        continue;
                    }
                }



            }


        }
        free(centroid);free(reflected);free(expanded);free(co);free(ci);
        free(order);
        free_2(Npar + 1, points1);
        free(chi2s1);
        for (i = 0;i < Npar;i++)
            P[i] = points[0][i];
        chi2 = chi2s[0];
        free_2(Npar + 1, points);
        free(chi2s);

    }
    else if (fit_info.noderiv) {
        double init_chi2 = 1;
        double loop_chi2 = 20;
        int iterations = 0;
        while (fabs(init_chi2 - loop_chi2) > acc) {
            init_chi2 = chi2;
            for (j = 0;j < Npar;j++) {
                int dir = 1;
                double lmax = 100;
                double scale = pow(2, iterations);
                double lam = lambda / scale;
                if (fit_info.Prange.size() == Npar)
                    lam = fit_info.Prange[j] / scale;
                // while (lam < lmax) {
                if (verbosity > 2) {
                    printf("current set: ");
                    for (int l = 0;l < Npar;l++) printf("%g\t", P_tmp[l]); printf("\t scanning par=%d\n", j);
                }
                for (int dir = -1; dir < 2;dir++) {
                    P_tmp[j] = P[j] + dir * lam;
                    chi2_tmp = chi2_fun(N, ensemble, x, y, P_tmp, Nvar, Npar, fun, fit_info);
                    while (chi2_tmp < chi2 && !isnan(chi2_tmp)) {
                        if (verbosity > 1) printf("found a better chi2: dir %d new=%.8g  old=%.8g  param=%d lambda=%g P=%g   Pnew=%g\n",
                            j, chi2_tmp, chi2, j, lam, P[j], P_tmp[j]);
                        chi2 = chi2_tmp;
                        P[j] = P_tmp[j];
                        P_tmp[j] = P[j] + dir * lam;
                        chi2_tmp = chi2_fun(N, ensemble, x, y, P_tmp, Nvar, Npar, fun, fit_info);
                    }
                    // dir *= -1;
                    // if does not find a better chi2 in both directions
                    // if (dir == 1) lam = 1e+6;
                    P_tmp[j] = P[j];
                }
            }
            loop_chi2 = chi2;
            iterations++;
            if (iterations == 10) { printf("\nmax iter reach!!\n\n"); break; }
        }
        if (verbosity > 2) {
            printf("final set: ");
            for (int l = 0;l < Npar;l++) printf("%g\t", P_tmp[l]); printf("\t chi2=%f\n", chi2);
        }
    }
    else {
        while (res > acc) {
            chi2_tmp = chi2 + 10;
            if (Niter > fit_info.maxiter) { if (verbosity >= 0) printf("Niter=%d of the Levenberg-Marquardt chi2 minimization: exeeds max number\n", Niter); break; }
            Niter++;
            nerror = 0;
            bool computed_alpha = false;

            while (chi2_tmp >= chi2) {  //do {} while() , at least one time is done. if chi is too big chi_tmp=chi+1 = chi 
                if (!computed_alpha) {
                    computed_alpha = alpha_fun(beta, alpha, N, ensemble,
                        x, y, Nvar, Npar, P_tmp, fun, der_fun_Nf_h, h, fit_info);
                }

                for (j = 0;j < Npar;j++) {
                    alpha_l[j][j] = (lambda + 1.) * alpha[j][j];
                    for (k = 0;k < j;k++)
                        alpha_l[j][k] = alpha[k][j];
                    for (k = j + 1;k < Npar;k++)
                        alpha_l[j][k] = alpha[j][k];
                }

                free(P_tmp);
                P_tmp = cholesky_solver_if_possible(Npar, alpha_l, beta);
                if (verbosity > 3) {
                    printf("lambda=%g\n", lambda);
                    for (j = 0;j < Npar;j++) {
                        printf("shift of par[%d]=%g  \n", j, P_tmp[j]);
                    }
                }
                //P_tmp=LU_decomposition_solver(Npar , alpha , beta);
                // check that the result is not nan

                if (Prange.size() == Npar) {
                    bool too_far = false;
                    for (j = 0;j < Npar;j++) {
                        if (fabs(P_tmp[j] + P[j] - guess[j]) > Prange[j]) {
                            too_far = true;
                        }

                    }
                    if (too_far) {
                        lambda *= 10;
                        if (verbosity > 0) {
                            printf(" non_linear_fit_Nf a parameter proposal was rejected because out of range Niter=%d\n", Niter);
                            printf("lambda=%f\n", lambda);
                            printf("chi2=%f chi2_tmp=%f  res=%f \n\n", chi2, chi2_tmp, res);
                            for (j = 0;j < Npar;j++) {
                                printf("P[%d]=%g   Ptmp[%d]=%g\n", j, P[j], j, P_tmp[j]);
                                P_tmp[j] = P[j];
                            }
                        }
                        continue;
                    }
                }

                for (j = 0;j < Npar;j++) {
                    P_tmp[j] += P[j];
                }
                for (j = 0;j < Npar;j++) {
                    if (isnan(P_tmp[j])) {
                        computed_alpha = false;
                        for (int ii = 0;ii < h.size();ii++) h[ii] /= 10.0;
                        P_tmp[j] = P[j] + (2. * rand() / RAND_MAX - 1.) * P[j] / 100.0;// recompute the hessian on this par
                        if (verbosity > 1) printf("parameter is nan, resetting: P[%d]=%g\n", j, P_tmp[j]);
                    }
                }
                chi2_tmp = chi2_fun(N, ensemble, x, y, P_tmp, Nvar, Npar, fun, fit_info);
                if (verbosity > 1) {
                    printf("%s:  new_chi2=%g   old_chi=%g  ", __func__, chi2_tmp, chi2);
                    for (int i = 0;i < Npar;i++) {
                        printf("  P[%d]=%g", i, P_tmp[i]);
                    }
                    printf("\n");
                }

                if (chi2_tmp != chi2_tmp) chi2_tmp = chi2 + 10;


                if (chi2_tmp >= chi2)
                    lambda *= 10;


                // free(C);
                 //error(lambda>1e+15,1,"non_linear_fit_Nf","lambda of the Levenberg-Marquardt chi2 minimization: exeeds 1e+15 lambda=%g",lambda);
                if (lambda > 1e+15) {
                    if (verbosity > 1) printf("lambda of the Levenberg-Marquardt chi2 minimization: exeeds 1e+15 lambda=%g\n RESET lambda=0.001\n", lambda);
                    for (j = 0;j < Npar;j++) {
                        P_tmp[j] = P[j] + (2. * rand() / ((double)RAND_MAX) - 1.) * P[j] / 1000.0;
                        computed_alpha = false;
                        if (verbosity > 0)  printf("  guess[%d]=%g\t", j, guess[j]);
                    }
                    lambda = 0.001;
                    nerror++;
                    if (nerror > 2) {
                        if (verbosity > 0)
                            printf("\n !!!!!!!!!!!!!!!! error:  Impossible to minimise the chi2 with Levenberg-Marquardt for this starting point   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                        if (verbosity > 0)  printf("\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
                        free(P_tmp);
                        free(beta);
                        free_2(Npar, alpha_l);
                        free_2(Npar, alpha);
                        non_linear_fit_result output;
                        output.P = P;
                        output.chi2 = chi2;
                        return output;
                    }
                }


            }
            res = chi2 - chi2_tmp;
            lambda /= 10;
            for (j = 0;j < Npar;j++) {
                P[j] = P_tmp[j];
            }
            chi2 = chi2_tmp;

            if (verbosity >= 2)  printf("\nres=%g  acc=%g  chi2=%g  \n", res, acc, chi2);
            error(res < 0, 2, "non_linear_fit", "The Levenberg-Marquardt accepted a configuration when the chi2 increased");

            computed_alpha = false;

        }
    }

    free_2(Npar, alpha_l);
    free_2(Npar, alpha);

    free(P_tmp);
    free(beta);
    non_linear_fit_result output;
    output.P = P;
    output.chi2 = chi2;
    return output;
}



// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(index_function,Nvariables,variables[], Nparameters,parameters[])
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
double* guess_for_non_linear_fit_Nf(int N, int* ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess, fit_type fit_info) {

    int i, j, jmax, k, e, n;
    double en_tot = 0;
    double chi2, chi2_tmp, chi2_tmp1;
    double* P, * P_tmp;
    double r, rm, gmax = 2;
    rm = (double)RAND_MAX;

    for (n = 0;n < N;n++)
        en_tot += ensemble[n];
    //for(i=0;i<Npar;i++)
      //      guess[i]=((r/rm)-0.5)*2;
    if (fit_info.unstable) for (i = 0;i < Npar;i++) P[i] = guess[i];
    else P = non_linear_fit_Nf(N, ensemble, x, y, Nvar, Npar, fun, guess, fit_info).P;
    chi2 = compute_chi_non_linear_Nf(N, ensemble, x, y, P, Nvar, Npar, fun) / (en_tot - Npar);
    if (!isnan(chi2))    jmax = 3 + ((int)(chi2 * 2));
    else jmax = 35;

    std::mt19937 mt_rand(123);

    if (jmax > 15 || chi2 != chi2 || chi2 < 0) jmax = 35;
    if (jmax <= 3) jmax = 15;

    chi2_tmp1 = chi2;

    for (j = 0;j < jmax;j++) {
        for (int ir = 0;ir < fit_info.repeat_start;ir++) {
            for (i = 0;i < Npar;i++) {
                r = mt_rand() / ((double)mt_rand.max());
                guess[i] = ((r / rm) - 0.5) * (j + 1);
                //  printf("%f\t",guess[i]);
            }
            // printf("\n");
            if (fit_info.unstable) for (i = 0;i < Npar;i++) P_tmp[i] = guess[i];
            else P_tmp = non_linear_fit_Nf(N, ensemble, x, y, Nvar, Npar, fun, guess, fit_info).P;
            chi2_tmp = compute_chi_non_linear_Nf(N, ensemble, x, y, P_tmp, Nvar, Npar, fun) / (en_tot - Npar);
            if (fit_info.verbosity > 2) printf("chi2=%.10f \tchi2_tmp=%.10f   dof=%f\n", chi2, chi2_tmp, en_tot - Npar);
            /*for(i=0;i<Npar;i++)
                printf("P[%d]=%g\t",i,P[i]);
            printf("\n");*/
            if (fabs(chi2 - chi2_tmp) < 1e-3) {
                gmax = gmax * 10;//printf("the chi2 didn't change\n\n");
                free(P_tmp);
            }
            else if (chi2_tmp - chi2 < -1e-3) {
                free(P); P = P_tmp;
                chi2 = chi2_tmp;
                //printf("chi2 smaller founded\n\n");
            }
            else {
                free(P_tmp);//printf("chi2 LARGER\n\n");
                chi2_tmp1 = chi2_tmp;
            }
        }
    }
    //do an other loop if chi2 still large
    if (chi2 > 10 || chi2 != chi2) {
        for (j = 0;j < jmax;j++) {
            for (int ir = 0;ir < fit_info.repeat_start;ir++) {

                for (i = 0;i < Npar;i++) {
                    r = mt_rand() / ((double)mt_rand.max());
                    guess[i] = ((r / rm) - 0.5) * exp(j - jmax / 2);
                    //  printf("%f\t",guess[i]);
                }
                // printf("\n");
                if (fit_info.unstable) for (i = 0;i < Npar;i++) P_tmp[i] = guess[i];
                else P_tmp = non_linear_fit_Nf(N, ensemble, x, y, Nvar, Npar, fun, guess, fit_info).P;
                chi2_tmp = compute_chi_non_linear_Nf(N, ensemble, x, y, P_tmp, Nvar, Npar, fun) / (en_tot - Npar);
                if (fit_info.verbosity > 2) printf("chi2=%.10f \tchi2_tmp=%.10f   dof=%f\n", chi2, chi2_tmp, en_tot - Npar);
                /*for(i=0;i<Npar;i++)
                    printf("P[%d]=%g\t",i,P[i]);
                printf("\n");*/
                if (fabs(chi2 - chi2_tmp) < 1e-3) {
                    gmax = gmax * 10;//printf("the chi2 didn't change\n\n");
                    free(P_tmp);
                }
                else if (chi2_tmp - chi2 < -1e-3) {
                    free(P); P = P_tmp;
                    chi2 = chi2_tmp;
                    //printf("chi2 smaller founded\n\n");
                }
                else {
                    free(P_tmp);//printf("chi2 LARGER\n\n");
                    chi2_tmp1 = chi2_tmp;
                }
            }
        }
    }
    //do an other loop if chi2 still large
    if (chi2 > 10 || chi2 != chi2) {
        for (j = 0;j < jmax;j++) {
            for (int ir = 0;ir < fit_info.repeat_start;ir++) {
                for (i = 0;i < Npar;i++) {
                    r = mt_rand() / ((double)mt_rand.max());
                    guess[i] = ((r / rm) - 0.5) / exp(j - jmax / 2);
                    //  printf("%f\t",guess[i]);
                }
                // printf("\n");
                if (fit_info.unstable) for (i = 0;i < Npar;i++) P_tmp[i] = guess[i];
                else P_tmp = non_linear_fit_Nf(N, ensemble, x, y, Nvar, Npar, fun, guess, fit_info).P;
                chi2_tmp = compute_chi_non_linear_Nf(N, ensemble, x, y, P_tmp, Nvar, Npar, fun) / (en_tot - Npar);
                if (fit_info.verbosity > 2) printf("chi2=%.10f \tchi2_tmp=%.10f   dof=%f\n", chi2, chi2_tmp, en_tot - Npar);
                /*for(i=0;i<Npar;i++)
                *       printf("P[%d]=%g\t",i,P[i]);
                *   printf("\n");*/
                if (fabs(chi2 - chi2_tmp) < 1e-3) {
                    gmax = gmax * 10;//printf("the chi2 didn't change\n\n");
                    free(P_tmp);
                }
                else if (chi2_tmp - chi2 < -1e-3) {
                    free(P); P = P_tmp;
                    chi2 = chi2_tmp;
                    //printf("chi2 smaller founded\n\n");
                }
                else {
                    free(P_tmp);//printf("chi2 LARGER\n\n");
                    chi2_tmp1 = chi2_tmp;
                }
            }
        }
    }
    if (fit_info.verbosity >= 0) printf("final chi2=%f\n", chi2);
    free(guess);
    return P;
}



double compute_chi_non_linear_Nf_cov1(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), double** cov1) {
    double chi2 = 0, f, f1;
    int e, n, count, e1, n1, count1;

    int en_tot = 0;
    for (n = 0;n < N;n++)
        for (e = 0;e < ensemble[n];e++)
            en_tot++;

    double* tmp = (double*)malloc(sizeof(double) * en_tot);
    count = 0;
    for (n = 0;n < N;n++) {
        for (e = 0;e < ensemble[n];e++) {
            tmp[count] = fun(n, Nvar, x[count], Npar, P) - y[count][0];// f1(n,e,N,N1,Nvar,x1[count],Npar,Npar1,P,fun)-y1[count][0];
            // printf("predicted=%g    latt=%g   n=%d  e=%d\n",fun(n, Nvar, x[count], Npar, P), y[count][0], n, e);
            count++;
        }
    }

    for (int i = 0;i < en_tot;i++)
        chi2 += tmp[i] * cov1[i][i] * tmp[i];

    for (int i = 0;i < en_tot;i++)
        for (int j = i + 1;j < en_tot;j++)
            chi2 += 2. * tmp[i] * cov1[i][j] * tmp[j];

    free(tmp);
    return chi2;
}

// double compute_chi_non_linear_Nf_cov(int N, int* ensemble, double** x, double** y, double* P, int Nvar, int Npar, double fun(int, int, double*, int, double*), double** cov) {
//     int n, e;
//     int en_tot = 0;
//     for (n = 0;n < N;n++)
//         for (e = 0;e < ensemble[n];e++)
//             en_tot += 1;

//     double** cov1 = symmetric_matrix_inverse(en_tot, cov);
//     double chi2 = compute_chi_non_linear_Nf_cov1(N, ensemble, x, y, P, Nvar, Npar, fun, cov1);
//     free_2(en_tot, cov1);
//     return chi2;
// }
// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(index_function,Nvariables,variables[], Nparameters,parameters[]), 
//cov[en_tot][en_tot]  is the covariance matrix, with en_tot=sum_i^N ensemble[i],
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
double* non_linear_fit_Nf_covariance(int N, int* ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess, double** cov1, int devorder) {


    double* (*der_fun_Nf_h)(int, int, double*, int, double*, double(int, int, double*, int, double*), std::vector<double>);
    if (devorder == 4) {
        der_fun_Nf_h = der_O4_fun_Nf_h;
    }
    else if (devorder == 2) {
        der_fun_Nf_h = der_O2_fun_Nf_h;
    }
    else {
        error(true, 1, "non_linear_fit_Nf", "order of the derivative must be 4 (default) or 2");
    }
    double** alpha, * X, * beta, ** a, ** C, * sigma;
    int i, j, k, e;
    double f, * fk, f1, * fk1;
    double chi2, chi2_tmp;
    double* P, * P_tmp, lambda, res;
    int n, count, n1, count1, e1, Niter = 0;
    std::vector<double> h(Npar);
    for (int i = 0; i < Npar;i++) h[i] = 1e-5;
    lambda = 0.001;
    res = 1;
    int en_tot = 0;
    for (n = 0;n < N;n++)
        for (e = 0;e < ensemble[n];e++)
            en_tot += 1;
    /*
  int yn=is_it_positive( cov,  en_tot);
  while(yn==1){
      printf("covariance matrix not positive defined adding 0.0001*cov[0][0]*I \n");
      for(i=0;i<en_tot;i++)
           cov[i][i]+=cov[0][0]*1e-12;
      yn=is_it_positive( cov,  en_tot);
      printf("now the matrix is positive defined.  %d\n",yn);
  }
  double **cov1=symmetric_matrix_inverse(en_tot, cov  );
  */

  /*   printf("covariance\n");
     for (j=0;j<en_tot;j++){
         for (i=0;i<en_tot;i++){
         printf("%g\t",cov[j][i]);
         }
         printf("\n");
     }
     printf("inverse\n");
     for (j=0;j<en_tot;j++){
         for (i=0;i<en_tot;i++){
         printf("%g\t",cov1[j][i]);
         }
         printf("\n");
     }*/

     //   double *fkk;
    P = (double*)malloc(Npar * sizeof(double));
    P_tmp = (double*)malloc(Npar * sizeof(double));


    for (j = 0;j < Npar;j++) {
        P[j] = guess[j];
        P_tmp[j] = P[j];
    }


    beta = (double*)calloc(Npar, sizeof(double));
    alpha = (double**)malloc(sizeof(double*) * Npar);
    for (j = 0;j < Npar;j++) {
        alpha[j] = (double*)calloc(Npar, sizeof(double));
    }

    chi2 = compute_chi_non_linear_Nf_cov1(N, ensemble, x, y, P_tmp, Nvar, Npar, fun, cov1);

    chi2_tmp = chi2 + 1;
    double* f_value = (double*)malloc(sizeof(double) * en_tot);
    double** df_value = (double**)malloc(sizeof(double*) * en_tot);

    while (res > 0.001) {
        chi2_tmp = chi2 + 1;
        if (Niter > 200) { printf("Niter=%d of the Levenberg-Marquardt chi2 minimization: exeeds max number\n", Niter); break; }
        Niter++;
        while (chi2_tmp > chi2) {

            count = 0;
            for (n = 0;n < N;n++) {
                for (e = 0;e < ensemble[n];e++) {
                    f_value[count] = fun(n, Nvar, x[e], Npar, P);
                    df_value[count] = der_fun_Nf_h(n, Nvar, x[e], Npar, P, fun, h);
                    count++;
                }
            }
            for (j = 0;j < Npar;j++) {
                /*for(i=0;i<en_tot;i++)
                    beta[j]+=(y[i][0]-f_value[i])*cov1[i][i]*df_value[i][j];*/
                for (i = 0;i < en_tot;i++)
                    for (int ii = 0;ii < en_tot;ii++)
                        beta[j] += (y[i][0] - f_value[i]) * cov1[i][ii] * df_value[ii][j];


                for (k = j;k < Npar;k++) {
                    //for(i=0;i<en_tot;i++)
                    //    alpha[j][k]+=df_value[i][j]*cov1[i][i]*df_value[i][k];
                    for (i = 0;i < en_tot;i++)
                        for (int ii = 0;ii < en_tot;ii++)
                            alpha[j][k] += df_value[i][j] * cov1[i][ii] * df_value[ii][k];
                }
            }
            for (i = 0;i < en_tot;i++)
                free(df_value[i]);

            double* beta1 = (double*)calloc(Npar, sizeof(double));
            double** alpha1 = (double**)malloc(Npar * sizeof(double));
            for (j = 0;j < Npar;j++)
                alpha1[j] = (double*)calloc(Npar, sizeof(double));

            count = 0;
            for (n = 0;n < N;n++) {
                for (e = 0;e < ensemble[n];e++) {

                    f = fun(n, Nvar, x[e + count], Npar, P);
                    fk = der_fun_Nf_h(n, Nvar, x[e + count], Npar, P, fun, h);
                    count1 = 0;
                    for (n1 = 0;n1 < N;n1++) {
                        for (e1 = 0;e1 < ensemble[n1];e1++) {
                            f1 = fun(n1, Nvar, x[e1 + count1], Npar, P);
                            fk1 = der_fun_Nf_h(n1, Nvar, x[e1 + count1], Npar, P, fun, h);

                            for (j = 0;j < Npar;j++) {
                                beta1[j] += (y[e1 + count1][0] - f1) * cov1[e1 + count1][e + count] * fk[j];
                                for (k = j;k < Npar;k++) {
                                    alpha1[j][k] += fk1[j] * cov1[e1 + count1][e + count] * fk[k];
                                }


                            }
                            free(fk1);
                        }
                        count1 += ensemble[n];
                    }
                    free(fk);
                }
                count += ensemble[n];
            }

            for (j = 0;j < Npar;j++)
                error(fabs(beta[j] - beta1[j]) > 1e-6, 1, "fit", "beta differs:   %d    %g    %g", j, beta[j], beta1[j]);
            for (j = 0;j < Npar;j++)
                for (k = j;k < Npar;k++)
                    error(fabs(alpha[j][k] - alpha1[j][k]) > 1e-6, 1, "fit", "alpha differs:   %d  %d   %g    %g", j, k, alpha[j][k], alpha1[j][k]);

            free(beta1);


            for (j = 0;j < Npar;j++) {
                alpha[j][j] *= (lambda + 1.);
                for (k = 0;k < j;k++)
                    alpha[j][k] = alpha[k][j];

            }
            free(P_tmp);
            P_tmp = cholesky_solver_if_possible(Npar, alpha, beta);

            for (j = 0;j < Npar;j++)
                P_tmp[j] += P[j];

            chi2_tmp = compute_chi_non_linear_Nf_cov1(N, ensemble, x, y, P_tmp, Nvar, Npar, fun, cov1);
            if (chi2_tmp != chi2_tmp) chi2_tmp = chi2 + 1;



            if (chi2_tmp > chi2)
                lambda *= 10;

            for (j = 0;j < Npar;j++) {

                beta[j] = 0;
                for (k = 0;k < Npar;k++) {
                    alpha[j][k] = 0;
                }
            }
            if (lambda > 1e+15) {
                printf("lambda of the Levenberg-Marquardt chi2 minimization: exeeds 1e+15 lambda=%g\n RESET lambda=0.001\n", lambda);
                lambda = 0.001;
            }


        }
        res = chi2 - chi2_tmp;

        chi2 = chi2_tmp;
        lambda /= 10;
        for (j = 0;j < Npar;j++) {
            P[j] = P_tmp[j];

        }
    }

    //free_2(en_tot,cov1);
    for (j = 0;j < Npar;j++) {
        free(alpha[j]);
    }
    free(f_value);free(df_value);

    free(P_tmp);
    free(alpha);free(beta);
    return P;
}


/*********************************************************************************
 * x[ensemble][variable number] ,   y[ensemble][0=mean,1=error],
 * fun(index_function,Nvariables,variables[], Nparameters,parameters[])
 * the function return an array[Nparameter]  with the value of the parameters that minimise the chi2
 * devorder can be negative, in that case each parameter has is own h[i]=Parameter[i] *h
 ***********************************************************************************/
 //double* non_linear_fit_Nf(int N, int* ensemble, double** x, double** y, int Nvar, int Npar,  double fun(int, int, double*, int, double*) , double* guess, double lambda, double acc, double h, std::vector< double > Prange, int devorder, int verbosity , int precision_sum )
double* non_linear_fit_Nf_cov(int N, int* ensemble, double** x, double** y, int Nvar, int Npar, double fun(int, int, double*, int, double*), double* guess, fit_type fit_info, double** cov1) {


    double lambda = fit_info.lambda;
    double acc = fit_info.acc;
    std::vector< double > h;
    if (fit_info.h.h.size() == 1) {
        h = std::vector< double >(Npar);
        for (int i = 0; i < Npar;i++)
            h[i] = fit_info.h.h[0];
    }
    else if (fit_info.h.h.size() == Npar) { h = fit_info.h; }
    else { printf("non_linear_fit_Nf: h must me a std::vector<double> of 1 o Npar lenght"); exit(1); }
    std::vector< double > Prange = fit_info.Prange;
    int devorder = fit_info.devorder;
    int verbosity = fit_info.verbosity;
    int precision_sum = fit_info.precision_sum;

    double* (*der_fun_Nf_h)(int, int, double*, int, double*, double(int, int, double*, int, double*), std::vector<double>);
    double (*chi2_fun)(int, int*, double**, double**, double*, int, int, double(int, int, double*, int, double*), double**);
    chi2_fun = compute_chi_non_linear_Nf_cov1;

    if (devorder == -2) {
        der_fun_Nf_h = der_O2_fun_Nf_h_adaptive;
    }
    else if (devorder == 4) {
        der_fun_Nf_h = der_O4_fun_Nf_h;
    }
    else if (devorder == 2) {
        der_fun_Nf_h = der_O2_fun_Nf_h;
    }
    else {
        error(true, 1, "non_linear_fit_Nf", "order of the derivative must be 4 (default) , 2,  -2 to set a step different for each parameter h[i]=P[i]*h    ");
    }

    double** alpha, * X, * beta, ** a, ** C, * sigma;
    int i, j, k, e;
    int en_tot = 0;
    for (int n = 0;n < N;n++)
        for (e = 0;e < ensemble[n];e++)
            en_tot += 1;
    double* f_value = (double*)malloc(sizeof(double) * en_tot);
    double** df_value = (double**)malloc(sizeof(double*) * en_tot);

    double chi2, chi2_tmp;
    double* P, * P_tmp, * P_tmp1, res;
    int n, count, Niter = 0;
    int nerror = 0;
    res = 1;

    P = (double*)malloc(Npar * sizeof(double));
    P_tmp = (double*)malloc(Npar * sizeof(double));

    for (j = 0;j < Npar;j++) {
        P[j] = guess[j];
        P_tmp[j] = P[j];

    }

    beta = (double*)calloc(Npar, sizeof(double));
    alpha = (double**)malloc(sizeof(double*) * Npar);
    for (j = 0;j < Npar;j++) {
        alpha[j] = (double*)calloc(Npar, sizeof(double));
    }
    double** alpha_l = (double**)malloc(sizeof(double*) * Npar);
    for (j = 0;j < Npar;j++) {
        alpha_l[j] = (double*)malloc(Npar * sizeof(double));
    }

    chi2 = chi2_fun(N, ensemble, x, y, P_tmp, Nvar, Npar, fun, cov1);//printf("chi2 in fit function=%g\n",chi2/(ensemble[0]*N-Npar));

    if (verbosity > 0) {
        printf("non_linear_fit_Nf:\ninitial chi2=%g\n", chi2);
        if (verbosity > 1) {
            for (j = 0;j < Npar;j++)
                printf("  P[%d]=%g", j, P_tmp[j]);
            printf("\n");
        }
    }
    if (fit_info.noderiv) {
        double init_chi2 = 1;
        double loop_chi2 = 20;

        for (j = 0;j < Npar;j++) {
            int dir = 1;
            double lmax = 100;

            // while (lam < lmax) {
            if (verbosity > 2) {
                printf("current set: ");
                for (int l = 0;l < Npar;l++) printf("%g\t", P_tmp[l]); printf("\t scanning par=%d\n", j);
            }
            double lam = lambda;
            if (fit_info.Prange.size() == Npar)
                lam = fit_info.Prange[j];
            while (lam > h[j]) {
                // for (int iterations = 0; fit_info.Prange[j]/pow(10,iterations) > h[j];iterations++) {
                    // lam=fit_info.Prange[j]/pow(10,iterations);
                for (int dir = -1; dir < 2;dir++) {
                    P_tmp[j] = P[j] + dir * lam;
                    chi2_tmp = chi2_fun(N, ensemble, x, y, P_tmp, Nvar, Npar, fun, cov1);
                    while (chi2_tmp < chi2 && !isnan(chi2_tmp)) {
                        if (verbosity > 1) printf("found a better chi2: dir %d new=%.8g  old=%.8g  param=%d lambda=%g P=%g   Pnew=%g\n",
                            j, chi2_tmp, chi2, j, lam, P[j], P_tmp[j]);
                        chi2 = chi2_tmp;
                        P[j] = P_tmp[j];
                        P_tmp[j] = P[j] + dir * lam;
                        chi2_tmp = chi2_fun(N, ensemble, x, y, P_tmp, Nvar, Npar, fun, cov1);
                    }

                    // dir *= -1;
                    // if does not find a better chi2 in both directions
                    // if (dir == 1) lam = 1e+6;
                    P_tmp[j] = P[j];
                }
                lam /= 2.;
            }

            loop_chi2 = chi2;
            // iterations++;
            // if (iterations == 9) break;
        }
        if (verbosity > 2) {
            printf("final set: ");
            for (int l = 0;l < Npar;l++) printf("%g\t", P_tmp[l]); printf("\t chi2=%f\n", chi2);
        }

    }
    else {
        while (res > acc) {
            chi2_tmp = chi2 + 10;
            if (Niter > 200) { if (verbosity > 0) printf("Niter=%d of the Levenberg-Marquardt chi2 minimization: exeeds max number\n", Niter); break; }
            Niter++;
            nerror = 0;
            bool computed_alpha = false;

            while (chi2_tmp >= chi2) {  //do {} while()   , at least one time is done. if chi is too big chi_tmp=chi+1 = chi 
                if (!computed_alpha) {
                    count = 0;
                    for (n = 0;n < N;n++) {
                        for (e = 0;e < ensemble[n];e++) {
                            f_value[count] = fun(n, Nvar, x[e], Npar, P);
                            df_value[count] = der_fun_Nf_h(n, Nvar, x[e], Npar, P, fun, h);
                            count++;
                        }
                    }
                    for (j = 0;j < Npar;j++) {
                        for (i = 0;i < en_tot;i++)
                            for (int ii = 0;ii < en_tot;ii++)
                                beta[j] += (y[i][0] - f_value[i]) * cov1[i][ii] * df_value[ii][j];


                        for (k = j;k < Npar;k++) {
                            for (i = 0;i < en_tot;i++)
                                for (int ii = 0;ii < en_tot;ii++)
                                    alpha[j][k] += df_value[j][i] * cov1[i][ii] * df_value[ii][k];
                        }
                    }
                    for (i = 0;i < en_tot;i++)
                        free(df_value[i]);
                    computed_alpha = true;
                }

                for (j = 0;j < Npar;j++) {
                    alpha_l[j][j] = (lambda + 1.) * alpha[j][j];
                    for (k = 0;k < j;k++)
                        alpha_l[j][k] = alpha[k][j];
                    for (k = j + 1;k < Npar;k++)
                        alpha_l[j][k] = alpha[j][k];
                }

                free(P_tmp);
                P_tmp = cholesky_solver_if_possible(Npar, alpha_l, beta);

                if (Prange.size() == Npar) {
                    bool too_far = false;
                    for (j = 0;j < Npar;j++) {
                        if (fabs(P_tmp[j] + P[j] - guess[j]) > Prange[j]) {
                            too_far = true;
                        }
                    }
                    if (too_far) {
                        lambda *= 10;
                        if (verbosity > 0) {
                            printf(" non_linear_fit_Nf a parameter proposal was rejected because out of range Niter=%d\n", Niter);
                            printf("lambda=%f\n", lambda);
                            printf("chi2=%f chi2_tmp=%f  res=%f \n\n", chi2, chi2_tmp, res);
                            for (j = 0;j < Npar;j++) {
                                printf("P[%d]=%g   Ptmp[%d]=%g\n", j, P[j], j, P_tmp[j]);
                                P_tmp[j] = P[j];
                            }
                        }
                        continue;
                    }
                }

                for (j = 0;j < Npar;j++)     P_tmp[j] += P[j];

                for (j = 0;j < Npar;j++) {
                    if (isnan(P_tmp[j])) {
                        computed_alpha = false;
                        for (int ii = 0; ii < Npar;ii++) h[ii] /= 10.0;
                        P_tmp[j] = P[j] + (2. * rand() / RAND_MAX - 1.) * P[j] / 100.0;// recompute the hessian on this par
                        if (verbosity > 1) printf("parameter is nan, resetting: P[%d]=%g\n", j, P_tmp[j]);
                    }
                }
                chi2_tmp = chi2_fun(N, ensemble, x, y, P_tmp, Nvar, Npar, fun, cov1);
                if (verbosity > 1) {
                    printf("%s:  new_chi2=%g   old_chi=%g  ", __func__, chi2_tmp, chi2);
                    for (int i = 0;i < Npar;i++) {
                        printf("  P[%d]=%g", i, P_tmp[i]);
                    }
                    printf("\n");
                }

                if (chi2_tmp != chi2_tmp) chi2_tmp = chi2 + 10;

                if (chi2_tmp >= chi2)
                    lambda *= 10;

                if (lambda > 1e+15) {
                    if (verbosity > 1) printf("lambda of the Levenberg-Marquardt chi2 minimization: exeeds 1e+15 lambda=%g\n RESET lambda=0.001\n", lambda);
                    for (j = 0;j < Npar;j++) {
                        P_tmp[j] = P[j] + (2. * rand() / ((double)RAND_MAX) - 1.) * P[j] / 1000.0;
                        computed_alpha = false;
                        if (verbosity > 0)  printf("  guess[%d]=%g\t", j, guess[j]);
                    }
                    lambda = 0.001;
                    nerror++;
                    if (nerror > 2) {
                        if (verbosity > 0)
                            printf("\n !!!!!!!!!!!!!!!! error:  Impossible to minimise the chi2 with Levenberg-Marquardt for this starting point   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                        if (verbosity > 0)  printf("\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
                        for (j = 0;j < Npar;j++)      free(alpha[j]);
                        free(P_tmp);
                        free(alpha);free(beta);
                        return P;
                    }
                }
            }
            res = chi2 - chi2_tmp;
            lambda /= 10;
            for (j = 0;j < Npar;j++) {
                P[j] = P_tmp[j];
                //             printf("P[%d]= %f    \t",j,P[j]);
            }
            if (verbosity >= 2)  printf("\nres=%g  acc=%g  chi2=%g  \n", res, acc, chi2);
            error(res < 0, 2, "non_linear_fit", "The Levenberg-Marquardt accepted a configuration when the chi2 increased");
            chi2 = chi2_tmp;
            for (j = 0;j < Npar;j++) {
                beta[j] = 0;
                for (k = 0;k < Npar;k++) {
                    alpha[j][k] = 0;
                }
            }
            computed_alpha = false;

        }
    }

    for (j = 0;j < Npar;j++) {
        free(alpha[j]);
        free(alpha_l[j]);
    }
    free(f_value);free(df_value);
    free(P_tmp);
    free(alpha);free(beta);free(alpha_l);
    return P;
}

double rtsafe(void (*funcd)(double, int, double*, double*, double*), int Npar, double* P, double x1, double x2,
    double xacc)
    //Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
    //between x1 and x2. The root, returned as the function value rtsafe, will be refined until
    //its accuracy is known within ±xacc. funcd is a user-supplied routine that returns both the
    //function value and the first derivative of the function.
{
    int j;
    double df, dx, dxold, f, fh, fl;
    double temp, xh, xl, rts;

    (*funcd)(x1, Npar, P, &fl, &df);
    (*funcd)(x2, Npar, P, &fh, &df);
    error(((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)), 1, "rtsafe", "Root must be bracketed in rtsafe");
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    if (fl < 0.0) {// Orient the search so that f (xl) < 0.
        xl = x1;
        xh = x2;
    }
    else {
        xh = x1;
        xl = x2;
    }
    rts = 0.5 * (x1 + x2);// Initialize the guess for root,
    dxold = fabs(x2 - x1);  //the “stepsize before last,”
    dx = dxold;// and the last step.

    (*funcd)(rts, Npar, P, &f, &df);
    for (j = 1;j <= MAXIT;j++) { // Loop over allowed iterations.
        if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0)//Bisect if Newton out of range,
            || (fabs(2.0 * f) > fabs(dxold * df))) {//or not decreasing fast enough.
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if (xl == rts) return rts;// Change in root is negligible.
        }
        else {//Newton step acceptable. Take it.
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if (temp == rts) return rts;
        }
        if (fabs(dx) < xacc) return rts; //Convergence criterion.

        (*funcd)(rts, Npar, P, &f, &df);
        if (f < 0.0)
            xl = rts;
        else
            xh = rts;
    }
    error(1 > 0, 1, "rtsafe", "Maximum number of iterations exceeded in rtsafe");
    return 0.0;
}

double rtbis(double (*func)(double, double, int, double*), double input, int Npar, double* P, double x1, double x2, double xacc)
//Using bisection, find the root of a function func known to lie between x1 and x2. The root,
//returned as rtbis, will be refined until its accuracy is ±xacc.
{

    int j;
    double dx, f, fmid, xmid, rtb;

    f = (*func)(input, x1, Npar, P);
    fmid = (*func)(input, x2, Npar, P);

    error(f * fmid >= 0.0, 1, "rtbis", "Root must be bracketed for bisection in rtbis f(x1)=%f   f(x2)=%f", f, fmid);
    rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);// Orient the search so that f>0
    for (j = 1;j <= MAXIT;j++) {// lies at x+dx.
        fmid = (*func)(input, xmid = rtb + (dx *= 0.5), Npar, P); //Bisection loop.
        if (fmid <= 0.0) rtb = xmid;
        if (fabs(dx) < xacc || fmid == 0.0) return rtb;
    }
    error(1 > 0, 1, "rtbis", "Too many bisections in rtbis");
    return 0.0;
}


/*****************************************************************************************************************
 * Using bisection, find the root of a function func-input known to lie between x1 and x2. The root,
 * returned as rtbis, will be refined until its accuracy is ±xacc.
 * /func return different values for different n
 * it solves function=input
 * if root not in range:
 *  Pedanticness==0 try other range if faliure return NaN (default)
 *  Pedanticness==1 return NaN
 *  Pedanticness==2 try other range if faliure exit
 *  Pedanticness==3  exit
 *****************************************************************************************************************/
double rtbis_func_eq_input(double (*func)(int, int, double*, int, double*), int n, int Nvar, double* x, int Npar, double* P, int ivar, double input, double x1, double x2, double xacc, int Pedanticness) {

    double* xt = (double*)malloc(sizeof(double) * Nvar);
    for (int i = 0;i < Nvar; i++)
        xt[i] = x[i];


    int j;
    double dx, f, fmid, xmid, rtb;
    xt[ivar] = x1;
    f = (*func)(n, Nvar, xt, Npar, P) - input;
    xt[ivar] = x2;
    fmid = (*func)(n, Nvar, xt, Npar, P) - input;

    error(f * fmid >= 0.0 && Pedanticness >= 2, 1, "rtbis", "Root must be bracketed for bisection in rtbis f(%g)=%f   f(%g)=%f", x1, f, x2, fmid);
    bool try_range = Pedanticness == 0 || Pedanticness == 2;
    bool return_nan = Pedanticness < 2;
    if (f * fmid >= 0.0 && try_range) {
        //printf("#f(x)  x\n");
        for (int i = 0;i < 100; i++) {
            xt[ivar] = x1 - 2 * (x2 - x1) * i;
            f = (*func)(n, Nvar, xt, Npar, P) - input;
            //printf("%f   %f\n",f, xt[ivar] );
            if (f * fmid <= 0.0) break;
        }
        if (f * fmid <= 0.0) {
            x1 = xt[ivar];
            f = (*func)(n, Nvar, xt, Npar, P) - input;
        }
        else {
            xt[ivar] = x1;
            f = (*func)(n, Nvar, xt, Npar, P) - input;
            for (int i = 0;i < 100; i++) {
                xt[ivar] = x2 + 2 * (x2 - x1) * i;
                fmid = (*func)(n, Nvar, xt, Npar, P) - input;
                //printf("%f   %f\n",f, xt[ivar] );
                if (f * fmid <= 0.0) break;
            }
            x2 = xt[ivar];
        }
        //xt[ivar]=x2;
        //fmid=(*func)(n,Nvar,xt,Npar,P)-input;
    }
    error(f * fmid >= 0.0 && !return_nan, 1, "rtbis", "Root must be bracketed for bisection in rtbis f(x1)=%f   f(x2)=%f", f, fmid);
    if (f * fmid >= 0.0 && return_nan) { printf("error: bisect f(%g)=%f   f(%g)=%f\n", x1, f, x2, fmid); return NAN; }

    rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);// Orient the search so that f>0
    for (j = 1;j <= MAXIT;j++) {// lies at x+dx.
        //printf(" x=%f \n",xt[ivar]   );
        xt[ivar] = rtb + (dx *= 0.5);
        fmid = (*func)(n, Nvar, xt, Npar, P) - input; //Bisection loop.
        //printf(" x+dx=%f fmind=%f   %d\n",xt[ivar], fmid ,MAXIT  );
        if (fmid <= 0.0) rtb = xt[ivar];
        if (fabs(dx) < xacc || fmid == 0.0) {
            free(xt);
            return rtb;
        }
    }
    free(xt);

    error(!return_nan, 1, "rtbis", "Too many bisections in rtbis");
    //printf("Too many bisections in rtbis\n");
    return NAN;
}


inline double derivative_O4_ivar(double (*func)(int, int, double*, int, double*), int n, int Nvar, double* x, int Npar, double* P, int ivar, double h) {

    x[ivar] = x[ivar] - 2. * h;      // x-2*h
    double df = func(n, Nvar, x, Npar, P);

    x[ivar] = x[ivar] + h;        // x-h
    df -= 8 * func(n, Nvar, x, Npar, P);

    x[ivar] = x[ivar] + 2 * h;     // x+h
    df += 8 * func(n, Nvar, x, Npar, P);

    x[ivar] = x[ivar] + h;      //h+2h
    df -= func(n, Nvar, x, Npar, P);

    x[ivar] = x[ivar] - 2. * h;//you need to leave the parameter as it was before you move to the next parameter
    df /= (12. * h);

    return df;
}

inline double derivative_O2_ivar(double (*func)(int, int, double*, int, double*), int n, int Nvar, double* x, int Npar, double* P, int ivar, double h) {

    x[ivar] = x[ivar] - h;      // x-h
    double df = -func(n, Nvar, x, Npar, P);

    x[ivar] = x[ivar] + 2 * h;        // x+h
    df += func(n, Nvar, x, Npar, P);


    x[ivar] = x[ivar] - h;//you need to leave the parameter as it was before you move to the next parameter
    df /= (2. * h);

    return df;
}

/*Using the Newton-Raphson method, find the root of a function known to lie in the interval  [ x1, x2]. The root rtnew*t will be refined until its accuracy is known within ±xacc. funcd
 i *s a user-supplied routine that returns both the function value and the first derivative of the
 function at the point x.
 *
 *
 */
double rtnewt(double (*func)(int, int, double*, int, double*), int n, int Nvar, double* x, int Npar, double* P, int ivar, double input, double xstart, double xmin, double xmax, float xacc, int JMAX, double h) {
    void nrerror(char error_text[]);
    int j;
    double df, dx, f, rtn;
    double* xt = (double*)malloc(sizeof(double) * Nvar);
    for (int i = 0;i < Nvar; i++)
        xt[i] = x[i];
    xt[ivar] = xstart;

    for (j = 1;j <= JMAX;j++) {

        f = (*func)(n, Nvar, xt, Npar, P) - input;
        df = derivative_O2_ivar(func, n, Nvar, xt, Npar, P, ivar, h);
        dx = f / df;
        xt[ivar] -= dx;
        printf("j=%d   k=%g    f=%g  df=%g  dx=%g\n", j, xt[ivar], f, df, dx);

        if ((xmin - xt[ivar]) * (xt[ivar] - xmax) < 0.0)
            error(0 == 0, 1, "rtnwt", "Jumped out of brackets [%g,%g] in rtnewt x=%g ", xmin, xmax, xt[ivar]);
        if (fabs(dx) < xacc) return xt[ivar];

    }

    error(0 == 0, 1, "rtnwt", "Maximum number of iterations exceeded in rtnewt");
    //nrerror("Maximum number of iterations exceeded in rtnewt");
    return 0.0;

}

/***********************************************************************
 * Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
 * b etween x1 and x2. The root, returned as the function va*lue rtsafe, will be refined until
 * its accuracy is known within ±xacc. funcd is a user-supplied routine that returns both the
 * function value and the first derivative of the function.
 ***********************************************************************/
 // float rtsafe(void (*funcd)(float, float *, float *), float x1, float x2,
 //              float xacc)
double  rtsafe(double (*func)(int, int, double*, int, double*), int n, int Nvar, double* x, int Npar, double* P, int ivar, double input, double x1, double x2, float xacc, int JMAX, double h) {
    int j;
    double df, dx, dxold, f, fh, fl;
    double temp, xh, xl, rts;
    double* xt = (double*)malloc(sizeof(double) * Nvar);
    for (int i = 0;i < Nvar; i++)
        xt[i] = x[i];


    xt[ivar] = x1;
    fl = (*func)(n, Nvar, xt, Npar, P) - input;
    //     df=derivative_ivar(func,n,Nvar,xt,Npar,P,ivar,h );
    xt[ivar] = x2;
    fh = (*func)(n, Nvar, xt, Npar, P) - input;
    //     df=derivative_ivar(func,n,Nvar,xt,Npar,P,ivar,h );


    error((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0), 2, "rtsafe", "Root must be bracketed in rtsafe  f(%g)=%g  f(%g)=%g", x1, fl, x2, fh);
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    if (fl < 0.0) {
        //         Orient the search so that f (xl) < 0.
        xl = x1;
        xh = x2;
    }
    else {
        xh = x1;
        xl = x2;
    }
    rts = 0.5 * (x1 + x2);           // Initialize the guess for root,
    dxold = fabs(x2 - x1);         // the “stepsize before last,”
    dx = dxold;                  //and the last step.

    xt[ivar] = rts;
    f = (*func)(n, Nvar, xt, Npar, P) - input;
    df = derivative_O2_ivar(func, n, Nvar, xt, Npar, P, ivar, h);
    for (j = 1;j <= MAXIT;j++) {//     Loop over allowed iterations.
        if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) || (fabs(2.0 * f) > fabs(dxold * df))) { //      Bisect if Newton out of range,   or not decreasing fast enough.
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if (xl == rts) return rts; // Change in root is negligible.
        }
        else {  //Newton step acceptable. Take it.
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if (temp == rts) return rts;
        }
        if (fabs(dx) < xacc) return rts; //    Convergence criterion.
        xt[ivar] = rts;   // The one new function evaluation per iteration.
        f = (*func)(n, Nvar, xt, Npar, P) - input;
        df = derivative_O2_ivar(func, n, Nvar, xt, Npar, P, ivar, h);
        if (f < 0.0)//            Maintain the bracket on the root.
            xl = rts;
        else
            xh = rts;
    }
    error(0 == 0, 1, "rtsafe:", "Maximum number of iterations exceeded in rtsafe");
    return 0.0;//Never get here.
}


/*

double f1( int Nvar, double *x,int Npar,double  *P)
{
    return P[0]*log(P[0]*x[0])+P[0]*log(P[1]*x[1]);
}
double *f1k( int Nvar, double *x,int Npar,double  *P)
{
    double *r;
    r=(double*) malloc(sizeof(double)*Npar);
    r[0]=log(P[0]*x[0])+1+log(P[1]*x[1]);
    r[1]=(P[0]/P[1]);
    return r;

}

int main(){

    double *P;
    double **x,**y;
    int i;

    x=(double**) malloc(sizeof(double*)*4);
    y=(double**) malloc(sizeof(double*)*4);
    for (i=0;i<4;i++){
       x[i]=(double*) malloc(sizeof(double)*2);
       y[i]=(double*) malloc(sizeof(double)*2);
    }
   x[0][0]=1;x[1][0]=4;x[2][0]=7;x[3][0]=10;
   x[0][1]=2;x[1][1]=3;x[2][1]=5;x[3][1]=9;

   y[0][0]=2;y[1][0]=5;y[2][0]=8;y[3][0]=11;
   y[0][1]=0.5;y[1][1]=2;y[2][1]=2;y[3][1]=2;
    int ensemble=4;
    int Nvar=2;
    int Npar=2;


//P=non_linear_fit( ensemble,x, y , Nvar,  Npar,  f1 , f1k);
P=non_linear_fit( ensemble,x, y , Nvar,  Npar,  f1 );
printf("P[0]=%f\n",P[0]);
printf("P[1]=%f\n",P[1]);

    return 0;
}
*/

/*Gnuplot result
gnuplot> f(x)=A*log(A*x)
gnuplot> fit f(x) 'tmp' i 0 u 1:2:3 yerr via A
iter      chisq       delta/lim  lambda   A
   0 1.0935357306e+02   0.00e+00  1.35e+01    4.343868e+00
   1 9.0384567941e+00  -1.11e+06  1.35e+00    3.012440e+00
   2 2.7981042656e+00  -2.23e+05  1.35e-01    2.553096e+00
   3 2.7944989111e+00  -1.29e+02  1.35e-02    2.541241e+00
   4 2.7944983222e+00  -2.11e-02  1.35e-03    2.541377e+00
iter      chisq       delta/lim  lambda   A

After 4 iterations the fit converged.
final sum of squares of residuals : 2.7945
rel. change during last iteration : -2.10747e-07

degrees of freedom    (FIT_NDF)                        : 3
rms of residuals      (FIT_STDFIT) = sqrt(WSSR/ndf)    : 0.965142
variance of residuals (reduced chisquare) = WSSR/ndf   : 0.931499
p-value of the Chisq distribution (FIT_P)              : 0.424406

Final set of parameters            Asymptotic Standard Error
=======================            ==========================
A               = 2.54138          +/- 0.1895       (7.455%)
*/
