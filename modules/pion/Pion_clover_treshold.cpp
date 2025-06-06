#define Pion_clover_treshold_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "global.hpp"
#include "resampling.hpp"
//#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"

// #include <omp.h> 
static void init_fit(int N, struct header* head, int Njack, struct data_jack* gJ, int Nvar, int Npar, int** en, int* en_tot, double**** x, double*** sigmax, double** chi2m, double** rm, double*** r, double*** fit, double**** y, double** chi2, double**** C, double threshold) {
    int imoms, imomt, imom0, iG, i, n, e, j;
    int count;
    *en_tot = 0;

    *en = (int*)calloc(N, sizeof(int));

    for (e = 0;e < ensembles;e++) {
        for (n = 0;n < N;n++) {
            if (gJ[e].M_PS_jack[0][Njack - 1] * gJ[e].w0[Njack - 1] < threshold)
                (*en)[n] += 1;
        }
    }



    for (n = 0;n < N;n++) {
        *en_tot += (*en)[n];
    }

    *x = double_malloc_3(Njack, *en_tot, Nvar);//(double**) malloc(sizeof(double*)*(*en_tot));
    *sigmax = double_malloc_2((*en)[0], Nvar);

    //*chi2m=(double*) malloc(sizeof(double)*(Npar));
    *rm = (double*)malloc(sizeof(double) * (Njack));

    *fit = (double**)malloc(sizeof(double*) * (*en_tot)); printf("tot  = %d\n", *en_tot);

    *r = (double**)malloc(sizeof(double*) * (Npar));
    for (i = 0;i < Npar;i++) {
        (*r)[i] = (double*)malloc(sizeof(double) * Njack);
    }

    *chi2 = (double*)malloc(sizeof(double) * Njack);
    (*y) = (double***)malloc(sizeof(double**) * Njack);
    (*C) = (double***)malloc(sizeof(double**) * Njack);

    for (j = 0;j < Njack;j++) {
        (*y)[j] = (double**)malloc(sizeof(double*) * (*en_tot));
        count = 0;
        for (n = 0;n < N;n++) {
            for (i = 0;i < ensembles;i++) {
                if (gJ[i].M_PS_jack[0][Njack - 1] * gJ[i].w0[Njack - 1] < threshold) {
                    (*y)[j][count] = (double*)malloc(sizeof(double) * 2);
                    count++;
                }
            }

        }
    }
}
static struct fit_result close_fit(int N, struct header* head, int Njack, struct data_jack* gJ, int Npar, int** en, int* en_tot, double**** x, double*** sigmax, double** chi2m, double** rm, double*** r, double*** fit, double**** y, double** chi2, double**** C, double threshold) {
    int imoms, imomt, imom0, iG, i, n, e, j;
    int count;

    free(*chi2m);
    free(*rm);


    count = 0;
    for (n = 0;n < N;n++) {
        for (i = 0;i < (*en)[n];i++) {
            if (gJ[i].M_PS_jack[0][Njack - 1] * gJ[i].w0[Njack - 1] < threshold)
                free((*fit)[i + count]);
            //free((*x)[i+count]);
        }
        count += (*en)[n];
    }
    free(*fit);
    free_3(Njack, *en_tot, *x);
    free_2((*en)[0], *sigmax);


    for (j = 0;j < Njack;j++) {

        count = 0;
        for (n = 0;n < N;n++) {
            for (i = 0;i < (*en)[n];i++) {
                free((*y)[j][i + count]);
            }
            count += (*en)[n];
        }
        free((*y)[j]);
    }
    free((*y));
    free(*en);
    struct fit_result fit_out;
    fit_out.Njack = Njack;
    fit_out.P = (double**)malloc(sizeof(double*) * Npar);
    fit_out.chi2 = (double*)malloc(sizeof(double*) * Njack);
    fit_out.C = (double***)malloc(sizeof(double**) * Npar);
    for (i = 0;i < Npar;i++) {
        fit_out.P[i] = (double*)malloc(sizeof(double*) * Njack);
        for (j = 0;j < Njack;j++) {
            fit_out.P[i][j] = (*r)[i][j];
        }
        free((*r)[i]);
        fit_out.C[i] = (double**)malloc(sizeof(double*) * Npar);
        for (n = 0;n < Npar;n++) {
            fit_out.C[i][n] = (double*)malloc(sizeof(double) * Njack);
            for (j = 0;j < Njack;j++) {
                fit_out.C[i][n][j] = (*C)[j][i][n];
            }
        }
    }
    for (j = 0;j < Njack;j++) {
        fit_out.chi2[j] = (*chi2)[j];
    }
    free(*r);
    free(*chi2);
    for (j = 0;j < Njack;j++) {
        for (n = 0;n < Npar;n++)
            free((*C)[j][n]);
        free((*C)[j]);
    }
    free(*C);
    return fit_out;
}




struct fit_result fit_Mpi_fw_chiral_FVE_clover_treshold(struct database_file_jack* jack_files, struct header* head, int Njack, int*** mass_index, struct data_jack* gJ, struct fit_type fit_info, double threshold) {
    double*** y, *** x, ** sigmax, ** r, * chi2, * tmp, * rm, * chi2m, ** fit, *** C;
    int i, j, e, im;
    int Npar = fit_info.Npar;
    int Nvar = 5;//fit_info.Nvar;//m_l, w0,M_PS^2,f_PS
    int ik1 = 0, ik2 = 0;
    int n, count, N = fit_info.N;
    int* en;

    int en_tot = 0;


    double* guess = (double*)malloc(sizeof(double) * Npar);
    for (i = 0;i < Npar;i++)
        guess[i] = rand();

    init_fit(N, head, Njack, gJ, Nvar, Npar, &en, &en_tot, &x, &sigmax, &chi2m, &rm, &r, &fit, &y, &chi2, &C, threshold);

    count = 0;
    for (n = 0;n < N;n++) {
        for (e = 0;e < ensembles;e++) {
            im = mass_index[e][ik2][ik1];
            if (gJ[e].M_PS_jack[0][Njack - 1] * gJ[e].w0[Njack - 1] < threshold) {
                if (n == 0) {
                    for (j = 0;j < Njack;j++) {
                        rm[j] = gJ[e].M_PS_jack[im][j] * gJ[e].w0[j];
                        rm[j] *= rm[j];
                    }
                    fit[count] = mean_and_error(jack_files[0].sampling, Njack, rm);
                }
                if (n == 1) {
                    for (j = 0;j < Njack;j++) {
                        rm[j] = gJ[e].f_PS_jack[im][j] * gJ[e].w0[j];
                    }
                    fit[count] = mean_and_error(jack_files[0].sampling, Njack, rm);
                }

                for (j = 0;j < jack_tot;j++) {
                    y[j][count][0] = rm[j];
                    y[j][count][1] = fit[count][1];
                }
                //x[e+count]=(double*) malloc(sizeof(double)*Nvar);

                count++;
            }
            //   printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }

    }

    double KM, Kf, K;




    //#pragma omp parallel for  private(tmp,i,count,n,e,im,x)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2)
    for (j = 0;j < Njack;j++) {
        count = 0;
        //x=(double**) malloc(sizeof(double*)*(en_tot));
        for (n = 0;n < N;n++) {
            for (e = 0;e < ensembles;e++) {
                im = mass_index[e][ik2][ik1];
                if (gJ[e].M_PS_jack[0][Njack - 1] * gJ[e].w0[Njack - 1] < threshold) {

                    x[j][count][0] = head[e].k[head[e].nk + ik2] * gJ[e].w0[j] / gJ[e].Zp[j];//ml*w0
                    x[j][count][1] = gJ[e].w0[j];//w0
                    x[j][count][2] = gJ[e].M_PS_jack[im][j] * gJ[e].M_PS_jack[im][j];//MPS^2
                    x[j][count][3] = gJ[e].f_PS_jack[im][j];//f_PS
                    x[j][count][4] = double(head[e].l1);//f_PS
                    count++;
                }
            }

        }
    }

    count = 0;
    for (n = 0;n < 1;n++) {
        for (e = 0;e < en[n];e++) {
            if (gJ[e].M_PS_jack[0][Njack - 1] * gJ[e].w0[Njack - 1] < threshold) {
                for (int v = 0;v < Nvar;v++) {
                    for (j = 0;j < Njack;j++)
                        rm[j] = x[j][count][v];
                    tmp = mean_and_error(jack_files[0].sampling, Njack, rm);
                    if (fabs(tmp[1]) < 1e-6) { printf("e=%d    v=%d   %g +- %g\n", e, v, tmp[0], tmp[1]); tmp[1] = tmp[0] / 1.0e+8; }
                    sigmax[count][v] = tmp[1];
                    free(tmp);
                }
                count++;
            }
        }

    }

    double** yy = double_malloc_2(en_tot, Njack);
    for (i = 0;i < en_tot;i++)
        for (j = 0;j < Njack;j++)
            yy[i][j] = y[j][i][0];

    double** cov = covariance(jack_files[0].sampling, en_tot, Njack, yy);
    free_2(en_tot, yy);

    double** yx = double_malloc_2(en_tot + Nvar * en[0], Njack);
    for (i = 0;i < en_tot;i++)
        for (j = 0;j < Njack;j++)
            yx[i][j] = y[j][i][0];

    for (i = en_tot;i < en_tot + Nvar * en[0];i++) {
        for (j = 0;j < Njack;j++) {
            yx[i][j] = x[j][(i - en_tot) / Nvar][(i - en_tot) % Nvar];
        }
    }

    double** cov_yx = covariance(jack_files[0].sampling, en_tot + Nvar * en[0], Njack, yx);
    free_2(en_tot + Nvar * en[0], yx);
    int yn = is_it_positive(cov_yx, en_tot + Nvar * en[0]);
    while (yn == 1) {
        printf("covariance matrix not positive defined adding eps*cov[0][0]*I \n");
        for (i = 0;i < en_tot + Nvar * en[0];i++)
            cov_yx[i][i] += cov_yx[0][0] * 1e-12;
        yn = is_it_positive(cov_yx, en_tot + Nvar * en[0]);
        printf("now the matrix is positive defined.  %d\n", yn);
    }
    double** cov_yx1 = symmetric_matrix_inverse(en_tot + Nvar * en[0], cov_yx);

    yn = is_it_positive(cov, en_tot);
    while (yn == 1) {
        printf("covariance matrix not positive defined adding eps*cov[0][0]*I \n");
        for (i = 0;i < en_tot;i++)
            cov[i][i] += cov[0][0] * 1e-12;
        yn = is_it_positive(cov, en_tot);
        printf("now the matrix is positive defined.  %d\n", yn);
    }
    double** cov1 = symmetric_matrix_inverse(en_tot, cov);

    guess = guess_for_non_linear_fit_Nf(N, en, x[0], y[0], Nvar, Npar, fit_info.function, guess);


    /*  double *guess1=(double*)  malloc(sizeof(double)*(Npar+Nvar*en[0]));
      for(i=0;i<Npar;i++)
          guess1[i]=guess[i];

      for(i=Npar;i<Npar+Nvar*en[0];i++)
          guess1[i]=x[0][ (i-Npar)/Nvar ][  (i-Npar)%Nvar  ];
      free(guess);
      guess=guess1;


      guess1=non_linear_fit_Nf_sigmax_covariance( N, en ,x[0], sigmax, y[0] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
      free(guess);
      guess=guess1;
      */
      //#pragma omp parallel for  private(tmp,i,count,n,e,im)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2,C,x,cov,cov_yx1,cov1)
    for (j = 0;j < Njack;j++) {
        //if (j==0){     }
        tmp = non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess).P;
        //tmp=non_linear_fit_Nf_sigmax( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_iterative( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_covariance( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
        //tmp=non_linear_fit_Nf_covariance(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess ,cov1);


        chi2[j] = compute_chi_non_linear_Nf(N, en, x[j], y[j], tmp, Nvar, Npar, fit_info.function) / (en_tot - Npar);
        //for (i=Npar;i<Npar+en[0]*Nvar;i++)
         //   printf("%.10f   ",fabs(tmp[i]-x[j][(i-Npar)/Nvar ][(i-Npar)%Nvar ] ) );
                //printf("\n");
        //printf("jacknife=%d of %d   chi2=%g\n",j, Njack,chi2[j]);

        C[j] = covariance_non_linear_fit_Nf(N, en, x[j], y[j], tmp, Nvar, Npar, fit_info.function);
        for (i = 0;i < Npar;i++) {
            r[i][j] = tmp[i];
        }
        if (j == Njack - 1) {
            printf("w0/a[fm]     mu*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf           (Mpi^2/fpi^2)* (Kf^2/KM^2)\n");
            count = 0;
            for (n = 0;n < N;n++) {
                printf("#function %d\n", n);
                for (e = 0;e < ensembles;e++) {
                    if (gJ[e].M_PS_jack[0][Njack - 1] * gJ[e].w0[Njack - 1] < threshold) {
                        KM = fit_info.function(2, Nvar, x[j][count], Npar, tmp);
                        Kf = fit_info.function(3, Nvar, x[j][count], Npar, tmp);
                        double* tmp1 = (double*)malloc(sizeof(double) * Njack);
                        for (int jj = 0;jj < Njack;jj++) {
                            tmp1[jj] = gJ[e].M_PS_jack[0][jj] / gJ[e].f_PS_jack[0][jj];
                            tmp1[jj] *= Kf / KM;
                            tmp1[jj] *= tmp1[jj];
                        }
                        double* tmp2 = mean_and_error(jack_files[0].sampling, Njack, tmp1);
                        if (n == 0)
                            K = KM * KM;
                        else if (n == 1)
                            K = Kf;
                        printf("%.5f     %.5f      %.5f   %.5f     %.5f \t\t %.5f  %.5f\n", x[j][count][1], x[j][count][0], fit[count][0] / K, fit[count][1] / K, K, tmp2[0], tmp2[1]);
                        free(tmp2);free(tmp1);
                        count++;
                    }
                }

            }

        }
        free(tmp);

    }
    free_2(en_tot, cov);
    free_2(en_tot, cov1);
    free_2(en_tot + Nvar * en[0], cov_yx);
    free_2(en_tot + Nvar * en[0], cov_yx1);

    chi2m = mean_and_error(jack_files[0].sampling, Njack, chi2);
    printf("$\\chi^2/dof=%f+-%f$  \n", chi2m[0], chi2m[1]);



    struct fit_result fit_out = close_fit(N, head, Njack, gJ, Npar, &en, &en_tot, &x, &sigmax, &chi2m, &rm, &r, &fit, &y, &chi2, &C, threshold);

    return fit_out;

}








struct fit_result fit_Mpi_fwMpi4_chiral_FVE_clover_threshold(struct database_file_jack* jack_files, struct header* head, int Njack, int*** mass_index, struct data_jack* gJ, struct fit_type fit_info, double threshold) {
    double*** y, *** x, ** sigmax, ** r, * chi2, * tmp, * rm, * chi2m, ** fit, *** C;
    int i, j, e, im;
    int Npar = fit_info.Npar;
    int Nvar = 5;//fit_info.Nvar;//m_l, w0,M_PS^2,f_PS
    int ik1 = 0, ik2 = 0;
    int n, count, N = fit_info.N;
    int* en;

    int en_tot = 0;


    double* guess = (double*)malloc(sizeof(double) * Npar);
    for (i = 0;i < Npar;i++)
        guess[i] = rand();
    init_fit(N, head, Njack, gJ, Nvar, Npar, &en, &en_tot, &x, &sigmax, &chi2m, &rm, &r, &fit, &y, &chi2, &C, threshold);
    printf("en_tot=%d  =%d +  %d\n", en_tot, en[0], en[1]);
    count = 0;
    for (n = 0;n < N;n++) {
        for (e = 0;e < ensembles;e++) {
            im = mass_index[e][ik2][ik1];
            if (gJ[e].M_PS_jack[0][Njack - 1] * gJ[e].w0[Njack - 1] < threshold) {
                printf("ensemble=%d  count=%d   \n", e, count);
                if (n == 0) {
                    for (j = 0;j < Njack;j++) {
                        rm[j] = gJ[e].M_PS_jack[im][j] * gJ[e].w0[j];
                        rm[j] *= rm[j];
                    }
                    fit[count] = mean_and_error(jack_files[0].sampling, Njack, rm);
                }
                if (n == 1) {
                    for (j = 0;j < Njack;j++) {
                        rm[j] = gJ[e].f_PS_jack[im][j];
                        rm[j] *= gJ[e].M_PS_jack[im][j] * gJ[e].M_PS_jack[im][j] * gJ[e].M_PS_jack[im][j] * gJ[e].M_PS_jack[im][j];
                        rm[j] *= gJ[e].w0[j] * gJ[e].w0[j] * gJ[e].w0[j] * gJ[e].w0[j] * gJ[e].w0[j];
                    }
                    fit[count] = mean_and_error(jack_files[0].sampling, Njack, rm);
                }

                for (j = 0;j < jack_tot;j++) {
                    y[j][count][0] = rm[j];
                    y[j][count][1] = fit[count][1];
                }
                count++;
            }

        }

    }

    double KM, Kf, K;




    //#pragma omp parallel for  private(tmp,i,count,n,e,im,x)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2)
    for (j = 0;j < Njack;j++) {
        count = 0;
        //x=(double**) malloc(sizeof(double*)*(en_tot));
        for (n = 0;n < N;n++) {
            for (e = 0;e < ensembles;e++) {
                im = mass_index[e][ik2][ik1];
                if (gJ[e].M_PS_jack[0][Njack - 1] * gJ[e].w0[Njack - 1] < threshold) {

                    x[j][count][0] = head[e].k[head[e].nk + ik2] * gJ[e].w0[j] / gJ[e].Zp[j];//ml*w0
                    x[j][count][1] = gJ[e].w0[j];//w0
                    x[j][count][2] = gJ[e].M_PS_jack[im][j] * gJ[e].M_PS_jack[im][j];//MPS^2
                    x[j][count][3] = gJ[e].f_PS_jack[im][j];//f_PS
                    x[j][count][4] = double(head[e].l1);//f_PS
                    count++;
                }
            }

        }
    }

    count = 0;
    for (n = 0;n < 1;n++) {
        for (e = 0;e < ensembles;e++) {
            if (gJ[e].M_PS_jack[0][Njack - 1] * gJ[e].w0[Njack - 1] < threshold) {

                for (int v = 0;v < Nvar;v++) {
                    for (j = 0;j < Njack;j++)
                        rm[j] = x[j][count][v];
                    tmp = mean_and_error(jack_files[0].sampling, Njack, rm);
                    //if (fabs(tmp[1])<1e-6) {printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] ); tmp[1]=tmp[0]/1.0e+8; }
                    sigmax[count][v] = tmp[1];
                    free(tmp);
                }
                count++;
            }
        }

    }

    double** yy = double_malloc_2(en_tot, Njack);
    for (i = 0;i < en_tot;i++)
        for (j = 0;j < Njack;j++)
            yy[i][j] = y[j][i][0];

    double** cov = covariance(jack_files[0].sampling, en_tot, Njack, yy);
    free_2(en_tot, yy);

    double** yx = double_malloc_2(en_tot + Nvar * en[0], Njack);
    for (i = 0;i < en_tot;i++)
        for (j = 0;j < Njack;j++)
            yx[i][j] = y[j][i][0];

    for (i = en_tot;i < en_tot + Nvar * en[0];i++) {
        for (j = 0;j < Njack;j++) {
            yx[i][j] = x[j][(i - en_tot) / Nvar][(i - en_tot) % Nvar];
        }
    }

    double** cov_yx = covariance(jack_files[0].sampling, en_tot + Nvar * en[0], Njack, yx);
    free_2(en_tot + Nvar * en[0], yx);
    int yn = is_it_positive(cov_yx, en_tot + Nvar * en[0]);
    while (yn == 1) {
        printf("covariance matrix not positive defined adding eps*cov[0][0]*I \n");
        for (i = 0;i < en_tot + Nvar * en[0];i++)
            cov_yx[i][i] += cov_yx[0][0] * 1e-12;
        yn = is_it_positive(cov_yx, en_tot + Nvar * en[0]);
        printf("now the matrix is positive defined.  %d\n", yn);
    }
    double** cov_yx1 = symmetric_matrix_inverse(en_tot + Nvar * en[0], cov_yx);

    yn = is_it_positive(cov, en_tot);
    while (yn == 1) {
        printf("covariance matrix not positive defined adding eps*cov[0][0]*I \n");
        for (i = 0;i < en_tot;i++)
            cov[i][i] += cov[0][0] * 1e-12;
        yn = is_it_positive(cov, en_tot);
        printf("now the matrix is positive defined.  %d\n", yn);
    }
    double** cov1 = symmetric_matrix_inverse(en_tot, cov);

    guess = guess_for_non_linear_fit_Nf(N, en, x[0], y[0], Nvar, Npar, fit_info.function, guess);


    /*  double *guess1=(double*)  malloc(sizeof(double)*(Npar+Nvar*en[0]));
      for(i=0;i<Npar;i++)
          guess1[i]=guess[i];

      for(i=Npar;i<Npar+Nvar*en[0];i++)
          guess1[i]=x[0][ (i-Npar)/Nvar ][  (i-Npar)%Nvar  ];
      free(guess);
      guess=guess1;


      guess1=non_linear_fit_Nf_sigmax_covariance( N, en ,x[0], sigmax, y[0] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
      free(guess);
      guess=guess1;
      */
      //#pragma omp parallel for  private(tmp,i,count,n,e,im)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2,C,x,cov,cov_yx1,cov1)
    for (j = 0;j < Njack;j++) {
        //if (j==0){     }
        tmp = non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess).P;
        //tmp=non_linear_fit_Nf_sigmax( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_iterative( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_covariance( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
        //tmp=non_linear_fit_Nf_covariance(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess ,cov1);


        chi2[j] = compute_chi_non_linear_Nf(N, en, x[j], y[j], tmp, Nvar, Npar, fit_info.function) / (en_tot - Npar);
        //for (i=Npar;i<Npar+en[0]*Nvar;i++)
         //   printf("%.10f   ",fabs(tmp[i]-x[j][(i-Npar)/Nvar ][(i-Npar)%Nvar ] ) );
                //printf("\n");
        //printf("jacknife=%d of %d   chi2=%g\n",j, Njack,chi2[j]);

        C[j] = covariance_non_linear_fit_Nf(N, en, x[j], y[j], tmp, Nvar, Npar, fit_info.function);
        for (i = 0;i < Npar;i++) {
            r[i][j] = tmp[i];
        }
        if (j == Njack - 1) {
            printf("w0/a[fm]     mu*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf           (Mpi^2/fpi^2)* (Kf^2/KM^2)\n");
            count = 0;
            for (n = 0;n < N;n++) {
                printf("#function %d\n", n);
                for (e = 0;e < ensembles;e++) {
                    if (gJ[e].M_PS_jack[0][Njack - 1] * gJ[e].w0[Njack - 1] < threshold) {

                        KM = fit_info.function(2, Nvar, x[j][count], Npar, tmp);
                        Kf = fit_info.function(3, Nvar, x[j][count], Npar, tmp);
                        double* tmp1 = (double*)malloc(sizeof(double) * Njack);
                        for (int jj = 0;jj < Njack;jj++) {
                            tmp1[jj] = gJ[e].M_PS_jack[0][jj] / gJ[e].f_PS_jack[0][jj];
                            tmp1[jj] *= Kf / KM;
                            tmp1[jj] *= tmp1[jj];
                        }
                        double* tmp2 = mean_and_error(jack_files[0].sampling, Njack, tmp1);
                        if (n == 0)
                            K = KM * KM;
                        else if (n == 1)
                            K = Kf;
                        printf("%.5f     %.5f      %.5f   %.5f     %.5f \t\t %.5f  %.5f\n", x[j][count][1], x[j][count][0], fit[count][0] / K, fit[count][1] / K, K, tmp2[0], tmp2[1]);
                        free(tmp2);free(tmp1);
                        count++;
                    }
                }

            }

        }
        free(tmp);

    }
    free_2(en_tot, cov);
    free_2(en_tot, cov1);
    free_2(en_tot + Nvar * en[0], cov_yx);
    free_2(en_tot + Nvar * en[0], cov_yx1);


    chi2m = mean_and_error(jack_files[0].sampling, Njack, chi2);
    printf("$\\chi^2/dof=%f+-%f$  \n", chi2m[0], chi2m[1]);



    struct fit_result fit_out = close_fit(N, head, Njack, gJ, Npar, &en, &en_tot, &x, &sigmax, &chi2m, &rm, &r, &fit, &y, &chi2, &C, threshold);

    return fit_out;

}
