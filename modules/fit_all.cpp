#define fit_all_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <iterator>
#include <random>

#include "tower.hpp"
#include "fit_all.hpp"
#include "non_linear_fit.hpp"
#include "mutils.hpp"
#include "resampling.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// print output
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_fit_output(char** argv, std::vector<data_all> gjack, struct fit_type fit_info,
    const char* label, struct fit_result fit_out, int* en, double*** x, double*** y, std::vector<int> myen) {
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;

    char namefile[NAMESIZE];
    FILE* f;

    mysprintf(namefile, NAMESIZE, "%s/%s_fit_data.txt", argv[3], label);
    f = open_file(namefile, "w+");
    int count = 0;
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < en[n];e++) {
            for (int v = 0; v < Nvar;v++) {
                fprintf(f, " %g\t ", x[Njack - 1][e + count][v]);
            }
            fprintf(f, " %g   %g  \t ", y[Njack - 1][e + count][0], y[Njack - 1][e + count][1]);
            fprintf(f, " %d   \n ", n);
        }
        count += en[n];
        fprintf(f, "\n\n");
    }

    fclose(f);
    ////////////// parameter print and correlation matrix
    mysprintf(namefile, NAMESIZE, "%s/%s_fit_P.tex", argv[3], label);
    f = open_file(namefile, "w+");
    fprintf(f, "\\begin{gather}\n");
    fprintf(f, "\\chi^2/d.o.f.=%g \\\\ \n", fit_out.chi2[Njack - 1]);//error_jackboot(argv[1], Njack, fit_out.chi2)
    for (int i = 0;i < Npar;i++) {
        fprintf(f, "P[%d]=%g\\pm (%.2g) \\\\ \n", i, fit_out.P[i][Njack - 1], error_jackboot(argv[1], Njack, fit_out.P[i]));
    }
    fprintf(f, "\\end{gather}\n");
    double** cov = covariance(argv[1], Npar, Njack, fit_out.P);
    fprintf(f, "{\\tiny\\begin{gather}\n C=\\begin{pmatrix}\n");
    for (int i = 0;i < fit_info.Npar;i++) {
        for (int k = 0;k < i;k++)
            cov[i][k] /= sqrt(cov[i][i] * cov[k][k]);
        for (int k = i + 1;k < fit_info.Npar;k++)
            cov[i][k] /= sqrt(cov[i][i] * cov[k][k]);
    }


    for (int i = 0;i < fit_info.Npar;i++) {
        for (int j = 0;j < fit_info.Npar;j++) {
            if (j == 0)  fprintf(f, "%.3g", cov[i][j]);
            else       fprintf(f, "& %.3g", cov[i][j]);
        }
        if (i != fit_info.Npar) fprintf(f, "\\\\ \n");
        else fprintf(f, "\n");
    }
    fprintf(f, "\\end{pmatrix}\n\\end{gather}}\n");
    free_2(Npar, cov);
    fclose(f);



}


struct fit_result fit_all_data(char** argv, std::vector<data_all> gjack, double lhs_fun(int, int, int, std::vector<data_all>, struct fit_type), struct fit_type fit_info, const char* label) {
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;
    auto myen = fit_info.myen;
    ////// allocation
    int* en = (int*)malloc(sizeof(int) * fit_info.N);// we need to init en and en_tot to allocate the other 
    for (int e = 0;e < fit_info.N; e++) { en[e] = myen.size(); }
    int en_tot = 0;      for (int n = 0;n < N;n++) { en_tot += en[n]; }// total data to fit

    double*** y = double_malloc_3(Njack, en_tot, 2);// 2 is mean/error
    double*** x = double_malloc_3(Njack, en_tot, Nvar);
    fit_result fit_out = malloc_fit(fit_info);
    double* guess = (double*)malloc(sizeof(double) * Npar);
    double** fit = (double**)malloc(sizeof(double*) * Njack);//result of the fit, the other dimension is allocated by the function non_linear_fit_Nf()

    printf("///// fit name:  %s \n", label);
    ////// allocation end
    /////////////// init
    std::mt19937 mt_rand(123);
    if (fit_info.guess.size() == 0) {
        for (int i = 0;i < Npar;i++)
            guess[i] = mt_rand() / ((double)mt_rand.max());//rand()/((double)RAND_MAX);
    }
    else {
        for (int i = 0;i < fit_info.guess.size();i++)
            guess[i] = fit_info.guess[i];
        for (int i = fit_info.guess.size();i < Npar;i++)
            guess[i] = 1;
    }



    //init x
    for (int j = 0;j < Njack;j++) {
        int count = 0;
        for (int n = 0;n < N;n++) {
            for (int e = 0;e < en[n];e++) {
                for (int i = 0; i < fit_info.Nvar; i++) {
                    x[j][count][i] = fit_info.x[i][e][j];
                }

                count++;
            }
        }
    }


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
            printf(" %g   %g   %g\n", x[Njack - 1][e + count][0], y[Njack - 1][e + count][0], y[Njack - 1][e + count][1]);
            free(tmpj);
        }
        count += en[n];
    }
    //////  init end

    ///////////////// the fit 
    // scan the parameter of the fit with the last jack
    if (fit_info.guess.size() == 0 || fit_info.repeat_start > 1) {
        guess = guess_for_non_linear_fit_Nf(N, en, x[Njack - 1], y[Njack - 1], Nvar, Npar, fit_info.function, guess);
    }
    if (fit_info.mean_only == false) {
        for (int j = Njack - 1;j >= 0;j--) {

            double a = timestamp();
            fit[j] = non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess, fit_info);
            fit_out.chi2[j] = compute_chi_non_linear_Nf(N, en, x[j], y[j], fit[j], Nvar, Npar, fit_info.function) / (en_tot - Npar);

            if (fit_info.verbosity > 0) {
                printf("jack =%d  chi2/dof=%g   chi2=%g   time=%g   \n", j, fit_out.chi2[j], fit_out.chi2[j] * (en_tot - Npar), timestamp() - a);
                if (fit_info.verbosity > 1) {
                    for (int i = 0;i < Npar;i++)
                        printf("P[%d]=%g \t", i, fit[j][i]);
                    printf("\n");
                }
            }

            // we do not need the covariance of the fit, it will be computed with jackboot
            //double **C=covariance_non_linear_fit_Nf(N, en,x[j], y[j],fit[j] , Nvar,  Npar, fit_info.function );            
            //for(int i=0;i<Npar;i++)
            //    for(int k=0;k<Npar;k++)
            //        fit_out.C[i][k][j]=C[i][k];
            //free_2(Npar, C);


        }
    }
    else if (fit_info.mean_only == true) {
        int j = Njack - 1;
        fit[j] = non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess, fit_info);
        fit_out.chi2[j] = compute_chi_non_linear_Nf(N, en, x[j], y[j], fit[j], Nvar, Npar, fit_info.function) / (en_tot - Npar);
        // for the other jackboot add a white noise to the mean
        for (j = Njack - 2;j >= 0;j--) {
            fit[j] = (double*)malloc(sizeof(double) * fit_info.Npar);

            for (int i = 0; i < fit_info.Npar;i++) {
                double noise = fit[Njack - 1][i] * mt_rand() / ((double)mt_rand.max() * 100);
                fit[j][i] = fit[Njack - 1][i] * noise;
                fit_out.chi2[j] = fit_out.chi2[Njack - 1] * noise;
            }
        }

    }
    for (int i = 0;i < Npar;i++)
        for (int j = 0;j < Njack;j++)
            fit_out.P[i][j] = fit[j][i];


    /////////////////////////////////////////////////////////////////////writing the result
    print_fit_output(argv, gjack, fit_info, label, fit_out, en, x, y, myen);

    ////// free
    free(en);
    //free(chi2j);
    free_3(Njack, en_tot, y);
    free_3(Njack, en_tot, x);
    free(guess);
    free_2(Njack, fit);

    ////// free end


    return fit_out;

}
