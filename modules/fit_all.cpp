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
#include "linear_fit.hpp"
#include "non_linear_fit.hpp"
#include "mutils.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"


// void data_single::cut_confs(int N) {
//     if (N >= Njack - 1) { printf("error: trying to cut the configuration up to N= %d > Njack = %d", N, Njack); exit(1); }

//     double** jack1 = (double**)malloc(sizeof(double*) * Nobs);//[obs][jack]
//     for (int i = 0; i < Nobs;i++) {
//         jack1[i] = (double*)malloc(sizeof(double) * (N + 1));
//         for (int j = 0; j < N ;j++) {
//             jack1[i][j] = jack[i][j];
//         }
//         jack1[i][N] = jack[i][Njack - 1];
//         free(jack[i]);
//     }

//     free(jack);

//     jack = jack1;
//     // for (int i = 0; i < Nobs;i++)
//     //     jack[i] = jack1[i];

//     Njack = N + 1;
// }

// void data_all::cut_confs(int N) {
//     for (int i = 0; i < ens;i++) {
//         en[i].cut_confs(N);
//     }
// }

void data_all::add_fit(struct fit_result fit_out) {

    int N = Nfits + 1;
    Nfits = N;

    fit_result* fit_tmp = (struct fit_result*)malloc(sizeof(struct fit_result) * N);
    for (int i = 0;i < N - 1;i++) {
        fit_tmp[i] = malloc_copy_fit_result(fits[i]);
        // fit_tmp[i].name=fits[i].name;
        for (int j = 0;j < fits[i].Npar;j++) {
            free(fits[i].P[j]);
            for (int k = 0;k < fits[i].Npar;k++) {
                free(fits[i].C[j][k]);
            }
            free(fits[i].C[j]);
        }
        free(fits[i].chi2);
        free(fits[i].P);
        free(fits[i].C);
        if (i == N - 2) { free(fits); }
    }

    fit_tmp[N - 1] = malloc_copy_fit_result(fit_out);
    fits = fit_tmp;


}

void data_single::init_error() {
    if (!init_err) {
        errors = (double*)malloc(sizeof(double) * Nobs);
        computed = (bool*)malloc(sizeof(bool) * Nobs);
        if (strcmp(resampling.c_str(), "jack") == 0)
            compute_mean_and_error = mean_and_error_jack_biased;
        else if (strcmp(resampling.c_str(), "boot") == 0)
            compute_mean_and_error = mean_and_error_boot;
        for (int j = 0;j < Nobs;j++) {
            computed[j] = false;
        }
        init_err = true;
    }
}

void data_single::reset_error() {
    if (init_err) { free(errors); free(computed); }
}

double data_single::error_jack(int i) {
    if (!init_err) init_error();
    if (!computed[i]) {
        double* r = compute_mean_and_error(Njack, jack[i]);
        errors[i] = r[1];
        free(r);
        computed[i] = true;
    }
    return errors[i];
}


void data_all::create_generalised_resampling() {
    // if the length is the same return dataj
    int same = 0;
    for (int i = 0; i < ens;i++) {
        printf("(jacks en %d)=%d\n", i, en[i].Njack);
        if (en[i].Njack == en[0].Njack)
            same++;
    }

    if (same == ens) {
        std::cout << "all the files have the same number of jack/boot , do nothing" << std::endl;
        return;
    }
    else {
        std::cout << "creating generalised jack" << std::endl;

        //jac_tot is the summ of all jackknife +1 
        //remember alle the dataj have one extra entry for the mean
        int jack_tot = 0;
        for (int e = 0; e < ens;e++)
            jack_tot += en[e].Njack;
        jack_tot = jack_tot - ens + 1;
        std::cout << "jack tot= " << jack_tot << std::endl;

        //get Nobs the minimum number of observable between the files
        int Nobs = en[0].Nobs;
        for (int i = 1; i < ens;i++) {
            if (Nobs > en[i].Nobs)
                Nobs = en[i].Nobs;
        }



        for (int e = 0;e < ens;e++) {

            double** tmp = double_malloc_2(Nobs, jack_tot);
            int counter = 0;
            for (int e1 = 0;e1 < ens;e1++) {
                for (int j = 0;j < (en[e1].Njack - 1);j++) {
                    for (int o = 0;o < Nobs;o++) {
                        if (e == e1) {
                            tmp[o][j + counter] = en[e].jack[o][j];
                        }
                        else {
                            tmp[o][j + counter] = en[e].jack[o][en[e].Njack - 1];
                        }
                    }
                }
                counter += en[e1].Njack - 1;
            }
            for (int o = 0;o < Nobs;o++)
                tmp[o][jack_tot - 1] = en[e].jack[o][en[e].Njack - 1];

            en[e].reset_error();
            free_2(en[e].Nobs, en[e].jack);
            en[e].reset_error();
            en[e].jack = tmp;
            en[e].Nobs = Nobs;
            en[e].Njack = jack_tot;
        }
        return;
    }

}

void data_all::add_space_for_n_observables(int n) {
    int Nobs = en[0].Nobs;
    int Njack = en[0].Njack;
    for (int e = 0;e < ens;e++) {

        double** tmp = double_malloc_2(Nobs + n, Njack);

        for (int j = 0; j < en[e].Njack; j++) {
            for (int o = 0; o < Nobs; o++) {
                tmp[o][j] = en[e].jack[o][j];
            }
            for (int o = Nobs; o < Nobs + n - 1; o++) {
                tmp[o][j] = 0;
            }
        }
        en[e].reset_error();
        free_2(en[e].Nobs, en[e].jack);
        en[e].jack = tmp;
        en[e].Nobs = Nobs + n;
        en[e].Njack = Njack;
    }
    return;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// print band
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// **argv=[   ???, resampling, ???, output_dir ]
void print_fit_band(char** argv, data_all gjack, struct fit_type fit_info,
    struct fit_type fit_info_m0, const char* label, const char* dir_name,
    struct fit_result fit_out, struct fit_result fit_out_m0, int var, int en, double h, std::vector<double> xval) {

    std::vector<std::vector<int>> myen = fit_info.Nxen;
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar;
    int Njack = gjack.en[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double** tif_m0 = swap_indices(fit_info_m0.Npar, Njack, fit_out_m0.P);
    double min, max;
    if (fit_info.band_range.size() != 2) {
        min = fit_info.x[var][0][Njack - 1];
        max = fit_info.x[var][0][Njack - 1];
        for (int e = 1; e < fit_info.entot; e++) {
            if (fit_info.x[var][e][Njack - 1] < min)
                min = fit_info.x[var][e][Njack - 1];
            if (fit_info.x[var][e][Njack - 1] > max)
                max = fit_info.x[var][e][Njack - 1];
        }
        if (min == max)  printf("no variation in this direction\n");
        min -= (max - min) / 10.0;
        max += (max - min) / 10.0;

    }
    else {
        min = fit_info.band_range[0];
        max = fit_info.band_range[1];
    }

    max += h / 2.0;

    for (int n = 0;n < N; n++) {

        mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_n%d_%s.txt", argv[3], label, n, dir_name);
        f = open_file(namefile, "w+");
        double* tmpx = (double*)malloc(sizeof(double) * Nvar);
        double* tmpy = (double*)malloc(sizeof(double) * Njack);
        printf("writing: %s\n", namefile);
        // band range
        double pos = min;
        while (pos < max) {
            for (int j = 0;j < Njack;j++) {
                for (int i = 0; i < xval.size(); i++) {
                    tmpx[i] = xval[i];
                }
                for (int i = xval.size(); i < Nvar; i++) {
                    tmpx[i] = fit_info.x[i][en ][j];
                }
                tmpx[var] = pos;
                tmpy[j] = fit_info.function(n, Nvar, tmpx, Npar, tif[j]);
            }
            fprintf(f, "%g  \t %g  %g\n", pos, tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));
            pos += h;
        }
        free(tmpy);free(tmpx);
        fclose(f);
    }
    free_2(Njack, tif);
    free_2(Njack, tif_m0);

}

void subMatrix(double** mat, double** temp, int p, int q, int N) {
    int i = 0, j = 0;
    // filling the sub matrix
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            // skipping if the current row or column is not equal to the current
            // element row and column
            if (row != p && col != q) {
                temp[i][j++] = mat[row][col];
                if (j == N - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
int determinantOfMatrix(double** matrix, int N) {
    int determinant = 0;
    if (N == 1) {
        return matrix[0][0];
    }
    if (N == 2) {
        return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
    }
    int sign = 1;
    double** temp = double_malloc_2(N, N);
    for (int i = 0; i < N; i++) {
        subMatrix(matrix, temp, 0, i, N);
        determinant += sign * matrix[0][i] * determinantOfMatrix(temp, N - 1);
        sign = -sign;
    }
    return determinant;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// print output
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// **argv=[   ???, resampling, ???, output_dir ]
void print_fit_output(char** argv, data_all gjack, struct fit_type fit_info,
    const char* label, struct fit_result fit_out, int* en, double*** x, double*** y) {
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar;
    int Njack = gjack.en[0].Njack;
    int N = fit_info.N;

    char namefile[NAMESIZE];
    FILE* f;

    mysprintf(namefile, NAMESIZE, "%s/%s_fit_data.txt", argv[3], label);
    f = open_file(namefile, "w+");
    int count = 0;
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < en[n];e++) {
            for (int v = 0; v < Nvar;v++) {
                fprintf(f, " %.12g\t ", x[Njack - 1][e + count][v]);
            }
            fprintf(f, " %.12g   %.12g  \t ", y[Njack - 1][e + count][0], y[Njack - 1][e + count][1]);
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
    double** cov;
    if (!fit_info.mean_only)
        cov = covariance(argv[1], Npar, Njack, fit_out.P);
    else if (fit_info.mean_only) {
        cov = double_malloc_2(Npar, Npar);
        for (int i = 0;i < fit_info.Npar;i++) {
            for (int k = 0;k < fit_info.Npar;k++) {
                cov[i][k] = fit_out.C[i][k][Njack - 1];
            }
        }
    }

    double det = determinantOfMatrix(cov, Npar);

    fprintf(f, "{\\tiny\\begin{gather}\n C=\\begin{pmatrix}\n");
    for (int i = 0;i < fit_info.Npar;i++) {
        for (int k = 0;k < i;k++)
            cov[i][k] /= sqrt(cov[i][i] * cov[k][k]);
        for (int k = i + 1;k < fit_info.Npar;k++)
            cov[i][k] /= sqrt(cov[i][i] * cov[k][k]);
    }
    for (int i = 0;i < fit_info.Npar;i++) {
        cov[i][i] /= sqrt(cov[i][i] * cov[i][i]);
    }



    for (int i = 0;i < fit_info.Npar;i++) {
        for (int j = 0;j < fit_info.Npar;j++) {
            if (j == 0)  fprintf(f, "%.3g", cov[i][j]);
            else       fprintf(f, "& %.3g", cov[i][j]);
        }
        if (i != fit_info.Npar) fprintf(f, "\\\\ \n");
        else fprintf(f, "\n");
    }
    fprintf(f, "\\end{pmatrix}\n\\\\det=%g\\\\ \\end{gather}}\n", det);
    fclose(f);
    ////////////////////////////////////////////////////////////////////////////////////////////
    mysprintf(namefile, NAMESIZE, "%s/%s_fit_P.dat", argv[3], label);
    f = open_file(namefile, "w+");
    fprintf(f, "npar   %d \n", fit_out.Npar);//error_jackboot(argv[1], Njack, fit_out.chi2)
    fprintf(f, "chi2dof   %g \n", fit_out.chi2[Njack - 1]);//error_jackboot(argv[1], Njack, fit_out.chi2)
    for (int i = 0;i < Npar;i++) {
        fprintf(f, "P[%d]   %-20.12g    %-20.12g  \n", i, fit_out.P[i][Njack - 1], myres->comp_error(fit_out.P[i]));
    }
    for (int i = 0;i < fit_info.Npar;i++) {
        for (int j = 0;j < fit_info.Npar;j++) {
            fprintf(f, "%-16.3g", cov[i][j]);

        }
        fprintf(f, "\n");
    }
    fclose(f);
    free_2(Npar, cov);


}

// **argv=[   ???, resampling, ???, output_dir ]
struct fit_result fit_all_data(char** argv, data_all gjack,
    double lhs_fun(int, int, int, data_all, struct fit_type),
    struct fit_type fit_info, const char* label) {
    int& Npar = fit_info.Npar;
    int& Nvar = fit_info.Nvar;
    int& Njack = gjack.en[0].Njack;
    int& N = fit_info.N;
    auto& myen = fit_info.myen;
    ////// allocation
    if (fit_info.Nxen.size() > 0) {
        fit_info.init_N_etot_form_Nxen();
        N = fit_info.N;
    }
    else {
        fit_info.init_Nxen_from_N_myen();
    }

    int* en = (int*)malloc(sizeof(int) * fit_info.N);// we need to init en and en_tot to allocate the other 
    for (int n = 0;n < fit_info.N; n++) { en[n] = fit_info.Nxen[n].size(); }
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
        error(fit_info.guess.size() > Npar, 1, "fit_all_data", "there are more guees than parameter in fit %s", label);
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
                    x[j][count][i] = fit_info.x[i][count][j];
                }
                count++;
            }
        }
    }

    ////////////////////////////////////////// y
    for (int i = 0; i < fit_info.Nvar; i++)  printf("x%-12d", i);
    printf("y            dy           n\n");
    int count = 0;
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < en[n];e++) {
            double* tmpj = (double*)malloc(sizeof(double) * Njack);
            for (int j = 0;j < Njack;j++) {
                tmpj[j] = lhs_fun(n, fit_info.Nxen[n][e], j, gjack, fit_info);
            }
            double err = error_jackboot(argv[1], Njack, tmpj);
            for (int j = 0;j < Njack;j++) {
                y[j][e + count][0] = tmpj[j];
                y[j][e + count][1] = err;
            }

            for (int i = 0; i < fit_info.Nvar; i++)  printf("%-13g", x[Njack - 1][e + count][i]);

            printf("%-13g%-13g%-4d\n", y[Njack - 1][e + count][0], y[Njack - 1][e + count][1], n);
            free(tmpj);
        }
        count += en[n];
    }
    //////  init end

    ///////////////// the fit 
    // scan the parameter of the fit with the last jack
    if (!fit_info.linear_fit) {
        if (fit_info.guess.size() == 0 || fit_info.repeat_start > 1) {
            guess = guess_for_non_linear_fit_Nf(N, en, x[Njack - 1], y[Njack - 1], Nvar, Npar, fit_info.function, guess, fit_info);
        }
    }
    if (fit_info.mean_only == false) {
        for (int j = Njack - 1;j >= 0;j--) {

            double a = timestamp();
            non_linear_fit_result single_jack_fit = non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess, fit_info);
            fit[j] = single_jack_fit.P;
            fit_out.chi2[j] = single_jack_fit.chi2 / (en_tot - Npar);

            if (fit_info.verbosity > 0) {
                printf("jack =%d  chi2/dof=%g   chi2=%g   time=%g   \nfinal set: ", j, fit_out.chi2[j], fit_out.chi2[j] * (en_tot - Npar), timestamp() - a);
                if (fit_info.verbosity > 1) {
                    for (int i = 0;i < Npar;i++)
                        printf("P[%d]=%g \t", i, fit[j][i]);
                    printf("\n");
                }
            }


        }
    }
    else if (fit_info.mean_only == true) {
        int j = Njack - 1;
        non_linear_fit_result single_jack_fit = non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess, fit_info);
        fit[j] = single_jack_fit.P;
        fit_out.chi2[j] = single_jack_fit.chi2 / (en_tot - Npar);

        // generate the other jackboot with the covariance matrix
        // it needs to be computed fully, with the second derivative term.
        printf("computing covariance\n");

        bool store_value_seconderiv = fit_info.second_deriv;
        fit_info.second_deriv = true;
        double** C = covariance_non_linear_fit_Nf(N, en, x[j], y[j], fit[j], Nvar, Npar, fit_info.function, fit_info);
        printf("HERE\n");
        fit_info.second_deriv = store_value_seconderiv;

        for (int ip = 0; ip < Npar;ip++)
            for (int ik = 0; ik < Npar;ik++)
                fit_out.C[ip][ik][j] = C[ip][ik];
        printf("err=%g\n", C[0][0]);

        double** tmp = fake_sampling_covariance(argv[1], fit[j], Njack, Npar, C, 1);
        for (int j1 = Njack - 2; j1 >= 0; j1--) {
            fit[j1] = (double*)malloc(sizeof(double) * Npar);
            for (int ip = 0; ip < Npar;ip++) {
                fit[j1][ip] = tmp[ip][j1];
                for (int ik = 0; ik < Npar; ik++) {
                    fit_out.C[ip][ik][j1] = C[ip][ik];
                }
            }
            fit_out.chi2[j1] = single_jack_fit.chi2 / (en_tot - Npar);
        }
        free_2(Npar, tmp);
        free_2(Npar, C);
        // for (j = Njack - 2;j >= 0;j--) {
        //     fit[j] = (double*)malloc(sizeof(double) * fit_info.Npar);

        //     for (int i = 0; i < fit_info.Npar;i++) {
        //         double noise = fit[Njack - 1][i] * (mt_rand()) / ((double)mt_rand.max() * 10000);
        //         fit[j][i] = fit[Njack - 1][i] + noise;
        //         printf("%g \t", fit[j][i]);
        //         fit_out.chi2[j] = fit_out.chi2[Njack - 1] * noise;
        //     }
        //     printf("\n");
        // }

    }
    if (fit_info.error_chi2 == true) {
        int Nboot = 10;
        // first we need to compute the boot of boot, from the jack
        double**** bby = create_boot_of_boot_from_jack(Njack, Nboot, en_tot, y);
        double* boot_chi2 = (double*)malloc(sizeof(double) * Nboot);

        for (int j = 0;j < Nboot - 1;j++) {
            // ari-calcola la correlazione
            if (fit_info.covariancey) {
                double** yy = double_malloc_2(en_tot, Nboot);
                for (int i = 0;i < en_tot;i++)
                    for (int ib = 0;ib < Nboot;ib++) {
                        yy[i][ib] = bby[j][ib][i][0];
                    }
                double** cov = covariance("boot", en_tot, Nboot, yy);

                for (int ie = 0;ie < en_tot;ie++)
                    for (int ie1 = 0;ie1 < en_tot;ie1++)
                        if (fit_info.cov[ie][ie1] != 0) {
                            // printf("%d  %d   %g   %g\n", ie, ie1, fit_info.cov[ie][ie1], cov[ie][ie1] * (Njack - 2) * (Njack - 2));
                            fit_info.cov[ie][ie1] = cov[ie][ie1] * (Njack - 2) * (Njack - 2);

                        }

                fit_info.compute_cov1_fit();

            }
            non_linear_fit_result single_jack_fit = non_linear_fit_Nf(N, en, x[j], bby[j][Nboot - 1], Nvar, Npar, fit_info.function, guess, fit_info);
            // fit[j] = single_jack_fit.P;  // we do not care at the value
            boot_chi2[j] = single_jack_fit.chi2 / (en_tot - Npar);
        }
        double dchi2 = error_jackboot("boot", Nboot, boot_chi2) * (Njack - 2);
        printf("boot chi2:\n");
        for (int j = 0;j < Nboot - 1;j++) {
            printf("%g\t", boot_chi2[j]);
            if (j < Njack) printf("%g\n", fit_out.chi2[j]);
            else printf("\n");
        }
        printf("chi2=  %g   +- %g   or %g\n", fit_out.chi2[Njack - 1], dchi2, myres->comp_error(fit_out.chi2));
        free_4(Nboot, Nboot, en_tot, bby);
        free(boot_chi2);

    }

    for (int i = 0;i < Npar;i++)
        for (int j = 0;j < Njack;j++)
            fit_out.P[i][j] = fit[j][i];

    fit_out.Npar = fit_info.Npar;
    fit_out.dof = (en_tot - Npar);
    // fit_out.name=label;
    mysprintf(fit_out.name, NAMESIZE, "%s", label);
    /////////////////////////////////////////////////////////////////////writing the result
    print_fit_output(argv, gjack, fit_info, label, fit_out, en, x, y);

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
