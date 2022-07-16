#define MINIMIZER_C

#include <vector>
#include "non_linear_fit.hpp"
#include "mutils.hpp"
#include <random>


// minimize a function respect to the parameters Npar
struct fit_result minimize_functions_Nf(struct fit_type fit_info) {
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar;
    int N = fit_info.N;
    int Njack = fit_info.Njack;
    std::vector<int> myen = { 0 };

    error(fit_info.covariancey == true, 1, "minimize_functions_Nf", " fit_info.covariancey must be false");
    error(fit_info.second_deriv == false, 1, "minimize_functions_Nf", " fit_info.second_deriv must be true");
    ////// allocation
    int* en = (int*)malloc(sizeof(int) * fit_info.N);// we need to init en and en_tot to allocate the other 
    for (int e = 0;e < fit_info.N; e++) { en[e] = myen.size(); }
    int en_tot = 0;      for (int n = 0;n < N;n++) { en_tot += en[n]; }// total data to fit

    double*** y = double_malloc_3(Njack, en_tot, 2);// 2 is mean/error
    double*** x = double_malloc_3(Njack, en_tot, Nvar);
    fit_result fit_out = malloc_fit(fit_info);
    double* guess = (double*)malloc(sizeof(double) * Npar);
    double** fit = (double**)malloc(sizeof(double*) * Njack);//result of the fit, the other dimension is allocated by the function non_linear_fit_Nf()

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
                    x[j][count][i] = fit_info.x[i][count][j];
                }
                count++;
            }
        }
    }

    ////////////////////////////////////////// y

    int count = 0;
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < en[n];e++) {
            for (int j = 0;j < Njack;j++) {
                y[j][e + count][0] = 0;
                y[j][e + count][1] = 1.;
            }
            for (int i = 0; i < fit_info.Nvar; i++)  printf("%-13g", x[Njack - 1][e + count][i]);
            printf("\n");
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
            non_linear_fit_result fitj = non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess, fit_info);
            fit[j] = fitj.P;
            fit_out.chi2[j] = fitj.chi2 ;

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

        non_linear_fit_result fitj = non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess, fit_info);
        fit[j] = fitj.P;
        fit_out.chi2[j] = fitj.chi2 ;

        // for the other jackboot add a white noise to the mean
        double** C = covariance_non_linear_fit_Nf(N, en, x[j], y[j], fit[j], Nvar, Npar, fit_info.function, fit_info);
        for (int ip = 0; ip < Npar;ip++)
            for (int ik = 0; ik < Npar;ik++)
                fit_out.C[ip][ik][j] = C[ip][ik];

        double** tmp = fake_sampling_covariance(fit_info.resampling.c_str(), fit[j], Njack, Npar, C, 1);
        for (int j1 = Njack - 2; j1 >= 0; j1--) {
            fit[j1] = (double*)malloc(sizeof(double) * Npar);
            for (int ip = 0; ip < Npar;ip++) {
                fit[j1][ip] = tmp[ip][j1];

            }
        }
        free_2(Npar, tmp);
        free_2(Npar, C);


    }

    for (int i = 0;i < Npar;i++)
        for (int j = 0;j < Njack;j++)
            fit_out.P[i][j] = fit[j][i];

    fit_out.Npar = fit_info.Npar;
    // fit_out.name=label;
    // mysprintf(fit_out.name, NAMESIZE, "%s", label);
    /////////////////////////////////////////////////////////////////////writing the result
    // print_fit_output(argv, gjack, fit_info, label, fit_out, en, x, y, myen);

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