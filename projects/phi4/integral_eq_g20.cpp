#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <iostream>
#include "integral_eq_QC3.hpp"
#include "mutils.hpp"
#include "non_linear_fit.hpp"
#include "linear_fit.hpp"
#include "minimizer.hpp"
#include "tower.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "global.hpp"

#include "common_integral_eq.hpp"

int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./integral_eq  whatever   whatever   output_dir");

    mysprintf(argv[1], NAMESIZE, "%s", "jack");
    mysprintf(argv[2], NAMESIZE, "%s", "./");


    int Npar = 1;
    
    int NPkiso = 3;
    

    int Njack = 15;
    int seed = 1;
    int tot_parK = NPkiso + Npar;
    double* mean = (double*)malloc(sizeof(double) * tot_parK);
    mean[0] = 210.135;
    mean[1] = 9.12063;
    mean[2] = 2227.37;
    mean[3] = -0.148414;
    double** cov = double_calloc_2(tot_parK, tot_parK);
    cov[0][0] = pow(23, 2);
    cov[1][1] = pow(0.0013, 2);
    cov[2][2] = pow(6e+2, 2);
    cov[3][3] = pow(0.0016, 2);
    cov[0][1] = 0.06;   cov[0][2] = 0.629;  cov[0][3] = -0.0076;
    ;                     cov[1][2] = 0.188;   cov[1][3] = -0.0878;
    ; ;                                       cov[2][3] = -0.254;

    for (int i = 0; i < tot_parK;i++) {
        for (int j = i + 1; j < tot_parK;j++) {
            cov[i][j] = cov[i][j] * sqrt(cov[i][i] * cov[j][j]);
            cov[j][i] = cov[i][j];
        }
    }
    make_the_matrix_positive(cov, tot_parK);


    double** tmp = fake_sampling_covariance("jack", mean, Njack, tot_parK, cov, seed);


    int N = 1000;
    // int N = 1000;
    double d = 0.005;
    double eps = 0.000005;
    printf("E3    M3_re   M3_im      Kdf_re  kdf_im  Finf_re  Finf_im\n");


    double Emin = 3.021;
    double Emax = 3.02125;
    
    // double Emin = 3.015;
    // double Emax = 3.0264;
    
    // double Emin = 3.02;
    // double Emax = 3.03;

    int NE = 20;
    double dE = (Emax - Emin) / ((double)NE);

    std::vector<double> E3(NE);
    // std::vector<std::complex<double>> M3(NE);
    printf("resampling parameters:\n");
    for (int i = 0; i < tot_parK;i++)
        printf("%g   %g\n", tmp[i][Njack - 1], error_jackboot("jack", Njack, tmp[i]));






    double** P = double_malloc_2(Njack, Npar);
    double** PKiso = double_malloc_2(Njack, NPkiso);
    for (int j = 0;j < Njack;j++) {
        P[j][0] = tmp[3][j];
        PKiso[j][0] = tmp[0][j];
        PKiso[j][1] = tmp[1][j];
        PKiso[j][2] = tmp[2][j];
    }
    double*** M3, *** F;

    // compute_M3(NE, Emin, dE, Njack, E3, M3, F, N, Npar, P, compute_kcot, PKiso, compute_kiso, eps);
    // write_M3(NE, Njack, E3, M3, F, "data_M3_kcot_1par_kiso_2par.txt", argv[3]);
    read_M3(NE, Njack, E3, M3, F, "data_M3_kcot_1par_kiso_3par.txt", argv[3]);

    data_all jackall = setup_data_for_fits(NE, Njack, M3, F);
    myres=new resampling_jack(Njack-1);

    {
        fit_type fit_info;
        fit_info.Njack = Njack;
        fit_info.N = 2;
        fit_info.myen = std::vector<int>(NE);
        for (int e = 0; e < fit_info.myen.size(); e++) fit_info.myen[e] = e;
        fit_info.Nvar = 2;
        fit_info.Npar = 4;
        fit_info.function = rhs_laurent_pole;

        fit_info.malloc_x();

        int scount = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < fit_info.myen.size();e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][scount][j] = E3[e];
                    fit_info.x[1][scount][j] = 0;
                }
                scount++;
            }
        }


        fit_info.acc = 0.01;
        fit_info.h = 1e-5;
        fit_info.devorder = -2;
        fit_info.verbosity = 0;
        fit_info.repeat_start = 10;
        fit_info.guess = { 3,   2.08327e-06,   254.085,   -10.0645 };
        fit_info.precision_sum = 0;
        char namefile[NAMESIZE];

        mysprintf(namefile, NAMESIZE, "int_eq_g20_npar%d", fit_info.Npar);
        struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_M3, fit_info, namefile);
        fit_info.band_range = { Emin , Emax };
        print_fit_band(argv, jackall, fit_info, fit_info, namefile, "E_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, fit_info.myen.size() - 1, 0.0002);
        fit_info.restore_default();

    }

    {
        fit_type fit_info;
        fit_info.Njack = Njack;
        fit_info.N = 2;
        fit_info.myen = std::vector<int>(NE - 10);
        for (int e = 0; e < fit_info.myen.size(); e++) fit_info.myen[e] = e + 5;
        fit_info.Nvar = 2;
        fit_info.Npar = 4;
        fit_info.function = rhs_laurent_pole;

        fit_info.malloc_x();

        int scount = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < fit_info.myen.size();e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][scount][j] = E3[fit_info.myen[e]];
                    fit_info.x[1][scount][j] = 0;
                }
                scount++;
            }
        }


        fit_info.acc = 0.01;
        fit_info.h = 1e-5;
        fit_info.devorder = -2;
        fit_info.verbosity = 0;
        fit_info.repeat_start = 10;
        fit_info.guess = { 9.1145,   2.08327e-06,   254.085,   -10.0645 };
        fit_info.precision_sum = 0;
        char namefile[NAMESIZE];

        mysprintf(namefile, NAMESIZE, "int_eq_g20_close_to_the_pole_npar%d", fit_info.Npar);
        struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_M3, fit_info, namefile);
        fit_info.band_range = { Emin , Emax };
        print_fit_band(argv, jackall, fit_info, fit_info, namefile, "E_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, fit_info.myen.size() - 1, 0.0002);
        fit_info.restore_default();

    }
    {
        fit_type fit_info;
        fit_info.Njack = Njack;
        fit_info.N = 2;
        fit_info.myen = std::vector<int>(NE);
        for (int e = 0; e < fit_info.myen.size(); e++) fit_info.myen[e] = e;
        fit_info.Nvar = 2;
        fit_info.Npar = 6;
        fit_info.function = rhs_laurent_pole;

        fit_info.malloc_x();

        int scount = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < fit_info.myen.size();e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][scount][j] = E3[e];
                    fit_info.x[1][scount][j] = 0;
                }
                scount++;
            }
        }


        fit_info.acc = 0.01;
        fit_info.h = 1e-5;
        fit_info.devorder = 2;
        fit_info.verbosity = 3;
        fit_info.repeat_start = 1;
        fit_info.guess = { 3.02112, -1.18582e-06, -252.892, 10.7525, -2951.33, 160.156 };
        fit_info.precision_sum = 0;
        // fit_info.noderiv=true;
        // fit_info.Prange={0.1, 1e-6,  1 ,1, 0.001,0.001};
        char namefile[NAMESIZE];

        mysprintf(namefile, NAMESIZE, "int_eq_g20_npar%d", fit_info.Npar);
        struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_M3, fit_info, namefile);
        fit_info.band_range = { Emin , Emax };
        print_fit_band(argv, jackall, fit_info, fit_info, namefile, "E_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, fit_info.myen.size() - 1, 0.0002);
        fit_info.restore_default();

    }

    {
        fit_type fit_info;
        fit_info.Njack = Njack;
        fit_info.N = 2;
        fit_info.myen = std::vector<int>(NE);
        for (int e = 0; e < fit_info.myen.size(); e++) fit_info.myen[e] = e;
        fit_info.Nvar = 2;
        fit_info.Npar = 4;
        fit_info.function = rhs_BW;

        fit_info.malloc_x();

        int scount = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < fit_info.myen.size();e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][scount][j] = E3[e];
                    fit_info.x[1][scount][j] = 0;
                }
                scount++;
            }
        }

        fit_info.acc = 0.01;
        fit_info.h = 1e-5;
        fit_info.devorder = 2;
        fit_info.verbosity = 0;
        fit_info.repeat_start = 1;
        fit_info.guess = { 3.02112, -1.18582e-06, -252.892, 10.7525, };
        fit_info.precision_sum = 0;
        // fit_info.noderiv=true;
        // fit_info.Prange={0.1, 1e-6,  1 ,1, 0.001,0.001};
        char namefile[NAMESIZE];

        mysprintf(namefile, NAMESIZE, "int_eq_g20_BW_npar%d", fit_info.Npar);
        struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_M3, fit_info, namefile);
        fit_info.band_range = { Emin , Emax };
        print_fit_band(argv, jackall, fit_info, fit_info, namefile, "E_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, fit_info.myen.size() - 1, 2e-5);
        fit_info.restore_default();

    }
    

    {//abs
        fit_type fit_info;
        fit_info.Njack = Njack;
        fit_info.N = 1;
        fit_info.myen = std::vector<int>(NE);
        for (int e = 0; e < fit_info.myen.size(); e++) fit_info.myen[e] = e;
        fit_info.Nvar = 2;
        fit_info.Npar = 4;
        fit_info.function = rhs_absBW;

        fit_info.malloc_x();

        int scount = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < fit_info.myen.size();e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][scount][j] = E3[e];
                    fit_info.x[1][scount][j] = 0;
                }
                scount++;
            }
        }

        fit_info.acc = 0.01;
        fit_info.h = 1e-5;
        fit_info.devorder = 2;
        fit_info.verbosity = 0;
        fit_info.repeat_start = 1;
        fit_info.guess = { 3.02112, -1.18582e-06, -252.892, 10.7525, };
        fit_info.precision_sum = 0;
        // fit_info.noderiv=true;
        // fit_info.Prange={0.1, 1e-6,  1 ,1, 0.001,0.001};
        char namefile[NAMESIZE];

        mysprintf(namefile, NAMESIZE, "int_eq_g20_absBW_npar%d", fit_info.Npar);
        struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_absM3, fit_info, namefile);
        fit_info.band_range = { Emin , Emax };
        print_fit_band(argv, jackall, fit_info, fit_info, namefile, "E_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, fit_info.myen.size() - 1, 2e-5);
        fit_info.restore_default();

    }


    {
        fit_type fit_info;
        fit_info.Njack = Njack;
        fit_info.N = 2;
        fit_info.myen = std::vector<int>(NE);
        for (int e = 0; e < fit_info.myen.size(); e++) fit_info.myen[e] = e;
        fit_info.Nvar = 2;
        fit_info.Npar = 4;
        fit_info.function = rhs_F;

        fit_info.malloc_x();

        int scount = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < fit_info.myen.size();e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][scount][j] = E3[e];
                    fit_info.x[1][scount][j] = 0;
                }
                scount++;
            }
        }

        fit_info.acc = 0.01;
        fit_info.h = 1e-5;
        fit_info.devorder = 2;
        fit_info.verbosity = 0;
        fit_info.repeat_start = 10;
        fit_info.guess = { 1,1 };
        fit_info.precision_sum = 0;
        char namefile[NAMESIZE];

        mysprintf(namefile, NAMESIZE, "F3_line%d", fit_info.Npar);
        struct fit_result fit_F = fit_all_data(argv, jackall, lhs_F, fit_info, namefile);
        fit_info.band_range = { Emin , Emax };
        print_fit_band(argv, jackall, fit_info, fit_info, namefile, "E_m", fit_F, fit_F, 0, fit_info.myen.size() - 1, 2e-5);
        fit_info.restore_default();

        // minimize |(1/K+F)|^2 to find the pole 
        fit_info.N = 1;
        fit_info.myen = { 1 };
        fit_info.Njack = Njack;
        fit_info.NM = true;
        fit_info.noderiv = false;
        fit_info.second_deriv = true;
        fit_info.mean_only = false;
        fit_info.Nvar = 7; // it is important that the value is correct since we need to pass all the x
        fit_info.Npar = 2; // what we are minimizing
        fit_info.function = denom_M;
        fit_info.verbosity = 0;
        fit_info.acc = 1e-12;
        fit_info.h = { 0.001, 1e-5 };
        fit_info.guess = { 3.0236, -2.77561e-08 };
        fit_info.malloc_x();

        scount = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < fit_info.myen.size();e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][scount][j] = fit_F.P[0][j];
                    fit_info.x[1][scount][j] = fit_F.P[1][j];
                    fit_info.x[2][scount][j] = fit_F.P[2][j];
                    fit_info.x[3][scount][j] = fit_F.P[3][j];
                    fit_info.x[4][scount][j] = PKiso[j][0];
                    fit_info.x[5][scount][j] = PKiso[j][1];
                    fit_info.x[6][scount][j] = PKiso[j][2];


                }
                scount++;
            }
        }
        fit_result min = minimize_functions_Nf(fit_info);
        printf("LM minimizer buildin functions\n min=%g   %g   chi2=%g\n", min.P[0][Njack - 1], min.P[1][Njack - 1], min.chi2[Njack - 1]);
        printf("min=%g (%g) +i  %g  (%g) chi2=%g\n", min.P[0][Njack - 1], error_jackboot("jack", Njack, min.P[0]), min.P[1][Njack - 1], error_jackboot("jack", Njack, min.P[1]), min.chi2[Njack - 1]);
        printf("Gamma= %g  (%g)\n", 2*min.P[1][Njack - 1], 2*error_jackboot("jack", Njack, min.P[1]));
        free_fit_result(fit_info, min);
        fit_info.restore_default();
    }




}