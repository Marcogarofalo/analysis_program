#define correlators_analysis_C
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "eigensystem.hpp"

#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <vector>
// #include <ranges>
#include <algorithm>
#include <cmath>
#include <memory>

#include "correlators_analysis.hpp"
#include "global.hpp"
#include "gnuplot.hpp"
#include "mutils.hpp"
#include "non_linear_fit.hpp"
#include "read.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "tower.hpp"
#include "eigensystem.hpp"

void check_correlatro_counter(int i) {
    if (corr_counter != i) {
        printf("correlator counter wrong\n");
        printf("corr_counter=%d     expected=%d\n", i, corr_counter);
        exit(-1);
    }
}

double constant_fit(int n, int Nvar, double* x, int Npar, double* P) {

    return P[0];
}

double sinh_mass(int n, int Nvar, double* x, int Npar, double* P) {

    double t = P[0];
    double T = P[1];
    double cp = sinh(x[0] * (t + 1. + (1 - T) / 2));
    double c = sinh(x[0] * (t + (1 - T) / 2));
    // printf("inside sinh = %f   %f \n", c,cp);
    return c / cp;
}

// in [varialbe] [t] [re_im]
double shift_and_M_eff_acosh(int t, int T, double** in) {
    double c[3];
    double P[2];
    double x[1];

    P[0] = t;
    P[1] = T;
    c[0] = in[t][0] - in[t + 1][0];
    c[1] = in[t + 1][0] - in[t + 2][0];
    double mass = log(c[0] / c[1]);
    // printf(" %d %g  %g   mass=%g\n", t,c[0],c[1],mass);
    /*double r=rtbis_func_eq_input(sinh_mass ,
                                     0 , 1, x ,//double n [n functions], Nvar, varibles
                                     0, P , 0, // fit_info.Npar, Parameters, ivar
                                     c[0]/c[1], mass-0.1, mass+0.1, mass*0.0001);
    */
    c[0] = in[t - 1][0] - in[t][0];
    c[1] = in[t][0] - in[t + 1][0];
    c[2] = in[t + 1][0] - in[t + 2][0];

    return acosh((c[0] + c[2]) / (2. * c[1]));
    // rp/r = sin(E2(t+1+(1-T/2)) /sin(E2(t+(1-T/2))
}

inline double second_der_time(int t, int T, double** in) {
    return (6 * in[t][0] - 4 * in[(t + 1) % T][0] + in[(t + 2) % T][0] - 4 * in[(t - 1 + T) % T][0] + in[(t - 2 + T) % T][0]);
}

double laplacian_M_eff_T(int t, int T, double** in) {
    double mass;
    double ct[1], ctp[1], res, tmp_mass, u, d;
    int i, L0;
    ct[0] = second_der_time(t, T, in);
    ctp[0] = second_der_time((t + 1) % T, T, in);
    mass = log(ct[0] / ctp[0]);
    res = 1;
    i = t;
    while (res > 1e-12) {
        u = 1. + exp(-mass * (T - 2. * i - 2.));
        d = 1. + exp(-mass * (T - 2. * i));
        tmp_mass = log((ct[0] / ctp[0]) * (u / d));
        res = fabs(tmp_mass - mass);
        mass = tmp_mass;
    }
    return mass;
}

inline double der2sym(int t, int T, double** in) {
    return (in[(t + 1) % T][0] - 2 * in[(t)][0] + in[(t - 1 + T) % T][0]);
}
double der2corr_M_eff_T(int t, int T, double** in) {
    double mass;
    double ct[1], ctp[1], res, tmp_mass, u, d;
    int i, L0;
    ct[0] = der2sym(t, T, in);
    ctp[0] = der2sym((t + 1) % T, T, in);
    mass = log(ct[0] / ctp[0]);
    res = 1;
    i = t;
    while (res > 1e-12) {
        u = 1. + exp(-mass * (T - 2. * i - 2.));
        d = 1. + exp(-mass * (T - 2. * i));
        tmp_mass = log((ct[0] / ctp[0]) * (u / d));
        res = fabs(tmp_mass - mass);
        mass = tmp_mass;
    }
    return mass;
}
double M_eff_T_ct_ctp1(int t, int T, double ct, double ctp1) {

    double  res, tmp_mass, u, d;
    int i;

    double mass = log(ct / ctp1);

    res = 1;
    i = t;
    while (res > 1e-12) {
        u = 1. + exp(-mass * (T - 2. * i - 2.));
        d = 1. + exp(-mass * (T - 2. * i));
        tmp_mass = log((ct / ctp1) * (u / d));
        res = fabs(tmp_mass - mass);
        mass = tmp_mass;
    }

    return mass;
}

template<int tau>
double M_eff_sinh_T_ct_ctptau(int t, int T, double ct, double ctp) {
    double mass;

    double res, tmp_mass, u, d;
    double dtau = (double)tau;
    mass = log(ct / ctp) / dtau;
    res = 1;
    // printf("mass  (t = %d) = %.12g\n", t, mass);
    while (res > 1e-12) {
        u = 1. + exp(-mass * (T - 2. * t - 2. * tau));
        d = 1. + exp(-mass * (T - 2. * t));
        tmp_mass = log((ct / ctp) * (u / d)) / dtau;
        res = fabs(tmp_mass - mass);
        // printf("mass iter = %.12g\n", mass);
        mass = tmp_mass;
    }
    // printf("mass end = %.12g\n#################################################\n", mass);
    return mass;
}

double M_eff_T(int t, int T, double** in) {

    double mass = M_eff_T_ct_ctp1(t, T, in[t][0], in[t + 1][0]);
    return mass;
}

template<int tau>
double M_eff_T_tau(int t, int T, double** in) {

    double mass = M_eff_sinh_T_ct_ctptau<tau>(t, T, in[t][0], in[(t + tau) % T][0]);
    return mass;
}
template double M_eff_T_tau<1>(int, int, double**);
template<> double M_eff_T_tau<2>(int t, int T, double** in) {
    double mass;
    if (t == T / 2 - 1)
        mass = 0;
    else
        mass = M_eff_sinh_T_ct_ctptau<2>(t, T, in[t][0], in[(t + 2) % T][0]);
    return mass;
}
template double M_eff_T_tau<3>(int, int, double**);
template double M_eff_T_tau<4>(int, int, double**);
template double M_eff_T_tau<5>(int, int, double**);

double M_eff_sinh_T(int t, int T, double** in) {
    double mass;

    double ct[1], ctp[1], res, tmp_mass, u, d;
    int i, L0;

    L0 = T;
    ct[0] = in[t][0];
    ctp[0] = in[t + 1][0];

    mass = log(ct[0] / ctp[0]);

    res = 1;
    i = t;
    while (res > 1e-12) {
        u = -1. + exp(-mass * (L0 - 1 - 2. * i - 2.));
        d = -1. + exp(-mass * (L0 - 1 - 2. * i));
        tmp_mass = log((ct[0] / ctp[0]) * (u / d));
        res = fabs(tmp_mass - mass);
        mass = tmp_mass;
    }

    return mass;
}

double M_eff_sinh_T_ct_ctp(int t, int T, double ct, double ctp) {
    double mass;

    double res, tmp_mass, u, d;

    mass = log(ct / ctp);
    res = 1;
    while (res > 1e-12) {
        u = -1. + exp(-mass * (T - 1 - 2. * t - 2.));
        d = -1. + exp(-mass * (T - 1 - 2. * t));
        tmp_mass = log((ct / ctp) * (u / d));
        res = fabs(tmp_mass - mass);
        mass = tmp_mass;
    }
    return mass;
}

double M_eff_t0_sinh_T(int t, int t0, int T, double** in) {
    double mass;

    double ct[1], ctp[1], res, tmp_mass, u, d;
    int i, L0;

    L0 = T;
    ct[0] = in[t][0];
    ctp[0] = in[t + 1][0];

    mass = log(ct[0] / ctp[0]);

    res = 1;
    i = t - t0;
    while (res > 1e-12) {
        u = -1. + exp(-mass * (L0 - 1 - 2. * i - 2.));
        d = -1. + exp(-mass * (L0 - 1 - 2. * i));
        tmp_mass = log((ct[0] / ctp[0]) * (u / d));
        res = fabs(tmp_mass - mass);
        mass = tmp_mass;
    }

    return mass;
}

double shift_and_M_eff_sinh_T(int t, int T, double** in) {
    double mass;

    double ct[1], ctp[1], res, tmp_mass, u, d;
    int i, L0;

    L0 = T;
    ct[0] = in[t][0] - in[t + 1][0];
    ctp[0] = in[t + 1][0] - in[t + 2][0];

    mass = log(ct[0] / ctp[0]);

    res = 1;
    i = t;
    while (res > 1e-12) {
        u = -1. + exp(-mass * (L0 - 1 - 2. * i - 2.));
        d = -1. + exp(-mass * (L0 - 1 - 2. * i));
        tmp_mass = log((ct[0] / ctp[0]) * (u / d));
        res = fabs(tmp_mass - mass);
        mass = tmp_mass;
    }

    return mass;
}


double shift_and_M_eff_log(int t, int T, double** in) {
    double mass;

    double ct[1], ctp[1], res, tmp_mass, u, d;
    int i, L0;

    L0 = T;
    ct[0] = in[t][0] - in[t + 1][0];
    ctp[0] = in[t + 1][0] - in[t + 2][0];

    mass = log(ct[0] / ctp[0]);


    return mass;
}


double M_eff_log(int t, int T, double** in) {
    double mass;

    double ct[1], ctp[1], res, tmp_mass, u, d;
    int i, L0;

    ct[0] = in[t][0];
    ctp[0] = in[t + 1][0];

    mass = log(ct[0] / ctp[0]);

    return mass;
}

double identity(int t, int T, double** in) {
    return in[t][0];
}
double identity_im(int t, int T, double** in) {
    return in[t][1];
}


double log_corr(int t, int T, double** in) {

    return -log(in[t][0]);
}

double shift_corr(int t, int T, double** in) {
    return in[t][0] - in[t + 1][0];
}

double M_eff_log_shift(int t, int T, double** in) {
    double mass;

    double ct[1], ctp[1], res, tmp_mass, u, d;
    int i, L0;

    ct[0] = in[t][0] - in[t + 1][0];
    ctp[0] = in[t + 1][0] - in[t + 2][0];

    mass = log(ct[0] / ctp[0]);

    return mass;
}


double lhs_function_f_PS(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    double corr = in[j][id][t][0];
    double M = fit_info.ext_P[0][j];
    double mu1 = fit_info.ext_P[1][j];
    double mu2 = fit_info.ext_P[2][j];

    double me = sqrt(corr * 2 * M / (exp(-t * M) + exp(-(fit_info.T - t) * M)));
    return (mu1 + mu2) * me / (M * sinh(M));
}


double lhs_function_f_PS_GEVP(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    double amp = 0;
    double M = fit_info.ext_P[0][j];
    double mu1 = fit_info.ext_P[1][j];
    double mu2 = fit_info.ext_P[2][j];


    double me = in[j][fit_info.corr_id[0]][t][0] * sqrt(2 * M) / sqrt((exp(-t * M) + exp(-(fit_info.T - t) * M)));
    return (mu1 + mu2) * me / (M * sinh(M));
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// return the line where to read the plateau
int line_read_plateaux(char** option, const char* corr, int& tmin, int& tmax, int& sep, const char* namefile_plateaux) {

    int line = 0;
    // option[3] path
    // option[6] file
    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], namefile_plateaux);
    std::fstream newfile;

    newfile.open(namefile, std::ios::in); // open a file to perform read operation using file object
    int match = 0;
    if (newfile.is_open()) { // checking whether the file is open
        std::string tp;
        while (getline(newfile, tp)) { // read data from file object and put it into string.
            line++;
            std::vector<std::string> x = split(tp, ' ');

            std::string name = option[6];
            std::string correlator = corr;
            if (x.empty() == 0) {
                if (x[0].compare(name) == 0 && x[1].compare(correlator) == 0 && x.size() >= 5) {
                    tmin = stoi(x[2]);
                    tmax = stoi(x[3]);
                    sep = stoi(x[4]);
                    printf("correlator %s  plateaux %d  %d %d\n", correlator.c_str(), tmin, tmax, sep);
                    match++;
                    // break;
                }
            }
        }
        newfile.close(); // close the file object.
    }
    else {
        error(0 == 0, 1, "correlators_analysis.cpp line_read_plateaux",
            "unable to open %s", namefile);
    }
    // error(match==0,1,"correlators_analysis.cpp line_read_plateaux",
    //       "no match for plateau %s   %s \n in the file %s ",option[6], corr,namefile);
    if (match == 0) {
        printf("no plateau found for %s in plateau file %s\n", corr, namefile);
        printf("looking for a line:\n %s  %s\n", option[6], corr);
        fprintf(stderr, "Please enter the plateau interval and separation:\n");
        myscanf(3, (char*)"%d %d %d", &tmin, &tmax, &sep);
        while (tmin > tmax && tmax >= file_head.l0 / 2 && tmin < 0) {
            fprintf(stderr, "please enter a valid number for tmin and tmax\n");
            myscanf(3, (char*)"%d %d %d", &tmin, &tmax, &sep);
        }
    }
    if (match > 1) {
        printf("multiple lines line:\n %s  %s\n", option[6], corr);
        exit(1);
    }
    if (tmin > tmax) {
        printf("\n\nerror tmin can not be more than tmax set tmin=tmax-1\n\n");
        tmin = tmax - 1;
    }
    if (tmax > file_head.l0) {
        printf("\n\n error tmax over T, set  tmax=T/2-1\n\n");
        tmin = file_head.l0 / 2 - 1;
    }

    return line;
}

struct fit_result try_fit(char** option, int tmin, int tmax, int sep, double** corr_ave, double** corr_J, int Njack, double** chi2, struct fit_type fit_info) {

    double*** y, *** x, * r, ** tmp;
    double* chi2j;
    int Npar = fit_info.Npar; // parameters
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int N = fit_info.N;                                // functions to fit
    int* en = (int*)malloc(sizeof(int) * fit_info.N); // data to fit for each function
    for (int i = 0; i < fit_info.N; i++)
        en[i] = (tmax - tmin) / sep + 1;
    int en_tot = 0; // total data to fit
    for (int n = 0; n < N; n++) {
        en_tot += en[n];
    }

    double* guess = (double*)malloc(sizeof(double) * fit_info.Npar); // initial guess for the parameter
    for (int i = 0; i < fit_info.Npar; i++)
        guess[i] = 1.;

    if (fit_info.guess.size() > 0) {
        error(fit_info.guess.size() > Npar, 1, "try_fit", " more guess than parametes");
        for (int i = 0; i < fit_info.guess.size(); i++)
            guess[i] = fit_info.guess[i];
        for (int i = fit_info.guess.size(); i < Npar; i++)
            guess[i] = 1;
    }

    chi2j = (double*)malloc(sizeof(double) * Njack);

    y = double_malloc_3(Njack, (tmax - tmin + 1) * fit_info.N, 2);
    x = double_malloc_3(Njack, en_tot, Nvar);

    for (int j = 0; j < Njack; j++) {
        int count = 0;
        for (int n = 0; n < N; n++) {
            for (int e = 0; e < en[n]; e++) {
                x[j][count][0] = e + tmin;
                for (int i = 0; i < fit_info.n_ext_P; i++) {
                    x[j][count][i + 1] = fit_info.ext_P[i][j];
                }

                count++;
            }
        }
    }

    double** fit = (double**)malloc(sizeof(double*) * Njack); // result of the fit, the other dimension is allocated by the function non_linear_fit_Nf()

    int count = 0;
    for (int n = 0; n < N; n++) {
        for (int i = tmin; i <= tmax; i += sep) {
            for (int j = 0; j < Njack; j++) {
                y[j][count][0] = corr_J[i + (file_head.l0 / 2) * n][j];
                // y[j][count][1] = error_jackboot(option[4], Njack, corr_J[i]);
                y[j][count][1] = corr_ave[i + (file_head.l0 / 2) * n][1];

            }
            count++;
        }
    }
    //     if(fit_info.guess.size()==0)
    guess = guess_for_non_linear_fit_Nf(N, en, x[Njack - 1], y[Njack - 1], Nvar, Npar, fit_info.function, guess, fit_info);
    // for (j=0;j<Njack;j++){
    for (int j = Njack - 1; j >= 0; j--) {
        // tmp=linear_fit( (tmax-tmin)/sep +1, x, y[j],  1,constant_fit_to_try );

        // if (j == Njack - 1) {
        //     printf("en_tot =%d\n", en_tot);
        //     count=0;
        //     for (int i = 0; i < N; i++) {
        //         for (int e= 0; e < en[i]; e++) {
        //             for (int v = 0; v < Nvar;v++) {
        //                 printf("%g  ", x[j][count][v]);
        //             }

        //             printf("\t %g   %g  \n", y[j][count][0], y[j][count][1]);
        //             count++;
        //         }
        //     }
        // }
        non_linear_fit_result output_fit = non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess, fit_info);
        fit[j] = output_fit.P;
        chi2j[j] = output_fit.chi2 / (en_tot - Npar);
        int max = 0;

        while (fabs(chi2j[j] - chi2j[Njack - 1]) / chi2j[Njack - 1] > fit_info.chi2_gap_jackboot && max < fit_info.guess_per_jack && !fit_info.linear_fit) {
            std::mt19937 mt_rand(max);
            printf("jack %d has a chi2= %g   while the mean has chi2=%g \n retry\n", j, chi2j[j], chi2j[Njack - 1]);

            double* guess1 = (double*)malloc(sizeof(double) * Npar);
            for (int i = 0; i < Npar; i++)
                guess1[i] = guess[i] + guess[i] * mt_rand() / ((double)10 * mt_rand.max());
            guess1 = guess_for_non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess1, fit_info);


            non_linear_fit_result tmp1 = non_linear_fit_Nf(N, en, x[j], y[j], Nvar, Npar, fit_info.function, guess1, fit_info);
            double* tmp_fit = tmp1.P;
            double tmp_chi2 = tmp1.chi2;

            if (tmp_chi2 < chi2j[j]) {
                chi2j[j] = tmp_chi2;
                for (int i = 0; i < Npar; i++)
                    fit[j][i] = tmp_fit[i];
            }

            printf("%d  chi= %g  P=\t", j, chi2j[j]);
            for (int i = 0; i < fit_info.Npar; i++) {
                printf("%g\t", fit[j][i]);
            }
            printf("\n");

            free(tmp_fit);
            free(guess1);

            max++;
        }
        // chi2j[j]=compute_chisqr((tmax-tmin)/sep +1, x, y[j],  1, tmp, constant_fit_to_try )/((tmax-tmin)/sep +1);
    }

    struct fit_result fit_out = malloc_fit(fit_info);
    for (int i = 0; i < Npar; i++) {
        for (int j = 0; j < Njack; j++) {
            fit_out.P[i][j] = fit[j][i];
        }

        /*fit_out.C[i]=(double**) malloc(sizeof(double*)*Npar);
        for(n=0;n<Npar;n++){
            fit_out.C[i][n]=(double*) malloc(sizeof(double)*Njack);
            for (j=0;j<Njack;j++){
                 fit_out.C[i][n][j]=(*C)[j][i][n];
            }
        }*/
    }
    for (int j = 0; j < Njack; j++) {
        fit_out.chi2[j] = chi2j[j];
    }

    free_2(Njack, fit);

    free_3(Njack, en_tot, x);
    free_3(Njack, (tmax - tmin + 1) * fit_info.N, y);
    free(chi2j);

    free(en);
    free(guess);
    fit_out.Npar = fit_info.Npar;
    return fit_out;
}

struct fit_result fit_fun_to_corr(char** option, struct kinematic kinematic_2pt, char* name, const char* description, double** mt, double** r, int Njack, const char* plateaux_masses, FILE* outfile, struct fit_type fit_info) {
    int yn;
    int tmin = 1, tmax = 1, sep = 1;
    double* m, * fit;
    struct fit_result fit_out;
    char jack_name[NAMESIZE];
    double* chi2;

    if (fit_info.codeplateaux) {
        tmin = fit_info.tmin;
        tmax = fit_info.tmax;
        if (tmin < 1) { printf("error: tmin=%d too small\n", tmin); exit(1); }
        if (tmax >= file_head.l0 / 2) { printf("error: tmax=%d too small\n", tmax); exit(1); }
        if (tmax < 0) { printf("error: tmax=%d not initilized\n", tmax); exit(1); }
    }
    else {
        yn = 1;
        if (strcmp(option[1], "read_plateaux") == 0) {
            int l = line_read_plateaux(option, description, tmin, tmax, sep, plateaux_masses);

        }
        else if (fit_info.plateaux_scan) {
            fprintf(fit_info.f_plateaux_scan, "\n\n#tmin tmax chi2 P1 P1err ...\n");
            for (tmin = 1; tmin < file_head.l0 / 2 - fit_info.Npar; tmin++) {
                for (tmax = tmin + fit_info.Npar; tmax < file_head.l0 / 2; tmax++) {
                    fit_result tmp = try_fit(option, tmin, tmax, sep, mt, r, Njack, &chi2, fit_info);
                    fprintf(fit_info.f_plateaux_scan, "%d   %d ", tmin, tmax);
                    fprintf(fit_info.f_plateaux_scan, " %.5g \t", tmp.chi2[Njack - 1]);
                    for (int i = 0; i < fit_info.Npar; i++) {
                        fprintf(fit_info.f_plateaux_scan, "%.15g    %.15g    \t", tmp.P[i][Njack - 1],
                            // error_jackboot(option[4], Njack, tmp.P[i])
                            myres->comp_error(tmp.P[i])
                        );
                    }
                    free_fit_result(fit_info, tmp);
                    fprintf(fit_info.f_plateaux_scan, "\n");
                }
            }
        }
        else if (strcmp(option[1], "see") == 0) {
            while (yn > 0) {
                printf("#%s 1) mu %g r %d theta %g 2) mu %g r %d theta %g \n",
                    description,
                    kinematic_2pt.k2, kinematic_2pt.r2, kinematic_2pt.mom2,
                    kinematic_2pt.k1, kinematic_2pt.r1, kinematic_2pt.mom1);
                plotting(file_head.l0, mt, &tmin, &tmax, &sep);
                fit_out = try_fit(option, tmin, tmax, sep, mt, r, Njack, &chi2, fit_info);
                double* m = mean_and_error(option[4], Njack, fit_out.P[0]);
                yn = plotting_fit(file_head.l0, mt, tmin, tmax, m, fit_out.chi2);
                // to_do:    add test of the fit
                // yn=0;
                free(m);
                free_fit_result(fit_info, fit_out);
            }
        }

    }

    fit_out = try_fit(option, tmin, tmax, sep, mt, r, Njack, &chi2, fit_info);
    // fit=give_jack_linear_fit( tmin,  tmax,sep , mt, r, Njack );
    // printing the correlator and the function
    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double* tmp = (double*)malloc(sizeof(double) * Njack);
    double* x = (double*)malloc(sizeof(double) * (fit_info.Nvar + fit_info.n_ext_P));
    int xs = fit_info.Nvar;
    for (int t = 1; t < file_head.l0 / 2; t++) {
        fprintf(outfile, "%d   %.15e    %.15e\t", t, mt[t][0], mt[t][1]);
        // variables and external parameters
        x[0] = t;
        for (int i = 0; i < fit_info.n_ext_P; i++)
            x[i + xs] = fit_info.ext_P[i][Njack - 1];

        for (int j = 0; j < Njack; j++)
            tmp[j] = fit_info.function(0, fit_info.Nvar + fit_info.n_ext_P, x, fit_info.Npar, tif[j]);
        double* f_val = mean_and_error(option[4], Njack, tmp);

        fprintf(outfile, "   %.15e    %.15e\t", f_val[0], f_val[1]);
        free(f_val);

        x[0] = t + 0.33;
        for (int j = 0; j < Njack; j++)
            tmp[j] = fit_info.function(0, fit_info.Nvar + fit_info.n_ext_P, x, fit_info.Npar, tif[j]);
        f_val = mean_and_error(option[4], Njack, tmp);
        fprintf(outfile, "%.15e   %.15e    %.15e\t", x[0], f_val[0], f_val[1]);
        free(f_val);
        x[0] = t + 0.66;
        for (int j = 0; j < Njack; j++)
            tmp[j] = fit_info.function(0, fit_info.Nvar + fit_info.n_ext_P, x, fit_info.Npar, tif[j]);
        f_val = mean_and_error(option[4], Njack, tmp);
        fprintf(outfile, "%.15e   %.15e    %.15e\n", x[0], f_val[0], f_val[1]);
        free(f_val);
    }
    free(tmp);
    free(x);
    free_2(Njack, tif);

    m = mean_and_error(option[4], Njack, fit_out.chi2);
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, tmin, tmax, m[0], m[1]);
    if (fit_info.plateaux_scan)
        fprintf(fit_info.f_plateaux_scan, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, tmin, tmax, m[0], m[1]);
    free(m);
    for (int i = 0; i < fit_info.Npar; i++) {
        m = mean_and_error(option[4], Njack, fit_out.P[i]);
        fprintf(outfile, "%.15g    %.15g    \t", m[0], m[1]);
        if (fit_info.plateaux_scan)
            fprintf(fit_info.f_plateaux_scan, "%.15g    %.15g    \t", m[0], m[1]);
        if (i == 0) {
            if (fit_info.verbosity >= 0)
                printf("#%s  fit in [%d,%d]:  %.15g    %.15g\n", description, tmin, tmax, m[0], m[1]);
        }
        free(m);
    }
    fprintf(outfile, "\n");
    if (fit_info.plateaux_scan)
        fprintf(fit_info.f_plateaux_scan, "\n");

    corr_counter++;

    fit_out.Npar = fit_info.Npar;
    return fit_out;
}

double* plateau_correlator_function(char** option, struct kinematic kinematic_2pt, char* name, double**** conf_jack, int Njack, const char* plateaux_masses, FILE* outfile, int index, const char* description, double (*fun)(int, int, double**), FILE* file_jack, struct fit_type fit_info) {
    // jackknife plateau_correlator_function(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,const char  *plateaux_masses,FILE *outfile,  int index , const char *description , double (*fun)(int ,int  , double ** ),  FILE * file_jack){
    /*int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
    if ( strcmp(option[1],"read_plateaux")==0 )
     go_to_line(*plateaux_masses,line);
    */
    double** r, * m, ** mt, * fit;
    int i, j, yn;

    r = double_malloc_2(file_head.l0, Njack);
    // mt = (double**)malloc(sizeof(double*) * file_head.l0);
    mt = double_malloc_2(file_head.l0 / 2, 2);
    fprintf(outfile, " \n\n");
    fprintf(outfile, "#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n", name,
        kinematic_2pt.k2, kinematic_2pt.r2, kinematic_2pt.mom2,
        kinematic_2pt.k1, kinematic_2pt.r1, kinematic_2pt.mom1);
    for (i = 1; i < file_head.l0 / 2; i++) {
        for (j = 0; j < Njack; j++) {
            // shift
            r[i][j] = fun(i, file_head.l0, conf_jack[j][index]);
        }
        mt[i][0] = r[i][Njack - 1];
        mt[i][1] = myres->comp_error(r[i]);
        // if (strcmp(option[4], "jack") == 0)
        //     mt[i] = mean_and_error_jack(Njack, r[i]);
        // if (strcmp(option[4], "boot") == 0)
        //     mt[i] = mean_and_error_boot(Njack, r[i]);
        // fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
    }

    // struct fit_type fit_info;
    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.function = constant_fit;
    fit_info.n_ext_P = 0;
    fit_info.linear_fit = true;

    // fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Njack,*plateaux_masses,outfile);
    struct fit_result fit_out = fit_fun_to_corr(option, kinematic_2pt, name, description, mt, r, Njack, plateaux_masses, outfile, fit_info);

    fit = (double*)malloc(sizeof(double) * Njack);
    for (j = 0; j < Njack; j++)
        fit[j] = fit_out.P[0][j];

    fwrite(fit_out.P[0], sizeof(double), Njack, file_jack);

    free_fit_result(fit_info, fit_out);

    for (i = 0; i < file_head.l0 / 2; i++)
        free(mt[i]);
    free(mt);
    free_2(file_head.l0, r);

    fflush(outfile);
    /*
     if ( strcmp(option[1],"read_plateaux")==0 ){
      fclose(*plateaux_masses);
      *plateaux_masses=open_file(kinematic_2pt.plateau_m_ll,"r");

     }*/
     //    jackknife tmp(option[4], Njack, fit);
     //    free(fit);
     //    return tmp;
    return fit;
}

struct fit_result fit_function_to_corr(char** option, struct kinematic kinematic_2pt, char* name, double**** conf_jack, const char* plateaux_masses, FILE* outfile, int index, int re_im, const char* description, struct fit_type fit_info, FILE* file_jack) {

    /*if ( strcmp(option[1],"read_plateaux")==0 )
     go_to_line(*plateaux_masses,line);
    */
    int Njack = fit_info.Njack;
    double** r, * m, ** mt, * fit;

    r = (double**)malloc(sizeof(double*) * file_head.l0);
    for (int i = 0; i < file_head.l0; i++)
        r[i] = (double*)malloc(sizeof(double) * Njack);
    // mt = (double**)malloc(sizeof(double*) * file_head.l0);
    mt = double_malloc_2(file_head.l0 / 2, 2);

    fprintf(outfile, " \n\n");
    fprintf(outfile, "#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n", name,
        kinematic_2pt.k2, kinematic_2pt.r2, kinematic_2pt.mom2,
        kinematic_2pt.k1, kinematic_2pt.r1, kinematic_2pt.mom1);
    for (int i = 1; i < file_head.l0 / 2; i++) {
        for (int j = 0; j < Njack; j++) {
            r[i][j] = conf_jack[j][index][i][re_im];
        }
        mt[i][0] = r[i][Njack - 1];
        mt[i][1] = myres->comp_error(r[i]);
        // if (strcmp(option[4], "jack") == 0)
        //     mt[i] = mean_and_error_jack(Njack, r[i]);
        // if (strcmp(option[4], "boot") == 0)
        //     mt[i] = mean_and_error_boot(Njack, r[i]);
        // fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
    }

    struct fit_result fit_out = fit_fun_to_corr(option, kinematic_2pt, name, description, mt, r, Njack, plateaux_masses, outfile, fit_info);
    // fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Njack,*plateaux_masses,outfile);

    fwrite(fit_out.P[0], sizeof(double), Njack, file_jack);

    for (int i = 0; i < file_head.l0 / 2; i++)
        free(mt[i]);
    free(mt);
    for (int i = 0; i < file_head.l0; i++)
        free(r[i]);
    free(r);

    fflush(outfile);

    /*if ( strcmp(option[1],"read_plateaux")==0 ){
     fclose(*plateaux_masses);
     *plateaux_masses=open_file(kinematic_2pt.plateau_m_ll,"r");

    }*/
    fit_out.Npar = fit_info.Npar;
    return fit_out;
}

// fit the function of fit info.function() to a combination of correlator  y= fun_of_corr( jack_index, data[jack][corr][time][reim], time, fit_info)
//
struct fit_result fit_fun_to_fun_of_corr(char** option, struct kinematic kinematic_2pt, char* name, double**** conf_jack, const char* plateaux_masses, FILE* outfile, double fun_of_corr(int, double****, int, struct fit_type), const char* description, struct fit_type fit_info, FILE* file_jack) {
    /*
    int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
       if ( strcmp(option[1],"read_plateaux")==0 )
        go_to_line(*plateaux_masses,line);
      */
    int Njack = fit_info.Njack;
    double** r, * m, ** mt, * fit;
    int i, j, yn;

    // error(fit_info.N != 1, 1, "fit_fun_to_fun_of_corr", "multiple correlator fit not implemented");
    error(file_head.l0 != fit_info.T, 1, "fit_fun_to_fun_of_corr", "T not equal to l0");
    int T = file_head.l0 * fit_info.N;
    r = (double**)malloc(sizeof(double*) * T / 2);
    for (i = 0; i < T / 2; i++)
        r[i] = (double*)malloc(sizeof(double) * Njack);
    // mt = (double**)malloc(sizeof(double*) * file_head.l0);
    mt = double_malloc_2(T / 2, 2);

    fprintf(outfile, " \n\n");
    fprintf(outfile, "#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n", name,
        kinematic_2pt.k2, kinematic_2pt.r2, kinematic_2pt.mom2,
        kinematic_2pt.k1, kinematic_2pt.r1, kinematic_2pt.mom1);
    for (i = 1; i < T / 2; i++) {
        for (j = 0; j < Njack; j++) {

            // r[i][j]= conf_jack[j][index][i][re_im];
            r[i][j] = fun_of_corr(j, conf_jack, i, fit_info);
        }
        mt[i][0] = r[i][Njack - 1];
        mt[i][1] = myres->comp_error(r[i]);
        // if (strcmp(option[4], "jack") == 0)
        //     mt[i] = mean_and_error_jack(Njack, r[i]);
        // if (strcmp(option[4], "boot") == 0)
        //     mt[i] = mean_and_error_boot(Njack, r[i]);
        // fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
    }

    struct fit_result fit_out = fit_fun_to_corr(option, kinematic_2pt, name, description, mt, r, Njack, plateaux_masses, outfile, fit_info);
    // fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Njack,*plateaux_masses,outfile);

    fwrite(fit_out.P[0], sizeof(double), Njack, file_jack);

    free_2(T / 2, mt);
    for (i = 0; i < T / 2; i++)
        free(r[i]);
    free(r);

    fflush(outfile);

    /* if ( strcmp(option[1],"read_plateaux")==0 ){
      fclose(*plateaux_masses);
      *plateaux_masses=open_file(kinematic_2pt.plateau_m_ll,"r");

     }*/
    fit_out.Npar = fit_info.Npar;
    return fit_out;
}



double** r_equal_value_or_vector(double** lambdat, double** vec, fit_type fit_info, int t, int t0, double** M) {
    int n = fit_info.n;
    int N = fit_info.N * fit_info.HENKEL_size; // if eigenvector N= component+ id_eigenvector * components
    double** r = double_malloc_2(N, 2);

    if (fit_info.value_or_vector == 0) {
        for (int n = 0; n < N; n++) {
            if ((t - t0) >= 0) {
                r[n][0] = lambdat[n][0];
                r[n][1] = lambdat[n][1];
            }
            else {
                r[n][0] = lambdat[N - 1 - n][0];
                r[n][1] = lambdat[N - 1 - n][1];
            }
        }
    }
    else if (fit_info.value_or_vector == 1) {
        for (int n = 0; n < N; n++) {
            if ((t - t0) >= 0) {
                r[n][0] = vec[n][0];
                r[n][1] = vec[n][1];
            }
            else {
                int comps = sqrt(N);
                error(comps * comps != N, 1, "r_equal_value_or_vector:", "N not a square!");
                int id = n / comps;
                int comp = n % comps;
                id = comp + (comps - 1 - id) * comps;
                r[n][0] = vec[id][0];
                r[n][1] = vec[id][1];
            }
        }
    }
    else if (fit_info.value_or_vector == 2) {
        int sqN = sqrt(N);
        double** amp = calloc_2<double>(N, 2);
        double* norm = (double*)calloc(sqN, sizeof(double));
        for (int n1 = 0; n1 < sqN; n1++) {
            for (int n2 = 0; n2 < sqN; n2++) {
                for (int ni = 0; ni < sqN; ni++) {
                    amp[n1 + n2 * sqN][0] += M[n1 + ni * sqN][0] * vec[ni + n2 * sqN][0];
                }
            }
        }

        for (int n1 = 0; n1 < sqN; n1++) {
            for (int ni = 0; ni < sqN; ni++) {
                norm[n1] += vec[ni + n1 * sqN][0] * amp[ni + n1 * sqN][0];
            }
        }
        for (int n1 = 0; n1 < sqN; n1++) {
            for (int n2 = 0; n2 < sqN; n2++) {
                amp[n1 + n2 * sqN][0] = amp[n1 + n2 * sqN][0] / sqrt((norm[n2]));
            }
        }
        for (int n = 0; n < N; n++) {
            if ((t - t0) >= 0) {
                r[n][0] = amp[n][0];
                r[n][1] = amp[n][1];
            }
            else {
                int comps = sqrt(N);
                error(comps * comps != N, 1, "r_equal_value_or_vector:", "N not a square!");
                int id = n / comps;
                int comp = n % comps;
                id = comp + (comps - 1 - id) * comps;
                r[n][0] = amp[id][0];
                r[n][1] = amp[id][1];
            }
        }
        free_2(N, amp);
        free(norm);
    }

    else {
        printf("you need to specify if you want the:\n\
                - eigenvalues (fit_info.value_or_vector=0) \n\
                - eigenvectors (fit_info.value_or_vector=1)");
        exit(1);
    }
    return r;
}

bool is_used(std::vector<int>& id_list, int i, int id) {
    bool id_used = false;
    for (int i1 = 0; i1 < i; i1++) {
        if (id_list[i1] == id) {
            id_used = true;
        }
    }
    return id_used;
}

double** GEVP_matrix(int j, double**** in, int t, struct fit_type fit_info) {
    double ct, ctp;
    int N = fit_info.N;
    if (fit_info.value_or_vector >= 1) {
        N = sqrt(fit_info.N);
        error(fit_info.N != (N * N), 1, "GEVP_matrix",
            "when you want the eigenvector N must be the square of the size of the matrix: fit_info.N=%d ", fit_info.N);
    }
    int ncorr = fit_info.corr_id.size();

    error(fit_info.HENKEL_size != 1, 1, __func__, "HENKEL_size need to be 1");
    int T = file_head.l0;
    if (t > T / 2 - 1) {
        double** r = double_calloc_2(fit_info.N * fit_info.HENKEL_size, 2);
        return r;
    }
    double** M = double_calloc_2(N * N, 2); // [NxN] [reim ]
    double** Mt0 = double_calloc_2(N * N, 2);

    double** lambdat = double_malloc_2(N, 2); // [N] [reim]
    double** vec = double_malloc_2(N * N, 2);
    int t0 = fit_info.t0_GEVP % T;
    if (fit_info.GEVP_swap_t_t0) {
        t0 = t;
        t = fit_info.t0_GEVP % T;
    }
    else if (fit_info.GEVP_tpt0)
        t0 = (fit_info.t0_GEVP + t) % T;

    if (ncorr == (N * N + N) / 2) {
        int count = 0;
        for (int i = 0; i < N; i++) {
            for (int k = i; k < N; k++) {
                int corr_ik = fit_info.corr_id[count];
                int ik = i + k * N;
                int ki = k + i * N;
                for (int reim = 0;reim < 2;reim++) {
                    M[ik][reim] = in[j][corr_ik][t][reim];
                    Mt0[ik][reim] = in[j][corr_ik][t0][reim];
                }
                M[ki][0] = M[ik][0];
                Mt0[ki][0] = Mt0[ik][0];

                M[ki][1] = -M[ik][1];
                Mt0[ki][1] = -Mt0[ik][1];
                count++;
            }
        }
    }
    else if (ncorr == N * N ) {
        for (int i = 0; i < N; i++) {// col
            for (int k = 0; k < N; k++) {// raw
                int ik = i + k * N;
                int corr_ik = fit_info.corr_id[ik];
                for (int reim = 0;reim < 2;reim++) {
                    M[ik][reim] = in[j][corr_ik][t][reim];
                    Mt0[ik][reim] = in[j][corr_ik][t0][reim];
                }
            }
        }
        for (int i = 0; i < N; i++) {
            for (int k = i ; k < N; k++) {
                int ik = i + k * N;
                int ki = k + i * N;
                M[ik][0] = (M[ik][0] + M[ki][0]) / 2.0;
                Mt0[ik][0] = (Mt0[ik][0] + Mt0[ki][0]) / 2.0;
                M[ki][0] = M[ik][0];
                Mt0[ki][0] = Mt0[ik][0];

                M[ik][1] = (M[ik][1] - M[ki][1]) / 2.0;
                Mt0[ik][1] = (Mt0[ik][1] - Mt0[ki][1]) / 2.0;
                M[ki][1] = -M[ik][1];
                Mt0[ki][1] = -Mt0[ik][1];
            }
        }
    }
    else {
        printf("GEVP_matrix\n");
        printf("you need to provide (N^2+N)/2 to populate the top triangular matrix NxN:\n  N=%d    ncorr=%d\n", N, ncorr);
        exit(1);
    }

    int verbosity = fit_info.verbosity;
    if (t > fit_info.GEVP_ignore_warning_after_t || j != 0)
        verbosity = -1;

    int err = generalysed_Eigenproblem(M, Mt0, N, &lambdat, &vec, verbosity);
    if (err > 0) {
        printf("above error at time t=%d\n", t);
    }

    int n = fit_info.n;
    double** r;
    if (fit_info.sort_by_vector) {
        fit_info.sort_by_vector = false;
        int tmp_vv = fit_info.value_or_vector;
        fit_info.value_or_vector = 1;
        int tmp_N = fit_info.N;
        if (tmp_vv == 0)
            fit_info.N *= fit_info.N;// don't worry when we ask for the values we define N = sqrt(fit_info.N)
        double** vec0 = GEVP_matrix(j, in, fit_info.t_sorting_vectors, fit_info);
        fit_info.sort_by_vector = true;
        fit_info.value_or_vector = tmp_vv;
        fit_info.N = tmp_N;
        std::vector<double> crossp(N);
        std::vector<int> id_list(N, -1);
        r = double_malloc_2(fit_info.N, 2);

        for (int i = 0; i < N; i++) {
            for (int j = 0;j < N;j++) {
                crossp[j] = 0;
                for (int k = 0;k < N;k++) {
                    crossp[j] += vec0[k + i * N][0] * vec[k + j * N][0] + vec0[k + i * N][1] * vec[k + j * N][1];
                    // crossp[j] += vec0[k + i * N][0] - vec[k + j * N][0] + vec0[k + i * N][1] - vec[k + j * N][1];
                }
            }
            
            // auto x = std::ranges::max_element(crossp);
            // int id = std::ranges::distance(crossp.begin(), x);
            auto x = std::max_element(crossp.begin(), crossp.end());
            int id = std::distance(crossp.begin(), x);
            while (is_used(id_list, i, id)) {
                x = std::max_element(crossp.begin(), crossp.end(), [x](auto const& e1, auto const& e2) {
                    if (e1 < *x)
                        return e2 < *x && e1 < e2;
                    else
                        return true; // true means that the new best is e2
                    });
                id = std::distance(crossp.begin(), x);
            }
            id_list[i] = id;
        }

        for (int i = 0; i < N; i++) {
            int id = id_list[i];
            if (fit_info.value_or_vector == 0) {
                r[i][0] = lambdat[id][0];
                r[i][1] = lambdat[id][1];
            }
            else if (fit_info.value_or_vector == 1) {
                for (int n = 0; n < N; n++) {
                    r[n + i * N][0] = vec[n + id * N][0];
                    r[n + i * N][1] = vec[n + id * N][1];
                }
            }
            else if (fit_info.value_or_vector == 2) {
                printf("fit_info.value_or_vector == 2 (ampitudes) not supported");
                exit(1);
            }

        }
    }
    else {
        r = r_equal_value_or_vector(lambdat, vec, fit_info, t, t0, M);
    }
    free_2(N * N, M);
    free_2(N * N, Mt0);
    free_2(N, lambdat);
    //     free_2(N,lambdatp1);
    free_2(N * N, vec);

    return r;
}
// fun_of_corr must return a double ** [N(correlators)][2(re/im)]
void add_correlators(char** option, int& ncorr_conf_jack, double****& conf_jack, double** fun_of_corr(int, double****, int, struct fit_type), struct fit_type fit_info) {

    int correlators_out = ncorr_conf_jack + fit_info.N;
    int Njack = fit_info.Njack;

    double**** corr_out = malloc_corr(Njack, correlators_out, file_head.l0);
    // copy the first part
    for (int j = 0; j < Njack; j++) {
        for (int v = 0; v < ncorr_conf_jack; v++) {
            for (int t = 0; t < file_head.l0; t++) {
                corr_out[j][v][t][0] = conf_jack[j][v][t][0];
                corr_out[j][v][t][1] = conf_jack[j][v][t][1];
            }
        }
    }



    for (int j = 0; j < Njack; j++) {
        for (int t = 0; t < file_head.l0; t++) {
            double** r = fun_of_corr(j, conf_jack, t, fit_info);
            for (int n = 0; n < fit_info.N; n++) {
                corr_out[j][ncorr_conf_jack + n][t][0] = r[n][0];
                corr_out[j][ncorr_conf_jack + n][t][1] = r[n][1];
            }
            free_2(fit_info.N, r);
        }
    }

    free_corr(Njack, ncorr_conf_jack, file_head.l0, conf_jack);
    conf_jack = corr_out;
    ncorr_conf_jack = correlators_out;
}

// fun_of_corr must return a double ** [N(correlators)][2(re/im)]
void add_correlators_no_alloc(char** option, int& ncorr_conf_jack, double****& conf_jack, double** fun_of_corr(int, double****, int, struct fit_type), struct fit_type fit_info, int Nmax) {

    int correlators_out = ncorr_conf_jack + fit_info.N;
    error(correlators_out > Nmax, 1, "add_correlators_no_alloc", "error conf_jack size: %d impossible to save %d correlators", Nmax, correlators_out);
    int Njack = fit_info.Njack;

    for (int j = 0; j < Njack; j++) {
        for (int t = 0; t < file_head.l0; t++) {
            double** r = fun_of_corr(j, conf_jack, t, fit_info);
            for (int n = 0; n < fit_info.N; n++) {
                conf_jack[j][ncorr_conf_jack + n][t][0] = r[n][0];
                conf_jack[j][ncorr_conf_jack + n][t][1] = r[n][1];
            }
            free_2(fit_info.N, r);
        }
    }

    ncorr_conf_jack = correlators_out;
}

// r[t][jack]
void add_one_correlators(char** option, int& ncorr_conf_jack, double****& conf_jack, struct fit_type fit_info, double** r) {

    int correlators_out = ncorr_conf_jack + 1;
    int Njack = fit_info.Njack;

    double**** corr_out = calloc_corr(Njack, correlators_out, file_head.l0);
    // copy the first part
    for (int j = 0; j < Njack; j++) {
        for (int v = 0; v < ncorr_conf_jack; v++) {
            for (int t = 0; t < file_head.l0; t++) {
                corr_out[j][v][t][0] = conf_jack[j][v][t][0];
                corr_out[j][v][t][1] = conf_jack[j][v][t][1];
            }
        }
    }


    for (int j = 0; j < Njack; j++) {
        for (int t = 0; t < file_head.l0; t++) {
            for (int n = 0; n < 1; n++) {
                corr_out[j][ncorr_conf_jack + n][t][0] = r[t][j];
                corr_out[j][ncorr_conf_jack + n][t][1] = 0;
            }
        }
    }

    free_corr(Njack, ncorr_conf_jack, file_head.l0, conf_jack);
    conf_jack = corr_out;
    ncorr_conf_jack = correlators_out;
}

void zero_corr(double* zeros, int Njack, FILE* jack_file) {
    fwrite(zeros, sizeof(double), Njack, jack_file);
    corr_counter++;
}


void write_jack(double* corr, int Njack, FILE* jack_file) {
    fwrite(corr, sizeof(double), Njack, jack_file);
    corr_counter++;
}

void print_result_in_file(FILE* outfile, double* res, const char* name, double chi2, int tmin, int tmax) {

    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < 2; t++) {
        fprintf(outfile, "%d   %.15g   %.15g\t", 0, 0.0, 0.0);
        fprintf(outfile, "%.15g   %.15g\t", 0.0, 0.0);
        fprintf(outfile, "%.15g   %.15g\t", 0.0, 0.0);
        fprintf(outfile, "%.15g   %.15g\n", 0.0, 0.0);
    }
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", name, tmin, tmax, chi2, 0.0);
    fprintf(outfile, "   %.15g   %.15g\n", res[myres->Njack - 1], myres->comp_error(res));

}