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

void write_M3(int NE, int Njack, std::vector<double>& E3, double*** M3, std::string namef, const char* prefix) {
    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", prefix, namef.c_str());
    FILE* f = open_file(namefile, "w+");
    fprintf(f, "%d\n", NE);
    fprintf(f, "%d\n", Njack);
    for (int i = 0;i < NE;i++) {
        for (int j = 0;j < Njack;j++) {
            printf("%.12g  %.12g  %.12g\n", E3[i], M3[i][0][j], M3[i][1][j]);
            fprintf(f, "%.12g  %.12g  %.12g\n", E3[i], M3[i][0][j], M3[i][1][j]);
        }
    }
    fclose(f);

}

void read_M3(int& NE, int& Njack, std::vector<double>& E3, double***& M3, std::string namef, const char* prefix) {
    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", prefix, namef.c_str());
    FILE* f = open_file(namefile, "r+");
    int is = 0;
    is += fscanf(f, "%d\n", &NE);
    is += fscanf(f, "%d\n", &Njack);
    printf("reading NE=%d  Njack=%d\n", NE, Njack);

    M3 = double_malloc_3(NE, 2, Njack);
    for (int i = 0;i < NE;i++) {
        for (int j = 0;j < Njack;j++) {
            double a, b;
            is += fscanf(f, "%lf %lf  %lf\n", &E3[i], &a, &b);
            M3[i][0][j] = a;
            M3[i][1][j] = b;
            // printf("%d %d  %g %g\n",i,j,M3[i][0][j], M3[i][1][j]);
        }
    }
    error(is != 3 * Njack * NE + 2, 1, "read_M3", "numer of elements read = %d not correct, expected=%d", is, 2 * Njack * NE + 3);
    fclose(f);

}

void compute_M3(int NE, double Emin, double dE, int Njack, std::vector<double>& E3, double*** &M3, int N, int Npar, double** P,
    double compute_kcot(int, double*, int, double*), double** PKiso, double compute_kiso(double, double*), double eps) {

    M3 = double_malloc_3(NE, 2, Njack);

    for (int i = 0;i < NE;i++) {
        for (int j = 0;j < Njack;j++) {
            E3[i] = Emin + i * dE;
            // E3[i] = 3.02;printf("\n\n MODIFY HERE \n\n");


            std::complex<double> m3 = compute_M3_sym(E3[i], N, Npar, P[j], compute_kcot, PKiso[j], compute_kiso, eps);
            M3[i][0][j] = m3.real(); M3[i][1][j] = m3.imag();

            // std::complex<double> Kdf = compute_kiso(E3[i], PKiso);
            // Eigen::MatrixXcd D = compute_D(E3[i], N, Npar, P, compute_kcot, eps);
            // std::complex<double> Finf = comput_Finf(E3[i], D, N, Npar, P, compute_kcot, eps);
            // printf("jack =%-4d%-18.8g%-14g%-18g%-14g%-18g||%-25g%-25g%-25g%-25g\n",
            //     j, E3[i], real(m3), imag(m3), real(Kdf), imag(Kdf), PKiso[0], PKiso[1], PKiso[2], P[3]);
            // printf("%-18.8g%-14g%-18g%-14g%-18g%-18.12g%-20.12g\n", E3[i], real(m3), imag(m3), real(Kdf), imag(Kdf), real(Finf), imag(Finf));
            // printf("%-18.8g%-14g%-18g%\n", E3[i], real(m3), imag(m3));
        }

        printf("%-18.8g%-14g%-18g%-14g%-18g\n", E3[i], M3[i][0][Njack - 1], error_jackboot("jack", Njack, M3[i][0]),
            M3[i][1][Njack - 1], error_jackboot("jack", Njack, M3[i][1]));
    }
}

double rhs_laurent_pole(int n, int Nvar, double* x, int Npar, double* P) {
    error(Npar % 2 != 0, 1, "rhs_laurent_pole:", "Npar=%d but it must be multiple of two since the parameters are complex", Npar);
    std::complex<double> z(x[0], x[1]);

    std::complex<double> p(P[0], P[1]);
    std::complex<double> am1(P[2], P[3]);

    std::complex<double> r = am1 / (z * z - p * p);

    if (Npar >= 6) {
        std::complex<double> a0(P[4], P[5]);
        r += a0;
    }
    r = 1. / r;

    if (n == 0)      return real(r);
    else if (n == 1) return imag(r);
    else { printf("%s\n", __func__);  exit(1); }
}

double lhs_M3(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;

    std::complex<double> M3(gjack.en[e].jack[0][j], gjack.en[e].jack[1][j]);
    M3 = 1.0 / M3;

    if (n == 0) {
        r = M3.real();
    }
    else if (n == 1) {
        r = M3.imag();
    }
    else {
        r = 0;  printf("lhs_M3 n=%d not implemented\n", n); exit(1);
    }
    return  r;
}

double compute_kiso(double E3_m, double* P) {
    return -P[0] / (E3_m * E3_m - P[1]) + P[2];
}

double compute_kcot(int Nvar, double* x, int Npar, double* P) {
    double r;
    r = 1. / P[0];

    return r;
}

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
    mean[0] = 96.629;
    mean[1] = 9.12853;
    mean[2] = 1773.18;
    mean[3] = -0.15631;
    double** cov = double_calloc_2(tot_parK, tot_parK);
    cov[0][0] = pow(16, 2);
    cov[1][1] = pow(0.0017, 2);
    cov[2][2] = pow(980, 2);
    cov[3][3] = pow(0.0027, 2);
    cov[0][1] = -0.463; cov[0][2] = 0.0674;     cov[0][3] = 0.0289;
     ;                  cov[1][2] = 0.105;      cov[1][3] = -0.123;
     ;   ;                                      cov[2][3] = -0.9;

    for (int i = 0; i < tot_parK;i++) {
        for (int j = i + 1; j < tot_parK;j++) {
            cov[i][j] = cov[i][j] * sqrt(cov[i][i] * cov[j][j]);
            cov[j][i] = cov[i][j];
        }
    }



    double** tmp = fake_sampling_covariance("jack", mean, Njack, tot_parK, cov, seed);


    int N = 600;
    // int N = 1000;
    double d = 0.005;
    double eps = 0.0005;
    printf("E3    M3_re   M3_im      Kdf_re  kdf_im  Finf_re  Finf_im\n");


    double Emin = 3.015;
    double Emax = 3.03;
    // double Emin = 3.02;
    // double Emax = 3.03;

    int NE = 20;
    double dE = (Emax - Emin) / ((double)NE);

    std::vector<double> E3(NE);
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
    double*** M3;

    compute_M3(NE, Emin, dE, Njack, E3, M3, N, Npar, P, compute_kcot, PKiso, compute_kiso, eps);
    write_M3(NE, Njack, E3, M3, "data_M3_kcot_1par_kiso_3par.txt", argv[3]);
    // read_M3(NE, Njack, E3, M3,  "data_M3_kcot_1par_kiso_3par.txt", argv[3]);

    data_all jackall;
    jackall.ens = NE;
    jackall.resampling = "jack";
    jackall.en = new data_single[jackall.ens];
    for (int i = 0;i < NE;i++) {
        jackall.en[i].Njack = Njack;
        jackall.en[i].Nobs = 2;
        jackall.en[i].resampling = jackall.resampling;
        jackall.en[i].jack = double_malloc_2(2, Njack);
        for (int j = 0;j < Njack;j++) {
            jackall.en[i].jack[0][j] = M3[i][0][j];
            jackall.en[i].jack[1][j] = M3[i][1][j];
        }
    }

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

        mysprintf(namefile, NAMESIZE, "int_eq_g10_npar%d", fit_info.Npar);
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

        mysprintf(namefile, NAMESIZE, "int_eq_g10_close_to_the_pole_npar%d", fit_info.Npar);
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

        mysprintf(namefile, NAMESIZE, "int_eq_g10_npar%d", fit_info.Npar);
        struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_M3, fit_info, namefile);
        fit_info.band_range = { Emin , Emax };
        print_fit_band(argv, jackall, fit_info, fit_info, namefile, "E_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, fit_info.myen.size() - 1, 0.0002);
        fit_info.restore_default();

    }



}