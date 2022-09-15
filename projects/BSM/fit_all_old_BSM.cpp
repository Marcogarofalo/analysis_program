#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"

#include "resampling.hpp"
#include "read.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
#include <string>
#include "correlators_analysis.hpp"
#include "eigensystem.hpp"
#include "resampling_new.hpp"
#include "global.hpp"

#include <cstring> 
#include <string>
#include <fstream>
#include <memory>


//local folder
#include "header_BSM.hpp"
#include "lhs_functions.hpp"
#include "fit_data.hpp"



using namespace std;
// int NeW=0;


double* interpol_etacr(int Nmus, int Njack, double* x, double** y, double tau) {
    //  y(s1) = ( 1  x(s1) )  (P[0])
    //  y(s1) = ( 1  x(s2) )  (P[1])

    double** matrix = double_malloc_2(Nmus, Nmus);
    double* Zj = (double*)malloc(sizeof(double) * Nmus);
    double* Zint = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++) {
        for (int i = 0;i < Nmus;i++) {
            Zj[i] = y[i][j];
            for (int k = 0;k < Nmus;k++) {
                matrix[i][k] = pow(x[i], k);
            }
        }
        double* P = LU_decomposition_solver(Nmus, matrix, Zj);
        Zint[j] = 0;
        for (int k = 0;k < Nmus;k++) {
            Zint[j] += pow(tau, k) * P[k];
        }
        free(P);
    }
    free(Zj);
    free_2(Nmus, matrix);
    printf("tau    eta_cr\n");
    for (int i = 0;i < Nmus;i++) {
        printf("%g  %g   %g\n", x[i], y[i][Njack - 1], error_jackboot("jack", Njack - 1, y[i]));
    }
    printf("extrapolated to %g  :%g  %g\n", tau, Zint[Njack - 1], error_jackboot("jack", Njack - 1, Zint));
    return(Zint);
}

void print_fit_band_eta(char** argv, vector<data_BSM> gjack, struct fit_type fit_info, const char* label, struct fit_result fit_out, vector<header_BSM> params, std::vector<int> myenW, std::vector<double> range = { -1.7,-1 }) {
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double* tmpx = (double*)malloc(sizeof(double*) * Nvar);
    double* tmpy = (double*)malloc(sizeof(double*) * Njack);
    int counting_e = 0;
    for (int e : myenW) {
        for (int n = 0;n < N; n++) {

            mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_n%d_en%d_eta.txt", argv[3], label, n, counting_e);
            f = open_file(namefile, "w+");
            double* tmpx = (double*)malloc(sizeof(double*) * Nvar);
            double* tmpy = (double*)malloc(sizeof(double*) * Njack);
            printf("writing: %s\n", namefile);

            for (int i = 0; i < 100; i++) {
                for (int j = 0;j < Njack;j++) {
                    double finalL = i;

                    tmpx[0] = (double)params[e].L;
                    tmpx[1] = (double)params[e].T;
                    tmpx[2] = (double)params[e].rho;
                    tmpx[3] = range[0] + (range[1] - range[0]) * i / 100.0;
                    tmpx[4] = (double)params[e].csw;
                    tmpx[5] = (double)params[e].mu03;
                    tmpx[6] = (double)params[e].m0;
                    //m0   put for each n the mass of the last ensemble


                    for (int i = fit_info.Nvar; i < fit_info.Nvar + fit_info.n_ext_P; i++)
                        tmpx[i] = fit_info.ext_P[i - fit_info.Nvar][j];


                    tmpy[j] = fit_info.function(n, Nvar, tmpx, Npar, tif[j]);
                }
                fprintf(f, "%g  \t %g  %g\n", tmpx[3], tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));
            }
            free(tmpy);free(tmpx);
            fclose(f);

        }
        counting_e++;
    }
    free_2(Njack, tif);

}






int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");

    vector<header_BSM> paramsj;
    vector<data_BSM> dataj;

    int NeW = 0;
    int NeNG = 0;
    header_BSM params;

    char namefile[NAMESIZE];

    //mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L16_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
    //emplace_back_par_data(namefile,paramsj,dataj);

     //0
    mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.560_M02_0000_mu03_0224_rho1/jackknife/%s_T40_L16_rho1.960000_eta-1.098250_csw0.000000_mu030.022400_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.560_M02_0000_mu03_0316_rho1/jackknife/%s_T40_L16_rho1.960000_eta-1.098250_csw0.000000_mu030.031600_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.560_M02_0000_mu03_0387_rho1/jackknife/%s_T40_L16_rho1.960000_eta-1.098250_csw0.000000_mu030.038700_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    ////////

    mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.580_M02_0000_mu03_0120_rho1/jackknife/%s_T40_L16_rho1.960000_eta-1.137474_csw0.000000_mu030.012000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    // mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.580_M02_0000_mu03_0172_rho1/jackknife/%s_T40_L16_rho1.960000_eta-1.137474_csw0.000000_mu030.017200_m00.000000", argv[2], argv[1]);
    // emplace_back_par_data(namefile, paramsj, dataj);

    // mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.580_M02_0000_mu03_0600_rho1/jackknife/%s_T40_L16_rho1.960000_eta-0.580000_csw0.000000_mu030.060000_m00.000000", argv[2], argv[1]);
    // emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.580_M02_0000_mu03_0224_rho1/jackknife/%s_T40_L16_rho1.960000_eta-1.137474_csw0.000000_mu030.022400_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.580_M02_0000_mu03_0316_rho1/jackknife/%s_T40_L16_rho1.960000_eta-1.137474_csw0.000000_mu030.031600_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.580_M02_0000_mu03_0387_rho1/jackknife/%s_T40_L16_rho1.960000_eta-1.137474_csw0.000000_mu030.038700_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    //////
    mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.660_M02_0000_mu03_0224_rho1/jackknife/%s_T40_L16_rho1.960000_eta-1.294366_csw0.000000_mu030.022400_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_W/etam0.660_M02_0000_mu03_0387_rho1/jackknife/%s_T40_L16_rho1.960000_eta-1.294366_csw0.000000_mu030.038700_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    NeW = dataj.size();
    printf("number of ensembles W = %d\n", NeW);
    vector<int> myenW(NeW);
    for (int i = 0;i < NeW; i++)  myenW[i] = i;

    ////////////////////// b595

    mysprintf(namefile, NAMESIZE, "%s/b595_W/etam0.4987_M02_0000_mu03_0186_rho1.002001/jackknife/%s_T48_L20_rho1.960000_eta-0.976077_csw0.000000_mu030.018600_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b595_W/etam0.4987_M02_0000_mu03_0321_rho1.002001/jackknife/%s_T48_L20_rho1.960000_eta-0.976077_csw0.000000_mu030.032100_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b595_W/etam0.5290_M02_0000_mu03_0186_rho1.002001/jackknife/%s_T48_L20_rho1.960000_eta-1.035382_csw0.000000_mu030.018600_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b595_W/etam0.5290_M02_0000_mu03_0321_rho1.002001/jackknife/%s_T48_L20_rho1.960000_eta-1.035382_csw0.000000_mu030.032100_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b595_W/etam0.5503_M02_0000_mu03_0100_rho1.002001/jackknife/%s_T48_L20_rho1.960000_eta-1.077071_csw0.000000_mu030.010000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b595_W/etam0.5503_M02_0000_mu03_0186_rho1.002001/jackknife/%s_T48_L20_rho1.960000_eta-1.077071_csw0.000000_mu030.018600_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b595_W/etam0.5503_M02_0000_mu03_0321_rho1.002001/jackknife/%s_T48_L20_rho1.960000_eta-1.077071_csw0.000000_mu030.032100_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    //////


    int NeW595 = dataj.size() - NeW;
    printf("number of ensembles W = %d\n", NeW595);
    vector<int> myenW595(NeW595);
    for (int i = 0;i < NeW595; i++)  myenW595[i] = i + NeW;



    ////////////////////// b575

    mysprintf(namefile, NAMESIZE, "%s/b575_W/etam0.585_M02_0000_mu03_0180_rho0.997226/jackknife/%s_T32_L16_rho1.960000_eta-1.150470_csw0.000000_mu030.018000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b575_W/etam0.585_M02_0000_mu03_0280_rho0.997226/jackknife/%s_T32_L16_rho1.960000_eta-1.150470_csw0.000000_mu030.028000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b575_W/etam0.585_M02_0000_mu03_0480_rho0.997226/jackknife/%s_T32_L16_rho1.960000_eta-1.150470_csw0.000000_mu030.048000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b575_W/etam0.605_M02_0000_mu03_0180_rho0.997226/jackknife/%s_T32_L16_rho1.960000_eta-1.189802_csw0.000000_mu030.018000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b575_W/etam0.605_M02_0000_mu03_0280_rho0.997226/jackknife/%s_T32_L16_rho1.960000_eta-1.189802_csw0.000000_mu030.028000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b575_W/etam0.605_M02_0000_mu03_0480_rho0.997226/jackknife/%s_T32_L16_rho1.960000_eta-1.189802_csw0.000000_mu030.048000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b575_W/etam0.695_M02_0000_mu03_0180_rho0.997226/jackknife/%s_T32_L16_rho1.960000_eta-1.366797_csw0.000000_mu030.018000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b575_W/etam0.695_M02_0000_mu03_0280_rho0.997226/jackknife/%s_T32_L16_rho1.960000_eta-1.366797_csw0.000000_mu030.028000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b575_W/etam0.695_M02_0000_mu03_0480_rho0.997226/jackknife/%s_T32_L16_rho1.960000_eta-1.366797_csw0.000000_mu030.048000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    int NeW575 = dataj.size() - NeW - NeW595;
    printf("number of ensembles W = %d\n", NeW575);
    vector<int> myenW575(NeW575);
    for (int i = 0;i < NeW575; i++)  myenW575[i] = i + NeW + NeW595;

    ////////////////////// b585 NG

    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.600_M02_0000_mu03_0316_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.167748_csw0.000000_mu030.031600_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.6051_M02_0000_mu03_0070_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.177674_csw0.000000_mu030.007000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.6051_M02_0000_mu03_0100_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.177674_csw0.000000_mu030.010000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.6051_M02_0000_mu03_0140_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.177674_csw0.000000_mu030.014000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.6051_M02_0000_mu03_0224_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.177674_csw0.000000_mu030.022400_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.614_M02_0000_mu03_0100_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.194996_csw0.000000_mu030.010000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.618_M02_0000_mu03_0040_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.202781_csw0.000000_mu030.004000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.618_M02_0000_mu03_0070_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.202781_csw0.000000_mu030.007000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.618_M02_0000_mu03_0100_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.202781_csw0.000000_mu030.010000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.618_M02_0000_mu03_0120_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.202781_csw0.000000_mu030.012000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.618_M02_0000_mu03_0224_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.202781_csw0.000000_mu030.022400_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.6201_M02_0000_mu03_0040_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.206868_csw0.000000_mu030.004000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.6201_M02_0000_mu03_0070_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.206868_csw0.000000_mu030.007000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.6201_M02_0000_mu03_0100_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.206868_csw0.000000_mu030.010000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.620_M02_0000_mu03_0224_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.206673_csw0.000000_mu030.022400_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.622_M02_0000_mu03_0040_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.210566_csw0.000000_mu030.004000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.622_M02_0000_mu03_0070_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.210566_csw0.000000_mu030.007000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.622_M02_0000_mu03_0100_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.210566_csw0.000000_mu030.010000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b585_NG/etam0.622_M02_0000_mu03_0120_rho1.00766/jackknife/%s_T40_L20_rho1.960000_eta-1.210566_csw0.000000_mu030.012000_m00.000000", argv[2], argv[1]);


    int NeNG585 = dataj.size() - NeW - NeW595 - NeW575;
    printf("number of ensembles NG 585 = %d\n", NeNG585);
    vector<int> myenNG585(NeNG585);
    for (int i = 0;i < NeNG585; i++)  myenNG585[i] = i + NeW + NeW595 + NeW575;


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////   NG
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    int size_before = dataj.size();

    ////////////////////// b595 NG

    // mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.5718_M02_0000_mu03_0077_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.113413_csw0.000000_mu030.007700_m00.000000", argv[2], argv[1]);
    // emplace_back_par_data(namefile, paramsj, dataj);
    // mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.5718_M02_0000_mu03_0108_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.113413_csw0.000000_mu030.010800_m00.000000", argv[2], argv[1]);
    // emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.5760_M02_0000_mu03_0077_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.121591_csw0.000000_mu030.007700_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.588_M02_0000_mu03_0060_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.144957_csw0.000000_mu030.006000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.588_M02_0000_mu03_0077_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.144957_csw0.000000_mu030.007700_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    // mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.588_M02_0000_mu03_0116_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.144957_csw0.000000_mu030.011600_m00.000000", argv[2], argv[1]);
    // emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.588_M02_0000_mu03_0145_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.144957_csw0.000000_mu030.014500_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.5893_M02_0000_mu03_0060_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.147489_csw0.000000_mu030.006000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.5893_M02_0000_mu03_0077_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.147489_csw0.000000_mu030.007700_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.5893_M02_0000_mu03_0116_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.147489_csw0.000000_mu030.011600_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.5893_M02_0000_mu03_0145_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.147489_csw0.000000_mu030.014500_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b595_NG/etam0.5893_M02_0000_mu03_0185_rho1.007170/jackknife/%s_T48_L24_rho1.960000_eta-1.147489_csw0.000000_mu030.018500_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    int NeNG595 = dataj.size() - size_before;
    printf("number of ensembles NG 595 = %d\n", NeNG595);
    vector<int> myenNG595(NeNG595);
    for (int i = 0;i < NeNG595; i++)  myenNG595[i] = i + size_before;
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////   NG
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    size_before = dataj.size();

    ////////////////////// b575 NG

    mysprintf(namefile, NAMESIZE, "%s/b575_NG/etam0.6315_M02_0000_mu03_0087_rho1.008742/jackknife/%s_T40_L16_rho1.960000_eta-1.227740_csw0.000000_mu030.008700_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b575_NG/etam0.6315_M02_0000_mu03_0131_rho1.008742/jackknife/%s_T40_L16_rho1.960000_eta-1.227740_csw0.000000_mu030.013100_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b575_NG/etam0.6315_M02_0000_mu03_0183_rho1.008742/jackknife/%s_T40_L16_rho1.960000_eta-1.227740_csw0.000000_mu030.018300_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b575_NG/etam0.645_M02_0000_mu03_0183_rho1.008742/jackknife/%s_T40_L16_rho1.960000_eta-1.253986_csw0.000000_mu030.018300_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b575_NG/etam0.651_M02_0000_mu03_0131_rho1.008742/jackknife/%s_T40_L16_rho1.960000_eta-1.265651_csw0.000000_mu030.013100_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b575_NG/etam0.654_M02_0000_mu03_0050_rho1.008742/jackknife/%s_T40_L16_rho1.960000_eta-1.271484_csw0.000000_mu030.005000_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b575_NG/etam0.654_M02_0000_mu03_0087_rho1.008742/jackknife/%s_T40_L16_rho1.960000_eta-1.271484_csw0.000000_mu030.008700_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b575_NG/etam0.654_M02_0000_mu03_0131_rho1.008742/jackknife/%s_T40_L16_rho1.960000_eta-1.271484_csw0.000000_mu030.013100_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b575_NG/etam0.654_M02_0000_mu03_0183_rho1.008742/jackknife/%s_T40_L16_rho1.960000_eta-1.271484_csw0.000000_mu030.018300_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    mysprintf(namefile, NAMESIZE, "%s/b575_NG/etam0.654_M02_0000_mu03_0277_rho1.008742/jackknife/%s_T40_L16_rho1.960000_eta-1.271484_csw0.000000_mu030.027700_m00.000000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    int NeNG575 = dataj.size() - size_before;
    printf("number of ensembles NG 575 = %d\n", NeNG575);
    vector<int> myenNG575(NeNG575);
    for (int i = 0;i < NeNG575; i++)  myenNG575[i] = i + size_before;

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////   NG
    ///////////////////////////////////////////////////////////////////////////////////////////////////






    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////   rho 3 W
    ///////////////////////////////////////////////////////////////////////////////////////////////////



    int NeW3 = dataj.size() - (NeW + NeNG);
    printf("number of ensembles W(rho3) = %d\n", NeW3);
    vector<int> myenW3(NeW3);
    for (int i = (NeW + NeNG);i < (NeW + NeNG) + NeW3; i++)  myenW3[i - (NeW + NeNG)] = i;

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////   rho 3 NG
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    int NeNG3 = dataj.size() - (NeW + NeNG + NeW3);
    printf("number of ensembles NG(rho3) = %d\n", NeNG3);
    vector<int> myenNG3(NeNG3);
    for (int i = (NeW + NeNG + NeW3);i < (NeW + NeNG + NeW3) + NeNG3; i++) {
        myenNG3[i - (NeW + NeNG + NeW3)] = i;
        printf("en[%d]=%d\n", i - (NeW + NeNG + NeW3), myenNG3[i - (NeW + NeNG + NeW3)]);
    }


    printf("%g   %g\n", dataj[0].jack[3][dataj[0].Njack - 3], error_jackboot(argv[1], dataj[0].Njack, dataj[0].jack[3]));
    /* for (int j=0; j<dataj[0].Njack; j++)
         printf("%d   %g\n",j,dataj[0].jack[3][j]);
    */
    vector<data_BSM> gjack = create_generalised_resampling(dataj);


    int Njack = gjack[0].Njack;

    if (strcmp(argv[1], "jack") == 0) {
        myres = new resampling_jack(Njack - 1);
    }
    else if (strcmp(argv[1], "boot") == 0) {
        myres = new resampling_boot(Njack - 1);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////////////////////////////////////////////////////////////////////////
    struct fit_type fit_info;
    struct fit_type fit_info_NG;
    struct fit_result  fit_critical;
    struct fit_result  fit_critical_NG;
    char namefit[NAMESIZE];
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;


    fit_critical = fit_data(argv, paramsj, gjack, lhs_critical_eta_mu_m0, fit_info, "eta_critical_tau5_b585", myenW);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau5_b585", fit_critical, paramsj, myenW);
    printf("eta_cr 585 tau5=%g  +-  %g\n", fit_critical.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[0]));

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();

    printf("/////////////////////////////////////  r_AWI_tau2///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 13,4 };


    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau2_b585", myenW);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau2_b585", fit_critical, paramsj, myenW);


    printf("/////////////////////////////////////  r_AWI_tau3///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 14,4 };


    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau3_b585", myenW);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau3_b585", fit_critical, paramsj, myenW);
    double** eta_cr_585_taus = (double**)malloc(sizeof(double*) * 2);
    eta_cr_585_taus[0] = (double*)malloc(sizeof(double) * Njack);
    eta_cr_585_taus[1] = (double*)malloc(sizeof(double) * Njack);
    // fit eta_cr to 3.75
    double* b585_taus = (double*)malloc(sizeof(double) * 2);
    b585_taus[0] = 3;
    for (int j = 0; j < Njack;j++) eta_cr_585_taus[0][j] = fit_critical.P[0][j];

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();
    printf("/////////////////////////////////////  r_AWI_tau4///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 15,4 };


    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau4_b585", myenW);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau4_b585", fit_critical, paramsj, myenW);

    b585_taus[1] = 4;
    for (int j = 0; j < Njack;j++) eta_cr_585_taus[1][j] = fit_critical.P[0][j];

    double* eta_cr_b585 = interpol_etacr(2, Njack, b585_taus, eta_cr_585_taus, 3.75);
    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();
    printf("//////////////////  fit GPS///////////////////////// ");
    fit_info_NG.Nvar = 7;
    fit_info_NG.Npar = 7;
    fit_info_NG.N = 1;
    fit_info_NG.Njack = gjack[0].Njack;
    fit_info_NG.n_ext_P = 1;
    fit_info_NG.ext_P = (double**)malloc(sizeof(double*) * fit_info_NG.n_ext_P);
    fit_info_NG.ext_P[0] = eta_cr_b585;
    fit_info_NG.function = rhs_GPS;


    mysprintf(namefit, NAMESIZE, "fit_GPS_b585_etacr");
    fit_result GPS = fit_data(argv, paramsj, gjack, lhs_GPS, fit_info_NG, namefit, myenW);
    print_fit_band_eta(argv, gjack, fit_info_NG, namefit, GPS, paramsj, myenW);

    fit_info_NG.restore_default();
    printf("//////////////////  NG_585///////////////////////// ");
    fit_info_NG.Nvar = 7;
    fit_info_NG.Npar = 8;
    fit_info_NG.N = 2;
    fit_info_NG.Njack = gjack[0].Njack;
    fit_info_NG.n_ext_P = 2;
    fit_info_NG.ext_P = (double**)malloc(sizeof(double*) * fit_info_NG.n_ext_P);
    fit_info_NG.ext_P[0] = eta_cr_b585;
    fit_info_NG.ext_P[1] = GPS.P[0];
    fit_info_NG.function = rhs_NG_mpcac_MPS2_mu_eta;


    mysprintf(namefit, NAMESIZE, "fit_NG_mpcac_MPS_b585");
    fit_result fit_critical_NG_585 = fit_data(argv, paramsj, gjack, lhs_mpcac_MPS2, fit_info_NG, namefit, myenNG585);
    print_fit_band_eta(argv, gjack, fit_info_NG, namefit, fit_critical_NG_585, paramsj, myenNG585);

    printf("m_pcac(eta_cr)=%g  +-  %g\n", fit_critical_NG_585.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical_NG_585.P[0]));
    printf("M_PS(eta_cr) =%g  +-  %g\n", fit_critical_NG_585.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_critical_NG_585.P[1]));


    // free_fit_result(fit_info_NG, fit_critical_NG);
    fit_info_NG.restore_default();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("///////////////////////////////////////////    b595   ///////////////////////////////////////////////////////");
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    printf("/////////////////////////////////////  r_AWI_tau2///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 13,4 };


    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau2_b595", myenW595);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau2_b595", fit_critical, paramsj, myenW595);

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();
    printf("/////////////////////////////////////  r_AWI_tau3///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 14,4 };


    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau3_b595", myenW595);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau3_b595", fit_critical, paramsj, myenW595);

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();
    printf("/////////////////////////////////////  r_AWI_tau4///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 15,4 };


    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau4_b595", myenW595);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau4_b595", fit_critical, paramsj, myenW595);
    double** eta_cr_595_taus = (double**)malloc(sizeof(double*) * 2);
    eta_cr_595_taus[0] = (double*)malloc(sizeof(double) * Njack);
    eta_cr_595_taus[1] = (double*)malloc(sizeof(double) * Njack);
    double* b595_taus = (double*)malloc(sizeof(double) * 2);
    b595_taus[0] = 4;
    b595_taus[1] = 5;
    for (int j = 0; j < Njack;j++) eta_cr_595_taus[0][j] = fit_critical.P[0][j];
    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();

    printf("/////////////////////////////////////  r_AWI_tau5///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 3,4 };



    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau5_b595", myenW595);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau5_b595", fit_critical, paramsj, myenW595);

    for (int j = 0; j < Njack;j++) eta_cr_595_taus[1][j] = fit_critical.P[0][j];
    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();

    printf("/////////////////////////////////////  r_AWI_tau6///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 16,4 };

    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau6_b595", myenW595);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau6_b595", fit_critical, paramsj, myenW595);


    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();
    double* eta_cr_b595 = interpol_etacr(2, Njack, b595_taus, eta_cr_595_taus, 4.5);
    printf("//////////////////  fit GPS///////////////////////// ");
    fit_info_NG.Nvar = 7;
    fit_info_NG.Npar = 7;
    fit_info_NG.N = 1;
    fit_info_NG.Njack = gjack[0].Njack;
    fit_info_NG.n_ext_P = 1;
    fit_info_NG.ext_P = (double**)malloc(sizeof(double*) * fit_info_NG.n_ext_P);
    fit_info_NG.ext_P[0] = eta_cr_b595;
    fit_info_NG.function = rhs_GPS;


    mysprintf(namefit, NAMESIZE, "fit_GPS_b595_etacr");
    fit_result GPS_b595 = fit_data(argv, paramsj, gjack, lhs_GPS, fit_info_NG, namefit, myenW595);
    print_fit_band_eta(argv, gjack, fit_info_NG, namefit, GPS_b595, paramsj, myenW595);

    fit_info_NG.restore_default();
    printf("//////////////////  NG 595///////////////////////// ");
    fit_info_NG.Nvar = 7;
    fit_info_NG.Npar = 8;
    fit_info_NG.N = 2;
    fit_info_NG.Njack = gjack[0].Njack;
    fit_info_NG.n_ext_P = 2;
    fit_info_NG.ext_P = (double**)malloc(sizeof(double*) * fit_info_NG.n_ext_P);
    fit_info_NG.ext_P[0] = eta_cr_b595;
    fit_info_NG.ext_P[1] = GPS_b595.P[0];
    fit_info_NG.function = rhs_NG_mpcac_MPS2_mu_eta;


    mysprintf(namefit, NAMESIZE, "fit_NG_mpcac_MPS_b595");
    fit_result fit_critical_NG_595 = fit_data(argv, paramsj, gjack, lhs_mpcac_MPS2, fit_info_NG, namefit, myenNG595);
    print_fit_band_eta(argv, gjack, fit_info_NG, namefit, fit_critical_NG_595, paramsj, myenNG595);

    printf("m_pcac(eta_cr)=%g  +-  %g\n", fit_critical_NG_595.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical_NG_595.P[0]));
    printf("M_PS(eta_cr) =%g  +-  %g\n", fit_critical_NG_595.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_critical_NG_595.P[1]));


    // free_fit_result(fit_info_NG, fit_critical_NG_595);
    fit_info_NG.restore_default();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("///////////////////////////////////////////    b575   ///////////////////////////////////////////////////////");
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    printf("/////////////////////////////////////  r_AWI_tau2///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 13,4 };


    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau2_b575", myenW575);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau2_b575", fit_critical, paramsj, myenW575);

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();
    printf("/////////////////////////////////////  r_AWI_tau3///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 14,4 };


    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau3_b575", myenW575);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau3_b575", fit_critical, paramsj, myenW575);

    double* eta_cr_b575 = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0; j < Njack;j++) eta_cr_b575[j] = fit_critical.P[0][j];

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();
    printf("/////////////////////////////////////  r_AWI_tau4///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 15,4 };


    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau4_b575", myenW575);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau4_b575", fit_critical, paramsj, myenW575);

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();

    printf("/////////////////////////////////////  r_AWI_tau5///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 3,4 };



    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau5_b575", myenW575);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau5_b575", fit_critical, paramsj, myenW575);

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();

    printf("/////////////////////////////////////  r_AWI_tau6///////////////////////////////////// ");
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;
    fit_info.corr_id = { 16,4 };

    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau6_b575", myenW575);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau6_b575", fit_critical, paramsj, myenW575);


    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();


    printf("//////////////////  fit GPS///////////////////////// ");
    fit_info_NG.Nvar = 7;
    fit_info_NG.Npar = 7;
    fit_info_NG.N = 1;
    fit_info_NG.Njack = gjack[0].Njack;
    fit_info_NG.n_ext_P = 1;
    fit_info_NG.ext_P = (double**)malloc(sizeof(double*) * fit_info_NG.n_ext_P);
    fit_info_NG.ext_P[0] = eta_cr_b575;
    fit_info_NG.function = rhs_GPS;


    mysprintf(namefit, NAMESIZE, "fit_GPS_b575_etacr");
    fit_result GPS_b575 = fit_data(argv, paramsj, gjack, lhs_GPS, fit_info_NG, namefit, myenW575);
    print_fit_band_eta(argv, gjack, fit_info_NG, namefit, GPS_b575, paramsj, myenW575);

    fit_info_NG.restore_default();
    printf("//////////////////  NG_575///////////////////////// ");
    fit_info_NG.Nvar = 7;
    fit_info_NG.Npar = 8;
    fit_info_NG.N = 2;
    fit_info_NG.Njack = gjack[0].Njack;
    fit_info_NG.n_ext_P = 2;
    fit_info_NG.ext_P = (double**)malloc(sizeof(double*) * fit_info_NG.n_ext_P);
    fit_info_NG.ext_P[0] = eta_cr_b575;
    fit_info_NG.ext_P[1] = GPS_b575.P[0];
    fit_info_NG.function = rhs_NG_mpcac_MPS2_mu_eta;


    mysprintf(namefit, NAMESIZE, "fit_NG_mpcac_MPS_b575");
    fit_result fit_critical_NG_575 = fit_data(argv, paramsj, gjack, lhs_mpcac_MPS2, fit_info_NG, namefit, myenNG575);
    print_fit_band_eta(argv, gjack, fit_info_NG, namefit, fit_critical_NG_575, paramsj, myenNG575);

    printf("m_pcac(eta_cr)=%g  +-  %g\n", fit_critical_NG_575.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical_NG_575.P[0]));
    printf("M_PS(eta_cr) =%g  +-  %g\n", fit_critical_NG_575.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_critical_NG_575.P[1]));


    // free_fit_result(fit_info_NG, fit_critical_NG_575);
    fit_info_NG.restore_default();

    printf("etacr\n");
    printf("5.75:   %g   %g\n", eta_cr_b575[Njack - 1], error_jackboot(argv[1], Njack, eta_cr_b575));
    printf("5.85:   %g   %g\n", eta_cr_b585[Njack - 1], error_jackboot(argv[1], Njack, eta_cr_b585));
    printf("5.95:   %g   %g\n", eta_cr_b595[Njack - 1], error_jackboot(argv[1], Njack, eta_cr_b595));
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////// fit continuum
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int Nbetas = 3;

    data_all jackall;
    jackall.ens = Nbetas;
    jackall.en = new data_single[3];
    for (int e = 0; e < Nbetas; e++) {
        jackall.en[e].Njack = gjack[0].Njack;
        jackall.en[e].Nobs = 3;
        jackall.en[e].jack = double_malloc_2(jackall.en[e].Nobs, Njack);
    }
    for (int j = 0;j < Njack;j++) {
        jackall.en[0].jack[0][j] = 3.289;
        jackall.en[1].jack[0][j] = 4.067;
        jackall.en[2].jack[0][j] = 4.917;

        jackall.en[0].jack[1][j] = fit_critical_NG_575.P[0][j];//m_pcac
        jackall.en[1].jack[1][j] = fit_critical_NG_585.P[0][j];//m_pcac
        jackall.en[2].jack[1][j] = fit_critical_NG_595.P[0][j];//m_pcac

        jackall.en[0].jack[2][j] = fit_critical_NG_575.P[1][j];//M_PS^2
        jackall.en[1].jack[2][j] = fit_critical_NG_585.P[1][j];//
        jackall.en[2].jack[2][j] = fit_critical_NG_595.P[1][j];//

    }
    jackall.en[0].header.beta = 5.75;
    jackall.en[1].header.beta = 5.85;
    jackall.en[2].header.beta = 5.95;

    fit_info.Npar = 2;
    fit_info.N = 1;
    fit_info.Nvar = 1;
    fit_info.Njack = Njack;
    fit_info.myen = { 0,1,2 };
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    int count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(1. / jackall.en[e].jack[0][j], 2);// a/r0
            }
            count++;
        }
    }
    fit_info.function = rhs_one_line;
    fit_result amu_SD_l_common = fit_all_data(argv, jackall, lhs_continuum_mpcac, fit_info, "mpcac_continuum");
    fit_info.band_range = { 0,0.1 };
    print_fit_band(argv, jackall, fit_info, fit_info, "mpcac_continuum", "a_r02", amu_SD_l_common, amu_SD_l_common, 0, fit_info.myen.size() - 1, 0.001);

    free_fit_result(fit_info, amu_SD_l_common);
    amu_SD_l_common = fit_all_data(argv, jackall, lhs_continuum_MPS, fit_info, "MPS_continuum");
    fit_info.band_range = { 0,0.1 };
    print_fit_band(argv, jackall, fit_info, fit_info, "MPS_continuum", "a_r02", amu_SD_l_common, amu_SD_l_common, 0, fit_info.myen.size() - 1, 0.001);



    fit_info.restore_default();
    // int Nbetas = 3;
    // int Njack = gjack[0].Njack;
    // double*** x = double_malloc_3(Njack, 3, 1);
    // double*** y = double_malloc_3(Njack, 3, 2);
    // double* mawi = (double*)malloc(sizeof(double) * Njack);
    // int en[1] = { Nbetas };
    // double guess[1] = { 1 };
    // for (int j = 0;j < Njack;j++) {
    //     x[j][0][0] = 3.289;
    //     x[j][1][0] = 4.067;
    //     x[j][2][0] = 4.917;
    //     y[j][0][0] = fit_critical_NG_575.P[0][j];
    //     y[j][1][0] = fit_critical_NG_585.P[0][j];
    //     y[j][2][0] = fit_critical_NG_595.P[0][j];
    //     y[j][0][1] = error_jackboot(argv[1], Njack, fit_critical_NG_575.P[0]);
    //     y[j][1][1] = error_jackboot(argv[1], Njack, fit_critical_NG_585.P[0]);
    //     y[j][2][1] = error_jackboot(argv[1], Njack, fit_critical_NG_595.P[0]);
    //     mawi[j] = non_linear_fit_Nf(1, en, x[j], y[j], 1, 1, , guess);
    // }
    return 0;
}
