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


    //////


    int NeW575 = dataj.size() - NeW - NeW595;
    printf("number of ensembles W = %d\n", NeW575);
    vector<int> myenW575(NeW575);
    for (int i = 0;i < NeW575; i++)  myenW575[i] = i + NeW + NeW595;

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////   NG
    ///////////////////////////////////////////////////////////////////////////////////////////////////


    NeNG = dataj.size() - NeW;
    printf("number of ensembles NG = %d\n", NeNG);
    vector<int> myenNG(NeNG);
    for (int i = NeW;i < NeW + NeNG; i++)  myenNG[i - NeW] = i;

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


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////////////////////////////////////////////////////////////////////////
    struct fit_type fit_info;
    struct fit_result  fit_critical;
    char namefit[NAMESIZE];
    fit_info.Nvar = 7;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu;


    fit_critical = fit_data(argv, paramsj, gjack, lhs_critical_eta_mu_m0, fit_info, "eta_critical_b585", myenW);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_b585", fit_critical, paramsj, myenW);

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


    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, "eta_critical_tau3_b585", myenW);
    print_fit_band_eta(argv, gjack, fit_info, "eta_critical_tau3_b585", fit_critical, paramsj, myenW);

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

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();


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


    return 0;
}
