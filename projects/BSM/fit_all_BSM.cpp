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


void print_fit_band_eta(char** argv, vector<data_BSM> gjack, struct fit_type fit_info, const char* label, struct fit_result fit_out, vector<header_BSM> params, std::vector<int> myenW, std::vector<double> range={-1.7,-1} ) {
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
                    tmpx[3] = range[0] + (range[1]-range[0])*i / 100.0;
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
    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m1.0983_M02_-0.010396_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.098300_csw1.000000_mu030.022400_m0-0.010396", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    //1
//      mysprintf(namefile,NAMESIZE,"%s/b5.85/L20T40/eta_m1.0983_M02_-0.024604_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.098300_csw1.000000_mu030.022400_m0-0.024604",argv[2],argv[1]);
//      emplace_back_par_data(namefile,paramsj,dataj);
    //  1
    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m1.0983_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.098300_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    //      mysprintf(namefile,NAMESIZE,"%s/b5.85/L20T40/eta_m1.2944_M02_-0.024604_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.294400_csw1.000000_mu030.022400_m0-0.024604",argv[2],argv[1]);
    //      emplace_back_par_data(namefile,paramsj,dataj);
    //      mysprintf(namefile,NAMESIZE,"%s/b5.85/L20T40/eta_m1.2944_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.294400_csw1.000000_mu030.022400_m0-0.040000",argv[2],argv[1]);
    //      emplace_back_par_data(namefile,paramsj,dataj);
    //      mysprintf(namefile,NAMESIZE,"%s/b5.85/L20T40/eta_m1.2944_M02_-0.024604_mu03_0.0120_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.294400_csw1.000000_mu030.012000_m0-0.024604",argv[2],argv[1]);
    //      emplace_back_par_data(namefile,paramsj,dataj);
    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m1.2944_M02_-0.010396_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.294400_csw1.000000_mu030.022400_m0-0.010396", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m1.2944_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.294400_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m1.2944_M02_-0.040000_mu03_0.0120_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.294400_csw1.000000_mu030.012000_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m1.4063_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.406300_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    NeW = dataj.size();
    printf("number of ensembles W = %d\n", NeW);
    vector<int> myenW(NeW);
    for (int i = 0;i < NeW; i++)  myenW[i] = i;


    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m1.2_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96_NG/jackknife/%s_T40_L20_rho1.960000_eta-1.200000_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m1.3_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96_NG/jackknife/%s_T40_L20_rho1.960000_eta-1.300000_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m1.4063_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96_NG/jackknife/%s_T40_L20_rho1.960000_eta-1.406300_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m1.4063_M02_-0.040000_mu03_0.0120_csw_1.0_rho1.96_NG/jackknife/%s_T40_L20_rho1.960000_eta-1.406300_csw1.000000_mu030.012000_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m1.45_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96_NG/jackknife/%s_T40_L20_rho1.960000_eta-1.450000_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m1.45_M02_-0.040000_mu03_0.0120_csw_1.0_rho1.96_NG/jackknife/%s_T40_L20_rho1.960000_eta-1.450000_csw1.000000_mu030.012000_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m1.5_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96_NG/jackknife/%s_T40_L20_rho1.960000_eta-1.500000_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m1.6_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96_NG/jackknife/%s_T40_L20_rho1.960000_eta-1.600000_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    NeNG = dataj.size() - NeW;
    printf("number of ensembles NG = %d\n", NeNG);
    vector<int> myenNG(NeNG);
    for (int i = NeW;i < NeW + NeNG; i++)  myenNG[i - NeW] = i;


    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m1.5_M02_-0.040000_mu03_0.0224_csw_1.0_rho3/jackknife/%s_T40_L20_rho3.000000_eta-1.500000_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m2.0_M02_-0.010396_mu03_0.0224_csw_1.0_rho3/jackknife/%s_T40_L20_rho3.000000_eta-2.000000_csw1.000000_mu030.022400_m0-0.010396", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m2.0_M02_-0.040000_mu03_0.0224_csw_1.0_rho3/jackknife/%s_T40_L20_rho3.000000_eta-2.000000_csw1.000000_mu030.022400_m0-0.040000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m2.0_M02_-0.045000_mu03_0.0224_csw_1.0_rho3/jackknife/%s_T40_L20_rho3.000000_eta-2.000000_csw1.000000_mu030.022400_m0-0.045000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m2.0_M02_-0.045000_mu03_0.0400_csw_1.0_rho3/jackknife/%s_T40_L20_rho3.000000_eta-2.000000_csw1.000000_mu030.040000_m0-0.045000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/eta_m2.1_M02_-0.045000_mu03_0.0224_csw_1.0_rho3/jackknife/%s_T40_L20_rho3.000000_eta-2.100000_csw1.000000_mu030.022400_m0-0.045000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    

    int NeW3 = dataj.size() - (NeW + NeNG);
    printf("number of ensembles W(rho3) = %d\n", NeW3);
    vector<int> myenW3(NeW3);
    for (int i = (NeW + NeNG);i < (NeW + NeNG) + NeW3; i++)  myenW3[i - (NeW + NeNG)] = i;

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m2.05_M02_-0.045000_mu03_0.0224_csw_1.0_rho3_NG/jackknife/%s_T40_L20_rho3.000000_eta-2.050000_csw1.000000_mu030.022400_m0-0.045000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);


    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m2.05_M02_-0.045000_mu03_0.0400_csw_1.0_rho3_NG/jackknife/%s_T40_L20_rho3.000000_eta-2.050000_csw1.000000_mu030.040000_m0-0.045000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m2.12_M02_-0.045000_mu03_0.0224_csw_1.0_rho3_NG/jackknife/%s_T40_L20_rho3.000000_eta-2.120000_csw1.000000_mu030.022400_m0-0.045000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    mysprintf(namefile, NAMESIZE, "%s/b5.85/L20T40/NG/eta_m2.12_M02_-0.045000_mu03_0.0400_csw_1.0_rho3_NG/jackknife/%s_T40_L20_rho3.000000_eta-2.120000_csw1.000000_mu030.040000_m0-0.045000", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);

    
    
    int NeNG3 = dataj.size() - (NeW + NeNG+ NeW3);
    printf("number of ensembles NG(rho3) = %d\n", NeNG3);
    vector<int> myenNG3(NeNG3);
    for (int i = (NeW + NeNG+ NeW3);i < (NeW + NeNG+ NeW3) + NeW3; i++)  myenNG3[i - (NeW + NeNG+ NeW3)] = i;


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
    fit_info.Npar = 8;
    fit_info.N = 2;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu_m0;


    fit_critical = fit_data(argv, paramsj, gjack, lhs_critical_eta_mu_m0, fit_info, "eta_m0_critical_b585", myenW);
    print_fit_band_eta(argv, gjack, fit_info, "eta_m0_critical_b585", fit_critical, paramsj, myenW);

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.Nvar = 7;
    fit_info.Npar = 8;
    fit_info.N = 2;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu_m0_shifted;


    mysprintf(namefit, NAMESIZE, "eta_m0_critical_b585_shifted");
    fit_critical = fit_data(argv, paramsj, gjack, lhs_critical_eta_mu_m0, fit_info, namefit, myenW);
    print_fit_band_eta(argv, gjack, fit_info, namefit, fit_critical, paramsj, myenW);

    printf("eta_cr=%g  +-  %g\n", fit_critical.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[0]));
    printf("m0_cr =%g  +-  %g\n", fit_critical.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[1]));

    free_fit_result(fit_info, fit_critical);
    fit_info.restore_default();
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.Nvar = 7;
    fit_info.Npar = 6;
    fit_info.N = 2;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu_m0_simple;


    mysprintf(namefit, NAMESIZE, "eta_m0_critical_b585_simple");
    fit_critical = fit_data(argv, paramsj, gjack, lhs_critical_eta_mu_m0, fit_info, namefit, myenW);
    print_fit_band_eta(argv, gjack, fit_info, namefit, fit_critical, paramsj, myenW);

    printf("eta_cr=%g  +-  %g\n", fit_critical.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[0]));
    printf("m0_cr =%g  +-  %g\n", fit_critical.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[1]));
    free_fit_result(fit_info, fit_critical);


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // local quantities fitting
    //////////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.Nvar = 7;
    fit_info.Npar = 6;
    fit_info.N = 2;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu_m0_simple;


    mysprintf(namefit, NAMESIZE, "eta_m0_loc_critical_b585_simple");
    fit_critical = fit_data(argv, paramsj, gjack, lhs_critical_eta_mu_m0_loc, fit_info, namefit, myenW);
    print_fit_band_eta(argv, gjack, fit_info, namefit, fit_critical, paramsj, myenW);

    printf("eta_cr=%g  +-  %g\n", fit_critical.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[0]));
    printf("m0_cr =%g  +-  %g\n", fit_critical.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[1]));
    free_fit_result(fit_info, fit_critical);

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // simple fit tau2 
    //////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n r_awi tau2\n\n");
    fit_info.Nvar = 7;
    fit_info.Npar = 6;
    fit_info.N = 2;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu_m0_simple;
    fit_info.corr_id = { 13,4 }; //r_awi_tau2 m_pcac

    mysprintf(namefit, NAMESIZE, "r_awi_tau2_critical_b585_simple");
    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, namefit, myenW);
    print_fit_band_eta(argv, gjack, fit_info, namefit, fit_critical, paramsj, myenW);

    printf("eta_cr=%g  +-  %g\n", fit_critical.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[0]));
    printf("m0_cr =%g  +-  %g\n", fit_critical.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[1]));
    free_fit_result(fit_info, fit_critical);


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // simple fit tau2  loc
    //////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n r_awi tau2 loc\n\n");
    fit_info.Nvar = 7;
    fit_info.Npar = 6;
    fit_info.N = 2;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu_m0_simple;
    fit_info.corr_id = { 17,4 }; //r_awi_tau2 m_pcac

    mysprintf(namefit, NAMESIZE, "r_awi_loc_tau2_critical_b585_simple");
    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, namefit, myenW);
    print_fit_band_eta(argv, gjack, fit_info, namefit, fit_critical, paramsj, myenW);

    printf("eta_cr=%g  +-  %g\n", fit_critical.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[0]));
    printf("m0_cr =%g  +-  %g\n", fit_critical.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[1]));
    // free_fit_result(fit_info, fit_critical);




    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // NG
    //////////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.Nvar = 7;
    fit_info.Npar = 7;
    fit_info.N = 2;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 2;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
    fit_info.ext_P[0] = fit_critical.P[0];
    fit_info.ext_P[1] = fit_critical.P[1];
    fit_info.function = rhs_NG_mpcac_MPS2;


    mysprintf(namefit, NAMESIZE, "fit_NG_mpcac_MPS_b585");
    fit_result fit_NG = fit_data(argv, paramsj, gjack, lhs_mpcac_MPS2, fit_info, namefit, myenNG);
    print_fit_band_eta(argv, gjack, fit_info, namefit, fit_NG, paramsj, myenNG);

    printf("m_pcac(eta_cr)=%g  +-  %g\n", fit_NG.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_NG.P[0]));
    printf("M_PS(eta_cr) =%g  +-  %g\n", fit_NG.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_NG.P[1]));

    fit_info.restore_default();



    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n rho=3 \n\n");
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // simple fit tau2  loc
    //////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n r_awi tau2 loc rho3\n\n");
    fit_info.Nvar = 7;
    fit_info.Npar = 6;
    fit_info.N = 2;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu_m0_simple;
    fit_info.corr_id = { 17,4 }; //r_awi_tau2 m_pcac

    mysprintf(namefit, NAMESIZE, "r_awi_loc_tau2_critical_b585_rho3_simple");
    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, namefit, myenW3);
    print_fit_band_eta(argv, gjack, fit_info, namefit, fit_critical, paramsj, myenW3, {-2.5,-1.5});

    printf("eta_cr=%g  +-  %g\n", fit_critical.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[0]));
    printf("m0_cr =%g  +-  %g\n", fit_critical.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[1]));
    free_fit_result(fit_info, fit_critical);

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // simple fit tau2  loc
    //////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n r_awi tau2 loc rho3\n\n");
    fit_info.Nvar = 7;
    fit_info.Npar = 7;
    fit_info.N = 2;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = rhs_critical_eta_mu_m0_7par;
    fit_info.corr_id = { 17,4 }; //r_awi_tau2 m_pcac

    mysprintf(namefit, NAMESIZE, "r_awi_loc_tau2_critical_b585_rho3_7par");
    fit_critical = fit_data(argv, paramsj, gjack, lhs_fit_two_func, fit_info, namefit, myenW3);
    print_fit_band_eta(argv, gjack, fit_info, namefit, fit_critical, paramsj, myenW3, {-2.5,-1.5});

    printf("eta_cr=%g  +-  %g\n", fit_critical.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[0]));
    printf("m0_cr =%g  +-  %g\n", fit_critical.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_critical.P[1]));
    //free_fit_result(fit_info, fit_critical);


    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n   NG rho 3  \n");
    //////////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.Nvar = 7;
    fit_info.Npar = 6;
    fit_info.N = 2;
    fit_info.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 2;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
    fit_info.ext_P[0] = fit_critical.P[0];
    fit_info.ext_P[1] = fit_critical.P[1];
    fit_info.function = rhs_NG_mpcac_MPS2_twoline;


    mysprintf(namefit, NAMESIZE, "fit_NG_mpcac_MPS_b585_rho3_lines");
    fit_result fit_NG3 = fit_data(argv, paramsj, gjack, lhs_mpcac_MPS2, fit_info, namefit, myenNG3);
    print_fit_band_eta(argv, gjack, fit_info, namefit, fit_NG3, paramsj, myenNG3, {-2.5,-1.5});

    printf("m_pcac(eta_cr)=%g  +-  %g\n", fit_NG3.P[0][Njack - 1], error_jackboot(argv[1], Njack, fit_NG3.P[0]));
    printf("M_PS(eta_cr)  =%g  +-  %g\n", fit_NG3.P[1][Njack - 1], error_jackboot(argv[1], Njack, fit_NG3.P[1]));

    fit_info.restore_default();





    return 0;
}
