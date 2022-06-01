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
#include "header_phi4.hpp"
#include "zeta_interpolation.hpp"
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>


//local folder
// #include "mass_phi4.hpp"
#include "fit_function.hpp"

#include <QC3_interface.hpp>
#include "fit_all.hpp"
#include "extra_func_phi4.hpp"


using namespace std;
int Ne = 0;


void print_fit_band_L_M(char** argv, vector<data_phi> gjack, struct fit_type fit_info, struct fit_type fit_info_m0, const char* label, struct fit_result fit_out, struct fit_result fit_out_m0, vector<cluster::IO_params> params, std::vector<int> myen, std::vector<int> Lrange = { 16,50 }) {
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double** tif_m0 = swap_indices(fit_info_m0.Npar, Njack, fit_out_m0.P);

    for (int n = 0;n < N; n++) {

        mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_n%d_L.txt", argv[3], label, n);
        f = open_file(namefile, "w+");
        double* tmpx = (double*)malloc(sizeof(double) * Nvar);
        double* tmpy = (double*)malloc(sizeof(double) * Njack);
        printf("writing: %s\n", namefile);

        for (int i = Lrange[0]; i < Lrange[1]; i++) {
            double finalL = i;
            tmpx[0] = finalL;
            double* E3_m = (double*)malloc(sizeof(double) * Njack);
            for (int j = 0;j < Njack;j++) {
                E3_m[j] = lhs_E3_m(n, 0, j, params, gjack, fit_info);
            }
            double E3_m_err = error_jackboot(argv[1], Njack, E3_m);
            for (int j = 0;j < Njack;j++) {
                int eL = -1;
                for (int e = 0; e < myen.size();e++) {
                    if (tmpx[0] == params[e].data.L[1])
                        eL = e;
                }
                if (eL == -1)
                    tmpx[1] = fit_info_m0.function(0, fit_info_m0.Nvar, tmpx, fit_info_m0.Npar, tif_m0[j]); //m0   put for each n the mass of the last ensemble
                else
                    tmpx[1] = gjack[eL].jack[1][j];
                // tmpx[1] = fit_out_m0.P[0][j];

                tmpx[2] = gjack[0].jack[2][j];//m1  //
                tmpx[3] = gjack[0].jack[4][j];//E20
                tmpx[4] = gjack[0].jack[5][j];//E21
                tmpx[5] = (double)params[0].data.L[0];//T
                tmpx[6] = 1;//k
                tmpx[7] = 1;//MpiL

                double* x = (double*)malloc(sizeof(double) * myen.size());
                double* y = (double*)malloc(sizeof(double) * myen.size());
                for (int e = 0; e < myen.size();e++) {
                    y[e] = lhs_E3_m(n, myen[e], j, params, gjack, fit_info);// E3( \vec{n} )/mass
                    if (params[0].data.gC > 0) {
                        y[e] = lhs_E3orE1_m_complex(n, myen[e], j, params, gjack, fit_info);
                    }
                    x[e] = params[myen[e]].data.L[1];
                }
                tmpx[8] = inter_spline(tmpx[0], myen.size(), x, y); // tmpx[0]=L 

                if (params[0].data.gC > 0) {
                    // tmpx[1] = gjack[0].jack[443][j];//m0
                    tmpx[10] = compute_k_m_g(n, 0, j, params, gjack, fit_info);//k_m
                    tmpx[6] = tmpx[10] * tmpx[1];//k
                    tmpx[3] = gjack[0].jack[594][j];//E20
                    tmpx[7] = gjack[0].jack[443][j] * (double)params[0].data.L[1] / (2. * pi_greco);//mL_2pi
                }

                free(x);free(y);
                //tmpx[8]=E3_m[j];// E3( \vec{n} )/mass

                // as error on E3/m  we take the error on the ensemble=0
                tmpx[9] = E3_m_err;


                // for (int i = fit_info.Nvar; i < fit_info.Nvar + fit_info.n_ext_P; i++)
                // 	tmpx[i] = fit_info.ext_P[i - fit_info.Nvar][j];
                for (int i = 0; i < fit_info.n_ext_P; i++) {
                    // printf("HERE %d  %d\n",i + fit_info.Nvar,Nvar);
                    tmpx[i + fit_info.Nvar] = fit_info.ext_P[i][j];
                }


                tmpy[j] = fit_info.function(n, Nvar, tmpx, Npar, tif[j]);
                //                 tmpy[j]=tmpx[8];
                //                 printf("L=%g  j=%d  E=%g\n",finalL,j,tmpx[8]);
                //                 if (fabs(finalL-36)<1e-5  && n==2)
                //                     printf("%g   j=%d \t k=%g  m=%g P0=%g  P1=%g\n",finalL,j,tmpy[j], tmpx[1], tif[j][0] ,tif[j][1]  );
                //                 
            }
            // if (fabs(finalL-20)<1e-5  && n==0){
            // 	printf("%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
            // 	printf(" %g  %g\n",fit_out.P[0][Njack-1],fit_out.P[1][Njack-1] );
            // }
            fprintf(f, "%g  \t %g  %g\n", finalL, tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));

            free(E3_m);
        }

        free(tmpy);free(tmpx);
        fclose(f);
    }
    free_2(Njack, tif);
    free_2(Njack, tif_m0);

}





void print_fit_band_E3_vs_L(char** argv, vector<data_phi> gjack, struct fit_type fit_info, struct fit_type fit_info_m0, const char* label, struct fit_result fit_out, struct fit_result fit_out_m0, vector<cluster::IO_params> params, std::vector<int> myen, struct fit_type fit_info_E3_poly, fit_result fit_E3_poly, std::vector<int> Lrange = { 16,50 }) {
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double** tif_m0 = swap_indices(fit_info_m0.Npar, Njack, fit_out_m0.P);
    double** tif_E3_poly = swap_indices(fit_info_E3_poly.Npar, Njack, fit_E3_poly.P);


    for (int n = 0;n < N; n++) {

        mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_n%d_L.txt", argv[3], label, n);
        f = open_file(namefile, "w+");
        double* tmpx = (double*)malloc(sizeof(double) * Nvar);
        double* tmpy = (double*)malloc(sizeof(double) * Njack);
        printf("writing: %s\n", namefile);

        for (int i = Lrange[0]; i < Lrange[1]; i++) {
            double finalL = i;
            tmpx[0] = finalL;
            double* E3_m = (double*)malloc(sizeof(double) * Njack);
            for (int j = 0;j < Njack;j++) {
                E3_m[j] = fit_info_E3_poly.function(n, fit_info_E3_poly.Nvar, tmpx, fit_info_E3_poly.Npar, tif_E3_poly[j]);
            }
            double E3_m_err = error_jackboot(argv[1], Njack, E3_m);



            for (int j = 0;j < Njack;j++) {
                tmpx[1] = fit_info_m0.function(0, fit_info_m0.Nvar, tmpx, fit_info_m0.Npar, tif_m0[j]); //m0   put for each n the mass of the last ensemble

                tmpx[2] = gjack[0].jack[2][j];//m1  //
                tmpx[3] = gjack[0].jack[4][j];//E20
                tmpx[4] = gjack[0].jack[5][j];//E21
                tmpx[5] = (double)params[0].data.L[0];//T
                tmpx[6] = 1;//k
                tmpx[7] = 1;//MpiL

                //                 double *x=(double*) malloc(sizeof(double)*myen.size());
                //                 double *y=(double*) malloc(sizeof(double)*myen.size());
                //                 for (int e=0; e<myen.size();e++){
                //                     y[e]= lhs_E3_m(n,myen[e],j,params,gjack,fit_info);// E3( \vec{n} )/mass
                //                     x[e]=params[myen[e]].data.L[1];
                //                 }
                //                 tmpx[8]=inter_spline( tmpx[0], myen.size(), x, y   ); // tmpx[0]=L 
                //                 free(x);free(y);

                tmpx[8] = E3_m[Njack - 1];// E3( \vec{n} )/mass
                // as error on E3/m  we take the error on the ensemble=0
                tmpx[9] = E3_m_err;

                tmpx[11] = E3_m[Njack - 1];// E3( \vec{n} )/mass
                // as error on E3/m  we take the error on the ensemble=0
                tmpx[12] = E3_m_err;

                if (params[0].data.gC > 0) {
                    tmpx[1] = gjack[0].jack[443][j];//m0
                    tmpx[10] = compute_k_m_g(n, 0, j, params, gjack, fit_info);//k_m
                    tmpx[6] = tmpx[10] * tmpx[1];//k
                    tmpx[3] = gjack[0].jack[594][j];//E20
                    tmpx[7] = gjack[0].jack[443][j] * (double)params[0].data.L[1] / (2. * pi_greco);//mL_2pi
                }

                for (int i = fit_info.Nvar; i < fit_info.Nvar + fit_info.n_ext_P; i++)
                    tmpx[i] = fit_info.ext_P[i - fit_info.Nvar][j];



                tmpy[j] = fit_info.function(n, Nvar, tmpx, Npar, tif[j]);
                //                 tmpy[j]=tmpx[8];
                //                 printf("L=%g  j=%d  E=%g\n",finalL,j,tmpx[8]);
                //                 if (fabs(finalL-36)<1e-5  && n==2)
                //                     printf("%g   j=%d \t k=%g  m=%g P0=%g  P1=%g\n",finalL,j,tmpy[j], tmpx[1], tif[j][0] ,tif[j][1]  );
                //                 
            }
            /*if (fabs(finalL-36)<1e-5  && n==2){
             *                printf("%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
        }*/
            fprintf(f, "%g  \t %g  %g\n", finalL, tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));

            free(E3_m);
        }

        free(tmpy);free(tmpx);
        fclose(f);
    }
    free_2(Njack, tif);
    free_2(Njack, tif_m0);
    free_2(Njack, tif_E3_poly);
}





void print_kiso_P0_inf_L_M(char** argv, vector<data_phi> gjack, struct fit_type fit_info, struct fit_type fit_info_m0, const char* label, struct fit_result fit_out, struct fit_result fit_out_m0, vector<cluster::IO_params> params, std::vector<int> myen, struct fit_type fit_info_E3_poly, fit_result fit_E3_poly, std::vector<int> Lrange = { 16,50 }) {
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    double** tif = double_malloc_2(Njack, Npar);
    for (int i = 0;i < Npar; i++)
        for (int j = 0;j < Njack;j++)
            tif[j][i] = 0;
    double** tif_m0 = swap_indices(fit_info_m0.Npar, Njack, fit_out_m0.P);
    double** tif_E3_poly = swap_indices(fit_info_E3_poly.Npar, Njack, fit_E3_poly.P);

    for (int n = 0;n < N; n++) {

        mysprintf(namefile, NAMESIZE, "%s/kiso_P0_n%d_L.txt", argv[3], n);
        f = open_file(namefile, "w+");
        double* tmpx = (double*)malloc(sizeof(double) * Nvar);
        double* tmpy = (double*)malloc(sizeof(double) * Njack);
        printf("writing: %s\n", namefile);

        for (int i = Lrange[0]; i < Lrange[1]; i++) {
            double finalL = i;
            tmpx[0] = finalL;
            double* E3_m = (double*)malloc(sizeof(double) * Njack);
            for (int j = 0;j < Njack;j++) {
                E3_m[j] = fit_info_E3_poly.function(n, fit_info_E3_poly.Nvar, tmpx, fit_info_E3_poly.Npar, tif_E3_poly[j]);
            }
            double E3_m_err = error_jackboot(argv[1], Njack, E3_m);
            for (int j = 0;j < Njack;j++) {
                tmpx[1] = fit_info_m0.function(0, fit_info_m0.Nvar, tmpx, fit_info_m0.Npar, tif_m0[j]); //m0   put for each n the mass of the last ensemble

                tmpx[2] = gjack[0].jack[2][j];//m1  //
                tmpx[3] = gjack[0].jack[4][j];//E20
                tmpx[4] = gjack[0].jack[5][j];//E21
                tmpx[5] = (double)params[0].data.L[0];//T
                tmpx[6] = 1;//k
                tmpx[7] = 1;//MpiL


                tmpx[8] = E3_m[Njack - 1];// E3( \vec{n} )/mass

                tmpx[9] = E3_m_err;


                for (int i = fit_info.Nvar; i < fit_info.Nvar + fit_info.n_ext_P; i++)
                    tmpx[i] = fit_info.ext_P[i - fit_info.Nvar][j];



                tmpy[j] = fit_info.function(n, Nvar, tmpx, Npar, tif[j]);

            }

            fprintf(f, "%g  \t %g  %g\n", finalL, tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));

            free(E3_m);
        }

        free(tmpy);free(tmpx);
        fclose(f);
    }
    for (int i = 0;i < Npar; i++)
        for (int j = 0;j < Njack;j++)
            tif[j][i] = -1e+3;
    for (int n = 0;n < N; n++) {

        mysprintf(namefile, NAMESIZE, "%s/kiso_P-1e+3_n%d_L.txt", argv[3], n);
        f = open_file(namefile, "w+");
        double* tmpx = (double*)malloc(sizeof(double) * Nvar);
        double* tmpy = (double*)malloc(sizeof(double) * Njack);
        printf("writing: %s\n", namefile);

        for (int i = Lrange[0]; i < Lrange[1]; i++) {
            double finalL = i;
            tmpx[0] = finalL;
            double* E3_m = (double*)malloc(sizeof(double) * Njack);
            for (int j = 0;j < Njack;j++) {
                E3_m[j] = fit_info_E3_poly.function(n, fit_info_E3_poly.Nvar, tmpx, fit_info_E3_poly.Npar, tif_E3_poly[j]);
            }
            double E3_m_err = error_jackboot(argv[1], Njack, E3_m);
            for (int j = 0;j < Njack;j++) {
                tmpx[1] = fit_info_m0.function(0, fit_info_m0.Nvar, tmpx, fit_info_m0.Npar, tif_m0[j]); //m0   put for each n the mass of the last ensemble

                tmpx[2] = gjack[0].jack[2][j];//m1  //
                tmpx[3] = gjack[0].jack[4][j];//E20
                tmpx[4] = gjack[0].jack[5][j];//E21
                tmpx[5] = (double)params[0].data.L[0];//T
                tmpx[6] = 1;//k
                tmpx[7] = 1;//MpiL


                tmpx[8] = E3_m[Njack - 1];// E3( \vec{n} )/mass

                tmpx[9] = E3_m_err;


                for (int i = fit_info.Nvar; i < fit_info.Nvar + fit_info.n_ext_P; i++)
                    tmpx[i] = fit_info.ext_P[i - fit_info.Nvar][j];



                tmpy[j] = fit_info.function(n, Nvar, tmpx, Npar, tif[j]);

            }

            fprintf(f, "%g  \t %g  %g\n", finalL, tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));

            free(E3_m);
        }

        free(tmpy);free(tmpx);
        fclose(f);
    }

    free_2(Njack, tif);
    free_2(Njack, tif_m0);
    free_2(Njack, tif_E3_poly);
}

void print_phase_shift(char** argv, vector<data_phi> gjack, struct fit_type fit_info, const char* label, struct fit_result fit_out) {

    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s_fit_phase_shift.txt", argv[3], label);
    FILE* f = open_file(namefile, "w+");
    int Npar = fit_info.Npar;
    int Njack = gjack[0].Njack;

    int ikmax = 100;
    double dk = 0.03;
    double* delta = (double*)malloc(sizeof(double) * Njack);
    for (int ik = 0;ik < ikmax;ik++) {

        double k_m = ik * dk;
        for (int j = 0; j < Njack;j++) {

            double a0m0 = fit_out.P[0][j];
            double r0m0 = fit_out.P[1][j];

            double kcotdelta_m = 1.0 / a0m0 + +r0m0 * k_m * k_m / 2.;  //   (k cot(d) )/ mass
            if (Npar >= 3) {
                kcotdelta_m += fit_out.P[2][j] * r0m0 * r0m0 * r0m0 * k_m * k_m * k_m * k_m;
            }
            delta[j] = std::atan(k_m / kcotdelta_m);

        }
        fprintf(f, "%g    %g   %g\n", k_m, delta[Njack - 1], error_jackboot(argv[1], Njack, delta));
    }
    free(delta);
    fclose(f);
}


void read_single_Njack_Nobs(FILE* stream, int header_size, int& Njack, int& Nobs) {

    long int tmp;
    int s = header_size;

    size_t i = fread(&Njack, sizeof(int), 1, stream);


    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    tmp -= header_size + sizeof(int);

    s = Njack;

    Nobs = (tmp) / ((s) * sizeof(double));

    fseek(stream, header_size + sizeof(int), SEEK_SET);



}

// void read_single_dataj(FILE* stream, cluster::IO_params params, data_single* dj) {

// 	int Njack;
// 	int Nobs;
// 	read_single_Njack_Nobs(stream, params.data.header_size, Njack, Nobs);

// 	dj->jack = double_malloc_2(Nobs, Njack);
// 	dj->Njack = Njack;
// 	dj->Nobs = Nobs;

// 	size_t i = 0;
// 	for (int obs = 0; obs < dj->Nobs; obs++) {
// 		i += fread(dj->jack[obs], sizeof(double), dj->Njack, stream);
// 	}
// 	dj->header.L = params.data.L[1];
// 	dj->header.T = params.data.L[0];

// }

data_single read_single_dataj(FILE* stream, cluster::IO_params params) {

    int Njack;
    int Nobs;
    read_single_Njack_Nobs(stream, params.data.header_size, Njack, Nobs);
    // data_single dj(Nobs,Njack);
    data_single dj;
    dj.jack = double_malloc_2(Nobs, Njack);
    dj.Nobs = Nobs;
    dj.Njack = Njack;
    //
    size_t i = 0;
    for (int obs = 0; obs < dj.Nobs; obs++) {
        i += fread(dj.jack[obs], sizeof(double), dj.Njack, stream);
    }
    dj.header.L = params.data.L[1];
    dj.header.T = params.data.L[0];
    return dj;

}

void read_all_the_files(std::vector<std::string> files, const char* resampling, data_all* jackall) {
    jackall->resampling = resampling;
    //jackall->en = (data_single*)malloc(sizeof(data_single) * files.size());
    jackall->en = new data_single[files.size()];
    jackall->ens = files.size();
    int count = 0;
    for (std::string s : files) {
        cluster::IO_params params;
        FILE* f = open_file(s.c_str(), "r");
        read_header_phi4(f, params);
        // read_single_dataj(f, params, &(jackall->en[count]));
        jackall->en[count] = read_single_dataj(f, params);
        jackall->en[count].resampling = resampling;
        count++;
        fclose(f);
    }

}
void print_test(data_all jackall) {

    printf("GEVP_E2_01 =%f   %f\n", jackall.en[1].jack[19][jackall.en[1].Njack - 1], error_jackboot("jack", jackall.en[1].Njack, jackall.en[1].jack[19]));
    for (int e = 0; e < jackall.ens; e++) {
        printf("%d  %d\n", jackall.en[e].header.T, jackall.en[e].header.L);
    }
}


int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");



    vector<cluster::IO_params> paramsj;
    vector<data_phi> dataj;

    int Ne = 0;
    cluster::IO_params params;

    char namefile[NAMESIZE];
    char jackboot[NAMESIZE];
    char resampling[NAMESIZE];
    mysprintf(jackboot, NAMESIZE, "%s", argv[1]);
    mysprintf(resampling, NAMESIZE, "%s", argv[1]);


    std::vector<std::string> files;
    // mysprintf(namefile, NAMESIZE, "%s/%s_G2t_T64_L20_msq0-1.230000_msq1-0.370000_l00.000000_l10.000000_mu0.000000_g10.000000_rep0", argv[2], argv[1]);
    // emplace_back_par_data(namefile, paramsj, dataj);
    // files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_G2t_T64_L20_msq0-1.241300_msq1-0.380000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_G2t_T64_L21_msq0-1.241300_msq1-0.380000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_G2t_T64_L22_msq0-1.241300_msq1-0.380000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_G2t_T64_L23_msq0-1.241300_msq1-0.380000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_G2t_T64_L24_msq0-1.241300_msq1-0.380000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_G2t_T64_L25_msq0-1.241300_msq1-0.380000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0", argv[2], argv[1]);
    emplace_back_par_data(namefile, paramsj, dataj);
    files.emplace_back(namefile);



    data_all jackall;
    read_all_the_files(files, argv[1], &jackall);
    jackall.init_error();


    printf("E1_0 =%f   %f\n", dataj[0].jack[1][dataj[0].Njack - 1], error_jackboot(argv[1], dataj[0].Njack, dataj[0].jack[1]));


    vector<data_phi> gjack = create_generalised_resampling_phi(dataj, paramsj);
    printf("GEVP_E2_01 =%f   %f\n", gjack[1].jack[19][gjack[1].Njack - 1], error_jackboot(argv[1], gjack[1].Njack, gjack[1].jack[19]));
    print_test(jackall);
    Ne = gjack.size();
    printf("number of ensembles = %d\n", Ne);

    std::vector<int> myen(Ne);
    for (int i = 0;i < Ne; i++)  myen[i] = i;

    int Njack = gjack[0].Njack;
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // init zeta func
    //////////////////////////////////////////////////////////////////////////////////////////////////



    std::vector<int> Ls;
    for (int i = 10;i < 42;i++) { Ls.emplace_back(i); }

    std::vector<double> masses;
    std::vector<double> err_mass;
    int e1;
    for (int L : Ls) {
        e1 = Ne - 1;
        for (int e = 0; e < Ne;e++) {
            if (L == paramsj[e].data.L[1]) {
                e1 = e;
                printf("L=%d found   mass=%g\n", L, gjack[e1].jack[443][Njack - 1]);
            }
        }
        masses.emplace_back(gjack[e1].jack[443][Njack - 1]);// 443 ->GEVP_0_3_1_001_011_meffl0
        err_mass.emplace_back(error_jackboot(argv[1], Njack, gjack[e1].jack[443]));
    }

    // zeta.Init_Lmq_g(Ls, masses, err_mass  );
    // zeta.write("zeta_complex_g5_largeL.dat");
    zeta.read("zeta_complex_g5_largeL.dat");

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // printing files
    //////////////////////////////////////////////////////////////////////////////////////////////////

    printing_file_for_maxim_and_fernando_complex(argv, paramsj, gjack, myen);

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////////////////////////////////////////////////////////////////////////
    fit_type fit_info, fit_info_m0;

    struct fit_result  fit_m1, fit_m0;
    fit_info.Nvar = 14;                  fit_info_m0.Nvar = 14;
    fit_info.Npar = 2;                   fit_info_m0.Npar = 2;
    fit_info.N = 1;                      fit_info_m0.N = 1;
    fit_info.Njack = gjack[0].Njack;     fit_info_m0.Njack = gjack[0].Njack;
    fit_info.n_ext_P = 0;                fit_info_m0.n_ext_P = 0;
    fit_info.function = M_finite_volume; fit_info_m0.function = M_finite_volume;

    fit_info.n_ext_P = 3;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
    fit_info.ext_P[0] = gjack[0].jack[1];
    fit_info.ext_P[1] = gjack[0].jack[1];
    fit_info.ext_P[2] = gjack[0].jack[1];


    fit_m0 = fit_data(argv, paramsj, gjack, M0_finite_volume_lhs, fit_info, "M0_finite_vol", myen);
    jackall.add_fit(fit_m0);

    fit_m1 = fit_data(argv, paramsj, gjack, M1_finite_volume_lhs, fit_info, "M1_finite_vol", myen);

    // corr 75 =E1_1_px
    fit_result fit_m1_p1 = fit_data(argv, paramsj, gjack, M1_p_finite_volume_lhs<75>, fit_info, "M1_p1_finite_vol", myen);
    free_fit_result(fit_info, fit_m1_p1);
    printf("\n/////////////////////////////////     E2_0//////////////////\n");
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // E20
    //////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // start fitting

    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = gjack[0].Njack;// E1_0
    fit_info.n_ext_P = 3;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
    fit_info.ext_P[0] = fit_m0.P[0];
    fit_info.ext_P[1] = fit_m0.P[0];
    fit_info.ext_P[2] = fit_m0.P[0];	 //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
    fit_info.function = muDE_rhs;

    struct fit_result fit_a_00 = fit_data(argv, paramsj, gjack, muDE_00_lhs, fit_info, "a_00_luscher", myen);

    printf("\n/////////////////////////////////     k cot delta    //////////////////\n");
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // kcot
    //////////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.Npar = 2;
    fit_info.N = 3;
    fit_info.Nvar = 2;
    fit_info.Njack = jackall.en[0].Njack;
    fit_info.myen = myen;
    fit_info.precision_sum = 2;
    fit_info.verbosity = 0;
    fit_info.guess = { -0.141739, -2.89287 };
    fit_info.lambda = 0.001;
    fit_info.acc = 0.001;
    fit_info.h = 1e-3;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    int count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = jackall.en[e].jack[1][j];
                fit_info.x[1][count][j] = compute_k_m_g_new(n, myen[e], j, jackall, fit_info);

            }
            count++;
        }
    }
    fit_info.function = rhs_kcotd_m_new;
    fit_result kcotd_m = fit_all_data(argv, jackall, lhs_kcotd_m_new, fit_info, "kcotd_m");
    fit_info.band_range = { };
    print_fit_band(argv, jackall, fit_info, fit_info_m0, "kcotd_m", "k_m", kcotd_m, fit_m0, 1, myen.size() - 1, 0.1);

    fit_info.restore_default();


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////     k cot delta deltaE    //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Npar = 2;
    fit_info.N = 3;
    fit_info.Nvar = 2;
    fit_info.Njack = jackall.en[0].Njack;
    fit_info.myen = myen;
    fit_info.precision_sum = 2;
    fit_info.verbosity = 0;
    fit_info.guess = { -0.141739, -2.89287 };
    fit_info.lambda = 0.001;
    fit_info.acc = 0.001;
    fit_info.h = 1e-3;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = jackall.en[e].jack[1][j];
                fit_info.x[1][count][j] = compute_k_m_g_new(n, myen[e], j, jackall, fit_info);

            }
            count++;
        }
    }
    fit_info.function = rhs_kcotd_m_new;
    fit_result fit_kcotd_DeltaE = fit_all_data(argv, jackall, lhs_kcotd_m_deltaE_g_new, fit_info, "kcotd_m_deltaE");
    fit_info.band_range = { };
    // print_fit_band_phi4(argv, jackall, fit_info, fit_info_m0, "kcotd_m_deltaE", "L", deltaE2_m_QC2, fit_m0, 0, 0, 1);
    print_fit_band(argv, jackall, fit_info, fit_info_m0, "kcotd_m_deltaE", "k_m", fit_kcotd_DeltaE, fit_m0, 1, myen.size() - 1, 0.1);

    fit_info.restore_default();
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////     k cot delta / m_inf deltaE    //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Npar = 2;
    fit_info.N = 3;
    fit_info.Nvar = 3;
    fit_info.Njack = jackall.en[0].Njack;
    fit_info.myen = myen;
    fit_info.precision_sum = 2;
    fit_info.verbosity = 0;
    fit_info.guess = { -0.141739, -2.89287 };
    fit_info.lambda = 0.001;
    fit_info.acc = 0.001;
    fit_info.h = 1e-3;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = jackall.en[e].jack[1][j];
                fit_info.x[1][count][j] = compute_k_m_g_new(n, myen[e], j, jackall, fit_info);
                fit_info.x[2][count][j] = fit_m0.P[0][j];

            }
            count++;
        }
    }
    fit_info.function = rhs_kcotd_m_new;
    fit_result fit_kcotd_minf_DeltaE = fit_all_data(argv, jackall, lhs_kcotd_minf_deltaE_g_new, fit_info, "kcotd_minf_deltaE");
    fit_info.band_range = { };
    // print_fit_band_phi4(argv, jackall, fit_info, fit_info_m0, "kcotd_m_deltaE", "L", deltaE2_m_QC2, fit_m0, 0, 0, 1);
    print_fit_band(argv, jackall, fit_info, fit_info_m0, "kcotd_minf_deltaE", "k_m", fit_kcotd_minf_DeltaE, fit_m0, 1, myen.size() - 1, 0.1);

    fit_info.restore_default();
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////     k cot delta deltaE_infm    //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Npar = 2;
    fit_info.N = 3;
    fit_info.Nvar = 2;
    fit_info.Njack = jackall.en[0].Njack;
    fit_info.myen = myen;
    fit_info.precision_sum = 2;
    fit_info.verbosity = 0;
    fit_info.guess = { -0.141739, -2.89287 };
    fit_info.lambda = 0.001;
    fit_info.acc = 0.001;
    fit_info.h = 1e-3;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = fit_m0.P[0][j];
                fit_info.x[1][count][j] = compute_k_m_g_new(n, myen[e], j, jackall, fit_info);

            }
            count++;
        }
    }
    fit_info.function = rhs_kcotd_m_new;
    fit_result fit_kcotd_DeltaE_infm = fit_all_data(argv, jackall, lhs_kcotd_m_deltaE_g_new, fit_info, "kcotd_m_deltaE_infm");
    fit_info.band_range = { };
    // print_fit_band_phi4(argv, jackall, fit_info, fit_info_m0, "kcotd_m_deltaE", "L", deltaE2_m_QC2, fit_m0, 0, 0, 1);
    print_fit_band(argv, jackall, fit_info, fit_info_m0, "kcotd_m_deltaE_infm", "k_m", fit_kcotd_DeltaE_infm, fit_m0, 1, myen.size() - 1, 0.1);

    fit_info.restore_default();





    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////   fit  deltaE2_m_QC3  //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Npar = 2;
    fit_info.N = 3;
    fit_info.Nvar = 3;
    fit_info.Njack = jackall.en[0].Njack;
    fit_info.myen = myen;
    fit_info.precision_sum = 2;
    fit_info.verbosity = 0;
    fit_info.guess = { -0.141739, -2.89287 };
    fit_info.lambda = 0.001;
    fit_info.acc = 0.001;
    fit_info.h = 1e-3;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = paramsj[e].data.L[1];
                fit_info.x[1][count][j] = jackall.en[fit_info.myen[e]].jack[1][j];
                fit_info.x[2][count][j] = fit_m0.P[0][j];
            }
            count++;
        }
    }
    fit_info.function = rhs_deltaE2_m_QC2;
    fit_result deltaE2_m_QC2 = fit_all_data(argv, jackall, lhs_deltaE2_m_latt_QC2, fit_info, "deltaE2_m_QC2");
    fit_info.band_range = { 19, 26 };
    print_fit_band_phi4(argv, jackall, fit_info, fit_info_m0, "deltaE2_m_QC2", "L", deltaE2_m_QC2, fit_m0, 0, 0, 1);

    jackall.add_fit(deltaE2_m_QC2);
    fit_info.restore_default();



    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////   fit  E3 poly  //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////

    printf("//////////////////// poly fit E3   ////////////////////////////////////\n");

    fit_type fit_info_E3_poly;
    fit_info_E3_poly.N = 2;
    fit_info_E3_poly.Nvar = 1;
    fit_info_E3_poly.function = rhs_poly_order_E3_m<4>;
    fit_info_E3_poly.myen = myen;
    fit_info_E3_poly.Njack = jackall.en[0].Njack;
    fit_info_E3_poly.Npar = fit_info_E3_poly.N * 4;
    fit_info_E3_poly.verbosity = 0;
    fit_info_E3_poly.x = double_malloc_3(fit_info_E3_poly.Nvar, fit_info_E3_poly.myen.size() * fit_info_E3_poly.N, fit_info_E3_poly.Njack);

    printf("size=%ld\n", fit_info_E3_poly.guess.size());
    count = 0;
    for (int n = 0;n < fit_info_E3_poly.N;n++) {
        for (int e = 0;e < fit_info_E3_poly.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info_E3_poly.x[0][count][j] = paramsj[e].data.L[1] * jackall.fits[0].P[0][j];
            }
            count++;
        }
    }

    fit_result fit_QC3_poly = fit_all_data(argv, jackall, lhs_E3_m_new, fit_info_E3_poly, "fit_QC3_poly");

    fit_info_E3_poly.band_range = { 5.5, 7.5 };
    print_fit_band(argv, jackall, fit_info_E3_poly, fit_info_m0, "fit_QC3_poly", "L", fit_QC3_poly, fit_m0, 0, 0, 0.1);

    //fit_info.restore_default();

#ifdef PYTHON
    //// we need python
    wchar_t* program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }
    Py_SetProgramName(program);  /* optional but recommended */

    Py_Initialize();


    printf("////////////////////  kiso pole fit   ////////////////////////////////////\n");
    init_python_detQC();
    init_python_detQC_kcot_kiso("kcot_2par", "kiso_pole", "find_2sol_QC3_pole");

    fit_info.Npar = 3;
    fit_info.N = 2;
    fit_info.Njack = jackall.en[0].Njack;

    fit_info.function = rhs_E3_m_QC3_pole_new;
    fit_info.n_ext_P = 0;
    fit_info.Nvar = 5;
    fit_info.myen = myen;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);

    int scount = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][scount][j] = jackall.en[e].header.L * jackall.fits[0].P[0][j];
                fit_info.x[1][scount][j] = jackall.fits[1].P[0][j];
                fit_info.x[2][scount][j] = jackall.fits[1].P[1][j];
                if (n == 0) {
                    fit_info.x[3][scount][j] = jackall.en[e].jack[354][j] / jackall.fits[0].P[0][j];
                    fit_info.x[4][scount][j] = jackall.en[e].error_jack(354) / jackall.fits[0].P[0][j];
                }
                else if (n == 1) {
                    fit_info.x[3][scount][j] = jackall.en[e].jack[355][j] / jackall.fits[0].P[0][j];
                    fit_info.x[4][scount][j] = jackall.en[e].error_jack(355) / jackall.fits[0].P[0][j];
                }
            }
            scount++;
        }
    }

    // fit_info.lambda = 0.001;
    fit_info.acc = 0.01;
    fit_info.h = 1e-2;
    // fit_info.Prange={100,100, 100};
    fit_info.devorder = -2;
    fit_info.verbosity = 3;
    fit_info.repeat_start = 1;
    // fit_info.guess = { -0.142262, -2.96471 };
    //fit_info.guess = { -150.299, 9.72572 };
    // fit_info.guess = { 318.061, 9.64639 };
    fit_info.guess = { 60.0257, 9.1358, -2859.6 };
    fit_info.guess = { 58.371, 9.13589, -2803.26 };
    // fit_info.guess = {-8.48047e+06, 10.66685e+06,-425.446 };
    fit_info.mean_only = true;
    fit_info.precision_sum = 2;

    // fit_info.noderiv=true;
    // fit_info.lambda=0.01;
    // fit_info.Prange={1,0.001, 10};


    mysprintf(namefile, NAMESIZE, "QC3_%dpar_pole", fit_info.Npar);
    // struct fit_result fit_QC3_2par = fit_all_data(argv, jackall, lhs_E3orE1_m_complex_new, fit_info, namefile);
    // fit_info.band_range = { 5.65,7.2 };
    // print_fit_band(argv, jackall, fit_info, fit_info, namefile, "L_m", fit_QC3_2par, fit_QC3_2par, 0, 0, 0.15);

    fit_info.restore_default();

    printf("//////////////////// kcot and kiso pole fit   ////////////////////////////////////\n");
    // init_python_detQC();
    init_python_detQC_kcot_kiso("kcot_2par", "kiso_pole", "find_2sol_QC3_pole");

    fit_info.Npar = 3 + 2;
    fit_info.N = 2 + 1;
    fit_info.Njack = jackall.en[0].Njack;

    fit_info.function = rhs_E3_m_QC3_pole_E2_QC2;
    fit_info.n_ext_P = 0;
    fit_info.Nvar = 2;
    fit_info.myen = myen;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);

    scount = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][scount][j] = jackall.en[e].header.L * jackall.en[e].jack[1][j];
                fit_info.x[1][scount][j] = jackall.en[e].jack[1][j];
            }
            scount++;
        }
    }

    // fit_info.lambda = 0.001;
    fit_info.acc = 0.1;
    fit_info.h = 1e-2;
    // //fit_info.Prange={1000,10000};
    fit_info.devorder = -2;
    fit_info.verbosity = 3;
    fit_info.repeat_start = 1;

    // fit_info.guess = { 74.37, 9.129, -710, -0.160, -32 };
    // fit_info.guess = {74.0805, 9.12889, -706.844, -0.159792, -31.5321 };
    fit_info.guess = { 73.7538, 9.12887, -706.766, -0.159647, -31.1049 };
    // fit_info.mean_only = true;
    fit_info.precision_sum = 2;

    // fit_info.noderiv=true;
    //fit_info.Prange={1,0.001, 10, 0.01, 1};


    mysprintf(namefile, NAMESIZE, "kcot_1lev_and_kiso_pole_3par", fit_info.Npar);
    // struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_E3_E1_E2_m_complex_new, fit_info, namefile);
    // fit_info.band_range = { 5.65,7.2 };
    // print_fit_band(argv, jackall, fit_info, fit_info, namefile, "L_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, myen.size()-1, 0.15);
    fit_info.restore_default();

    printf("//////////////////// kcot 1par and kiso pole fit   ////////////////////////////////////\n");
    // init_python_detQC();
    init_python_detQC_kcot_kiso("kcot_1par", "kiso_pole", "find_2sol_QC3_pole");

    fit_info.Npar = 3 + 1;
    fit_info.N = 2 + 1;
    fit_info.Njack = jackall.en[0].Njack;

    fit_info.function = rhs_E3_m_QC3_pole_E2_QC2_1par;
    fit_info.n_ext_P = 0;
    fit_info.Nvar = 2;
    fit_info.myen = myen;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);

    scount = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][scount][j] = jackall.en[e].header.L * jackall.en[e].jack[1][j];
                fit_info.x[1][scount][j] = jackall.en[e].jack[1][j];
            }
            scount++;
        }
    }

    // fit_info.lambda = 0.001;
    fit_info.acc = 0.1;
    fit_info.h = 1e-3;
    // //fit_info.Prange={1000,10000};
    fit_info.devorder = -2;
    fit_info.verbosity = 3;
    fit_info.repeat_start = 1;
    // fit_info.guess = {22.2822 ,        9.19729   ,      1317.02  ,       -0.155816 };
    fit_info.guess = { 23.6528,    9.14173, 2622.92, -0.156144 };
    // fit_info.mean_only = true;
    fit_info.precision_sum = 2;

    //  fit_info.noderiv=true;
    //  fit_info.Prange={1,0.0004, 10,  0.001};


    mysprintf(namefile, NAMESIZE, "kcot_1par_1lev_and_kiso_pole_3par", fit_info.Npar);
    // struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_E3_E1_E2_m_complex_new, fit_info, namefile);
    // fit_info.band_range = { 5.0,7.0 };
    // print_fit_band(argv, jackall, fit_info, fit_info, namefile, "L_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, myen.size()-1, 0.15);
    fit_info.restore_default();


    printf("//////////////////// kcot 1par and kiso pole 2par fit   ////////////////////////////////////\n");
    // init_python_detQC();
    init_python_detQC_kcot_kiso("kcot_1par", "kiso_pole", "find_2sol_QC3_pole");

    fit_info.Npar = 2 + 1;
    fit_info.N = 2 + 1;
    fit_info.Njack = jackall.en[0].Njack;

    fit_info.function = rhs_E3_m_QC3_pole_E2_QC2_1par;
    fit_info.n_ext_P = 0;
    fit_info.Nvar = 2;
    fit_info.myen = myen;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);

    scount = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][scount][j] = jackall.en[e].header.L * jackall.en[e].jack[1][j];
                fit_info.x[1][scount][j] = jackall.en[e].jack[1][j];
            }
            scount++;
        }
    }

    // fit_info.lambda = 0.001;
    fit_info.acc = 0.01;
    // fit_info.h = 1e-2;
    // //fit_info.Prange={1000,10000};
    fit_info.devorder = 2;
    fit_info.verbosity = 3;
    fit_info.repeat_start = 1;
    // fit_info.guess = { 177.049 , 9.12102, -0.14583 };
    fit_info.guess = { 15.688,     9.14136, -0.150566 };
    fit_info.h = {0.01,0.00001, 0.0001};
    fit_info.mean_only = false;
    fit_info.precision_sum = 2;
    fit_info.covariancey = false;
    // fit_info.noderiv = true;
    // fit_info.Prange={0.01, 0.00001, 0.0001};


    mysprintf(namefile, NAMESIZE, "kcot_1par_1lev_and_kiso_pole_2par", fit_info.Npar);
    struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_E3_E1_E2_m_complex_new, fit_info, namefile);
    fit_info.band_range = { 5, 7 };
    print_fit_band(argv, jackall, fit_info, fit_info, namefile, "L_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, myen.size() - 1, 0.15);
    fit_info.restore_default();



    printf("//////////////////// kcot 1par and kiso no coupling fit   ////////////////////////////////////\n");
    // init_python_detQC();
    init_python_detQC_kcot_kiso("kcot_1par", "kiso_1par", "find_2sol_QC3_pole");

    fit_info.Npar = 2 + 1;
    fit_info.N = 2 + 1;
    fit_info.Njack = jackall.en[0].Njack;

    fit_info.function = rhs_E3_m_QC3_const_E2_QC2_1par;
    fit_info.n_ext_P = 0;
    fit_info.Nvar = 2;
    fit_info.myen = myen;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);

    scount = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][scount][j] = jackall.en[e].header.L * jackall.en[e].jack[1][j];
                fit_info.x[1][scount][j] = jackall.en[e].jack[1][j];
            }
            scount++;
        }
    }

    // fit_info.lambda = 0.001;
    fit_info.acc = 0.05;
    fit_info.h = 1e-3;
    // //fit_info.Prange={1000,10000};
    fit_info.devorder = -2;
    fit_info.verbosity = 3;
    fit_info.repeat_start = 1;
    //  fit_info.guess = { 85.28, 9.1294,  -0.155809 };
    fit_info.guess = { 1239.9, 3.03258, -0.157269 };
    // fit_info.mean_only = true;
    fit_info.precision_sum = 2;

    // fit_info.noderiv=true;
    // fit_info.Prange={10,0.001, 0.01};


    mysprintf(namefile, NAMESIZE, "kcot_1par_1lev_and_kiso_const", fit_info.Npar);
    // struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_E3_E1_E2_m_complex_new, fit_info, namefile);
    // fit_info.band_range = { 5.0,7.0  };
    // print_fit_band(argv, jackall, fit_info, fit_info, namefile, "L_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, myen.size() - 1, 0.15);
    fit_info.restore_default();


    printf("//////////////////// kcot and kiso no coupling  fit   ////////////////////////////////////\n");
    // init_python_detQC();
    init_python_detQC_kcot_kiso("kcot_2par", "kiso_1par", "find_2sol_QC3_pole");

    fit_info.Npar = 2 + 2;
    fit_info.N = 2 + 1;
    fit_info.Njack = jackall.en[0].Njack;

    fit_info.function = rhs_E3_m_QC3_const_E2_QC2;
    fit_info.n_ext_P = 0;
    fit_info.Nvar = 2;
    fit_info.myen = myen;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);

    scount = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][scount][j] = jackall.en[e].header.L * jackall.en[e].jack[1][j];
                fit_info.x[1][scount][j] = jackall.en[e].jack[1][j];
            }
            scount++;
        }
    }

    // fit_info.lambda = 0.001;
    fit_info.acc = 0.6;
    fit_info.h = 1e-2;
    // //fit_info.Prange={1000,10000};
    fit_info.devorder = -2;
    fit_info.verbosity = 3;
    fit_info.repeat_start = 1;

    fit_info.guess = { -744.738,  9.12657, -0.159187, -23.5928 };

    // fit_info.mean_only = true;
    fit_info.precision_sum = 2;

    // fit_info.noderiv=true;
    //fit_info.Prange={1,0.001, 10, 0.01, 1};


    mysprintf(namefile, NAMESIZE, "kcot_1lev_and_kiso_const_3par", fit_info.Npar);
    // struct fit_result kcot_1lev_and_kiso_pole_3par = fit_all_data(argv, jackall, lhs_E3_E1_E2_m_complex_new, fit_info, namefile);
    // fit_info.band_range = { 5.65,7.2 };
    // print_fit_band(argv, jackall, fit_info, fit_info, namefile, "L_m", kcot_1lev_and_kiso_pole_3par, kcot_1lev_and_kiso_pole_3par, 0, myen.size()-1, 0.05);
    fit_info.restore_default();





    ///// close python
    python_detQC_write_database();
    python_detQC_free();
    if (Py_FinalizeEx() < 0) {
        exit(120);
    }
    PyMem_RawFree(program);

#endif
    return 0;
}
