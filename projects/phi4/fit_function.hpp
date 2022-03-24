#ifndef fit_function_H
#define fit_function_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <random>

#include "linear_fit.hpp"
#include "non_linear_fit.hpp"
#include "mutils.hpp"
#include "resampling.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "global.hpp"
#include "eigensystem.hpp"
#include "gamma_analysis.hpp"
#include "tower.hpp"

#include "header_phi4.hpp"
#include "mass_phi4.hpp"
#include "zeta_interpolation.hpp"
#include "mutils.hpp"
#include "linear_fit.hpp"
#ifdef PYTHON
#include <QC3_interface.hpp>
#endif
//using namespace std;

void printing_file_for_maxim_and_fernando(char** argv, std::vector<cluster::IO_params> paramsj, std::vector<data_all> gjack, std::vector<int> myen) {
    int Njack = gjack[0].Njack;
    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, "%s", argv[1]);
    for (int j = 0; j < Njack; j++) {
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "data_j%d.txt", j);
        FILE* f = open_file(name, "w+");
        fprintf(f, "L	T	g          m02       m12        l0        E1         	E2         	E2p1      	 E2A1        	E3_0     	     	E3_1  \n");
        for (int e : myen) {
            fprintf(f, "%-4d%-4d%-10g%-10g%-10g%-10g%-14.10g%-14.10g%-14.10g%-14.10g%-14.10g%-14.10g\n",
                paramsj[e].data.L[1], paramsj[e].data.L[0], paramsj[e].data.gC,
                paramsj[e].data.msq0, paramsj[e].data.msq1, paramsj[e].data.lambdaC0,
                gjack[e].jack[443][Njack - 1],
                gjack[e].jack[594][Njack - 1],
                gjack[e].jack[597][Njack - 1],
                gjack[e].jack[601][Njack - 1],
                gjack[e].jack[444][Njack - 1],
                gjack[e].jack[445][Njack - 1]);
        }
        fclose(f);
    }
    {
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "data.txt");
        FILE* f = open_file(name, "w+");
        fprintf(f, "L	T	g          m02       m12        l0        E1    	  E1err     	E2      	E2err       	E2p1     	E2p1err 	 E2A1       	E2A1err     	E3_0     	E3_0err     	E3_1    	E3_1err  \n");
        for (int e : myen) {
            fprintf(f, "%-4d%-4d%-10g%-10g%-10g%-10g%-14g%-18.10g%-14.10g%-18.10g%-14.10g%-18.10g%-14.10g%-18.10g%-14.10g%-18.10g%-14.10g%-18.10g\n",
                paramsj[e].data.L[1], paramsj[e].data.L[0], paramsj[e].data.gC,
                paramsj[e].data.msq0, paramsj[e].data.msq1, paramsj[e].data.lambdaC0,
                gjack[e].jack[443][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[443]),
                gjack[e].jack[594][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[594]),
                gjack[e].jack[597][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[597]),
                gjack[e].jack[601][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[601]),
                gjack[e].jack[444][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[444]),
                gjack[e].jack[445][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[445]));
        }
        fclose(f);
    }
    {
        for (int e : myen) {
            char name[NAMESIZE];
            mysprintf(name, NAMESIZE, "covariance_spectrum_L%d.txt", paramsj[e].data.L[1]);
            FILE* f = open_file(name, "w+");
            double** to_cov = double_malloc_2(6, Njack);
            for (int j = 0;j < Njack;j++) {
                to_cov[0][j] = gjack[e].jack[443][j];
                to_cov[1][j] = gjack[e].jack[594][j];
                to_cov[2][j] = gjack[e].jack[597][j];
                to_cov[3][j] = gjack[e].jack[601][j];
                to_cov[4][j] = gjack[e].jack[444][j];
                to_cov[5][j] = gjack[e].jack[445][j];
            }
            double** cov = covariance(resampling, 6, Njack, to_cov);
            for (int i = 0;i < 6;i++) {
                for (int j = 0;j < 6;j++) {
                    fprintf(f, "%-18.10g", cov[i][j]);
                }
                fprintf(f, "\n");
            }
            fclose(f);
            free(to_cov);
            free(cov);
        }
    }
}

void printing_file_for_maxim_and_fernando_complex(char** argv, std::vector<cluster::IO_params> paramsj, std::vector<data_all> gjack, std::vector<int> myen) {
    int Njack = gjack[0].Njack;
    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, "%s", argv[1]);
    for (int j = 0; j < Njack; j++) {
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "data_j%d.txt", j);
        FILE* f = open_file(name, "w+");
        fprintf(f, "L	T	g          m02       m12        l0        E1         	E2         	E2p1      	 E2A1        	E3_0     	     	E3_1  \n");
        for (int e : myen) {
            fprintf(f, "%-4d%-4d%-10g%-10g%-10g%-10g%-14.10g%-14.10g%-14.10g%-14.10g%-14.10g%-14.10g\n",
                paramsj[e].data.L[1], paramsj[e].data.L[0], paramsj[e].data.gC,
                paramsj[e].data.msq0, paramsj[e].data.msq1, paramsj[e].data.lambdaC0,
                gjack[e].jack[1][Njack - 1],
                gjack[e].jack[4][Njack - 1],
                gjack[e].jack[100][Njack - 1],
                gjack[e].jack[80][Njack - 1],
                gjack[e].jack[355][Njack - 1],
                gjack[e].jack[354][Njack - 1]);
        }
        fclose(f);
    }
    {
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "data.txt");
        FILE* f = open_file(name, "w+");
        fprintf(f, "L	T	g          m02       m12        l0        E1    	  E1err     	E2      	E2err       	E2p1     	E2p1err 	 E2A1       	E2A1err     	E3_0     	E3_0err     	E3_1    	E3_1err  \n");
        for (int e : myen) {
            fprintf(f, "%-4d%-4d%-10g%-10g%-10g%-10g%-14g%-18.10g%-14.10g%-18.10g%-14.10g%-18.10g%-14.10g%-18.10g%-14.10g%-18.10g%-14.10g%-18.10g\n",
                paramsj[e].data.L[1], paramsj[e].data.L[0], paramsj[e].data.gC,
                paramsj[e].data.msq0, paramsj[e].data.msq1, paramsj[e].data.lambdaC0,
                gjack[e].jack[1][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[1]),
                gjack[e].jack[4][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[4]),
                gjack[e].jack[100][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[100]),
                gjack[e].jack[80][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[80]),
                gjack[e].jack[355][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[355]),
                gjack[e].jack[354][Njack - 1], error_jackboot(resampling, Njack, gjack[e].jack[354]));
        }
        fclose(f);
    }
    {
        for (int e : myen) {
            char name[NAMESIZE];
            mysprintf(name, NAMESIZE, "covariance_spectrum_L%d.txt", paramsj[e].data.L[1]);
            FILE* f = open_file(name, "w+");
            double** to_cov = double_malloc_2(6, Njack);
            for (int j = 0;j < Njack;j++) {
                to_cov[0][j] = gjack[e].jack[1][j];
                to_cov[1][j] = gjack[e].jack[4][j];
                to_cov[2][j] = gjack[e].jack[100][j];
                to_cov[3][j] = gjack[e].jack[80][j];
                to_cov[4][j] = gjack[e].jack[355][j];
                to_cov[5][j] = gjack[e].jack[354][j];
            }
            double** cov = covariance(resampling, 6, Njack, to_cov);
            for (int i = 0;i < 6;i++) {
                for (int j = 0;j < 6;j++) {
                    fprintf(f, "%-18.10g", cov[i][j]);
                }
                fprintf(f, "\n");
            }
            fclose(f);
            free(to_cov);
            free(cov);
        }
    }
}


void init_dvec(int n, int* dvec, int* dvec1, int* dvec2, int* dmax1, int* dmax2) {

    if (n == 0) {//E2_0
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
        dvec1[0] = 0; dvec1[1] = 0; dvec1[2] = 0;
        dvec2[0] = 0; dvec2[1] = 0; dvec2[2] = 0;
        dmax1[0] = 1; dmax1[1] = 0; dmax1[2] = 0;
        dmax2[0] = -1; dmax2[1] = 0; dmax2[2] = 0;
    }
    else if (n == 1) {//E2_0_p1
        dvec[0] = 1; dvec[1] = 0; dvec[2] = 0;
        dvec1[0] = 1; dvec1[1] = 0; dvec1[2] = 0;
        dvec2[0] = 0; dvec2[1] = 0; dvec2[2] = 0;
        dmax1[0] = 1; dmax1[1] = 1; dmax1[2] = 0;
        dmax2[0] = 0; dmax2[1] = -1; dmax2[2] = 0;
    }
    else if (n == 2) {//E2_0_p11
        dvec[0] = 1; dvec[1] = 1; dvec[2] = 0;
        dvec1[0] = 1; dvec1[1] = 1; dvec1[2] = 0;
        dvec2[0] = 0; dvec2[1] = 0; dvec2[2] = 0;
        dmax1[0] = 1; dmax1[1] = 0; dmax1[2] = 0;
        dmax2[0] = 0; dmax2[1] = 1; dmax2[2] = 0;
    }
    else if (n == 4) {//E2_0_p111
        dvec[0] = 1; dvec[1] = 1; dvec[2] = 1;
        dvec1[0] = 1; dvec1[1] = 1; dvec1[2] = 1;
        dvec2[0] = 0; dvec2[1] = 0; dvec2[2] = 0;
        dmax1[0] = 1; dmax1[1] = 1; dmax1[2] = 0;
        dmax2[0] = 0; dmax2[1] = 0; dmax2[2] = 1;
    }
    else if (n == 3) {//E2_0_A1
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
        dvec1[0] = 1; dvec1[1] = 0; dvec1[2] = 0;
        dvec2[0] = -1; dvec2[1] = 0; dvec2[2] = 0;
        dmax1[0] = 1; dmax1[1] = 0; dmax1[2] = 1;
        dmax2[0] = -1; dmax2[1] = 0; dmax2[2] = -1;
    }
    else if (n == 5) {//E2_0_p111
        dvec[0] = 1; dvec[1] = 1; dvec[2] = 1;
        dvec1[0] = 1; dvec1[1] = 1; dvec1[2] = 1;
        dvec2[0] = 0; dvec2[1] = 0; dvec2[2] = 0;
        dmax1[0] = 1; dmax1[1] = 1; dmax1[2] = 1;
        dmax2[0] = 0; dmax2[1] = 0; dmax2[2] = 0;
    }
    else {
        printf("init_dvec n=%d not implemented\n", n); exit(1);
    }
}


void init_dvec_E2_g(int n, int* dvec, int* dvec1, int* dvec2, int* dmax1, int* dmax2) {

    if (n == 0) {//E2_0
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
        dvec1[0] = 0; dvec1[1] = 0; dvec1[2] = 0;
        dvec2[0] = 0; dvec2[1] = 0; dvec2[2] = 0;
        dmax1[0] = 1; dmax1[1] = 0; dmax1[2] = 0;
        dmax2[0] = -1; dmax2[1] = 0; dmax2[2] = 0;
    }
    else if (n == 1) {//E2_0_p1
        dvec[0] = 1; dvec[1] = 0; dvec[2] = 0;
        dvec1[0] = 1; dvec1[1] = 0; dvec1[2] = 0;
        dvec2[0] = 0; dvec2[1] = 0; dvec2[2] = 0;
        dmax1[0] = 1; dmax1[1] = 1; dmax1[2] = 0;
        dmax2[0] = 0; dmax2[1] = -1; dmax2[2] = 0;
    }
    else if (n == 2) {//E2_0_A1
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
        dvec1[0] = 1; dvec1[1] = 0; dvec1[2] = 0;
        dvec2[0] = -1; dvec2[1] = 0; dvec2[2] = 0;
        dmax1[0] = 1; dmax1[1] = 0; dmax1[2] = 1;
        dmax2[0] = -1; dmax2[1] = 0; dmax2[2] = -1;
    }
    else {
        printf("%s n=%d not implemented\n", __func__, n); exit(1);
    }
}

void init_dvec_g_E3(int n, int* dvec, int* dvec1, int* dvec2, int* dmax1, int* dmax2) {

    if (n == 0 && n == 1) {//E3_0
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
        dvec1[0] = 0; dvec1[1] = 0; dvec1[2] = 0;
        dvec2[0] = 0; dvec2[1] = 0; dvec2[2] = 0;
        dmax1[0] = 1; dmax1[1] = 0; dmax1[2] = 0;
        dmax2[0] = -1; dmax2[1] = 0; dmax2[2] = 0;
    }

    else {
        printf("init_dvec n=%d not implemented\n", n); exit(1);
    }
}

void init_dvec_QC3_pole(int n, int* dvec) {
    if (n == 0) {//E2_0
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
    }
    else if (n == 1) {//E2_0_p1
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
    }
    else if (n == 2) {//E2_0_p11
        dvec[0] = 1; dvec[1] = 0; dvec[2] = 0;
    }
    else if (n == 3) {//E2_0_p1
        dvec[0] = 1; dvec[1] = 0; dvec[2] = 0;
    }
    else if (n == 4) {//E2_0_p1
        dvec[0] = 1; dvec[1] = 1; dvec[2] = 0;
    }
    else if (n == 5) {//E2_0_p1
        dvec[0] = 1; dvec[1] = 1; dvec[2] = 0;
    }

    else {
        printf("init_dvec n=%d not implemented\n", n); exit(1);
    }
}

double get_E2_n(int n, int e, int j, std::vector<data_all> gjack) {
    double E2;
    if (n == 0) {//E2_0
        E2 = gjack[e].jack[4][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 1) {//E2_0_p1
        E2 = gjack[e].jack[100][j];
        //         dvec[0]=1; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 2) {//E2_0_p11
        E2 = gjack[e].jack[102][j];
        //         dvec[0]=1; dvec[1]=1; dvec[2]=0;
    }
    else if (n == 4) {//E2_0_p111
        E2 = gjack[e].jack[104][j];
        //         dvec[0]=1; dvec[1]=1; dvec[2]=1;
    }
    else if (n == 3) {//E2_0_A1
        E2 = gjack[e].jack[80][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else { printf("%s n=%d not implemented\n", __func__, n);  exit(1); }
    return E2;

}

double get_E2_g_n(int n, int e, int j, vector<data_all> gjack) {
    double E2;
    if (n == 0) {//E2_0
        // if(gjack[e].jack[594][j]!=0) 
        //     E2=gjack[e].jack[594][j];
        // else 
        E2 = gjack[e].jack[4][j];

    }
    else if (n == 1) {//E2_p1
        // if(gjack[e].jack[597][j]!=0) 
        //     E2=gjack[e].jack[597][j];
        // else
        E2 = gjack[e].jack[100][j];
    }
    else if (n == 2) {//E2_A1
        // if(gjack[e].jack[601][j]!=0) 
        //     E2=gjack[e].jack[601][j];
        // else
        E2 = gjack[e].jack[80][j];
    }

    else { printf("%s n=%d not implemented\n", __func__, n);  exit(1); }
    return E2;

}
///////////////////////////////////////////////////////////////////////////////////////////
///////kcotd
//////////////////////////////////////////////////////////////////////////////////////////
inline double kcotd(double E2, double mass, int* dvec, int L) {
    double E2_CM = energy_CM(E2, dvec, L);
    double k = sqrt(E2_CM * E2_CM / 4. - mass * mass);

    double q = k * L / (2. * pi_greco);
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];
    dzeta_function(z, q * q, 0, 0, dvec, gamma, A, 1.e-3, 1.e+6, 5);

    return  z[0] * 2 * pi_greco / (pow(pi_greco, 3. / 2.) * gamma * L);
}

double lhs_kcotd(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double r;
    if (n == 0) {//E2_0
        int dvec[3] = { 0,0,0 };
        r = kcotd(gjack[e].jack[4][j], gjack[e].jack[1][j], dvec, params[e].data.L[1]);
    }
    else if (n == 1) {//E2_0_p1
        int dvec[3] = { 1,0,0 };
        r = kcotd(gjack[e].jack[100][j], gjack[e].jack[1][j], dvec, params[e].data.L[1]);
    }
    else if (n == 2) {//E2_0_p1
        int dvec[3] = { 1,1,0 };
        r = kcotd(gjack[e].jack[102][j], gjack[e].jack[1][j], dvec, params[e].data.L[1]);
    }
    else if (n == 4) {//E2_0_p1
        int dvec[3] = { 1,1,1 };
        r = kcotd(gjack[e].jack[104][j], gjack[e].jack[1][j], dvec, params[e].data.L[1]);
    }
    else if (n == 3) {//E2_0_A1
        int dvec[3] = { 0,0,0 };
        r = kcotd(gjack[e].jack[80][j], gjack[e].jack[1][j], dvec, params[e].data.L[1]);
    }
    else { r = 1; exit(1); }

    return r;
}



double lhs_kcotd_g(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double r;
    if (n == 0) {//E2_0
        int dvec[3] = { 0,0,0 };
        if (gjack[e].jack[594][j] != 0) {
            r = kcotd(gjack[e].jack[594][j], gjack[e].jack[443][j], dvec, params[e].data.L[1]);
        }
        else {
            r = kcotd(gjack[e].jack[4][j], gjack[e].jack[443][j], dvec, params[e].data.L[1]);
        }
    }
    else if (n == 1) {
        int dvec[3] = { 1,0,0 };
        r = kcotd(gjack[e].jack[601][j], gjack[e].jack[443][j], dvec, params[e].data.L[1]);
    }
    else if (n == 2) {
        int dvec[3] = { 0,0,0 };
        r = kcotd(gjack[e].jack[595][j], gjack[e].jack[443][j], dvec, params[e].data.L[1]);
    }
    else { printf("%s", __func__); exit(1); }

    return r;
}

double lhs_kcotd_m_g(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double r;

    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2 = get_E2_g_n(n, e, j, gjack);
    r = kcotd(E2, gjack[e].jack[443][j], dvec, params[e].data.L[1]);
    return r / gjack[e].jack[443][j];
}


double lhs_kcotd_m_deltaE_g(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double r;

    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2 = get_E2_g_n(n, e, j, gjack);
    double mass = gjack[e].jack[1][j];


    double hatp2 = 4. * sin(dvec1[0] * pi_greco / params[e].data.L[1]) * sin(dvec1[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec1[1] * pi_greco / params[e].data.L[2]) * sin(dvec1[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec1[2] * pi_greco / params[e].data.L[3]) * sin(dvec1[2] * pi_greco / params[e].data.L[3]);

    double E2fL = acosh(cosh(mass) + 0.5 * (hatp2));
    hatp2 = 4. * sin(dvec2[0] * pi_greco / params[e].data.L[1]) * sin(dvec2[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec2[1] * pi_greco / params[e].data.L[2]) * sin(dvec2[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec2[2] * pi_greco / params[e].data.L[3]) * sin(dvec2[2] * pi_greco / params[e].data.L[3]);
    E2fL += acosh(cosh(mass) + 0.5 * (+hatp2));

    double L = params[e].data.L[1];
    double Ef1 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double Ef2 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double E2f = Ef1 + Ef2;

    E2 = E2 - E2fL + E2f;

    r = kcotd(E2, mass, dvec, params[e].data.L[1]);
    double mass_inf = fit_info.ext_P[2][j];
    return r / mass_inf;
}


double lhs_kcotd_Elatt(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2;
    double mass = gjack[e].jack[1][j];

    E2 = get_E2_n(n, e, j, gjack);
    //     double E20=gjack[e].jack[4][j];
    double hatp2 = 4. * sin(dvec1[0] * pi_greco / params[e].data.L[1]) * sin(dvec1[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec1[1] * pi_greco / params[e].data.L[2]) * sin(dvec1[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec1[2] * pi_greco / params[e].data.L[3]) * sin(dvec1[2] * pi_greco / params[e].data.L[3]);
    //     double E20=gjack[e].jack[4][j];
    double E2fL = acosh(cosh(mass) + 0.5 * (hatp2));
    hatp2 = 4. * sin(dvec2[0] * pi_greco / params[e].data.L[1]) * sin(dvec2[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec2[1] * pi_greco / params[e].data.L[2]) * sin(dvec2[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec2[2] * pi_greco / params[e].data.L[3]) * sin(dvec2[2] * pi_greco / params[e].data.L[3]);
    E2fL += acosh(cosh(mass) + 0.5 * (+hatp2));

    double L = params[e].data.L[1];
    double Ef1 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double Ef2 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double E2f = Ef1 + Ef2;

    E2 = E2 - E2fL + E2f;
    //     printf("%g  %g  %g \n",E2,E2fL,E2f);
    double E2_CM = energy_CM(E2, dvec, L);
    double k = sqrt(E2_CM * E2_CM / 4. - mass * mass);

    double q = k * L / (2. * pi_greco);
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];
    dzeta_function(z, q * q, 0, 0, dvec, gamma, A, 1.e-3, 1.e6, 3);

    return  z[0] * 2 * pi_greco / (pow(pi_greco, 3. / 2.) * gamma * L);


}


double lhs_delta(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2;
    double mass = gjack[e].jack[1][j];

    E2 = get_E2_n(n, e, j, gjack);

    double L = params[e].data.L[1];
    //     printf("%g  %g  %g \n",E2,E2fL,E2f);
    double E2_CM = energy_CM(E2, dvec, L);
    double k = sqrt(E2_CM * E2_CM / 4. - mass * mass);

    double q = k * L / (2. * pi_greco);
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];
    dzeta_function(z, q * q, 0, 0, dvec, gamma, A, 1.e-3, 1.e6, 3);

    return  atan((pow(pi_greco, 3. / 2.) * gamma * L * k) / z[0] * 2 * pi_greco);


}

double lhs_delta_g(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2;
    double mass = gjack[e].jack[443][j];

    E2 = get_E2_g_n(n, e, j, gjack);

    double L = params[e].data.L[1];
    //     printf("%g  %g  %g \n",E2,E2fL,E2f);
    double E2_CM = energy_CM(E2, dvec, L);
    double k = sqrt(E2_CM * E2_CM / 4. - mass * mass);

    double q = k * L / (2. * pi_greco);
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];
    dzeta_function(z, q * q, 0, 0, dvec, gamma, A, 1.e-3, 1.e6, 3);

    return  atan((pow(pi_greco, 3. / 2.) * gamma * L * k) / z[0] * 2 * pi_greco);

}
double lhs_delta_Elatt(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2;
    double mass = gjack[e].jack[1][j];

    E2 = get_E2_n(n, e, j, gjack);
    //     double E20=gjack[e].jack[4][j];
    double hatp2 = 4. * sin(dvec1[0] * pi_greco / params[e].data.L[1]) * sin(dvec1[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec1[1] * pi_greco / params[e].data.L[2]) * sin(dvec1[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec1[2] * pi_greco / params[e].data.L[3]) * sin(dvec1[2] * pi_greco / params[e].data.L[3]);
    //     double E20=gjack[e].jack[4][j];
    double E2fL = acosh(cosh(mass) + 0.5 * (hatp2));
    hatp2 = 4. * sin(dvec2[0] * pi_greco / params[e].data.L[1]) * sin(dvec2[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec2[1] * pi_greco / params[e].data.L[2]) * sin(dvec2[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec2[2] * pi_greco / params[e].data.L[3]) * sin(dvec2[2] * pi_greco / params[e].data.L[3]);
    E2fL += acosh(cosh(mass) + 0.5 * (+hatp2));

    double L = params[e].data.L[1];
    double Ef1 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double Ef2 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double E2f = Ef1 + Ef2;

    E2 = E2 - E2fL + E2f;
    //     printf("%g  %g  %g \n",E2,E2fL,E2f);
    double E2_CM = energy_CM(E2, dvec, L);
    double k = sqrt(E2_CM * E2_CM / 4. - mass * mass);

    double q = k * L / (2. * pi_greco);
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];
    dzeta_function(z, q * q, 0, 0, dvec, gamma, A, 1.e-3, 1.e6, 3);

    return  atan((pow(pi_greco, 3. / 2.) * gamma * L * k) / z[0] * 2 * pi_greco);


}

double compute_hatp2(int* dvec, double L) {
    double hatp2 = sin(dvec[0] * pi_greco / L) * sin(dvec[0] * pi_greco / L);
    hatp2 += sin(dvec[1] * pi_greco / L) * sin(dvec[1] * pi_greco / L);
    hatp2 += sin(dvec[2] * pi_greco / L) * sin(dvec[2] * pi_greco / L);
    hatp2 *= 4.;
    return hatp2;
}


double lhs_kcotd_ECM_latt(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2;
    double mass = gjack[e].jack[1][j];
    double L = params[e].data.L[1];


    E2 = get_E2_n(n, e, j, gjack);
    double normp = compute_hatp2(dvec, L);
    double E2_CM = acosh(cosh(E2) - 0.5 * (normp));

    double k = sqrt(E2_CM * E2_CM / 4. - mass * mass);
    //     double k=    acos(-cosh(E2_CM/2.) + cosh(mass)  +1) ;

    double q = k * L / (2. * pi_greco);
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];
    dzeta_function(z, q * q, 0, 0, dvec, gamma, A, 1.e-3, 1.e6, 3);

    return  z[0] * 2 * pi_greco / (pow(pi_greco, 3. / 2.) * gamma * L);


}

double lhs_delta_ECM_latt(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2;
    double mass = gjack[e].jack[1][j];
    double L = params[e].data.L[1];


    E2 = get_E2_n(n, e, j, gjack);
    double normp = compute_hatp2(dvec, L);
    double E2_CM = acosh(cosh(E2) - 0.5 * (normp));

    //     double k=sqrt(E2_CM*E2_CM/4. -mass*mass);
    double k = acos(-cosh(E2_CM / 2.) + cosh(mass) + 1);

    double q = k * L / (2. * pi_greco);
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];
    dzeta_function(z, q * q, 0, 0, dvec, gamma, A, 1.e-3, 1.e6, 3);

    return  atan((pow(pi_greco, 3. / 2.) * gamma * L * k) / z[0] * 2 * pi_greco);


}
double compute_k(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double L = params[e].data.L[1];
    double mass = gjack[e].jack[1][j];
    if (gjack[e].Nobs >= 443 && gjack[e].jack[443][j] != 0) {
        mass = gjack[e].jack[443][j];
    }
    double E2 = gjack[e].jack[4][j];
    double E2_0_p1 = gjack[e].jack[100][j];
    double E2_0_p11 = gjack[e].jack[102][j];
    double E2_0_p111 = gjack[e].jack[104][j];
    double E2_0_A1 = gjack[e].jack[80][j];
    double E2_CM;
    if (n == 0) {//E2_0
        int dvec[3] = { 0,0,0 };
        E2_CM = energy_CM(E2, dvec, L);
    }
    else if (n == 1) {//E2_0_p1
        int dvec[3] = { 1,0,0 };
        E2_CM = energy_CM(E2_0_p1, dvec, L);
    }
    else if (n == 2) {//E2_0_p11
        int dvec[3] = { 1,1,0 };
        E2_CM = energy_CM(E2_0_p11, dvec, L);
    }
    else if (n == 4) {//E2_0_p111
        int dvec[3] = { 1,1,1 };
        E2_CM = energy_CM(E2_0_p111, dvec, L);
    }
    else if (n == 3) {//E2_0_A1
        int dvec[3] = { 0,0,0 };
        E2_CM = energy_CM(E2_0_A1, dvec, L);
    }
    else {
        printf("%s n=%d not implemented\n", __func__, n);
        exit(1);
    }
    return sqrt(E2_CM * E2_CM / 4. - mass * mass);

}

double compute_k_m_g(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double L = params[e].data.L[1];
    double mass = gjack[e].jack[1][j];
    // if (gjack[e].Nobs>=443 && gjack[e].jack[443][j]!=0){
    //     mass=gjack[e].jack[443][j];
    // }
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2 = get_E2_g_n(n, e, j, gjack);
    double E2_CM = energy_CM(E2, dvec, L);

    return sqrt(E2_CM * E2_CM / 4. - mass * mass) / mass;

}

double rhs_delta(int n, int Nvar, double* x, int Npar, double* P) {
    double a0m = P[0], r0m = P[1];
    double k_m = x[10];
    double kcot = 1.0 / (a0m)+r0m * k_m * k_m / (2.);

    if (Npar > 2) {
        kcot = -P[2] * r0m * r0m * r0m * k_m * k_m * k_m * k_m;
    }

    return  atan(k_m / kcot);
}


double rhs_delta_g(int n, int Nvar, double* x, int Npar, double* P) {
    double a0m = P[0], r0m = P[1];
    double k_m = x[10];
    double kcot = 1.0 / (a0m)+r0m * k_m * k_m / (2.);

    if (Npar > 2) {
        kcot = -P[2] * r0m * r0m * r0m * k_m * k_m * k_m * k_m;
    }

    return  atan(k_m / kcot);

}
double rhs_kcotd(int n, int Nvar, double* x, int Npar, double* P) {
    double a0 = P[0], r0 = P[1];
    double k = x[6];
    double kcot = 1.0 / a0 + r0 * k * k / 2.;

    if (Npar > 2) {
        kcot = -P[2] * r0 * r0 * r0 * k * k * k * k;
    }

    return  kcot;

}

double rhs_kcotd_m(int n, int Nvar, double* x, int Npar, double* P) {
    double a0m = P[0], r0m = P[1];
    double k_m = x[10];
    double kcot_m = 1.0 / a0m + r0m * k_m * k_m / 2.;

    if (Npar > 2) {
        kcot_m = -P[2] * r0m * r0m * r0m * k_m * k_m * k_m * k_m;
    }

    return  kcot_m;

}

double rhs_kcotd_g(int n, int Nvar, double* x, int Npar, double* P) {
    double a0 = P[0], r0 = P[1];
    double k = x[13] * x[1]; // k/m   * m
    double kcot = 1.0 / a0 + r0 * k * k / 2.;

    if (Npar > 2) {
        kcot = -P[2] * r0 * r0 * r0 * k * k * k * k;
    }

    return  kcot;

}

double to_invert_k_from_phase_shift(int n, int Nvar, double* x, int Npar, double* P) {
    double a0m0 = P[0], r0m0 = P[1];//, P2=P[2];
    double L = x[0];
    double mass = x[1];
    double k = x[2];
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];

    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);

    double E2_CM = sqrt(k * k + mass * mass) * 2.;
    double E2 = sqrt((k * k + mass * mass) * 4. + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));

    double qsq = k * k * (L / (2 * pi_greco)) * (L / (2 * pi_greco));
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];
    //       printf("a0=%.12g   r0=%.12g    P2=%.12g  dvec=(%d,%d,%d)\n",a0,r0,P2,dvec[0] ,dvec[1] ,dvec[2]);

    double zinter = zeta.compute(L, n, mass/*(L/(2*pi_greco))*/, qsq);
    //     if(n==2 && fabs(L-36)<1e-3){
    //         dzeta_function(z,  qsq,0 , 0, dvec, gamma, A, 1.e-3, 1.e6 ,3);
    //         printf("zeta=%g  zint=%g    qsq=%g  gamma=%g  L=%g   mass=%g\n",z[0],zinter,qsq,gamma,L,mass);
    //     }
    //     double zinter=zeta_qsqg.compute(qsq,n,gamma);
    z[0] = zinter; z[1] = 0;
    //      if (Npar>=3){
    //          dzeta_function(z,  qsq,0 , 0, dvec, gamma, A, 1.e-3, 1.e6 ,2);
    //      }
    std::complex<double>  zc(z[0], z[1]);

    //     std::complex<double>   q=std::sqrt(std::complex<double> (qsq,0));

    double r = real(zc * 2. * pi_greco / (pow(pi_greco, 3. / 2.) * L * gamma * mass));
    //     double kcotdelta=1.0/a0  + r0*k*k/2. - P2*r0*r0*r0*k*k*k*k;
    double k_m = k / mass;
    double kcotdelta_m = 1.0 / a0m0;  //   (k cot(d) )/ mass
    if (Npar >= 2) {
        kcotdelta_m += r0m0 * k_m * k_m / 2.;
    }
    if (Npar >= 3) {
        kcotdelta_m += P[2] * r0m0 * r0m0 * r0m0 * k_m * k_m * k_m * k_m;
    }
    //      if(n==2 && fabs(L-36)<1e-3)  printf("kcotdelta_m =%g  k=%g   r=%g  qsq=%g  gamma=%g  z=%g, %g   fun=%g  dvec=(%d,%d,%d)\n",kcotdelta_m,k,r,qsq,gamma,z[0],z[1], kcotdelta_m-r,dvec[0] ,dvec[1] ,dvec[2]);
    //    std::cout<< "ZC " << zc << "       q "<< q<<endl;
        //   printf(" k=%g   f(k)=%g  kcotdelta_m=%g   r=%g  z=%g  gamma=%g   mass=%g  L=%g   dvec=(%d,%d,%d) E2CM=%g\n ",k, kcotdelta_m-r ,kcotdelta_m, r , z[0] ,gamma,mass, L, dvec[0] ,dvec[1] ,dvec[2], E2_CM);
    return         kcotdelta_m - r;


}

double to_invert_k_from_phase_shift_g(int n, int Nvar, double* x, int Npar, double* P) {
    double a0m0 = P[0], r0m0 = P[1];//, P2=P[2];
    double L = x[0];
    double mass = x[1];
    double k = x[2];
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];

    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);

    double E2_CM = sqrt(k * k + mass * mass) * 2.;
    double E2 = sqrt((k * k + mass * mass) * 4. + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));

    double qsq = k * k * (L / (2 * pi_greco)) * (L / (2 * pi_greco));
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];

    double zinter = zeta.compute(L, n, mass/*(L/(2*pi_greco))*/, qsq);

    z[0] = zinter; z[1] = 0;

    std::complex<double>  zc(z[0], z[1]);

    double r = real(zc * 2. * pi_greco / (pow(pi_greco, 3. / 2.) * L * gamma * mass));
    double k_m = k / mass;
    double kcotdelta_m = 1.0 / a0m0;  //   (k cot(d) )/ mass
    if (Npar >= 2) {
        kcotdelta_m += r0m0 * k_m * k_m / 2.;
    }
    if (Npar >= 3) {
        kcotdelta_m += P[2] * r0m0 * r0m0 * r0m0 * k_m * k_m * k_m * k_m;
    }
    // printf(" k=%g   f(k)=%g  kcotdelta_m=%g   r=%g  z=%g  gamma=%g   mass=%g  L=%g   dvec=(%d,%d,%d) E2CM=%g\n ",k, kcotdelta_m-r ,kcotdelta_m, r , z[0] ,gamma,mass, L, dvec[0] ,dvec[1] ,dvec[2], E2_CM);
    // printf("P[0]=%g \t",P[0]);
    // if (Npar>=2) printf("P[1]=%g \t",P[1]);
    // if (Npar>=3) printf("P[2]=%g \t",P[2]);
    // printf("\n");
    return         kcotdelta_m - r;


}
double rhs_k_from_phase_shift(int n, int Nvar, double* x, int Npar, double* P) {
    double E2;
    double xx[3] = { x[0],x[1],1 }; //{L,mass,   k to be find by the bisection}
    //double p[5]={0 ,2.*pi_greco/x[0], sqrt(2)*2.*pi_greco/x[0] , sqrt(3)*2.*pi_greco/x[0],0 };

//      printf("before   n=%d   L=%g   m=%g  a=%g  r=%g   P=%g\n",n,x[0],x[1],P[0],P[1],P[2]);
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double mass = x[1];
    double L = x[0];
    // xmin have to be k in free theory

    double E1f = sqrt(mass * mass + (2. * pi_greco / L) * (2. * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double E2f = sqrt(mass * mass + (2. * pi_greco / L) * (2. * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double Ef = E1f + E2f;
    double ECMfsq = Ef * Ef - (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]);
    // double gamma=E2f/sqrt(ECMfsq);

    double kf = sqrt(ECMfsq / 4. - mass * mass);

    E1f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax1[0]) * (dmax1[0]) + dmax1[1] * dmax1[1] + (dmax1[2]) * (dmax1[2])));
    E2f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax2[0]) * (dmax2[0]) + dmax2[1] * dmax2[1] + (dmax2[2]) * (dmax2[2])));
    Ef = E1f + E2f;
    ECMfsq = Ef * Ef - (2 * pi_greco / L) * (2 * pi_greco / L) * ((dvec[0]) * (dvec[0]) + dvec[1] * dvec[1] + dvec[2] * dvec[2]);
    double kf1 = sqrt(ECMfsq / 4. - mass * mass);


    double xmin = kf + 1e-8;
    double xmax = kf1 - 1e-6;//- (kf1-kf)/1e+3;
//     if(n==2 && fabs(L-36)<1e-3) 
    E2 = rtsafe(to_invert_k_from_phase_shift, n, 3, xx, Npar, P, 2/*ivar*/, 0. /*input*/, xmin, xmax, 1e-5, 100, 4e-5);


    return E2 / mass;
}


double rhs_k_from_phase_shift_g(int n, int Nvar, double* x, int Npar, double* P) {
    double E2;
    double xx[3] = { x[0],x[1],1 }; //{L,mass,   k to be find by the bisection}
    //double p[5]={0 ,2.*pi_greco/x[0], sqrt(2)*2.*pi_greco/x[0] , sqrt(3)*2.*pi_greco/x[0],0 };

//      printf("before   n=%d   L=%g   m=%g  a=%g  r=%g   P=%g\n",n,x[0],x[1],P[0],P[1],P[2]);
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double mass = x[1];
    double L = x[0];
    // xmin have to be k in free theory

    double E1f = sqrt(mass * mass + (2. * pi_greco / L) * (2. * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double E2f = sqrt(mass * mass + (2. * pi_greco / L) * (2. * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double Ef = E1f + E2f;
    double ECMfsq = Ef * Ef - (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]);
    // double gamma=E2f/sqrt(ECMfsq);

    double kf = sqrt(ECMfsq / 4. - mass * mass);

    E1f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax1[0]) * (dmax1[0]) + dmax1[1] * dmax1[1] + (dmax1[2]) * (dmax1[2])));
    E2f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax2[0]) * (dmax2[0]) + dmax2[1] * dmax2[1] + (dmax2[2]) * (dmax2[2])));
    Ef = E1f + E2f;
    ECMfsq = Ef * Ef - (2 * pi_greco / L) * (2 * pi_greco / L) * ((dvec[0]) * (dvec[0]) + dvec[1] * dvec[1] + dvec[2] * dvec[2]);
    double kf1 = sqrt(ECMfsq / 4. - mass * mass);


    double xmin = kf + 1e-6;
    double xmax = kf1 - 1e-6;//- (kf1-kf)/1e+3;
//     if(n==2 && fabs(L-36)<1e-3) 
    E2 = rtsafe(to_invert_k_from_phase_shift_g, n, 3, xx, Npar, P, 2/*ivar*/, 0. /*input*/, xmin, xmax, 1e-5, 100, 4e-5);


    return E2 / mass;
}



double rhs_deltaE2_m_quant_cond(int n, int Nvar, double* x, int Npar, double* P) {

    double xx[3] = { x[0],x[1],1 }; //{L,mass,   k to be find by the bisection}
    //double p[5]={0 ,2.*pi_greco/x[0], sqrt(2)*2.*pi_greco/x[0] , sqrt(3)*2.*pi_greco/x[0],0 };

    //      printf("before   n=%d   L=%g   m=%g  a=%g  r=%g   P=%g\n",n,x[0],x[1],P[0],P[1],P[2]);
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double mass = x[1];
    double L = x[0];
    // xmin have to be k in free theory
    double E1f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double E2f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double Ef = E1f + E2f;
    double ECMfsq = Ef * Ef - (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]);

    double kf = sqrt(ECMfsq / 4. - mass * mass);


    E1f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax1[0]) * (dmax1[0]) + dmax1[1] * dmax1[1] + (dmax1[2]) * (dmax1[2])));
    E2f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax2[0]) * (dmax2[0]) + dmax2[1] * dmax2[1] + (dmax2[2]) * (dmax2[2])));
    double Ef_p = E1f + E2f;
    double ECMfsq_p = Ef_p * Ef_p - (2 * pi_greco / L) * (2 * pi_greco / L) * ((dvec[0]) * (dvec[0]) + dvec[1] * dvec[1] + dvec[2] * dvec[2]);
    double kf1 = sqrt(ECMfsq_p / 4. - mass * mass);


    double xmin = kf + 1e-6;
    double xmax = kf1 - 1e-6;//- (kf1-kf)/1e+3;
//     printf("L=%g  n=%d      P0=%g  P1=%g  kmin=%g  kmax=%g   d=(%d,%d,%d)=(%d,%d,%d)+(%d,%d,%d)\n",L,n,P[0],P[1],xmin,xmax,
//            dvec[0],dvec[1],dvec[2], dvec1[0],dvec1[1],dvec1[2], dvec2[0],dvec2[1],dvec2[2]    );
    double k = rtsafe(to_invert_k_from_phase_shift, n, 3, xx, Npar, P, 2/*ivar*/, 0. /*input*/, xmin, xmax, 1e-5, 100, 1e-4);
    double E2 = sqrt((k * k + mass * mass) * 4. + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));

    //     return (E2-Ef)/mass;
    return (E2) / mass;
    //     double E2CM2=(k*k+mass*mass)*4.;
    //     double hatp2=4.*sin(dvec[0]*pi_greco/L) *sin(dvec[0]*pi_greco/L) ;
    //     hatp2+=4.*sin(dvec[1]*pi_greco/L) *sin(dvec[1]*pi_greco/L) ;
    //     hatp2+=4.*sin(dvec[2]*pi_greco/L) *sin(dvec[2]*pi_greco/L) ;
    //     
    //     return acosh(1.+0.5*(E2CM2 +hatp2) );
    //     return (E2)/mass;

}



double rhs_deltaE2_m_quant_cond_g(int n, int Nvar, double* x, int Npar, double* P) {

    double mass = x[1];
    // double mass=x[Nvar-1];
    double xx[3] = { x[0],mass,1 }; //{L,mass,   k to be find by the bisection}
    //double p[5]={0 ,2.*pi_greco/x[0], sqrt(2)*2.*pi_greco/x[0] , sqrt(3)*2.*pi_greco/x[0],0 };

    //      printf("before   n=%d   L=%g   m=%g  a=%g  r=%g   P=%g\n",n,x[0],x[1],P[0],P[1],P[2]);
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);

    double L = x[0];
    // xmin have to be k in free theory
    double E1f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double E2f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double Ef = E1f + E2f;
    double ECMfsq = Ef * Ef - (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]);

    double kf = sqrt(ECMfsq / 4. - mass * mass);


    E1f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax1[0]) * (dmax1[0]) + dmax1[1] * dmax1[1] + (dmax1[2]) * (dmax1[2])));
    E2f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax2[0]) * (dmax2[0]) + dmax2[1] * dmax2[1] + (dmax2[2]) * (dmax2[2])));
    double Ef_p = E1f + E2f;
    double ECMfsq_p = Ef_p * Ef_p - (2 * pi_greco / L) * (2 * pi_greco / L) * ((dvec[0]) * (dvec[0]) + dvec[1] * dvec[1] + dvec[2] * dvec[2]);
    double kf1 = sqrt(ECMfsq_p / 4. - mass * mass);


    double xmin = kf + 1e-6;
    double xmax = kf1 - 1e-6;

    double k = rtsafe(to_invert_k_from_phase_shift_g, n, 3, xx, Npar, P, 2/*ivar*/, 0. /*input*/, xmin, xmax, 1e-5, 100, 1e-4);
    double E2 = sqrt((k * k + mass * mass) * 4. + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));

    double mass_inf = x[Nvar - 1];
    // printf("m  %g  %g    %d\n",mass,mass_inf,Nvar-1);
    return (E2 / mass_inf);


}

double rhs_deltaE2_m_QC2(int n, int Nvar, double* x, int Npar, double* P) {

    double mass = x[1];
    double xx[3] = { x[0],mass,1 }; //{L,mass,   k to be find by the bisection}


    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);

    double L = x[0];
    // xmin have to be k in free theory
    double E1f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double E2f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double Ef = E1f + E2f;
    double ECMfsq = Ef * Ef - (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]);

    double kf = sqrt(ECMfsq / 4. - mass * mass);


    E1f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax1[0]) * (dmax1[0]) + dmax1[1] * dmax1[1] + (dmax1[2]) * (dmax1[2])));
    E2f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax2[0]) * (dmax2[0]) + dmax2[1] * dmax2[1] + (dmax2[2]) * (dmax2[2])));
    double Ef_p = E1f + E2f;
    double ECMfsq_p = Ef_p * Ef_p - (2 * pi_greco / L) * (2 * pi_greco / L) * ((dvec[0]) * (dvec[0]) + dvec[1] * dvec[1] + dvec[2] * dvec[2]);
    double kf1 = sqrt(ECMfsq_p / 4. - mass * mass);


    double xmin = kf + 1e-6;
    double xmax = kf1 - 1e-6;

    double k = rtsafe(to_invert_k_from_phase_shift_g, n, 3, xx, Npar, P, 2/*ivar*/, 0. /*input*/, xmin, xmax, 1e-5, 100, 1e-4);
    double E2 = sqrt((k * k + mass * mass) * 4. + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));

    double mass_inf = x[2];
    return (E2 / mass_inf);


}

#ifdef PYTHON
double rhs_E3_m_QC3_pole(int n, int Nvar, double* x, int Npar, double* P) {

    double Pkcot[2];
    Pkcot[0] = x[Nvar - 3];
    Pkcot[1] = x[Nvar - 2];

    int Nkcot = 2;
    int Nkiso = Npar;
    int dvec[3];//,dvec1[3],dvec2[3],dmax1[3],dmax2[3];
    init_dvec_QC3_pole(n, dvec);
    //init_dvec(n,dvec,dvec1,dvec2,dmax1,dmax2);
    double nnP[3];
    nnP[0] = (double)dvec[0]; nnP[1] = (double)dvec[1]; nnP[2] = (double)dvec[2];
    double mass = x[Nvar - 1];//x[1];
    double L = x[0];
    int steps = 4;

    L = L * mass;

    double Estart = x[11] - x[12];
    double Eend = x[11] + x[12];

    // printf("E=[%g , %g]   E1f=%g   E2f=%g  m=%g\n",Estart,Eend,Estart,Eend,mass);
    double r = python_detQC_call(Estart, Eend, steps, L, nnP, Nkcot, Pkcot, Nkiso, P);
    //     printf("res=%g\n",r);
    return r;
}


double rhs_E3_m_QC3_2sol(int n, int Nvar, double* x, int Npar, double* P) {

    double Pkcot[2];
    Pkcot[0] = x[Nvar - 3];
    Pkcot[1] = x[Nvar - 2];

    int Nkcot = 2;
    int Nkiso = Npar;
    int dvec[3];//,dvec1[3],dvec2[3],dmax1[3],dmax2[3];
    init_dvec_QC3_pole(n, dvec);
    //init_dvec(n,dvec,dvec1,dvec2,dmax1,dmax2);
    double nnP[3];
    nnP[0] = (double)dvec[0]; nnP[1] = (double)dvec[1]; nnP[2] = (double)dvec[2];
    double mass = x[Nvar - 1];
    double L = x[0];
    int steps = n;
    L = L * mass;

    double Estart = 3.05 + 1e-5;
    double DeltaE = 0.005;
    double r = python_detQC_call(Estart, DeltaE, n, L, nnP, Nkcot, Pkcot, Nkiso, P);
    // printf("r=%g   n=%d\n",r,n);
    return r;
}

double rhs_E3_m_QC3(int n, int Nvar, double* x, int Npar, double* P) {

    double Pkcot[2];
    Pkcot[0] = x[Nvar - 2];
    Pkcot[1] = x[Nvar - 1];

    double Pkiso[2];
    Pkiso[0] = P[0];
    if (Npar > 1)    Pkiso[1] = P[1];
    int Nkcot = 2;
    int Nkiso = Npar;
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double nnP[3];
    nnP[0] = (double)dvec[0]; nnP[1] = (double)dvec[1]; nnP[2] = (double)dvec[2];
    double mass = x[1];
    double L = x[0];
    int steps = 4;
    double E1f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double E2f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double E3f = mass;
    double Estart = (E1f + E2f + E3f) / mass + 1e-6;
    //     double Estart=(E1f+E2f+E3f)/mass+0.05;
    //     double Estart=x[8]-x[8]*1e-3 ;// E3/mass

    E1f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax1[0]) * (dmax1[0]) + dmax1[1] * dmax1[1] + (dmax1[2]) * (dmax1[2])));
    E2f = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * ((dmax2[0]) * (dmax2[0]) + dmax2[1] * dmax2[1] + (dmax2[2]) * (dmax2[2])));
    double Eend = (E1f + E2f + E3f) / mass - 1e-6;

    //     double Eend=x[8]+x[8]*1e-3 ;
    //     double Eend=Estart+0.10;

    L = L * mass;
    //if (n==0){ Estart=x[8]-x[9];  Eend=x[8]+x[9];}
    Estart = x[8] - x[9];  Eend = x[8] + x[9];

    //   printf("E=[%g , %g]   E1f=%g   E2f=%g  m=%g\n",Estart,Eend,E1f,E2f,mass);
    double r = python_detQC(Estart, Eend, steps, L, nnP, Nkcot, Pkcot, Nkiso, Pkiso);
    //    printf("res=%g\n",r);
    return r;
}


double rhs_E3_m_QC3_latt(int n, int Nvar, double* x, int Npar, double* P) {

    double Pkcot[2];
    Pkcot[0] = x[Nvar - 2];
    Pkcot[1] = x[Nvar - 1];

    double Pkiso[2];
    Pkiso[0] = P[0];
    if (Npar > 1) Pkiso[1] = P[1];

    int Nkcot = 2;
    int Nkiso = Npar;
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double nnP[3];
    nnP[0] = (double)dvec[0]; nnP[1] = (double)dvec[1]; nnP[2] = (double)dvec[2];
    double mass = x[1];
    double L = x[0];
    int steps = 4;


    double hatp2 = 4. * sin(dvec1[0] * pi_greco / L) * sin(dvec1[0] * pi_greco / L);
    hatp2 += 4. * sin(dvec1[1] * pi_greco / L) * sin(dvec1[1] * pi_greco / L);
    hatp2 += 4. * sin(dvec1[2] * pi_greco / L) * sin(dvec1[2] * pi_greco / L);
    //     double E20=gjack[e].jack[4][j];
    double E2fL = acosh(cosh(mass) + 0.5 * (hatp2));
    hatp2 = 4. * sin(dvec2[0] * pi_greco / L) * sin(dvec2[0] * pi_greco / L);
    hatp2 += 4. * sin(dvec2[1] * pi_greco / L) * sin(dvec2[1] * pi_greco / L);
    hatp2 += 4. * sin(dvec2[2] * pi_greco / L) * sin(dvec2[2] * pi_greco / L);
    E2fL += acosh(cosh(mass) + 0.5 * (+hatp2));

    double Ef1 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double Ef2 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double E2f = Ef1 + Ef2;



    double dE = (-E2fL + E2f) / mass;
    //if (n==0){ Estart=x[8]-x[9];  Eend=x[8]+x[9];}
    double Estart = x[8] - x[9] + dE;
    double Eend = x[8] + x[9] + dE;

    L = L * mass;
    //   printf("E=[%g , %g]   E1f=%g   E2f=%g  m=%g\n",Estart,Eend,E1f,E2f,mass);
    double r = python_detQC(Estart, Eend, steps, L, nnP, Nkcot, Pkcot, Nkiso, Pkiso);
    //    printf("res=%g\n",r);
    return r;
}
#endif


template<int order>
double rhs_poly_order_E3_m(int n, int Nvar, double* x, int Npar, double* P) {

    double L = x[0];
    int on = order * n;
    double r = 0;
    double oL = 1;
    for (int o = 0; o < order;o++) {
        r += P[on + o] * oL;
        oL *= L;
    }

    //double r=P[on]+P[on+1]*L+P[3*n+2]*L*L;

    return r;
}

double rhs_poly_E3_m(int n, int Nvar, double* x, int Npar, double* P) {

    double L = x[0];

    //     double r=P[2*n]+P[2*n+1]*L;
    double r = P[3 * n] + P[3 * n + 1] * L + P[3 * n + 2] * L * L;

    return r;
}


double rhs_poly2_E3_m(int n, int Nvar, double* x, int Npar, double* P) {

    double L = x[0];

    //     double r=P[2*n]+P[2*n+1]*L;
    double r = P[2 * n] + P[2 * n + 1] * L;

    return r;
}

double rhs_const_E3_m(int n, int Nvar, double* x, int Npar, double* P) {

    // double L=x[0];

    //     double r=P[2*n]+P[2*n+1]*L;
    double r = P[n];

    return r;
}
double to_invert_q_from_phase_shift(int n, int Nvar, double* x, int Npar, double* P) {
    double L_a0 = P[0], r0L = P[1];//, P2=P[2];

    double mL_2pi = x[0];
    double q = x[1];
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);

    double E2_CM = sqrt(q * q + mL_2pi * mL_2pi) * 2.;
    double E2 = sqrt((q * q + mL_2pi * mL_2pi) * 4. + (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));

    double qsq = q * q;
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];
    double L = x[0];
    double zinter = zeta.compute(L, n, mL_2pi, qsq);
    z[0] = zinter; z[1] = 0;

    std::complex<double>  zc(z[0], z[1]);

    double r = real(zc / (pow(pi_greco, 3. / 2.) * gamma));
    double qcotdelta = L_a0 + +r0L * q * q / 2.;  //   (k cot(d) )/ mass

    return         qcotdelta - r;


}

double rhs_q_from_phase_shift(int n, int Nvar, double* x, int Npar, double* P) {
    double q;
    double xx[3];
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double mL_2pi = x[7];

    double E1f = sqrt(mL_2pi * mL_2pi + (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double E2f = sqrt(mL_2pi * mL_2pi + (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double Ef = E1f + E2f;
    double ECMfsq = Ef * Ef - (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]);
    // double gamma=E2f/sqrt(ECMfsq);

    double qf = sqrt(ECMfsq / 4. - mL_2pi * mL_2pi);

    E1f = sqrt(mL_2pi * mL_2pi + ((dvec1[0]) * (dvec1[0]) + dvec1[1] * dvec1[1] + (dvec1[2] + 1) * (dvec1[2] + 1)));
    E2f = sqrt(mL_2pi * mL_2pi + ((dvec2[0]) * (dvec2[0]) + dvec2[1] * dvec2[1] + (dvec2[2] - 1) * (dvec2[2] - 1)));
    Ef = E1f + E2f;
    ECMfsq = Ef * Ef - ((dvec[0]) * (dvec[0]) + dvec[1] * dvec[1] + dvec[2] * dvec[2]);
    double qf1 = sqrt(ECMfsq / 4. - mL_2pi * mL_2pi);


    double xmin = qf + 1e-10;
    double xmax = qf1 - 1e-10;//- (kf1-kf)/1e+3;
//     printf("kf=%g   kf1=%g      x=%g   %g\n",qf,qf1,xmin,xmax);
    xx[0] = mL_2pi; xx[1] = 1;  xx[2] = x[0];
    q = rtsafe(to_invert_k_from_phase_shift, n, 2, xx, Npar, P, 1/*ivar*/, 0. /*input*/, xmin, xmax, 1e-5, 100, 1e-4);

    return q;
}

double lhs_k(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E2;
    int dvec[3];
    double mass = gjack[e].jack[1][j];
    if (n == 0) {//E2_0
        E2 = gjack[e].jack[4][j];
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
    }
    else if (n == 1) {//E2_0_p1
        E2 = gjack[e].jack[100][j];
        dvec[0] = 1; dvec[1] = 0; dvec[2] = 0;
    }
    else if (n == 2) {//E2_0_p11
        E2 = gjack[e].jack[102][j];
        dvec[0] = 1; dvec[1] = 1; dvec[2] = 0;
    }
    else if (n == 4) {//E2_0_p111
        E2 = gjack[e].jack[104][j];
        dvec[0] = 1; dvec[1] = 1; dvec[2] = 1;
    }
    else if (n == 3) {//E2_0_A1
        E2 = gjack[e].jack[80][j];
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
    }
    else { E2 = 0; dvec[0] = 0; dvec[1] = 0; dvec[2] = 0; exit(1); }
    double E2_CM = sqrt(E2 * E2 - (2. * pi_greco / params[e].data.L[1]) * (2. * pi_greco / params[e].data.L[1]) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));
    double k = sqrt((E2_CM * E2_CM / 4. - mass * mass)) / mass;
    return k;

}

double lhs_k_g(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E2;
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    double mass = gjack[e].jack[443][j];

    E2 = get_E2_g_n(n, e, j, gjack);

    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);

    double E2_CM = sqrt(E2 * E2 - (2. * pi_greco / params[e].data.L[1]) * (2. * pi_greco / params[e].data.L[1]) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));
    double k = sqrt((E2_CM * E2_CM / 4. - mass * mass)) / mass;
    return k;

}

double lhs_deltaE2_m_latt(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E2;
    double mass = gjack[e].jack[1][j];
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);

    if (n == 0) {//E2_0
        E2 = gjack[e].jack[4][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 1) {//E2_0_p1
        E2 = gjack[e].jack[100][j];
        //         dvec[0]=1; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 2) {//E2_0_p11
        E2 = gjack[e].jack[102][j];
        //         dvec[0]=1; dvec[1]=1; dvec[2]=0;
    }
    else if (n == 4) {//E2_0_p111
        E2 = gjack[e].jack[104][j];
        //         dvec[0]=1; dvec[1]=1; dvec[2]=1;
    }
    else if (n == 3) {//E2_0_A1
        E2 = gjack[e].jack[80][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else { E2 = 0; dvec[0] = 0; dvec[1] = 0; dvec[2] = 0; exit(1); }
    //     double E20=gjack[e].jack[4][j];
    double hatp2 = 4. * sin(dvec1[0] * pi_greco / params[e].data.L[1]) * sin(dvec1[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec1[1] * pi_greco / params[e].data.L[2]) * sin(dvec1[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec1[2] * pi_greco / params[e].data.L[3]) * sin(dvec1[2] * pi_greco / params[e].data.L[3]);
    //     double E20=gjack[e].jack[4][j];
    double E2fL = acosh(cosh(mass) + 0.5 * (hatp2));
    hatp2 = 4. * sin(dvec2[0] * pi_greco / params[e].data.L[1]) * sin(dvec2[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec2[1] * pi_greco / params[e].data.L[2]) * sin(dvec2[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec2[2] * pi_greco / params[e].data.L[3]) * sin(dvec2[2] * pi_greco / params[e].data.L[3]);
    E2fL += acosh(cosh(mass) + 0.5 * (+hatp2));

    double L = params[e].data.L[1];
    double Ef1 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double Ef2 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double E2f = Ef1 + Ef2;



    return (E2 - E2fL + E2f) / mass;
    //     return (E2)/mass;
    //     double E2_CM=sqrt(E2*E2-(2.*pi_greco/params[e].data.L[1])*(2.*pi_greco/params[e].data.L[1])*( dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]   ));
    //     double k=sqrt((E2_CM*E2_CM/4. - mass*mass))/mass;
    //     return k;
}



double lhs_deltaE2_m_latt_g(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E2;
    double mass = gjack[e].jack[1][j];
    // double mass=fit_info.ext_P[2][j];

    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);
    E2 = get_E2_g_n(n, e, j, gjack);
    double hatp2 = 4. * sin(dvec1[0] * pi_greco / params[e].data.L[1]) * sin(dvec1[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec1[1] * pi_greco / params[e].data.L[2]) * sin(dvec1[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec1[2] * pi_greco / params[e].data.L[3]) * sin(dvec1[2] * pi_greco / params[e].data.L[3]);
    //     double E20=gjack[e].jack[4][j];
    double E2fL = acosh(cosh(mass) + 0.5 * (hatp2));
    hatp2 = 4. * sin(dvec2[0] * pi_greco / params[e].data.L[1]) * sin(dvec2[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec2[1] * pi_greco / params[e].data.L[2]) * sin(dvec2[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec2[2] * pi_greco / params[e].data.L[3]) * sin(dvec2[2] * pi_greco / params[e].data.L[3]);
    E2fL += acosh(cosh(mass) + 0.5 * (+hatp2));

    double L = params[e].data.L[1];
    double Ef1 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double Ef2 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double E2f = Ef1 + Ef2;

    double mass_inf = fit_info.ext_P[2][j];
    return (E2 - E2fL + E2f) / mass_inf;

}



double lhs_deltaE2_m_latt_QC2(int n, int e, int j, vector<data_all> gjack, struct fit_type fit_info) {
    double E2;
    double mass = gjack[e].jack[1][j];

    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g(n, dvec, dvec1, dvec2, dmax1, dmax2);
    E2 = get_E2_g_n(n, e, j, gjack);
    double L = (double)gjack[e].header.L;
    double hatp2 = 4. * sin(dvec1[0] * pi_greco / L) * sin(dvec1[0] * pi_greco / L);
    hatp2 += 4. * sin(dvec1[1] * pi_greco / L) * sin(dvec1[1] * pi_greco / L);
    hatp2 += 4. * sin(dvec1[2] * pi_greco / L) * sin(dvec1[2] * pi_greco / L);
    //     double E20=gjack[e].jack[4][j];
    double E2fL = acosh(cosh(mass) + 0.5 * (hatp2));
    hatp2 = 4. * sin(dvec2[0] * pi_greco / L) * sin(dvec2[0] * pi_greco / L);
    hatp2 += 4. * sin(dvec2[1] * pi_greco / L) * sin(dvec2[1] * pi_greco / L);
    hatp2 += 4. * sin(dvec2[2] * pi_greco / L) * sin(dvec2[2] * pi_greco / L);



    double Ef1 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double Ef2 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double E2f = Ef1 + Ef2;

    double mass_inf = fit_info.x[2][0][j];

    return (E2 - E2fL + E2f) / mass_inf;

}

double lhs_k_p111(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E2;
    int dvec[3];
    double mass = gjack[e].jack[1][j];
    //E2_0_p111
    E2 = gjack[e].jack[104][j];
    dvec[0] = 1; dvec[1] = 1; dvec[2] = 1;

    double E2_CM = sqrt(E2 * E2 - (2. * pi_greco / params[e].data.L[1]) * (2. * pi_greco / params[e].data.L[1]) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));
    double k = sqrt((E2_CM * E2_CM / 4. - mass * mass)) / mass;
    return k;

}


double lhs_q(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E2;
    int dvec[3];
    double mass = gjack[e].jack[1][j];
    if (n == 0) {//E2_0
        E2 = gjack[e].jack[4][j];
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
    }
    else if (n == 1) {//E2_0_p1
        E2 = gjack[e].jack[100][j];
        dvec[0] = 1; dvec[1] = 0; dvec[2] = 0;
    }
    else if (n == 2) {//E2_0_p11
        E2 = gjack[e].jack[102][j];
        dvec[0] = 1; dvec[1] = 1; dvec[2] = 0;
    }
    else if (n == 4) {//E2_0_p111
        E2 = gjack[e].jack[104][j];
        dvec[0] = 1; dvec[1] = 1; dvec[2] = 1;
    }
    else if (n == 3) {//E2_0_A1
        E2 = gjack[e].jack[80][j];
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
    }
    else { E2 = 0; dvec[0] = 0; dvec[1] = 0; dvec[2] = 0; exit(1); }
    double E2_CM = sqrt(E2 * E2 - (2. * pi_greco / params[e].data.L[1]) * (2. * pi_greco / params[e].data.L[1]) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));
    double k = sqrt((E2_CM * E2_CM / 4. - mass * mass)) * params[e].data.L[1] / (2. * pi_greco);
    return k;

}


///////////////////////////////////////////////////////////////////////////////////////////
///////  lhs E3/m
//////////////////////////////////////////////////////////////////////////////////////////

double lhs_E3orE1_g_m(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E3;
    double mass = gjack[e].jack[443][j];


    if (n == 0) {//GEVP 1
        E3 = gjack[e].jack[444][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 1) {//GEVP 2
        E3 = gjack[e].jack[445][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else { E3 = 0; printf("lhs_E3orE1_m n=%d not implemented\n", n); exit(1); }
    return E3 / mass;
}



double lhs_E3orE1_m_complex(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E3;
    //double mass=gjack[e].jack[1][j];
    double mass = fit_info.ext_P[2][j];


    if (n == 0) {//GEVP 1
        E3 = gjack[e].jack[354][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 1) {//GEVP 2
        E3 = gjack[e].jack[355][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else {
        E3 = 0;
        // printf("lhs_E3orE1_m n=%d not implemented\n",n); exit(1);
    }

    return E3 / mass;
}

double lhs_E3orE1_m(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E3;
    double mass = gjack[e].jack[1][j];


    if (n == 0) {//GEVP 1
        E3 = gjack[e].jack[131][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 1) {//GEVP 2
        E3 = gjack[e].jack[132][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 2) {//GEVP 1 p1
        E3 = gjack[e].jack[142][j];
        //         dvec[0]=1; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 3) {//GEVP 1 p1
        E3 = gjack[e].jack[143][j];
        //         dvec[0]=1; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 4) {//GEVP 1 p11
        E3 = gjack[e].jack[156][j];
        //         dvec[0]=1; dvec[1]=1; dvec[2]=0;
    }
    else if (n == 5) {//GEVP 1 p11
        E3 = gjack[e].jack[157][j];
        //         dvec[0]=1; dvec[1]=1; dvec[2]=0;
    }
    else { E3 = 0; printf("lhs_E3orE1_m n=%d not implemented\n", n); exit(1); }
    return E3 / mass;
}

double lhs_E3_m(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E3;
    double mass = gjack[e].jack[1][j];
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);

    if (n == 0) {//E2_0
        E3 = gjack[e].jack[116][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 1) {//E2_0_p1
        E3 = gjack[e].jack[117][j];
        //         dvec[0]=1; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 2) {//E2_0_p11
        E3 = gjack[e].jack[118][j];
        //         dvec[0]=1; dvec[1]=1; dvec[2]=0;
    }
    else if (n == 4) {//E2_0_p111
        E3 = gjack[e].jack[119][j];
        //         dvec[0]=1; dvec[1]=1; dvec[2]=1;
    }
    else if (n == 3) {//E2_0_A1
        E3 = gjack[e].jack[120][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else { E3 = 0; printf("%s n=%d not implemented\n", __func__, n); }

    //     double hatp2=4.*sin(dvec1[0]*pi_greco/params[e].data.L[1]) *sin(dvec1[0]*pi_greco/params[e].data.L[1]) ;
    //     hatp2+=4.*sin(dvec1[1]*pi_greco/params[e].data.L[2]) *sin(dvec1[1]*pi_greco/params[e].data.L[2]) ;
    //     hatp2+=4.*sin(dvec1[2]*pi_greco/params[e].data.L[3]) *sin(dvec1[2]*pi_greco/params[e].data.L[3]) ;
    // 
    //     double E2fL=acosh(  cosh(mass) +0.5*( hatp2));
    //     hatp2=4.*sin(dvec2[0]*pi_greco/params[e].data.L[1]) *sin(dvec2[0]*pi_greco/params[e].data.L[1]) ;
    //     hatp2+=4.*sin(dvec2[1]*pi_greco/params[e].data.L[2]) *sin(dvec2[1]*pi_greco/params[e].data.L[2]) ;
    //     hatp2+=4.*sin(dvec2[2]*pi_greco/params[e].data.L[3]) *sin(dvec2[2]*pi_greco/params[e].data.L[3]) ;
    //     E2fL+=acosh(cosh(mass) +0.5*( + hatp2));
    //     
    //     double L=params[e].data.L[1];
    //     double Ef1=sqrt(mass*mass+(2*pi_greco/L)*(2*pi_greco/L)*(dvec1[0]*dvec1[0]+dvec1[1]*dvec1[1]+dvec1[2]*dvec1[2])   );
    //     double Ef2=sqrt(mass*mass+(2*pi_greco/L)*(2*pi_greco/L)*(dvec2[0]*dvec2[0]+dvec2[1]*dvec2[1]+dvec2[2]*dvec2[2])   );
    //     double E2f=Ef1+Ef2;
    //     return (E2-E2fL+E2f)/mass;
    return E3 / mass;
}

double lhs_E3_m_g(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double E3;
    double mass = gjack[e].jack[443][j];
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];

    if (n == 0) {//E2_0
        E3 = gjack[e].jack[444][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 1) {//E2_0_p1
        E3 = gjack[e].jack[445][j];
        //         dvec[0]=1; dvec[1]=0; dvec[2]=0;
    }

    else {
        E3 = 0;
        // printf("%s n=%d not implemented\n",__func__,n);
         //exit(1);
    }

    return E3 / mass;
}




double lhs_E3_m_latt(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {

    double E3 = lhs_E3_m(n, e, j, params, gjack, fit_info);

    double mass = gjack[e].jack[1][j];
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec(n, dvec, dvec1, dvec2, dmax1, dmax2);


    double hatp2 = 4. * sin(dvec1[0] * pi_greco / params[e].data.L[1]) * sin(dvec1[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec1[1] * pi_greco / params[e].data.L[2]) * sin(dvec1[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec1[2] * pi_greco / params[e].data.L[3]) * sin(dvec1[2] * pi_greco / params[e].data.L[3]);

    double E2fL = acosh(cosh(mass) + 0.5 * (hatp2));
    hatp2 = 4. * sin(dvec2[0] * pi_greco / params[e].data.L[1]) * sin(dvec2[0] * pi_greco / params[e].data.L[1]);
    hatp2 += 4. * sin(dvec2[1] * pi_greco / params[e].data.L[2]) * sin(dvec2[1] * pi_greco / params[e].data.L[2]);
    hatp2 += 4. * sin(dvec2[2] * pi_greco / params[e].data.L[3]) * sin(dvec2[2] * pi_greco / params[e].data.L[3]);
    E2fL += acosh(cosh(mass) + 0.5 * (+hatp2));

    double L = params[e].data.L[1];
    double Ef1 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double Ef2 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double E2f = Ef1 + Ef2;
    return E3 + (-E2fL + E2f) / mass;
}


///////////////////////////////////////////////////////////////////////////////////////////
///////kcotd
//////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// functions lhs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double M0_finite_volume_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    return gjack[e].jack[1][j];
}
double M0_g_finite_volume_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    return gjack[e].jack[443][j];
}


double M1_finite_volume_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    return gjack[e].jack[2][j];
}
template<int id>
double M1_p_finite_volume_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    return gjack[e].jack[id][j];
}


double DE_00_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    return gjack[e].jack[4][j] - 2 * gjack[e].jack[1][j];
}

double muDE_00_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double DE = gjack[e].jack[4][j] - 2 * gjack[e].jack[1][j];
    return DE * gjack[e].jack[1][j] / 2.;
}

double muDE_00_g_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double DE = gjack[e].jack[594][j] - 2 * gjack[e].jack[443][j];
    return DE * gjack[e].jack[443][j] / 2.;
}
double muDE_00_infm_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double DE = gjack[e].jack[4][j] - 2 * gjack[e].jack[1][j];
    double mu = fit_info.ext_P[0][j] * fit_info.ext_P[1][j] / (fit_info.ext_P[0][j] + fit_info.ext_P[1][j]);
    return DE * mu;
}


double muDE_01_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double DE = gjack[e].jack[19][j] - gjack[e].jack[1][j] - gjack[e].jack[2][j];
    double mu = gjack[e].jack[1][j] * gjack[e].jack[2][j] / (gjack[e].jack[1][j] + gjack[e].jack[2][j]);
    return DE * mu;
}

double muDE_01_div_shift_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double DE = gjack[e].jack[85][j] - gjack[e].jack[1][j] - gjack[e].jack[2][j];
    double mu = gjack[e].jack[1][j] * gjack[e].jack[2][j] / (gjack[e].jack[1][j] + gjack[e].jack[2][j]);
    return DE * mu;

}


double a_01_luescher_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {

    double* a = scattering_len_luscher(fit_info.Njack, gjack[e].jack[1], gjack[e].jack[2], gjack[e].jack[19], params[e].data.L[1]);
    double r = a[j];
    free(a);
    return r;//  
}


template<int id>
double lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    return gjack[e].jack[id][j];
}


template<int id1, int id2>
double lhs_diff(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    return gjack[e].jack[id1][j] - gjack[e].jack[id2][j];
}

template<int idBH>
double lhs_LminusBH(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double* a = scattering_len_luscher(fit_info.Njack, gjack[e].jack[1], gjack[e].jack[2], gjack[e].jack[19], params[e].data.L[1]);
    double r = a[j];
    free(a);
    return r - gjack[e].jack[idBH][j];
}
template<int id1, int id2>
double lhs_L_C2_m_C2shifted(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {

    double* a = scattering_len_luscher(fit_info.Njack, gjack[e].jack[1], gjack[e].jack[2], gjack[e].jack[id1], params[e].data.L[1]);
    double r = a[j];
    free(a);

    double* a1 = scattering_len_luscher(fit_info.Njack, gjack[e].jack[1], gjack[e].jack[2], gjack[e].jack[id2], params[e].data.L[1]);
    double r1 = a1[j];
    free(a1);

    return r - r1;
}

template<int id>
double lhs_Luescher_C(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    double* a = scattering_len_luscher(fit_info.Njack, gjack[e].jack[1], gjack[e].jack[2], gjack[e].jack[id], params[e].data.L[1]);
    double r = a[j];
    free(a);
    return r;
}


double a_01_BH_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    //return gjack[e].jack[12][j];// 0t_8tT_2
    //return gjack[e].jack[23][j];//   03t16
    return gjack[e].jack[53][j];//   03t16_shifted1
}
double a_00_BH_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    //return gjack[e].jack[10][j];
    return gjack[e].jack[21][j];//   03t16
}
double a_11_BH_lhs(int n, int e, int j, vector<cluster::IO_params> params, vector<data_all> gjack, struct fit_type fit_info) {
    return gjack[e].jack[11][j];
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// functions rhs
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double M_finite_volume(int n, int Nvar, double* x, int Npar, double* P) {
    double M = P[0];// mass at inif volume
    double L = x[0];
    //int K1=dbesk1(M*L);
    double K1 = exp(-M * L) / sqrt(M * L);
    double r = M;
    if (Npar > 1)
        r += P[1] * K1 / (M * L);
    return r;
}


double muDE_rhs(int n, int Nvar, double* x, int Npar, double* P) {
    double a = P[0];// mass at inif volume
    double L = x[0];
    //double mu=x[Nvar]*x[Nvar+1]/(x[Nvar]+x[Nvar+1]);
    //int K1=dbesk1(M*L);
    double r = -(2. * pi_greco * a) / (L * L * L);
    r *= 1 - 2.837297 * (a / L) + 6.375183 * (a / L) * (a / L);
    return r;

}

double a_luscher_infm(int n, int Nvar, double* x, int Npar, double* P) {
    double a = P[0];// mass at inif volume
    double L = x[0];
    double mu = x[Nvar] * x[Nvar + 1] / (x[Nvar] + x[Nvar + 1]);
    //int K1=dbesk1(M*L);
    double r = -(2. * pi_greco * a) / (mu * L * L * L);
    r *= 1 - 2.837297 * (a / L) + 6.375183 * (a / L) * (a / L);
    return r;

}
double a_luscher(int n, int Nvar, double* x, int Npar, double* P) {
    double a = P[0];// mass at inif volume
    double L = x[0];
    double mu = x[1] * x[1] / (x[1] + x[1]);
    //int K1=dbesk1(M*L);
    double r = -(2. * pi_greco * a) / (mu * L * L * L);
    r *= 1 - 2.837297 * (a / L) + 6.375183 * (a / L) * (a / L);
    return r;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// print fit band
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_fit_band_L(char** argv, vector<data_all> gjack, struct fit_type fit_info, const char* label, struct fit_result fit_out, int* en, double*** x, double*** y, vector<cluster::IO_params> params, std::vector<int> myen) {

    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_L.txt", argv[3], label);
    f = open_file(namefile, "w+");
    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double* tmpx = (double*)malloc(sizeof(double) * Nvar);
    double* tmpy = (double*)malloc(sizeof(double) * Njack);
    printf("writing: %s\n", namefile);
    for (int i = 0; i < 100; i++) {

        for (int j = 0;j < Njack;j++) {
            tmpx[0] = 14. + i * 0.5;
            tmpx[1] = x[j][0][1];//m0
            tmpx[2] = x[j][0][2];//m1
            tmpx[3] = x[j][0][3];//E20
            tmpx[4] = x[j][0][4];//E21
            tmpx[5] = (double)params[0].data.L[0];//T
            tmpx[6] = x[j][0][6];
            tmpx[7] = x[j][0][7];
            for (int i = 8; i < fit_info.Nvar; i++)
                tmpx[i] = x[j][0][i];

            for (int i = fit_info.Nvar; i < fit_info.Nvar + fit_info.n_ext_P; i++)
                tmpx[i] = fit_info.ext_P[i - fit_info.Nvar][j];



            tmpy[j] = fit_info.function(0, Nvar, tmpx, Npar, tif[j]);//N, Nvar, x ,Npar,P

        }
        fprintf(f, "%g  \t %g  %g\n", tmpx[0], tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));

    }

    free(tmpy);free(tmpx);
    fclose(f);
    ////////// end fit band k
    free_2(Njack, tif);
}

void print_fit_band_T(char** argv, vector<data_all> gjack, struct fit_type fit_info, const char* label, struct fit_result fit_out, int* en, double*** x, double*** y, vector<cluster::IO_params> params, std::vector<int> myen) {

    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_T.txt", argv[3], label);
    f = open_file(namefile, "w+");
    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double* tmpx = (double*)malloc(sizeof(double) * Nvar);
    double* tmpy = (double*)malloc(sizeof(double) * Njack);
    printf("writing: %s\n", namefile);
    for (int i = 0; i < 100; i++) {

        for (int j = 0;j < Njack;j++) {
            tmpx[0] = (double)params[0].data.L[1];//L
            tmpx[1] = x[j][0][1];//m0
            tmpx[2] = x[j][0][2];//m1
            tmpx[3] = x[j][0][3];//E20
            tmpx[4] = x[j][0][4];//E21
            tmpx[5] = 16 + i * 1;//T
            tmpx[6] = x[j][0][6];
            tmpx[7] = x[j][0][7];
            for (int i = 8; i < fit_info.Nvar; i++)
                tmpx[i] = x[j][0][i];

            for (int i = fit_info.Nvar; i < fit_info.Nvar + fit_info.n_ext_P; i++)
                tmpx[i] = fit_info.ext_P[i - fit_info.Nvar][j];

            tmpy[j] = fit_info.function(0, Nvar, tmpx, Npar, tif[j]);//N, Nvar, x ,Npar,P

        }
        fprintf(f, "%g  \t %g  %g\n", tmpx[5], tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));

    }

    free(tmpy);free(tmpx);
    fclose(f);
    ////////// end fit band k
    free_2(Njack, tif);
}

void print_fit_band_k(char** argv, vector<data_all> gjack, struct fit_type fit_info, const char* label, struct fit_result fit_out, int* en, double*** x, double*** y, vector<cluster::IO_params> params, std::vector<int> myen) {

    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_k.txt", argv[3], label);
    f = open_file(namefile, "w+");
    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double* tmpx = (double*)malloc(sizeof(double) * Nvar);
    double* tmpy = (double*)malloc(sizeof(double) * Njack);
    printf("writing: %s\n", namefile);
    for (int i = 0; i < 100; i++) {

        for (int j = 0;j < Njack;j++) {
            tmpx[0] = (double)params[0].data.L[1];//L
            tmpx[1] = x[j][0][1];//m0
            tmpx[2] = x[j][0][2];//m1
            tmpx[3] = x[j][0][3];//E20
            tmpx[4] = x[j][0][4];//E21
            tmpx[5] = params[0].data.L[1];//T
            tmpx[6] = 0 + i * 0.004;//k
            tmpx[7] = x[j][0][7];
            for (int i = 8; i < fit_info.Nvar; i++)
                tmpx[i] = x[j][0][i];

            for (int i = fit_info.Nvar; i < fit_info.Nvar + fit_info.n_ext_P; i++)
                tmpx[i] = fit_info.ext_P[i - fit_info.Nvar][j];

            tmpy[j] = fit_info.function(0, Nvar, tmpx, Npar, tif[j]);//N, Nvar, x ,Npar,P

        }
        fprintf(f, "%g  \t %g  %g   \n", tmpx[6], tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));

    }

    free(tmpy);free(tmpx);
    fclose(f);
    ////////// end fit band k
    free_2(Njack, tif);
}


void print_fit_band_k_m(char** argv, vector<data_all> gjack, struct fit_type fit_info, const char* label, struct fit_result fit_out, int* en, double*** x, double*** y, vector<cluster::IO_params> params, std::vector<int> myen) {

    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_k_m.txt", argv[3], label);
    f = open_file(namefile, "w+");
    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double* tmpx = (double*)malloc(sizeof(double) * Nvar);
    double* tmpy = (double*)malloc(sizeof(double) * Njack);
    printf("writing: %s\n", namefile);
    for (int i = 0; i < 100; i++) {

        for (int j = 0;j < Njack;j++) {
            tmpx[0] = (double)params[0].data.L[1];//L
            tmpx[1] = x[j][0][1];//m0
            tmpx[2] = x[j][0][2];//m1
            tmpx[3] = x[j][0][3];//E20
            tmpx[4] = x[j][0][4];//E21
            tmpx[5] = params[0].data.L[1];//T
            tmpx[6] = 0 + i * 0.004;//k
            tmpx[7] = x[j][0][7];
            tmpx[8] = x[j][0][8];
            tmpx[9] = x[j][0][9];

            tmpx[10] = 0 + i * 0.02;//k/m
            for (int i = 11; i < fit_info.Nvar; i++)
                tmpx[i] = x[j][0][i];

            for (int i = fit_info.Nvar; i < fit_info.Nvar + fit_info.n_ext_P; i++)
                tmpx[i] = fit_info.ext_P[i - fit_info.Nvar][j];

            tmpy[j] = fit_info.function(0, Nvar, tmpx, Npar, tif[j]);//N, Nvar, x ,Npar,P

        }
        fprintf(f, "%g  \t %g  %g   \n", tmpx[10], tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy));

    }

    free(tmpy);free(tmpx);
    fclose(f);
    ////////// end fit band k
    free_2(Njack, tif);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// print output
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_fit_output(char** argv, vector<data_all> gjack, struct fit_type fit_info, const char* label, struct fit_result fit_out, int* en, double*** x, double*** y, vector<cluster::IO_params> params, std::vector<int> myen) {
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;

    char namefile[NAMESIZE];
    FILE* f;

    mysprintf(namefile, NAMESIZE, "%s/%s_fit_data.txt", argv[3], label);
    printf("writing data in : %s\n", namefile);
    f = open_file(namefile, "w+");
    int count = 0;
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < en[n];e++) {
            fprintf(f, " %g   %g   %g   %d\t ", x[Njack - 1][e + count][0], y[Njack - 1][e + count][0], y[Njack - 1][e + count][1], params[myen[e]].data.L[0]);
            fprintf(f, " %g   \t ", x[Njack - 1][e + count][6]);
            fprintf(f, " %d   \t ", n);
            fprintf(f, " %g   \t ", x[Njack - 1][e + count][10]);
            fprintf(f, " %g   \n ", x[Njack - 1][e + count][1]);
        }
        count += en[n];
        fprintf(f, "\n\n");
    }

    fclose(f);
    ////////////// parameter print and correlation matrix
    mysprintf(namefile, NAMESIZE, "%s/%s_fit_P.tex", argv[3], label);
    f = open_file(namefile, "w+");
    fprintf(f, "\\begin{gather}\n");
    printf("chi^2/d.o.f.=%g  \n", fit_out.chi2[Njack - 1]);//error_jackboot(argv[1], Njack, fit_out.chi2)
    fprintf(f, "\\chi^2/d.o.f.=%g  \n", fit_out.chi2[Njack - 1]);//error_jackboot(argv[1], Njack, fit_out.chi2)
    for (int i = 0;i < Npar;i++) {
        fprintf(f, "P[%d]=%g\\pm (%.2g) \\\\ \n", i, fit_out.P[i][Njack - 1], error_jackboot(argv[1], Njack, fit_out.P[i]));
        printf("P[%d]=%g (%.2g) \n", i, fit_out.P[i][Njack - 1], error_jackboot(argv[1], Njack, fit_out.P[i]));
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

    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);

    //     if (strcmp(label,"k_from_phase_shift")!=0 && strcmp(label,"k_from_phase_shift_n4")!=0 && strcmp(label,"k_from_phase_shift_n5")!=0 && strcmp(label,"k_from_phase_shift_acotZ")!=0 && strcmp(label,"QC3")!=0 && strcmp(label,"deltaE2_m_quant_cond")!=0 ){
    if (strcmp(label, "M0_finite_vol") == 0 || strcmp(label, "M1_finite_vol") == 0 || strcmp(label, "M1_p1_finite_vol") == 0 || strcmp(label, "a_00_luscher") == 0
        || strcmp(label, "kcotd") == 0
        || strcmp(label, "kcotd_m") == 0 || strcmp(label, "kcotd_m_deltaE") == 0 || strcmp(label, "kcotd_Elatt") == 0
        || strcmp(label, "kcotd_2par") == 0 || strcmp(label, "kcotd_Elatt_2par") == 0 || strcmp(label, "kcotd_ECM_latt_2par") == 0
        || strcmp(label, "delta_2par") == 0 || strcmp(label, "delta_2par_g") == 0 || strcmp(label, "delta_Elatt_2par") == 0 || strcmp(label, "delta_ECM_latt_2par") == 0
        || strcmp(label, "a_00_luscher_infm") == 0 || strcmp(label, "a_01_luscher") == 0
        || strcmp(label, "a_01_luscher_div_shift") == 0 || strcmp(label, "a_00_BH") == 0 || strcmp(label, "a_01_BH_03t16_shifted") == 0
        || strcmp(label, "a_01_BH_04t16_shifted") == 0 || strcmp(label, "a_01_BH_02t10_shifted") == 0 || strcmp(label, "a_01_lusher_const") == 0 || strcmp(label, "fit_diff_L_BH02t16") == 0 || strcmp(label, "fit_diff_L_BH03t16") == 0
        || strcmp(label, "a_01_luscher_const") == 0
        || strcmp(label, "fit_diff_L_C2_m_C2shifted") == 0 || strcmp(label, "fit_L_C2") == 0 || strcmp(label, "fit_L_C2shifted") == 0
        || strcmp(label, "delta_2par_m1-4.705") == 0
        )
        //             || strcmp(label,"M1_finite_vol")==0  || strcmp(label,"a_00_luscher")==0  || strcmp(label,"kcotd")==0 
        //             || strcmp(label,"kcotd_Elatt")==0 || strcmp(label,"a_00_luscher_infm")==0 || strcmp(label,"a_01_luscher")==0 
        //             || strcmp(label,"a_01_luscher_div_shift")==0  || strcmp(label,"a_00_BH")==0 || strcmp(label,"a_01_BH_03t16_shifted")==0 || strcmp(label,"a_01_BH_04t16_shifted")==0 ,  strcmp(label,"a_01_BH_02t10_shifted")==0  || strcmp(label,"a_01_lusher_const")==0  )
    {

        /////////fit band L
        print_fit_band_L(argv, gjack, fit_info, label, fit_out, en, x, y, params, myen);
        print_fit_band_T(argv, gjack, fit_info, label, fit_out, en, x, y, params, myen);
        print_fit_band_k(argv, gjack, fit_info, label, fit_out, en, x, y, params, myen);
        print_fit_band_k_m(argv, gjack, fit_info, label, fit_out, en, x, y, params, myen);


    }

    free_2(Njack, tif);

}


struct fit_result fit_data(char** argv, vector<cluster::IO_params> params, vector<data_all> gjack,
    double lhs_fun(int, int, int, std::vector<cluster::IO_params>, std::vector<data_all>, struct fit_type),
    struct fit_type fit_info, const char* label, std::vector<int> myen) {
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar + fit_info.n_ext_P;
    int Njack = gjack[0].Njack;
    int N = fit_info.N;
    ////// allocation
    int* en = (int*)malloc(sizeof(int) * fit_info.N);// we need to init en and en_tot to allocate the other 
    for (int e = 0;e < fit_info.N; e++) { en[e] = myen.size(); }
    int en_tot = 0;      for (int n = 0;n < N;n++) { en_tot += en[n]; }// total data to fit

    double*** y = double_malloc_3(Njack, en_tot, 2);// 2 is mean/error
    double*** x = double_malloc_3(Njack, en_tot, Nvar);
    struct fit_result fit_out = malloc_fit(fit_info);
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
                x[j][count][0] = (double)params[myen[e]].data.L[1];//L
                x[j][count][1] = gjack[myen[e]].jack[1][j];//m0
                // if (params[myen[e]].data.gC>0 && gjack[0].Nobs>=443 && gjack[myen[e]].jack[443][j]!=0){
                //         x[j][count][1]=gjack[myen[e]].jack[443][j];//m0
                // }
                x[j][count][2] = gjack[myen[e]].jack[2][j];//m1
                x[j][count][3] = gjack[myen[e]].jack[4][j];//E20
                // if (params[myen[e]].data.gC>0 && gjack[0].Nobs>=594 && gjack[myen[e]].jack[594][j]!=0){
                //     x[j][count][3]=gjack[myen[e]].jack[594][j];//E20
                // }
                x[j][count][4] = gjack[myen[e]].jack[5][j];//E21
                x[j][count][5] = (double)params[myen[e]].data.L[0];//T

                x[j][count][10] = compute_k_m_g(n, myen[e], j, params, gjack, fit_info);//k/m
                x[j][count][6] = x[j][count][10] * x[j][count][1]; //k/m

                x[j][count][7] = gjack[myen[e]].jack[1][j] * (double)params[myen[e]].data.L[1] / (2. * pi_greco);//mL_2pi
                // if (params[myen[e]].data.gC>0 && gjack[0].Nobs>=443 && gjack[myen[e]].jack[443][j]!=0){
                //     x[j][count][7]=gjack[myen[e]].jack[443][j]*(double) params[myen[e]].data.L[1]/(2.*pi_greco);//mL_2pi
                // }

                // if (params[myen[e]].data.gC>0 && gjack[0].Nobs>=594 && gjack[myen[e]].jack[443][j]!=0){
                //     x[j][count][10]=compute_k_m_g(n,myen[e],j,params,gjack,fit_info);//k/m
                //     x[j][count][6]=x[j][count][10]*x[j][count][1];// k
                // }

                for (int i = 0; i < fit_info.n_ext_P; i++) {
                    x[j][count][i + fit_info.Nvar] = fit_info.ext_P[i][j];
                }
                //copy the other variables after 
                 //other var
                //x[j][count][fit_info.n_ext_P]=(double) params[e].data.L[1];

                count++;
            }
        }
    }
    int count = 0;
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < en[n];e++) {
            double* E3_m = (double*)malloc(sizeof(double) * Njack);
            for (int j = 0;j < Njack;j++) {
                E3_m[j] = lhs_E3_m(n, myen[e], j, params, gjack, fit_info);//E3( \vec{n} )/mass
                if (params[myen[e]].data.gC > 0)
                    E3_m[j] = lhs_E3orE1_m_complex(n, myen[e], j, params, gjack, fit_info);

            }

            double err = error_jackboot(argv[1], Njack, E3_m);
            for (int j = 0;j < Njack;j++) {
                x[j][count][8] = E3_m[Njack - 1];
                x[j][count][9] = err;

            }

            for (int j = 0;j < Njack;j++) {
                E3_m[j] = lhs_E3orE1_m_complex(n, myen[e], j, params, gjack, fit_info);//E3( \vec{n} )/mass
                if (params[myen[e]].data.gC > 0)
                    E3_m[j] = lhs_E3orE1_m_complex(n, myen[e], j, params, gjack, fit_info);

            }

            err = error_jackboot(argv[1], Njack, E3_m);
            for (int j = 0;j < Njack;j++) {
                x[j][count][11] = E3_m[Njack - 1];// E3 or E1 ???
                x[j][count][12] = err;


            }
            free(E3_m);


            count++;

        }
    }

    ////////////////////////////////////////// y
    count = 0;
    for (int n = 0;n < N;n++) {
        for (int e = 0;e < en[n];e++) {
            double* tmpj = (double*)malloc(sizeof(double) * Njack);
            for (int j = 0;j < Njack;j++) {
                //y[j][e+count][0]=gjack[e].jack[1][j];
                tmpj[j] = lhs_fun(n, myen[e], j, params, gjack, fit_info);
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
    print_fit_output(argv, gjack, fit_info, label, fit_out, en, x, y, params, myen);

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


#endif
