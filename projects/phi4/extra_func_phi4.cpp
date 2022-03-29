#define extra_func_phi4_C

#include "extra_func_phi4.hpp"
#include "mutils.hpp"
#include "tower.hpp"
#include "resampling.hpp"
extern "C" { 
    #include "../external/rzeta/src/dzeta_function.h"
}


#include "fit_all.hpp"


void print_fit_band_phi4(char** argv, data_all gjack, struct fit_type fit_info,
    struct fit_type fit_info_m0, const char* label, const char* dir_name,
    struct fit_result fit_out, struct fit_result fit_out_m0, int var, int en, double h) {

    std::vector<int> myen = fit_info.myen;
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
        for (int e = 1; e < myen.size() * N;e++) {
            if (fit_info.x[var][e][Njack - 1] < min)
                min = fit_info.x[var][e][Njack - 1];
            if (fit_info.x[var][e][Njack - 1] > max)
                max = fit_info.x[var][e][Njack - 1];
        }
        if (min = max)  printf("no variation in this direction\n");
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
                for (int i = 0; i < Nvar; i++) {
                    tmpx[i] = fit_info.x[i][en + n * myen.size()][j];
                }
                tmpx[var] = pos;
                tmpx[1] = fit_info_m0.function(n, Nvar, tmpx, Npar, tif_m0[j]);
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



void init_dvec_E2_g_extra(int n, int* dvec, int* dvec1, int* dvec2, int* dmax1, int* dmax2) {

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

double get_E2_g_n_extra(int n, int e, int j, data_all gjack) {
    double E2;
    if (n == 0) {//E2_0

        E2 = gjack.en[e].jack[4][j];

    }
    else if (n == 1) {//E2_p1

        E2 = gjack.en[e].jack[100][j];
    }
    else if (n == 2) {//E2_A1

        E2 = gjack.en[e].jack[80][j];
    }

    else { printf("%s n=%d not implemented\n", __func__, n);  exit(1); }
    return E2;

}

double energy_CM_extra(  double E  , int  *p,int L){
    
    double normp=p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
    normp*=2*pi_greco/((double) L);
    normp*=2*pi_greco/((double) L);
    return sqrt(E*E- normp);
}

inline double kcotd_extra(double E2, double mass, int* dvec, int L) {
    double E2_CM = energy_CM_extra(E2, dvec, L);
    double k = sqrt(E2_CM * E2_CM / 4. - mass * mass);

    double q = k * L / (2. * pi_greco);
    double gamma = E2 / E2_CM;
    double A = 1.;
    double z[2];
    dzeta_function(z, q * q, 0, 0, dvec, gamma, A, 1.e-3, 1.e+6, 5);

    return  z[0] * 2 * pi_greco / (pow(pi_greco, 3. / 2.) * gamma * L);
}



double lhs_kcotd_m_deltaE_g_new(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;

    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g_extra(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2 = get_E2_g_n_extra(n, e, j, gjack);
    // double mass = gjack.en[e].jack[1][j];
    double mass =fit_info.x[0][e][j];
    double L = gjack.en[e].header.L;

    double hatp2 = 4. * sin(dvec1[0] * pi_greco / L) * sin(dvec1[0] * pi_greco / L);
    hatp2 += 4. * sin(dvec1[1] * pi_greco / L) * sin(dvec1[1] * pi_greco / L);
    hatp2 += 4. * sin(dvec1[2] * pi_greco / L) * sin(dvec1[2] * pi_greco / L);

    double E2fL = acosh(cosh(mass) + 0.5 * (hatp2));
    hatp2 = 4. * sin(dvec2[0] * pi_greco / L) * sin(dvec2[0] * pi_greco / L);
    hatp2 += 4. * sin(dvec2[1] * pi_greco / L) * sin(dvec2[1] * pi_greco / L);
    hatp2 += 4. * sin(dvec2[2] * pi_greco / L) * sin(dvec2[2] * pi_greco / L);
    E2fL += acosh(cosh(mass) + 0.5 * (+hatp2));


    double Ef1 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec1[0] * dvec1[0] + dvec1[1] * dvec1[1] + dvec1[2] * dvec1[2]));
    double Ef2 = sqrt(mass * mass + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec2[0] * dvec2[0] + dvec2[1] * dvec2[1] + dvec2[2] * dvec2[2]));
    double E2f = Ef1 + Ef2;

    E2 = E2 - E2fL + E2f;

    r = kcotd_extra(E2, mass, dvec, L);
    
    return r / mass;
}


double rhs_kcotd_m_new(int n, int Nvar, double* x, int Npar, double* P) {
    double a0m = P[0], r0m = P[1];
    double k_m = x[1];
    double kcot_m = 1.0 / a0m + r0m * k_m * k_m / 2.;

    if (Npar > 2) {
        kcot_m = -P[2] * r0m * r0m * r0m * k_m * k_m * k_m * k_m;
    }

    return  kcot_m;

}


double compute_k_m_g_new(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double L = gjack.en[e].header.L;
    // double mass = gjack.en[e].jack[1][j];
    double mass = fit_info.x[0][e][j];

    
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g_extra(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2 = get_E2_g_n_extra(n, e, j, gjack);
    double E2_CM = energy_CM_extra(E2, dvec, L);

    return sqrt(E2_CM * E2_CM / 4. - mass * mass) / mass;

}


/////////////////////////////////////////////
//// QC3
////////////////////////////////////////////////


// void init_dvec_QC3_pole_new(int n, int* dvec) {
//     if (n == 0) {//E1_1
//         dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
//     }
//     else if (n == 1) {// E3_0
//         dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
//     }
//     else if (n == 2) {//E2_0_p11
//         dvec[0] = 1; dvec[1] = 0; dvec[2] = 0;
//     }
//     else if (n == 3) {//E2_0_p1
//         dvec[0] = 1; dvec[1] = 0; dvec[2] = 0;
//     }
//     else if (n == 4) {//E2_0_p1
//         dvec[0] = 1; dvec[1] = 1; dvec[2] = 0;
//     }
//     else if (n == 5) {//E2_0_p1
//         dvec[0] = 1; dvec[1] = 1; dvec[2] = 0;
//     }

//     else {
//         printf("init_dvec n=%d not implemented\n", n); exit(1);
//     }
// }

// #ifdef PYTHON
// double rhs_E3_m_QC3_pole_new(int n, int Nvar, double* x, int Npar, double* P) {

//     double Pkcot[2];
//         Pkcot[0] = x[1];
//         Pkcot[1] = x[2];

//     int Nkcot = 2;
//     int Nkiso = Npar;
//     int dvec[3];//,dvec1[3],dvec2[3],dmax1[3],dmax2[3];
//     init_dvec_QC3_pole_new(n, dvec);
//     //init_dvec(n,dvec,dvec1,dvec2,dmax1,dmax2);
//     double nnP[3];
//     nnP[0] = (double)dvec[0]; nnP[1] = (double)dvec[1]; nnP[2] = (double)dvec[2];
    
//     double Lm = x[0];// L* M
//     int steps = 4;

    

//     double Estart = x[3] - x[4]; // E3/m_inf +- error
//     double Eend = x[3] + x[4];

//     // printf("E=[%g , %g]   E1f=%g   E2f=%g  m=%g\n",Estart,Eend,Estart,Eend,mass);
//     double r = python_detQC_call(Estart, Eend, steps, Lm, nnP, Nkcot, Pkcot, Nkiso, P);
//     //     printf("res=%g\n",r);
//     return r;
// }
// #endif
