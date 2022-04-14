#define extra_func_phi4_C

#include "extra_func_phi4.hpp"
#include "mutils.hpp"
#include "tower.hpp"
#include "resampling.hpp"
#include "zeta_interpolation.hpp"
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

#ifdef PYTHON
void print_fit_band_QC3_phi4(char** argv, data_all gjack, struct fit_type fit_info,
    struct fit_type fit_info_E3_poly, const char* label, const char* dir_name,
    struct fit_result fit_out, struct fit_result fit_out_poly, int var, int en, double h) {

    std::vector<int> myen = fit_info.myen;
    int Npar = fit_info.Npar;
    int Nvar = fit_info.Nvar;
    int Njack = gjack.en[0].Njack;
    int N = fit_info.N;
    char namefile[NAMESIZE];
    FILE* f;

    double** tif = swap_indices(fit_info.Npar, Njack, fit_out.P);
    double** tif_poly = swap_indices(fit_info_E3_poly.Npar, Njack, fit_out_poly.P);
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

    double* tmpy = (double*)malloc(sizeof(double) * Njack);
    double* tmpx = (double*)malloc(sizeof(double) * Nvar);
    double* tmpx_poly = (double*)malloc(sizeof(double) * fit_info_E3_poly.Nvar);
    double* E3 = (double*)malloc(sizeof(double) * Njack);

    for (int n = 0;n < N; n++) {

        mysprintf(namefile, NAMESIZE, "%s/%s_fit_out_n%d_%s.txt", argv[3], label, n, dir_name);
        f = open_file(namefile, "w+");

        printf("writing: %s\n", namefile);
        // band range
        double pos = min;
        while (pos < max) {

            for (int j = 0;j < Njack;j++) {
                tmpx_poly[0] = pos;// gjack.fits[0].P[0][j];
                E3[j] = fit_info_E3_poly.function(n, fit_info_E3_poly.Nvar, tmpx_poly, fit_info_E3_poly.Npar, tif_poly[j]);
            }
            for (int j = 0;j < Njack;j++) {
                for (int i = 0; i < Nvar; i++) {
                    tmpx[i] = fit_info.x[i][en + n * myen.size()][j];
                }
                tmpx[var] = pos;
                tmpx[3] = E3[j];
                tmpx[4] = error_jackboot(argv[1], Njack, E3);
                // tmpx[1] = fit_info_m0.function(n, Nvar, tmpx, Npar, tif_m0[j]);

                tmpy[j] = fit_info.function(n, Nvar, tmpx, Npar, tif[j]);
            }
            fprintf(f, "%g  \t %g  %g   %g   %g\n", pos, tmpy[Njack - 1], error_jackboot(argv[1], Njack, tmpy), tmpx[3], tmpx[4]);
            pos += h;
        }

        fclose(f);
    }
    free(tmpx);free(tmpy);
    free(E3);free(tmpx_poly);
    free_2(Njack, tif);
    free_2(Njack, tif_poly);
}
#endif

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


void init_dvec_E2_with_E3(int n, int* dvec, int* dvec1, int* dvec2, int* dmax1, int* dmax2) {

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

double energy_CM_extra(double E, int* p, int L) {

    double normp = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
    normp *= 2 * pi_greco / ((double)L);
    normp *= 2 * pi_greco / ((double)L);
    return sqrt(E * E - normp);
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
    double mass = fit_info.x[0][e][j];
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


double lhs_kcotd_m_new(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;

    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g_extra(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2 = get_E2_g_n_extra(n, e, j, gjack);
    // double mass = gjack.en[e].jack[1][j];
    double mass = fit_info.x[0][e][j];
    double L = gjack.en[e].header.L;

    r = kcotd_extra(E2, mass, dvec, L);
    return r / mass;
}

double lhs_kcotd_minf_deltaE_g_new(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;

    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
    init_dvec_E2_g_extra(n, dvec, dvec1, dvec2, dmax1, dmax2);
    double E2 = get_E2_g_n_extra(n, e, j, gjack);
    // double mass = gjack.en[e].jack[1][j];
    double mass = fit_info.x[0][e][j];
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

    double mass_inf = fit_info.x[2][e][j];

    return r / mass_inf;
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


double to_invert_k_from_phase_shift_with_E3(int n, int Nvar, double* x, int Npar, double* P) {
    double a0m0 = P[0], r0m0 = P[1];//, P2=P[2];
    double L = x[0];
    double mass = x[1];
    double k = x[2];
    int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];

    init_dvec_E2_with_E3(n, dvec, dvec1, dvec2, dmax1, dmax2);

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
    // printf(" k=%g   f(k)=%g  kcotdelta_m=%g   r=%g  z=%g  gamma=%g   mass=%g  L=%g   dvec=(%d,%d,%d) E2CM=%g\n ", k, kcotdelta_m - r, kcotdelta_m, r, z[0], gamma, mass, L, dvec[0], dvec[1], dvec[2], E2_CM);
    // printf("P[0]=%g \t",P[0]);
    // if (Npar>=2) printf("P[1]=%g \t",P[1]);
    // if (Npar>=3) printf("P[2]=%g \t",P[2]);
    // printf("\n");
    return         kcotdelta_m - r;


}
/////////////////////////////////////////////
//// QC3
////////////////////////////////////////////////


void init_dvec_QC3_pole_new(int n, int* dvec) {
    if (n == 0) {//E1_1
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
    }
    else if (n == 1) {// E3_0
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
    }
    else if (n == 2) {//E2_0
        dvec[0] = 0; dvec[1] = 0; dvec[2] = 0;
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



double lhs_E3orE1_m_complex_new(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double E3;
    //double mass=gjack[e].jack[1][j];
    double mass = gjack.fits[0].P[0][j];


    if (n == 0) {//GEVP 1
        E3 = gjack.en[e].jack[354][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 1) {//GEVP 2
        E3 = gjack.en[e].jack[355][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else {
        E3 = 0;
        // printf("lhs_E3orE1_m n=%d not implemented\n",n); exit(1);
    }

    return E3 / mass;
}


double lhs_E3_E1_E2_m_complex_new(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double E3;
    double mass = gjack.en[e].jack[1][j];
    // double mass = gjack.fits[0].P[0][j];


    if (n == 0) {//GEVP 1
        E3 = gjack.en[e].jack[354][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 1) {//GEVP 2
        E3 = gjack.en[e].jack[355][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if (n == 2) {//GEVP 2
        E3 = gjack.en[e].jack[4][j];
        //         dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else {
        E3 = 0;
        // printf("lhs_E3orE1_m n=%d not implemented\n",n); exit(1);
    }

    return E3 / mass;
}


#ifdef PYTHON
double rhs_E3_m_QC3_pole_new(int n, int Nvar, double* x, int Npar, double* P) {

    double Pkcot[2];
    Pkcot[0] = x[1];
    Pkcot[1] = x[2];

    int Nkcot = 2;
    int Nkiso = Npar;
    int dvec[3];//,dvec1[3],dvec2[3],dmax1[3],dmax2[3];
    init_dvec_QC3_pole_new(n, dvec);
    //init_dvec(n,dvec,dvec1,dvec2,dmax1,dmax2);
    double nnP[3];
    nnP[0] = (double)dvec[0]; nnP[1] = (double)dvec[1]; nnP[2] = (double)dvec[2];

    double Lm = x[0];// L* M
    int steps = 4;
    // // when we load find_2sol
    double r = python_detQC_call(3.05, 2e-4, n, Lm, nnP, Nkcot, Pkcot, Nkiso, P);


    return r;
}

double rhs_E3_m_QC3_pole_E2_QC2(int n, int Nvar, double* x, int Npar, double* P) {

    double Pkcot[2];
    Pkcot[0] = P[3];
    Pkcot[1] = P[4];

    int Nkcot = 2;
    int Nkiso = Npar - Nkcot;
    int dvec[3];//,dvec1[3],dvec2[3],dmax1[3],dmax2[3];
    init_dvec_QC3_pole_new(n, dvec);
    //init_dvec(n,dvec,dvec1,dvec2,dmax1,dmax2);
    double nnP[3];
    nnP[0] = (double)dvec[0]; nnP[1] = (double)dvec[1]; nnP[2] = (double)dvec[2];

    double Lm = x[0];// L* M
    int steps = 4;
    double r;

    if (n < 2) {
        // // when we load find_2sol
        r = python_detQC_call(3.05, 2e-4, n, Lm, nnP, Nkcot, Pkcot, Nkiso, P);
    }
    else {
        int nr=n-2;
        double mass = x[1];
        double L = x[2];

        double xx[3] = { L,mass,1 }; //{L,mass,   k to be find by the bisection}


        int dvec[3], dvec1[3], dvec2[3], dmax1[3], dmax2[3];
        init_dvec_E2_with_E3(nr, dvec, dvec1, dvec2, dmax1, dmax2);
        

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

        double k = rtsafe(to_invert_k_from_phase_shift_with_E3, nr, 3, xx, Nkcot/*Npar*/, Pkcot, 2/*ivar*/, 0. /*input*/, xmin, xmax, 1e-5, 100, 1e-4);
        double E2 = sqrt((k * k + mass * mass) * 4. + (2 * pi_greco / L) * (2 * pi_greco / L) * (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]));

        // double mass_inf = x[0]/x[2];
        r = (E2 / mass);

    }

    return r;
}

#endif
