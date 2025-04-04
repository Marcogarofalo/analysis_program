#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <string>
#include <iostream>
#include <vector>

#include "global.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "read.hpp"
#include "mutils.hpp"
#include "mass_phi4.hpp"
#include "header_phi4.hpp"
#include "correlators_analysis.hpp"
#include "lhs_functions.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"

using namespace std;


int id_GEVP_031_p1;

struct  kinematic kinematic_2pt;

int r_value(int r) {
    int vr;
    if (r == 0) vr = 1;
    else if (r == 1) vr = -1;
    else { vr = -10; error(0 == 0, 0, "r_value", "r value is not 0 neither 1\n"); }
    return vr;
}



double C2(int n, int Nvar, double* x, int Npar, double* P) {

    double C2t;

    double E2 = P[0];
    double A2 = P[1];
    double A0 = P[2];
    double t = x[0];

    double T = (double)file_head.l0;
    C2t = A2 * A2 * exp(-E2 * T / 2.) * cosh(E2 * (t - T / 2.));
    C2t += A0 * A0;

    return C2t;

}


double C3(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double A3 = P[1];
    double A12 = P[2];
    double t = x[0];
    double M = x[1];
    double E2 = x[2];
    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C3 = A3 * A3 * exp(-E3 * T / 2.) * cosh(E3 * (t - T / 2.));
    C3 += A12 * A12 * exp(-(E2 + M) * T / 2.) * cosh((E2 - M) * (t - T / 2.));

    return C3;

}

double fun_l1_GEVP(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double E2 = P[1];
    //     double A10=P[3];
    double t = x[0];
    double M = x[1];

    if (x[0] == 3) return 1.0;

    double T = (double)file_head.l0;
    C3 = exp(-E3 * (t - 3.)) + exp(-E3 * (T - t - 3.));
    C3 += exp(-(E2 - M) * (t - 3.)) + exp(-(E2 - M) * (T - t - 3.));
    //     C3+=A10*A10 *  exp(-(M)*T/2.) * cosh( (M) *(t -T/2.));

    return C3;

}


double C3_vev(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double A3 = P[1];
    double A12 = P[2];
    double A10 = P[3];
    double t = x[0];
    double M = x[1];
    double E2 = x[2];
    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C3 = A3 * A3 * exp(-E3 * T / 2.) * cosh(E3 * (t - T / 2.));
    C3 += A12 * A12 * exp(-(E2 + M) * T / 2.) * cosh((E2 - M) * (t - T / 2.));
    C3 += A10 * A10 * exp(-(M)*T / 2.) * cosh((M) * (t - T / 2.));

    return C3;

}


double C3_2exp(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double A3 = P[1];
    double M = P[2];
    double A10 = P[3];
    double t = x[0];


    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C3 = A3 * A3 * exp(-E3 * T / 2.) * cosh(E3 * (t - T / 2.));
    C3 += A10 * A10 * exp(-(M)*T / 2.) * cosh((M) * (t - T / 2.));

    return C3;

}

double C3_vev_2par(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double A3 = P[1];
    double A10 = P[2];
    double t = x[0];
    double M = x[1];

    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C3 = A3 * A3 * exp(-E3 * T / 2.) * cosh(E3 * (t - T / 2.));
    C3 += A10 * A10 * exp(-(M)*T / 2.) * cosh((M) * (t - T / 2.));

    return C3;

}

double C3p(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double A3 = P[1];
    double A10_2p = P[2];
    double A1p_20 = P[3];
    double t = x[0];
    double E1 = x[1];
    double E2 = x[2];
    double E1p = x[3];
    double E2p = x[4];

    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C3 = A3 * A3 * exp(-E3 * T / 2.) * cosh(E3 * (t - T / 2.));
    C3 += A10_2p * A10_2p * exp(-(E2p + E1) * T / 2.) * cosh((E2p - E1) * (t - T / 2.));
    C3 += A1p_20 * A1p_20 * exp(-(E2 + E1p) * T / 2.) * cosh((E2 - E1p) * (t - T / 2.));

    return C3;

}

double C3p_vev(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double A3 = P[1];
    double A10_2p = P[2];
    double A1p_20 = P[3];
    double A1p_0 = P[4];
    double t = x[0];
    double E1 = x[1];
    double E2 = x[2];
    double E1p = x[3];
    double E2p = x[4];
    double e1 = x[5]; // energy of one particle if the when the other is zero

    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C3 = A3 * A3 * exp(-E3 * T / 2.) * cosh(E3 * (t - T / 2.));
    C3 += A10_2p * A10_2p * exp(-(E2p + E1) * T / 2.) * cosh((E2p - E1) * (t - T / 2.));
    C3 += A1p_20 * A1p_20 * exp(-(E2 + E1p) * T / 2.) * cosh((E2 - E1p) * (t - T / 2.));
    C3 += A1p_0 * A1p_0 * exp(-(e1)*T / 2.) * cosh((e1) * (t - T / 2.));
    return C3;

}



double C3_A1(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double A3 = P[1];
    double A10_2p = P[2];
    double A1p_20 = P[3];
    double A10_20 = P[4];
    //     double DE=P[5];
    double t = x[0];
    double E1 = x[1];
    double E2p = x[2];
    double E1p = x[3];
    double E2A1 = x[4];
    double E2 = x[5];

    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C3 = A3 * A3 * exp(-E3 * T / 2.) * cosh(E3 * (t - T / 2.));
    C3 += A10_2p * A10_2p * exp(-(E2A1 + E1) * T / 2.) * cosh((E2A1 - E1) * (t - T / 2.));
    C3 += A1p_20 * A1p_20 * exp(-(E2p + E1p) * T / 2.) * cosh((E2p - E1p) * (t - T / 2.));
    C3 += A10_20 * A10_20 * exp(-(E1)*T / 2.) * cosh((E1) * (t - T / 2.));
    //    C3+=A10_20*A10_20 *  exp(-(E2+E1)*T/2.) *cosh( (E2-E1) *(t -T/2.));//* cosh( (0.12591+0.004) *(t -T/2.));//

    //     C3+=A10_20*A10_20 *  exp(-(E2+E1)*T/2.) * cosh( (DE) *(t -T/2.));
    //        printf(" E3=%g   E1=%g   E1p=%g   E2=%g   E2p=%g   E2_A1=%g  E2-E1=%g\n", E3,E1,E1p,E2,E2p,E2A1,  E2-E1);  
    return C3;

}


double C3_A1_4par(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double A3 = P[1];
    double A1p_20 = P[2];
    double A10_20 = P[3];
    //     double DE=P[5];
    double t = x[0];
    double E1 = x[1];
    double E2p = x[2];
    double E1p = x[3];
    double E2A1 = x[4];
    double E2 = x[5];

    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C3 = A3 * A3 * exp(-E3 * T / 2.) * cosh(E3 * (t - T / 2.));
    C3 += A1p_20 * A1p_20 * exp(-(E2p + E1p) * T / 2.) * cosh((E2p - E1p) * (t - T / 2.));
    C3 += A10_20 * A10_20 * exp(-(E1)*T / 2.) * cosh((E1) * (t - T / 2.));

    return C3;

}


double C3_A1_3par(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double A3 = P[1];
    double A10_20 = P[2];
    //     double DE=P[5];
    double t = x[0];
    double E1 = x[1];
    double E2p = x[2];
    double E1p = x[3];
    double E2A1 = x[4];
    double E2 = x[5];

    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C3 = A3 * A3 * exp(-E3 * T / 2.) * cosh(E3 * (t - T / 2.));
    C3 += A10_20 * A10_20 * exp(-(E1)*T / 2.) * cosh((E1) * (t - T / 2.));

    return C3;

}


double C3_A1_old(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double E3 = P[0];
    double A3 = P[1];
    double A10_2p = P[2];
    double A1p_20 = P[3];
    //     double A10_20=P[4];
        //     double DE=P[5];
    double t = x[0];
    double E1 = x[1];
    double E2p = x[2];
    double E1p = x[3];
    double E2A1 = x[4];
    double E2 = x[5];

    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C3 = A3 * A3 * exp(-E3 * T / 2.) * cosh(E3 * (t - T / 2.));
    C3 += A10_2p * A10_2p * exp(-(E2A1 + E1) * T / 2.) * cosh((E2A1 - E1) * (t - T / 2.));
    C3 += A1p_20 * A1p_20 * exp(-(E2p + E1p) * T / 2.) * cosh((E2p - E1p) * (t - T / 2.));
    //     C3+=A10_20*A10_20 *  exp(-(E1)*T/2.) *cosh( (E1) *(t -T/2.));
        //    C3+=A10_20*A10_20 *  exp(-(E2+E1)*T/2.) *cosh( (E2-E1) *(t -T/2.));//* cosh( (0.12591+0.004) *(t -T/2.));//

        //     C3+=A10_20*A10_20 *  exp(-(E2+E1)*T/2.) * cosh( (DE) *(t -T/2.));
        //        printf(" E3=%g   E1=%g   E1p=%g   E2=%g   E2p=%g   E2_A1=%g  E2-E1=%g\n", E3,E1,E1p,E2,E2p,E2A1,  E2-E1);  
    return C3;

}

double me_k3pi_rhs_T_2(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double A = P[0];
    double t = x[0];
    double E1_1 = x[1];
    double E3_0 = x[2];

    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    double tf = T / 2;
    C3 = A * exp(-E1_1 * (tf - t) / 2. - E3_0 * t / 2);

    return C3;


}
template<int tf>
double me_k3pi_rhs(int n, int Nvar, double* x, int Npar, double* P) {

    double C3;

    double A = P[0];
    double t = x[0];
    double E1_1 = x[1];
    double E3_0 = x[2];

    //check if file_head.l0 arrives here
    //     double T=(double)file_head.l0;

    C3 = A * exp(-E1_1 * (tf - t) / 2. - E3_0 * t / 2.);
    return C3;


}

double me_3pik_rhs_T_2(int n, int Nvar, double* x, int Npar, double* P) {


    double C3;

    double A = P[0];
    double A1 = P[1];
    double A2 = P[2];
    double t = x[0];
    double E1_1 = x[1];
    double E3_0 = x[2];
    double E1_0 = x[3];
    double E2_01 = x[4];

    double me_phi1 = x[5];
    double me_phi3 = x[6];

    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    double tf = T / 2.;

    C3 = A * ((4 * E3_0 * E1_1) / (me_phi1 * me_phi3)) * exp(-E3_0 * (tf - t) - E1_1 * t);
    C3 += A1 * exp(-E1_0 * tf + E1_0 * t - E1_1 * t);
    C3 += A2 * exp(-E1_0 * T + E1_0 * tf - E2_01 * t);

    return C3;

}
template<int tf>
double me_3pik_rhs(int n, int Nvar, double* x, int Npar, double* P) {


    double C3;

    double A = P[0];
    double A1 = P[1];
    double A2 = P[2];
    double t = x[0];
    double E1_1 = x[1];
    double E3_0 = x[2];
    double E1_0 = x[3];
    double E2_01 = x[4];

    double me_phi1 = x[5];
    double me_phi3 = x[6];

    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;

    C3 = A * ((4 * E3_0 * E1_1) / (me_phi1 * me_phi3)) * exp(-E3_0 * (tf - t) - E1_1 * t);
    C3 += A1 * exp(-E1_0 * tf + E1_0 * t - E1_1 * t);
    C3 += A2 * exp(-E1_0 * T + E1_0 * tf - E2_01 * t);


    return C3;
}

double C2_diff_masses(int n, int Nvar, double* x, int Npar, double* P) {

    double C2;

    double E2 = P[0];
    double A2 = P[1];
    double A12 = P[2];

    double t = x[0];
    double M0 = x[1];
    double M1 = x[2];
    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;
    C2 = A2 * A2 * (exp(-t * E2) + exp(-(T - t) * E2));
    C2 += A12 * A12 * (exp(-T * M0 - t * M1 + t * M0) + exp(-T * M1 - t * M0 + t * M1));

    return C2;

}

double C2_diff_masses_weight_shift(int n, int Nvar, double* x, int Npar, double* P) {

    double C2;

    double E2 = P[0];
    double A2 = P[1];
    //double A12=P[2];

    double t = x[0];
    double tp = t + 1;
    double M0 = x[1];
    double M1 = x[2];
    //check if file_head.l0 arrives here
    double T = (double)file_head.l0;

    C2 = A2 * A2 * (exp(-tp * E2) + exp(-(T - tp) * E2)) / (exp(-T * M0 - tp * M1 + tp * M0) + exp(-T * M1 - tp * M0 + tp * M1));
    C2 -= A2 * A2 * (exp(-t * E2) + exp(-(T - t) * E2)) / (exp(-T * M0 - t * M1 + t * M0) + exp(-T * M1 - t * M0 + t * M1));



    return C2;

}


double four_pt_BH_t_T_8(int n, int Nvar, double* x, int Npar, double* P) {

    double C4;

    double aN = P[0];
    //double norm=P[1];


    double T = (double)file_head.l0;
    double t = x[0] - T / 8.;
    double M0 = x[1];
    double M1 = x[2];

    C4 = 8. * pi_greco * (M0 + M1) * aN * t;
    C4 -= 16 * aN * aN * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1 * t);
    //C4 *= norm*exp(M1*t) ;
    C4 /= 8 * M0 * M1 * t;

    return C4;

}
double four_pt_BH_t_t3(int n, int Nvar, double* x, int Npar, double* P) {
    double C4;
    double aN = P[0];
    double T = (double)file_head.l0;
    double t = x[0] - 3.;
    double M0 = x[1];
    double M1 = x[2];

    C4 = 8. * pi_greco * (M0 + M1) * aN * t;
    C4 -= 16 * aN * aN * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1 * t);
    C4 /= 8 * M0 * M1 * t;

    return C4;
}

double four_pt_BH_t_t4(int n, int Nvar, double* x, int Npar, double* P) {
    double C4;
    double aN = P[0];
    double T = (double)file_head.l0;
    double t = x[0] - 4.;
    double M0 = x[1];
    double M1 = x[2];

    C4 = 8. * pi_greco * (M0 + M1) * aN * t;
    C4 -= 16 * aN * aN * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1 * t);
    C4 /= 8 * M0 * M1 * t;

    return C4;
}
double four_pt_BH_t_t5(int n, int Nvar, double* x, int Npar, double* P) {
    double C4;
    double aN = P[0];
    double T = (double)file_head.l0;
    double t = x[0] - 5.;
    double M0 = x[1];
    double M1 = x[2];

    C4 = 8. * pi_greco * (M0 + M1) * aN * t;
    C4 -= 16 * aN * aN * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1 * t);
    C4 /= 8 * M0 * M1 * t;

    return C4;
}
double four_pt_BH_00_t_t3(int n, int Nvar, double* x, int Npar, double* P) {
    double C4;
    double aN = P[0];
    double T = (double)file_head.l0;
    double t = x[0] - 3.;
    double M0 = x[1];
    double M1 = x[2];

    C4 = 16. * pi_greco * (M0 + M1) * aN * t;
    C4 -= 16 * aN * aN * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1 * t);
    C4 /= 8 * M0 * M1 * t;

    return C4;
}

double four_pt_BH_00_t_t4(int n, int Nvar, double* x, int Npar, double* P) {
    double C4;
    double aN = P[0];
    double T = (double)file_head.l0;
    double t = x[0] - 4.;
    double M0 = x[1];
    double M1 = x[2];

    C4 = 16. * pi_greco * (M0 + M1) * aN * t;
    C4 -= 16 * aN * aN * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1 * t);
    C4 /= 8 * M0 * M1 * t;

    return C4;
}
double four_pt_BH_00_t_t5(int n, int Nvar, double* x, int Npar, double* P) {
    double C4;
    double aN = P[0];
    double T = (double)file_head.l0;
    double t = x[0] - 5.;
    double M0 = x[1];
    double M1 = x[2];

    C4 = 16. * pi_greco * (M0 + M1) * aN * t;
    C4 -= 16 * aN * aN * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1 * t);
    C4 /= 8 * M0 * M1 * t;

    return C4;
}
double four_pt_BH_t_T_8_const(int n, int Nvar, double* x, int Npar, double* P) {
    double r = four_pt_BH_t_T_8(n, Nvar, x, Npar, P);
    return r + P[1] / (x[0] - file_head.l0 / 8.);
}

double four_pt_BH_t_t3_const(int n, int Nvar, double* x, int Npar, double* P) {
    double r = four_pt_BH_t_t3(n, Nvar, x, Npar, P);
    return r + P[1] / (x[0] - 3.);
}

double four_pt_BH_t_t4_const(int n, int Nvar, double* x, int Npar, double* P) {
    double r = four_pt_BH_t_t4(n, Nvar, x, Npar, P);
    return r + P[1] / (x[0] - 4.);
}
double four_pt_BH_t_t5_const(int n, int Nvar, double* x, int Npar, double* P) {
    double r = four_pt_BH_t_t5(n, Nvar, x, Npar, P);
    return r + P[1] / (x[0] - 5.);
}

double four_pt_BH_00_t_t3_const(int n, int Nvar, double* x, int Npar, double* P) {
    double r = four_pt_BH_00_t_t3(n, Nvar, x, Npar, P);
    return r + P[1] / (x[0] - 3.);
}
double four_pt_BH_00_t_t4_const(int n, int Nvar, double* x, int Npar, double* P) {
    double r = four_pt_BH_00_t_t4(n, Nvar, x, Npar, P);
    return r + P[1] / (x[0] - 4.);
}
double four_pt_BH_00_t_t5_const(int n, int Nvar, double* x, int Npar, double* P) {
    double r = four_pt_BH_00_t_t5(n, Nvar, x, Npar, P);
    return r + P[1] / (x[0] - 5.);
}


template<int tx, int delta  >
double four_pt_BH_t_tx_shifted(int n, int Nvar, double* x, int Npar, double* P) {
    double C4;
    double aN = P[0];
    double T = (double)file_head.l0;
    double t = x[0] - tx;
    double M0 = x[1];
    double M1 = x[2];

    C4 = 8. * pi_greco * (M0 + M1) * aN;
    C4 -= 16 * aN * aN * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1) * (sqrt(t + delta) - sqrt(t)) / delta;
    C4 /= 8 * M0 * M1;

    return C4;
}

template<int tx, int delta  >
double four_pt_BH_00_t_tx_shifted(int n, int Nvar, double* x, int Npar, double* P) {
    double C4;
    double aN = P[0];
    double T = (double)file_head.l0;
    double t = x[0] - tx;
    double M0 = x[1];
    double M1 = x[2];

    C4 = 16. * pi_greco * (M0 + M1) * aN;
    C4 -= 16 * aN * aN * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1) * (sqrt(t + delta) - sqrt(t)) / delta;
    C4 /= 8 * M0 * M1;

    return C4;
}


double four_pt_BH_line(int n, int Nvar, double* x, int Npar, double* P) {

    double C4;

    double aN = P[0];
    //double norm=P[1];


    double T = (double)file_head.l0;
    double t = x[0] - T / 8.;
    double M0 = x[1];
    double M1 = x[2];

    C4 = 8. * pi_greco * (M0 + M1) * aN * t;
    C4 /= 8 * M0 * M1 * t;
    return C4;

}

double four_pt_BH_2par(int n, int Nvar, double* x, int Npar, double* P) {

    double C4;

    double A = P[0];
    double B = P[1];


    double T = (double)file_head.l0;
    double t = x[0] - T / 8.;
    double M0 = x[1];
    double M1 = x[2];

    C4 = A * 8. * pi_greco * (M0 + M1) * t;
    C4 += B * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1 * t);
    //C4 *= norm*exp(M1*t) ;
    //C4 *= norm ;
    C4 /= 8 * M0 * M1 * t;
    return C4;

}


double four_pt_BH_3par(int n, int Nvar, double* x, int Npar, double* P) {

    double C4;

    double A = P[0];
    double B = P[1];
    double C = P[2];


    double T = (double)file_head.l0;
    double t = x[0] - T / 8.;
    double M0 = x[1];
    double M1 = x[2];

    C4 = A * t;
    C4 += B * sqrt(2. * pi_greco * (M0 + M1) * M0 * M1 * t);
    C4 += C;
    //C4 *= norm*exp(M1*t) ;
    //C4 *= norm ;
    C4 /= 8 * M0 * M1 * t;
    return C4;

}



void get_kinematic(int ik2, int r2, int ik1, int r1, int imom2, int imom1) {
    kinematic_2pt.ik2 = ik2;
    kinematic_2pt.ik1 = ik1;
    kinematic_2pt.k2 = file_head.k[ik2 + file_head.nk];
    kinematic_2pt.k1 = file_head.k[ik1 + file_head.nk];
    kinematic_2pt.r2 = -r_value(r2);
    kinematic_2pt.r1 = r_value(r1);
    kinematic_2pt.mom2 = -file_head.mom[imom2][1];
    kinematic_2pt.mom1 = file_head.mom[imom1][1];

    kinematic_2pt.mom02 = file_head.mom[imom2][0];
    kinematic_2pt.mom01 = file_head.mom[imom1][0];

    kinematic_2pt.line = 1;
    printf("%g %g  %d  %d\n", kinematic_2pt.k1, kinematic_2pt.k2, ik1 + file_head.nk, ik2 + file_head.nk);

}


static int**** mass_index;

void init_mass_index() {
    int k1, k2, r1, r2, i;
    int nk = file_head.nk;

    mass_index = (int****)malloc(sizeof(int***) * nk);
    for (k1 = 0;k1 < nk;k1++) {
        mass_index[k1] = (int***)malloc(sizeof(int**) * 2);
        for (r1 = 0;r1 < 2;r1++) {
            mass_index[k1][r1] = (int**)malloc(sizeof(int*) * (k1 + 1));
            for (k2 = 0;k2 <= k1;k2++) {
                mass_index[k1][r1][k2] = (int*)malloc(sizeof(int) * 2);
            }
        }
    }

    i = 0;
    for (k1 = 0;k1 < nk;k1++)
        for (r1 = 0;r1 < 2;r1++)
            for (k2 = 0;k2 <= k1;k2++)
                for (r2 = 0;r2 < 2;r2++) {
                    mass_index[k1][r1][k2][r2] = i;
                    i++;
                }


}

static void  print_file_head(FILE* stream) {
    int i;

    fprintf(stream, "twist= %d\n", file_head.twist);
    fprintf(stream, "nf=%d\n", file_head.nf);
    fprintf(stream, "nsrc=%d\n", file_head.nsrc);
    fprintf(stream, "L0=%d\n", file_head.l0);
    fprintf(stream, "L1=%d\n", file_head.l1);
    fprintf(stream, "L2=%d\n", file_head.l2);
    fprintf(stream, "L3=%d\n", file_head.l3);
    fprintf(stream, "mus=%d\n", file_head.nk);
    fprintf(stream, "moms=%d\n", file_head.nmoms);

    fprintf(stream, "beta=%f\n", file_head.beta);
    fprintf(stream, "ksea=%f\n", file_head.ksea);
    fprintf(stream, "musea=%f\n", file_head.musea);
    fprintf(stream, "csw=%f\n", file_head.csw);

    fprintf(stream, "masses=");
    for (i = 0;i < 2 * file_head.nk;++i)
        fprintf(stream, "%f\t", file_head.k[i]);
    fprintf(stream, "\n");

    fprintf(stream, "momenta=");
    for (i = 0;i < file_head.nmoms;++i)
        fprintf(stream, "%f  %f   %f   %f\n", file_head.mom[i][0], file_head.mom[i][1], file_head.mom[i][2], file_head.mom[i][3]);
}



int read_nconfs(FILE* stream, cluster::IO_params params) {

    long int tmp;
    int s;


    std::cout << "header size=" << params.data.header_size << std::endl;
    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    error(tmp == -1 || tmp == 0, 1, "read_nconfs", "ftell returned -1");
    std::cout << "ftell= " << ftell(stream) << std::endl;
    std::cout << "stored= " << tmp << std::endl;
    tmp -= params.data.header_size;

    s = params.data.size;
    std::cout << "size=" << s << std::endl;

    int c = (tmp) / (sizeof(int) + (s) * sizeof(double));


    std::cout << "confs=" << c << std::endl;
    fseek(stream, params.data.header_size, SEEK_SET);

    return c;

}
static void  write_file_head(FILE* stream) {
    int i, dsize;
    double* dstd;

    fwrite(&file_head.twist, sizeof(int), 1, stream);
    fwrite(&file_head.nf, sizeof(int), 1, stream);
    fwrite(&file_head.nsrc, sizeof(int), 1, stream);
    fwrite(&file_head.l0, sizeof(int), 1, stream);
    fwrite(&file_head.l1, sizeof(int), 1, stream);
    fwrite(&file_head.l2, sizeof(int), 1, stream);
    fwrite(&file_head.l3, sizeof(int), 1, stream);
    fwrite(&file_head.nk, sizeof(int), 1, stream);
    fwrite(&file_head.nmoms, sizeof(int), 1, stream);

    fwrite(&file_head.beta, sizeof(double), 1, stream);
    fwrite(&file_head.ksea, sizeof(double), 1, stream);
    fwrite(&file_head.musea, sizeof(double), 1, stream);
    fwrite(&file_head.csw, sizeof(double), 1, stream);

    fwrite(file_head.k, sizeof(double), 2 * file_head.nk, stream);

    for (i = 0;i < file_head.nmoms;i++)
        fwrite(file_head.mom[i], sizeof(double), 4, stream);
}



double matrix_element_GEVP(int t, double** cor, double mass) {
    double me;

    me = cor[t][0] / sqrt(exp(-mass * t) + exp(-(file_head.l0 - t) * mass));
    me *= 2 * mass;

    return  me;
}




void read_twopt(FILE* stream, int iconf, double*** to_write, cluster::IO_params params, int index) {

    int tmp = params.data.header_size;// 
    tmp += sizeof(double) * iconf * params.data.size + sizeof(int) * (iconf + 1);


    double* obs = (double*)malloc(params.data.size * sizeof(double));

    fseek(stream, tmp, SEEK_SET);
    size_t i = fread(obs, sizeof(double), params.data.size, stream);

    for (int t = 0;t < params.data.L[0];t++) {
        size_t  id = index + t * params.data.ncorr;
        (*to_write)[t][0] = obs[id];

    }
    free(obs);


}



void setup_single_file_jack(char* save_name, char** argv, const char* name, int Njack) {
    FILE* f;
    mysprintf(save_name, NAMESIZE, "/dev/null");
    f = fopen(save_name, "w+");
    error(f == NULL, 1, "setup_file_jack ",
        "Unable to open output file /dev/null");
    write_file_head(f);
    fwrite(&Njack, sizeof(int), 1, f);
    fclose(f);
}

void setup_single_file_jack_ASCI(char* save_name, char** argv, const char* name, int Njack) {
    FILE* f;
    mysprintf(save_name, NAMESIZE, "/dev/null");
    f = fopen(save_name, "w+");
    error(f == NULL, 1, "setup_file_jack ",
        "Unable to open output file /dev/null");
    fclose(f);
}

void setup_file_jack(char** argv, int Njack) {
    if (strcmp(argv[4], "jack") == 0) {
        setup_single_file_jack(file_jack.M_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.Zf_PS, argv, "/dev/null", Njack);

        setup_single_file_jack(file_jack.M_PS_GEVP, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS_ls_ss, argv, "/dev/null", Njack);

    }

    if (strcmp(argv[4], "boot") == 0) {

        setup_single_file_jack(file_jack.M_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.Zf_PS, argv, "/dev/null", Njack);

        setup_single_file_jack(file_jack.M_PS_GEVP, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS_ls_ss, argv, "/dev/null", Njack);
    }
}

int main(int argc, char** argv) {
    int size;
    int i, j, t;

    int* iconf, confs;
    double**** data, **** data_bin, ** out, ** tmp;
    char c;
    double* in;

    double**** M, **** vec, **** projected_O;
    double**** lambda, **** lambda0;

    double* fit, *** y, * x, * m, * me;


    double**** conf_jack, ** r, ** mt, ** met;
    int Ncorr = 1;
    int t0 = 2;

    FILE* f_ll = NULL, * f_sl = NULL, * f_ls = NULL, * f_ss = NULL;

    FILE* plateaux_masses = NULL, * plateaux_masses_GEVP = NULL;
    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");
    FILE* plateaux_f = NULL;
    char namefile[NAMESIZE];
    srand(1);



    error(argc != 8, 1, "main ",
        "usage:./phi4  blind/see/read_plateaux -p path file -bin $bin  jack/boot \n separate path and file please");
    error(strcmp(argv[1], "blind") != 0 && strcmp(argv[1], "see") != 0 && strcmp(argv[1], "read_plateaux") != 0, 1, "main ",
        "argv[1] only options:  blind/see/read_plateaux ");

    cluster::IO_params params;
    mysprintf(namefile, NAMESIZE, "%s/%s", argv[3], argv[4]);
    FILE* infile = open_file(namefile, "r+");
    read_header_phi4(infile, params);




    error(strcmp(argv[5], "-bin") != 0, 1, "main", "argv[4] must be: -bin");
    error(strcmp(argv[7], "jack") != 0 && strcmp(argv[7], "boot") != 0, 1, "main",
        "argv[6] only options: jack/boot");


    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, argv[1]); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p"); // -p
    mysprintf(option[3], NAMESIZE, argv[3]); // path
    mysprintf(option[4], NAMESIZE, argv[7]); //resampling
    mysprintf(option[5], NAMESIZE, "no"); // pdf
    mysprintf(option[6], NAMESIZE, argv[4]); // infile

    printf("resampling %s\n", option[4]);
    int T = params.data.L[0];

    double mu1 = params.data.msq0;
    double mu2 = params.data.msq1;
    printf("mu=%g  %g\n", mu1, mu2);
    //mysprintf(argv[4],NAMESIZE,"jack");

    file_head.l0 = T;
    file_head.l1 = params.data.L[1];file_head.l2 = params.data.L[2];file_head.l3 = params.data.L[3];
    file_head.nk = 2;
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;file_head.k[1] = 0;
    file_head.k[2] = mu1;
    file_head.k[3] = mu2;

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (i = 0;i < file_head.nmoms;i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }


    mysprintf(namefile, NAMESIZE, "%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_output",
        argv[3], T, params.data.L[1], params.data.msq0, params.data.msq1,
        params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC, params.data.replica);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_meff_correlators",
        argv[3], T, params.data.L[1], params.data.msq0, params.data.msq1,
        params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC, params.data.replica);
    FILE* outfile_meff_corr = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_raw_correlators",
        argv[3], T, params.data.L[1], params.data.msq0, params.data.msq1,
        params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC, params.data.replica);
    FILE* outfile_raw_corr = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_shifted_correlators",
        argv[3], T, params.data.L[1], params.data.msq0, params.data.msq1,
        params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC, params.data.replica);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");


    mysprintf(namefile, NAMESIZE, "%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_log_meff_shifted",
        argv[3], T, params.data.L[1], params.data.msq0, params.data.msq1,
        params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC, params.data.replica);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_gamma",
        argv[3], T, params.data.L[1], params.data.msq0, params.data.msq1,
        params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC, params.data.replica);

    FILE* out_gamma = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d",
        argv[3], option[4],
        T, params.data.L[1], params.data.msq0, params.data.msq1,
        params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC, params.data.replica);

    FILE* jack_file = open_file(namefile, "w+");
    write_header_phi4(jack_file, params);

    // open infile and count the lines
    //


    int count = 0;
    confs = read_nconfs(infile, params);
    //confs=confs/10;
    cout << "correlators =" << params.data.ncorr << endl;
    // compute what will be the neff after the binning 
    int bin = atoi(argv[6]);
    int Neff = confs / bin;
    cout << "effective configurations after binning (" << bin << "):  " << Neff << endl;

    int Njack;
    if (strcmp(argv[7], "jack") == 0) {
        Njack = Neff +1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[7], "boot") == 0) {
        Njack = (Neff*2 + 1);
        myres = new resampling_boot(Neff*2);
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[7]= %s is not jack or boot", argv[7]);
    }
    fwrite(&Njack, sizeof(int), 1, jack_file);

    int var = params.data.ncorr;
    data = calloc_corr(confs, var, file_head.l0);

    setup_file_jack(option, Njack);

    get_kinematic(0, 0, 1, 0, 0, 0);
    printf("option[4]=%s\n", option[4]);

    for (int iconf = 0; iconf < confs;iconf++) {

        for (int i = 0; i < params.data.ncorr; i++)
            read_twopt(infile, iconf, &data[iconf][i], params, i);



        //read_twopt(infile, iconf, &data[iconf][0], params,0);//2pt 0
        //read_twopt(infile, iconf, &data[iconf][1], params,1);//2pt 1

        //read_twopt(infile, iconf, &data[iconf][2], params,2);//C2t0
        //read_twopt(infile, iconf, &data[iconf][3], params,3);//C2t1
        //read_twopt(infile, iconf, &data[iconf][4], params,4);//C2t

        //read_twopt(infile, iconf, &data[iconf][5], params,5);//C3t0
        //read_twopt(infile, iconf, &data[iconf][6], params,6);//C3t1
        //read_twopt(infile, iconf, &data[iconf][7], params,7);//C3t

        //read_twopt(infile, iconf, &data[iconf][8], params,8);//C4t0
        //read_twopt(infile, iconf, &data[iconf][9], params,9);//C4t1
        //read_twopt(infile, iconf, &data[iconf][10], params,10);//C4t


        //read_twopt(infile, iconf, &data[iconf][11], params,11);//C201

    }
    /*
        for(int i =0 ; i< 2; i++){
            std::string myname="one_to_one_"+ std::to_string(i) + ".txt";
            FILE *fi=open_file(myname.c_str() ,"w+" );
            for (int iconf=0; iconf< confs ;iconf++){
                for(int t=0;t<params.data.L[0];t++)
                    fprintf(fi,"%d   %d  %.12g\n",iconf,t,data[iconf][i][t][0]);
            }
            fclose(fi);
        }
      */

      // symmetrise_corr(confs, 0, file_head.l0,data);
      // symmetrise_corr(confs, 1, file_head.l0,data);

      // symmetrise_corr(confs, 2, file_head.l0,data);
      // symmetrise_corr(confs, 3, file_head.l0,data);
      // symmetrise_corr(confs, 4, file_head.l0,data);

      // symmetrise_corr(confs, 5, file_head.l0,data);
      // symmetrise_corr(confs, 6, file_head.l0,data);
      // symmetrise_corr(confs, 7, file_head.l0,data);

      // if (params.data.ncorr>11){
      //     symmetrise_corr(confs, 11, file_head.l0,data);
      // }
      // if(params.data.ncorr>33){
      //     for(int i =33 ;i< 48; i++)
      //         symmetrise_corr(confs, i, file_head.l0,data);
      // }
      // if(params.data.ncorr>105){
      //     for(int i =49 ;i< 116; i++)
      //         symmetrise_corr(confs, i, file_head.l0,data);
      // }   


  //    FILE *f3t16=fopen("E1_1.txt","w+");
  //    for (int iconf=0; iconf< confs ;iconf++){
  //        for (int t =0; t< T ;t++){
  //            fprintf(f3t16,"%.12g  ",data[iconf][1][t][0]);
  //        }
  //        fprintf(f3t16,"\n");
  //    }
  //    fclose(f3t16);

    data_bin = binning(confs, var, file_head.l0, data, bin);
    //if you want to do the gamma analysis you need to do before freeing the raw data
    effective_mass_phi4_gamma(option, kinematic_2pt, (char*)"P5P5", data_bin, Neff, namefile_plateaux, out_gamma, 0, "M_{PS}^{ll}");
    effective_mass_phi4_gamma(option, kinematic_2pt, (char*)"P5P5", data_bin, Neff, namefile_plateaux, out_gamma, 1, "M_{PS1}^{ll}");
    //effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data,  confs ,namefile_plateaux,out_gamma,3,"M_{PS}^{ll}");

    free_corr(confs, var, file_head.l0, data);

    // conf_jack = create_resampling(option[4], Neff, var, file_head.l0, data_bin);
    conf_jack = myres->create(Neff, var, file_head.l0, data_bin);
    free_corr(Neff, var, file_head.l0, data_bin);

    ////////////////// symmetrization/////////////////////////////////////////////
    for (int i = 0;i <= 7;i++) { symmetrise_jackboot(Njack, i, file_head.l0, conf_jack); }
    if (params.data.ncorr > 11) { symmetrise_jackboot(Njack, 11, file_head.l0, conf_jack); }
    if (params.data.ncorr > 33) { for (int i = 33;i < 48;i++)  symmetrise_jackboot(Njack, i, file_head.l0, conf_jack); }
    if (params.data.ncorr > 105) { for (int i = 49;i < 116;i++)  symmetrise_jackboot(Njack, i, file_head.l0, conf_jack); }
    ////////////////// symmetrization/////////////////////////////////////////////







    double* zeros = (double*)calloc(Njack, sizeof(double));

    double** mass = (double**)malloc(sizeof(double*) * 3);
    fprintf(outfile, "#correlator\n");
    for (int t = 1; t < T / 2;t++) {
        for (int v = 0;v < 2;v++) {
            double* mj0 = (double*)malloc(sizeof(double) * Njack);
            for (int j = 0;j < Njack;j++) {
                mj0[j] = conf_jack[j][v][t][0] - conf_jack[j][v][t + 1][0];
            }
            double* m = mean_and_error(option[4], Njack, mj0);
            if (v == 0) { fprintf(outfile, "%d  %g  %g\t", t, m[0], m[1]);   free(m); }
            else { fprintf(outfile, "%g  %g\t", m[0], m[1]);   free(m); }
            free(mj0);
        }
        fprintf(outfile, "\n");
    }



    /////////////////// some declarations///////////
    double* E1_0_px, * E1_1_px, * E1_0_py, * E1_1_py, * E1_0_pz, * E1_1_pz;
    double* E1_0_p1, * E1_0_p11, * E1_0_p111;
    double* E2_0_p1, * E2_0_p11, * E2_0_p111, * E2_01;
    double* E2_0_A1;
    double* E3_0;
    /////!!!!!!!!!!!!!!!!! write 0 in the first jackknife!!!!!!!!!!!!!
    double* mj0 = (double*)calloc(Njack, sizeof(double));
    fwrite(mj0, sizeof(double), Njack, jack_file);
    free(mj0);
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5f\n  %.15g    %.15g    %d   %d\n\n\n", "#need_for_gnuplot", 0, 0, 0.0, 0.0, 0.0, 0, 0);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //print all the effective masses correlators
    //set the option to not read for a plateaux
    double t1 = timestamp();
    char  save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");

    FILE* dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    for (int icorr = 0; icorr < params.data.ncorr; icorr++) {
        //log effective mass
        double* tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack
            , namefile_plateaux, outfile_meff_corr, icorr, "meff_corr", M_eff_log, dev_null, fit_info_silent);
        free(tmp_meff_corr);

        //raw correlator
        int save_T = file_head.l0;
        file_head.l0 *= 2;
        tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "raw_corr", identity, dev_null, fit_info_silent);
        file_head.l0 = save_T;
        free(tmp_meff_corr);


        // shifted correlator
        tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_shifted_corr, icorr, "shifted_corr", shift_corr, dev_null, fit_info_silent);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack
            , namefile_plateaux, outfile_log_meff_shifted, icorr, "log_meff_shifted", M_eff_log_shift, dev_null, fit_info_silent);
        free(tmp_meff_corr);



    }
    fit_info_silent.restore_default();
    sprintf(option[1], "%s", save_option);// restore option
    corr_counter = 0;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    file_head.k[2] = mu1;
    file_head.k[3] = mu1;
    struct fit_type fit_info;
    mysprintf(namefile, NAMESIZE, "%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_scan_plateaux",
        argv[3], T, params.data.L[1], params.data.msq0, params.data.msq1,
        params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC, params.data.replica);
    fit_info.name_plateaux_scan = namefile;
    fit_info.f_plateaux_scan = open_file(fit_info.name_plateaux_scan.c_str(), "w+");

    if (strcmp(argv[4], "G2t_T64_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
        fit_info.plateaux_scan = true;
    }

    //c++ 1 || r 2
    mass[0] = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 0, "E1_0", M_eff_T, jack_file, fit_info);
    //mass=compute_effective_mass(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,0,"M_{PS}^{ll}");
    fit_info.restore_default();
    file_head.k[2] = mu2;
    file_head.k[3] = mu2;

    //c++ 2 || r 3
    mass[1] = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1, "E1_1", M_eff_T, jack_file);

    //!!!!!
    //there is not this correlation function  <phi0 phi1>
    //!!!!//c++ 3 || r 4
    mass[2] = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1, "E1", M_eff_T, jack_file);
    //compute_effective_mass(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,1,"M_{PS}^{ll}");


    //c++ 4 || r 5
    double** E2 = (double**)malloc(sizeof(double*) * 3);
    file_head.k[2] = mu1;    file_head.k[3] = mu1;
    E2[0] = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 2, "E2_0", shift_and_M_eff_sinh_T, jack_file);
    double* a_0 = scattering_len_luscher(Njack, mass[0], mass[0], E2[0], params.data.L[3]);
    double* tmpj = (double*)malloc(sizeof(double) * Njack);
    sub_jackboot(Njack, tmpj, E2[0], mass[0]);
    sub_jackboot(Njack, tmpj, tmpj, mass[0]);
    fprintf(outfile, "#scattering length  a  err deltaE2 err    mu00 err    deltaE2*mu00  err   a_00*(m0)=-3lambda/4pi  err\n %g  %g     %g  %g\t",
        a_0[Njack - 1], error_jackboot(option[4], Njack, a_0), tmpj[Njack - 1], error_jackboot(option[4], Njack, tmpj));


    double* tmp_muj = (double*)malloc(sizeof(double) * Njack);

    //reduced mass
    for (j = 0; j < Njack;j++)
        tmp_muj[j] = mass[0][j] * mass[0][j] / (mass[0][j] + mass[0][j]);
    fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

    //reduced mass time DeltaE2
    for (j = 0; j < Njack;j++)
        tmp_muj[j] = tmp_muj[j] * tmpj[j];
    fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

    double* a0m0 = (double*)malloc(sizeof(double) * Njack);
    for (j = 0; j < Njack;j++)
        a0m0[j] = a_0[j] * (mass[0][j]);
    fprintf(outfile, "%g   %g\n", a0m0[Njack - 1], error_jackboot(option[4], Njack, a0m0));
    free(tmpj); free(tmp_muj);


    int dvec[3] = { 0,0,0 };
    if (params.data.lambdaC0 != 0 && E2[0][Njack - 1] > 2 * mass[0][Njack - 1]) phase_shift(E2[0], mass[0], dvec, params.data.L[1], outfile, Njack, option[4]);

    /*
    double *delta=(double*) malloc(sizeof(double)*Njack);
    double *k=(double*) malloc(sizeof(double)*Njack);
    double *kcotd=(double*) malloc(sizeof(double)*Njack);
    int dvec[3]= {0,0,0};
    phase_shift(fit_out.P[0],mass[0],dvec, params.data.L[1], outfile,  Njack, option[4] );

    fprintf(outfile,"#re(delta)   err   k  err   kcotd\n");
    for (int j=0;j< Njack;j++){
        delta[j]=phase_shift( E2[0][j], mass[0][j],dvec, params.data.L[1] );
        k[j]=sqrt(E2[0][j]*E2[0][j]/4.-mass[0][j]*mass[0][j]);
        kcotd[j]=k[j]/std::tan(delta[j]);
    }
    fprintf(outfile,"%.12g  %.12g %.12g  %.12g    %.12g  %.12g\n ", delta[Njack-1],error_jackboot(option[4],Njack,delta )  , k[Njack-1],  error_jackboot(option[4],Njack,k ),       kcotd[Njack-1],  error_jackboot(option[4],Njack,kcotd ));
    free(delta);free(k);*/

    ///////////////////////////////////////////////////////////////////////////////////////

    //c++ 5 || r 6
    file_head.k[2] = mu2;    file_head.k[3] = mu2;
    E2[1] = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 3, "E2_1", shift_and_M_eff_sinh_T, jack_file);
    double* a_1 = scattering_len_luscher(Njack, mass[1], mass[1], E2[1], params.data.L[1]);
    tmpj = (double*)malloc(sizeof(double) * Njack);
    sub_jackboot(Njack, tmpj, E2[1], mass[1]);
    sub_jackboot(Njack, tmpj, tmpj, mass[1]);
    fprintf(outfile, "#scattering length  a  err deltaE2 err\n %g  %g     %g  %g\n",
        a_1[Njack - 1], error_jackboot(option[4], Njack, a_1), tmpj[Njack - 1], error_jackboot(option[4], Njack, tmpj));
    free(tmpj);
    free(a_1);

    //c++ 6 || r 7
    file_head.k[2] = mu1;    file_head.k[3] = mu2;
    E2[2] = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 4, "E2", shift_and_M_eff_sinh_T, jack_file);


    //     struct fit_type fit_info;
    struct fit_result  fit_out;
    fit_info.Nvar = 1;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.function = C3;
    fit_info.n_ext_P = 2;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);

    //c++ 7 || r 8
    fit_info.ext_P[0] = mass[0];
    fit_info.ext_P[1] = E2[0];
    file_head.k[2] = mu1;    file_head.k[3] = mu1;
    fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 5, 0/*reim*/, "E3_0", fit_info, jack_file);
    E3_0 = malloc_copy_jackboot(Njack, fit_out.P[0]);
    int dvec_A1[3] = { 0,0,0 };
    E3_print_extra(E3_0, mass[0], dvec_A1, params.data.L[1], outfile, Njack, option[4]);
    free_fit_result(fit_info, fit_out);

    //c++ 8 || r 9
    fit_info.ext_P[0] = mass[1];
    fit_info.ext_P[1] = E2[1];
    file_head.k[2] = mu2;    file_head.k[3] = mu2;
    fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 6, 0/*reim*/, "E3_1", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);

    //c++ 9 || r 10
    file_head.k[2] = mu1;    file_head.k[3] = mu2;
    fit_info.ext_P[0] = mass[1];
    fit_info.ext_P[1] = E2[2];
    fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 7, 0/*reim*/, "E3", fit_info
        , jack_file);
    free_fit_result(fit_info, fit_out);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 2;
    fit_info.function = four_pt_BH_t_T_8;
    //c++ 10 || r 11
    fit_info.ext_P[0] = mass[0];
    fit_info.ext_P[1] = mass[0];
    file_head.k[2] = mu1;    file_head.k[3] = mu1;
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0, "E4_0", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);

    //c++ 11 || r 12
    file_head.k[2] = mu2;    file_head.k[3] = mu2;
    fit_info.ext_P[0] = mass[1];
    fit_info.ext_P[1] = mass[1];
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1, "E4_1", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);


    file_head.k[2] = mu1;    file_head.k[3] = mu2;
    fit_info.ext_P[0] = mass[0];
    fit_info.ext_P[1] = mass[1];
    //c++ 12 || r 13
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH, "E4", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);
    //c++ 13 || r 14
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH, "E4_plat1", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);
    //c++ 14 || r 15
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH, "E4_plat2", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);
    //c++ 15 || r 16
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH, "E4_line", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    fit_info.Nvar = 1;
    fit_info.Npar = 2;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 2;
    fit_info.function = four_pt_BH_2par;

    file_head.k[2] = mu1;    file_head.k[3] = mu2;
    fit_info.ext_P[0] = mass[0];
    fit_info.ext_P[1] = mass[1];
    //c++ 16 || r 17
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH, "E4_2p", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    fit_info.Nvar = 1;
    fit_info.Npar = 2;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 2;
    fit_info.function = four_pt_BH_t_T_8_const;

    file_head.k[2] = mu1;    file_head.k[3] = mu2;
    fit_info.ext_P[0] = mass[0];
    fit_info.ext_P[1] = mass[1];
    //c++ 17 || r 18
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH, "E4_const", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);



    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    fit_info.Nvar = 1;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 2;
    fit_info.function = four_pt_BH_3par;

    file_head.k[2] = mu1;    file_head.k[3] = mu2;
    fit_info.ext_P[0] = mass[0];
    fit_info.ext_P[1] = mass[1];
    //c++ 18 || r 19
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH, "E4_3p", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

    fit_info.Nvar = 1;
    fit_info.Npar = 3;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 2;
    fit_info.function = C2_diff_masses;

    file_head.k[2] = mu1;    file_head.k[3] = mu2;
    fit_info.ext_P[0] = mass[0];
    fit_info.ext_P[1] = mass[1];

    fit_info.repeat_start = 10;
    fit_info.acc = 1e-5;

    if (strcmp(argv[4], "G2t_T64_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
        fit_info.plateaux_scan = true;
    }

    //c++ 19 || r 20
    if (params.data.L[0] < 100) {
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 11, 0/*reim*/, "E2_01", fit_info, jack_file);

    }
    else {
        double* tmp = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 11, "E2_01", M_eff_T, jack_file);
        fit_out.P = (double**)malloc(sizeof(double*) * Njack);
        fit_out.P[0] = tmp;
        fit_out.P[1] = (double*)malloc(sizeof(double) * Njack);
        fit_out.P[2] = (double*)malloc(sizeof(double) * Njack);
        fit_out.chi2 = (double*)malloc(sizeof(double) * Njack);
        fit_out.C = double_malloc_3(Njack, fit_info.Npar, fit_info.Npar);
    }
    {
        E2_01 = malloc_copy_jackboot(Njack, fit_out.P[0]);
        double* a = scattering_len_luscher(Njack, mass[0], mass[1], fit_out.P[0], params.data.L[1]);
        tmpj = (double*)malloc(sizeof(double) * Njack);
        tmp_muj = (double*)malloc(sizeof(double) * Njack);
        sub_jackboot(Njack, tmpj, fit_out.P[0], mass[0]);
        sub_jackboot(Njack, tmpj, tmpj, mass[1]);
        fprintf(outfile, "#scattering length  a  err deltaE2 err    mu01 err    deltaE2*mu01  err   a_01*(m0+m1)=-mu/2pi  err      a_01*(m0+m1)/(a0m0)=4/3 a_01*pi/mu   err    \n %g  %g     %g  %g\t",
            a[Njack - 1], error_jackboot(option[4], Njack, a), tmpj[Njack - 1], error_jackboot(option[4], Njack, tmpj));

        /*
        FILE *fm0=open_file("mass_0.txt","w+");
        FILE *fm1=open_file("mass_1.txt","w+");
        FILE *fm2=open_file("mass_2.txt","w+");
        FILE *fE2_01=open_file("E2_01.txt","w+");
        for (j=0; j<Njack;j++){
            fprintf(fm0,"%.12g\n",mass[0][j]);
            fprintf(fm1,"%.12g\n",mass[1][j]);
            fprintf(fm2,"%.12g\n",mass[2][j]);
            fprintf(fE2_01,"%.12g\n",fit_out.P[0][j]);
        }
        */
        //reduced mass
        for (j = 0; j < Njack;j++)
            tmp_muj[j] = mass[0][j] * mass[1][j] / (mass[0][j] + mass[1][j]);
        fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

        //reduced mass time DeltaE2
        for (j = 0; j < Njack;j++)
            tmp_muj[j] = tmp_muj[j] * tmpj[j];
        fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));


        for (j = 0; j < Njack;j++)
            tmp_muj[j] = a[j] * (mass[0][j] + mass[1][j]);
        fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

        for (j = 0; j < Njack;j++)
            tmp_muj[j] = (a[j] * (mass[0][j] + mass[1][j])) / a0m0[j];
        fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

        // a_01*pi/mu
        for (j = 0; j < Njack;j++)
            tmp_muj[j] = a[j] * pi_greco * (mass[0][j] + mass[1][j]) / (mass[0][j] * mass[1][j]);
        fprintf(outfile, "%g   %g\n", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

        free(a);
        free(tmpj); free(tmp_muj);
    }
    free_fit_result(fit_info, fit_out);
    fflush(outfile);
    fit_info.restore_default();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //GEVP two particle
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 0;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
    fit_info.function = constant_fit;

    file_head.k[2] = mu1;    file_head.k[3] = mu2;
    fit_info.ext_P[0] = mass[0];
    fit_info.ext_P[1] = mass[1];
    //c++ 20 || r 21
    if (params.data.muC > 0) {
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, GEVP_shift_matrix, "GEVP_E2_01", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);
    }
    else { zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(20);



    if (params.data.ncorr > 15) {
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// C4_BH_03t16
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_t3;
        //c++ 21 || r 22
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;

        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0_03t16, "E4_0_03t16", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);


        //c++ 22 || r 23
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1_03t16, "E4_1_03t16", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 23 || r 24
        fit_info.function = four_pt_BH_t_t3;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        if (params.data.L[0] < 50 || strcmp(argv[4], "G2t_T96_L22_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0)
            fit_info.plateaux_scan = true;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_03t16, "E4_03t16", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();
        ///// const fit

        fit_info.Nvar = 1;
        fit_info.Npar = 2;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
        fit_info.function = four_pt_BH_00_t_t3_const;

        //c++ 24 || r 25
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0_03t16, "E4_0_03t16_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 25 || r 26
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1_03t16, "E4_1_03t16_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 26 || r 27
        fit_info.function = four_pt_BH_t_t3_const;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_03t16, "E4_03t16_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);



        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// C4_BH_04t16
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_t4;
        //c++ 27 || r 22
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0_04t16, "E4_0_04t16", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 28 || r 23
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1_04t16, "E4_1_04t16", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 29 || r 24
        fit_info.function = four_pt_BH_t_t4;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_04t16, "E4_04t16", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);
        ///// const fit
        fit_info.Nvar = 1;
        fit_info.Npar = 2;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_t4_const;

        //c++ 30 || r 25
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0_04t16, "E4_0_04t16_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 31 || r 26
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1_04t16, "E4_1_04t16_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 32 || r 27
        fit_info.function = four_pt_BH_t_t4_const;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_04t16, "E4_04t16_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// C4_BH_03t20
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_t3;
        //c++ 33 || r 22
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0_03t20, "E4_0_03t20", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 34 || r 23
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1_03t20, "E4_1_03t20", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 35 || r 24
        fit_info.function = four_pt_BH_t_t3;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_03t20, "E4_03t20", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);
        ///// const fit
        fit_info.Nvar = 1;
        fit_info.Npar = 2;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_t3_const;

        //c++ 36 || r 25
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0_03t20, "E4_0_03t20_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 37 || r 26
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1_03t20, "E4_1_03t20_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 38 || r 27
        fit_info.function = four_pt_BH_t_t3_const;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_03t20, "E4_03t20_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// C4_BH_04t20
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_t4;
        //c++ 39 || r 22
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0_04t20, "E4_0_04t20", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 40 || r 23
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1_04t20, "E4_1_04t20", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 41 || r 24
        fit_info.function = four_pt_BH_t_t4;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_04t20, "E4_04t20", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);
        ///// const fit
        fit_info.Nvar = 1;
        fit_info.Npar = 2;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_t4_const;

        //c++ 42 || r 25
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0_04t20, "E4_0_04t20_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 43 || r 26
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1_04t20, "E4_1_04t20_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 44 || r 27
        fit_info.function = four_pt_BH_t_t4_const;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_04t20, "E4_04t20_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);


        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       /// C4_BH_05t20
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_t5;
        //c++ 45 || r 22
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0_05t20, "E4_0_05t20", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 46 || r 23
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1_05t20, "E4_1_05t20", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 47 || r 24
        fit_info.function = four_pt_BH_t_t5;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_05t20, "E4_05t20", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);
        ///// const fit
        fit_info.Nvar = 1;
        fit_info.Npar = 2;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_t5_const;

        //c++ 48 || r 25
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_0_05t20, "E4_0_05t20_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 49 || r 26
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_1_05t20, "E4_1_05t20_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 50 || r 27
        fit_info.function = four_pt_BH_t_t5_const;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_05t20, "E4_05t20_const", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);



        ///////////////////////togrep:shifted fit of BH
        //03t16
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_tx_shifted<3, 1>;

        //c++ 51 || r 52
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_00_tx_tf_shifetd<1, 3, 16, 0, 0, 15>, "E4_0_03t16_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 52 || r 53
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<1, 3, 16, 1, 1, 16>, "E4_1_03t16_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        fit_info.function = four_pt_BH_t_tx_shifted<3, 1>;
        //c++ 53 || r 54
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<1, 3, 16, 0, 1, 17>, "E4_03t16_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);
        //////////////////////
        //04t16
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_tx_shifted<4, 1>;

        //c++ 54 || r 55
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_00_tx_tf_shifetd<1, 4, 16, 0, 0, 18>, "E4_0_04t16_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 55 || r 56
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<1, 4, 16, 1, 1, 19>, "E4_1_04t16_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 56 || r 57
        fit_info.function = four_pt_BH_t_tx_shifted<4, 1>;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<1, 4, 16, 0, 1, 20>, "E4_04t16_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //////////////////////
        //03t20
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_tx_shifted<3, 1>;

        //c++ 57 || r 58
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_00_tx_tf_shifetd<1, 3, 20, 0, 0, 21>, "E4_0_03t20_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 58 || r 59
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<1, 3, 20, 1, 1, 22>, "E4_1_03t20_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 59 || r 60
        fit_info.function = four_pt_BH_t_tx_shifted<3, 1>;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<1, 3, 20, 0, 1, 23>, "E4_03t20_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);


        //////////////////////
        //04t20
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_tx_shifted<4, 1>;

        //c++ 60 || r 61
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_00_tx_tf_shifetd<1, 4, 20, 0, 0, 24>, "E4_0_04t20_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 61 || r 62
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<1, 4, 20, 1, 1, 25>, "E4_1_04t20_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 62 || r 63
        fit_info.function = four_pt_BH_t_tx_shifted<4, 1>;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<1, 4, 20, 0, 1, 26>, "E4_04t20_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //////////////////////
        //05t20
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_00_t_tx_shifted<5, 1>;

        //c++ 63 || r 64
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_00_tx_tf_shifetd<1, 5, 20, 0, 0, 27>, "E4_0_05t20_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 64 || r 65
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<1, 5, 20, 1, 1, 28>, "E4_1_05t20_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        //c++ 65 || r 66
        fit_info.function = four_pt_BH_t_tx_shifted<5, 1>;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<1, 5, 20, 0, 1, 29>, "E4_05t20_shifted", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);



        /////HERE we try different shifts
        //03t16
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = four_pt_BH_t_tx_shifted<3, 2>;

        //c++ 66 || r 67
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<2, 3, 16, 0, 1, 17>, "E4_03t16_shifted_2", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        fit_info.function = four_pt_BH_t_tx_shifted<3, 3>;

        //c++ 67 || r 68
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<3, 3, 16, 0, 1, 17>, "E4_03t16_shifted_3", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        fit_info.function = four_pt_BH_t_tx_shifted<3, 4>;

        //c++ 68 || r 69
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<4, 3, 16, 0, 1, 17>, "E4_03t16_shifted_4", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

        fit_info.function = four_pt_BH_t_tx_shifted<3, 5>;

        //c++ 69 || r 70
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            lhs_four_BH_01_tx_tf_shifetd<4, 3, 16, 0, 1, 17>, "E4_03t16_shifted_5", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);


        if (params.data.ncorr > 30) {
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// C4_BH_10_03t16
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Nvar = 1;
            fit_info.Npar = 1;
            fit_info.N = 1;
            fit_info.Njack = Njack;
            fit_info.n_ext_P = 2;
            fit_info.function = four_pt_BH_t_t3;

            //c++ 70 || r 71
            file_head.k[2] = mu1;    file_head.k[3] = mu2;
            fit_info.ext_P[0] = mass[1];
            fit_info.ext_P[1] = mass[0];
            fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_10_03t16, "E4_10_03t16", fit_info, jack_file);
            free_fit_result(fit_info, fit_out);
            ///// const fit

            fit_info.Nvar = 1;
            fit_info.Npar = 2;
            fit_info.N = 1;
            fit_info.Njack = Njack;
            fit_info.n_ext_P = 2;
            fit_info.function = four_pt_BH_t_t3_const;


            //c++ 71 || r 72
            file_head.k[2] = mu1;    file_head.k[3] = mu2;
            fit_info.ext_P[0] = mass[1];
            fit_info.ext_P[1] = mass[0];
            fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, lhs_four_BH_10_03t16, "E4_10_03t16_const", fit_info, jack_file);
            free_fit_result(fit_info, fit_out);

            if (params.data.ncorr > 32) {

                //02t10
                fit_info.Nvar = 1;
                fit_info.Npar = 1;
                fit_info.N = 1;
                fit_info.Njack = Njack;
                fit_info.n_ext_P = 2;
                fit_info.function = four_pt_BH_t_tx_shifted<2, 1>;

                //c++ 72 || r 73
                fit_info.ext_P[0] = mass[0];
                fit_info.ext_P[1] = mass[1];
                file_head.k[2] = mu1;    file_head.k[3] = mu1;
                fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
                    lhs_four_BH_01_tx_tf_shifetd<1, 2, 10, 0, 1, 31>, "E4_02t10_shifted", fit_info, jack_file);
                free_fit_result(fit_info, fit_out);


                //02t12
                fit_info.Nvar = 1;
                fit_info.Npar = 1;
                fit_info.N = 1;
                fit_info.Njack = Njack;
                fit_info.n_ext_P = 2;
                fit_info.function = four_pt_BH_t_tx_shifted<2, 1>;

                //c++ 73 || r 74
                fit_info.ext_P[0] = mass[0];
                fit_info.ext_P[1] = mass[1];
                file_head.k[2] = mu1;    file_head.k[3] = mu1;
                fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
                    lhs_four_BH_01_tx_tf_shifetd<1, 2, 12, 0, 1, 32>, "E4_02t12_shifted", fit_info, jack_file);
                free_fit_result(fit_info, fit_out);

                if (params.data.ncorr > 33) {
                    file_head.k[2] = mu1;
                    file_head.k[3] = mu1;
                    //c++ 74 || r 75
                    E1_0_px = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 33, "E1_0_px", M_eff_T, jack_file);
                    //free(E1_0_px);

                    file_head.k[2] = mu2;
                    file_head.k[3] = mu2;
                    //c++ 75 || r 76
                    E1_1_px = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 34, "E1_1_px", M_eff_T, jack_file);
                    //free(E1_1_px);

                    file_head.k[2] = mu1;
                    file_head.k[3] = mu1;
                    //c++ 76 || r 77
                    E1_0_py = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 35, "E1_0_py", M_eff_T, jack_file);
                    //free(E1_0_py);

                    file_head.k[2] = mu2;
                    file_head.k[3] = mu2;
                    //c++ 77 || r 78
                    E1_1_py = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 36, "E1_1_py", M_eff_T, jack_file);
                    //free(E1_1_py);



                    file_head.k[2] = mu1;
                    file_head.k[3] = mu1;
                    //c++ 78 || r 79
                    E1_0_pz = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 37, "E1_0_pz", M_eff_T, jack_file);
                    //free(E1_0_pz);

                    file_head.k[2] = mu2;
                    file_head.k[3] = mu2;
                    //c++ 79 || r 80
                    E1_1_pz = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 38, "E1_1_pz", M_eff_T, jack_file);
                    //free(E1_1_pz);


                    //c++ 80 || r 81
                    file_head.k[2] = mu1;    file_head.k[3] = mu1;
                    E2_0_A1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 39, "E2_0_A1", shift_and_M_eff_sinh_T, jack_file);
                    int dvec[3] = { 0,0,0 };
                    if (params.data.lambdaC0 != 0) phase_shift(E2_0_A1, mass[0], dvec, params.data.L[1], outfile, Njack, option[4]);

                    //     fprintf(outfile,"#E2_CM   err  E2_CM/M   err   re(delta)   err\n " );
                    //         double *E2_CM=(double*) malloc(sizeof(double)*Njack);
                    //         double *k=(double*) malloc(sizeof(double)*Njack);
                    //         for (int j=0;j< Njack;j++){
                    //             E2_CM[j]=energy_CM(E2_0_A1[j],dvec,params.data.L[1]);
                    //             k[j]=sqrt(E2_CM[j]*E2_CM[j]/4. -mass[0][j]*mass[0][j]);
                    //         }
                    //         fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
                    //         for (int j=0;j< Njack;j++)
                    //             E2_CM[j]/=mass[0][j];
                    //         fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );fflush(outfile);
                    //         double *delta=(double*) malloc(sizeof(double)*Njack);
                    //         for (int j=0;j< Njack;j++)
                    //             delta[j]=phase_shift( E2_0_A1[j], mass[0][j],dvec, params.data.L[1] );
                    //         fprintf(outfile,"%.12g  %.12g  %.12g  %.12g \n ", delta[Njack-1],error_jackboot(option[4],Njack,delta ),k[Njack-1], error_jackboot(option[4],Njack,delta ));
                    //     free(delta);free(k);


                        //c++ 81 || r 82
                    file_head.k[2] = mu2;    file_head.k[3] = mu2;
                    double* E2_1_A1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 40, "E2_1_A1", shift_and_M_eff_sinh_T, jack_file);
                    free(E2_1_A1);

                    fit_info.Nvar = 1;
                    fit_info.Npar = 3;
                    fit_info.N = 1;
                    fit_info.Njack = Njack;
                    fit_info.n_ext_P = 2;
                    fit_info.function = C2_diff_masses;

                    file_head.k[2] = mu1;    file_head.k[3] = mu2;
                    fit_info.ext_P[0] = mass[0];
                    fit_info.ext_P[1] = mass[1];
                    //c++ 82 || r 83
                    fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 41, 0/*reim*/, "E2_01_A1", fit_info, jack_file);
                    free_fit_result(fit_info, fit_out);

                    //c++ 83 || r 84
                    file_head.k[2] = mu1;    file_head.k[3] = mu1;
                    double* E2_0_E1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 42, "E2_0_E1", shift_and_M_eff_sinh_T, jack_file);
                    free(E2_0_E1);

                    //c++ 84 || r 85
                    file_head.k[2] = mu2;    file_head.k[3] = mu2;
                    double* E2_1_E1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 43, "E2_1_E1", shift_and_M_eff_sinh_T, jack_file);
                    free(E2_1_E1);

                    fit_info.Nvar = 1;
                    fit_info.Npar = 3;
                    fit_info.N = 1;
                    fit_info.Njack = Njack;
                    fit_info.n_ext_P = 2;
                    fit_info.function = C2_diff_masses;

                    file_head.k[2] = mu1;    file_head.k[3] = mu2;
                    fit_info.ext_P[0] = mass[0];
                    fit_info.ext_P[1] = mass[1];
                    //c++ 85 || r 86
                    fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 44, 0/*reim*/, "E2_01_E1", fit_info, jack_file);

                    //c++ 86 || r 87
                    file_head.k[2] = mu1;    file_head.k[3] = mu1;
                    double* E2_0_E2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 45, "E2_0_E2", shift_and_M_eff_sinh_T, jack_file);
                    free(E2_0_E2);

                    //c++ 87 || r 88
                    file_head.k[2] = mu2;    file_head.k[3] = mu2;
                    double* E2_1_E2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 46, "E2_1_E2", shift_and_M_eff_sinh_T, jack_file);
                    free(E2_1_E2);

                    fit_info.Nvar = 1;
                    fit_info.Npar = 3;
                    fit_info.N = 1;
                    fit_info.Njack = Njack;
                    fit_info.n_ext_P = 2;
                    fit_info.function = C2_diff_masses;

                    file_head.k[2] = mu1;    file_head.k[3] = mu2;
                    fit_info.ext_P[0] = mass[0];
                    fit_info.ext_P[1] = mass[1];
                    //c++ 88 || r 89
                    fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 47, 0/*reim*/, "E2_01_E2", fit_info, jack_file);


                } //if ncor>33   
                else { for (int i = 74;i < 89;i++)  zero_corr(zeros, Njack, jack_file); }

            } //if ncorr>32
            else { for (int i = 72;i < 89;i++)  zero_corr(zeros, Njack, jack_file); }

        } //if ncorr>30  
        else { for (int i = 70;i < 89;i++)  zero_corr(zeros, Njack, jack_file); }
    }//if ncorr>15  
    else { for (int i = 21;i < 89;i++)  zero_corr(zeros, Njack, jack_file); }

    //E2_01_div_shift
    fit_info.Nvar = 1;
    fit_info.Npar = 2;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 2;
    fit_info.function = C2_diff_masses_weight_shift;

    file_head.k[2] = mu1;    file_head.k[3] = mu2;
    fit_info.ext_P[0] = mass[0];
    fit_info.ext_P[1] = mass[1];
    //c++ 89 || r 90
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
        lhs_E2_div_shift<11>, "E2_01_div_shift", fit_info, jack_file);

    //a=scattering_len_luscher(  Njack,  mass[0], mass[1], fit_out.P[0] ,params.data.L[1]);
    double* a = (double*)calloc(Njack, sizeof(double));
    tmpj = (double*)malloc(sizeof(double) * Njack);
    tmp_muj = (double*)malloc(sizeof(double) * Njack);
    sub_jackboot(Njack, tmpj, fit_out.P[0], mass[0]);
    sub_jackboot(Njack, tmpj, tmpj, mass[1]);
    fprintf(outfile, "#scattering length  a  err deltaE2 err    mu01 err    deltaE2*mu01  err   a_01*(m0+m1)=-mu/2pi  err      a_01*(m0+m1)/(a0m0)=4/3 a_01*pi/mu   err    \n %g  %g     %g  %g\t",
        a[Njack - 1], error_jackboot(option[4], Njack, a), tmpj[Njack - 1], error_jackboot(option[4], Njack, tmpj));

    //reduced mass
    for (j = 0; j < Njack;j++)
        tmp_muj[j] = mass[0][j] * mass[1][j] / (mass[0][j] + mass[1][j]);
    fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

    //reduced mass time DeltaE2
    for (j = 0; j < Njack;j++)
        tmp_muj[j] = tmp_muj[j] * tmpj[j];
    fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));


    for (j = 0; j < Njack;j++)
        tmp_muj[j] = a[j] * (mass[0][j] + mass[1][j]);
    fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

    for (j = 0; j < Njack;j++)
        tmp_muj[j] = (a[j] * (mass[0][j] + mass[1][j])) / a0m0[j];
    fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

    // a_01*pi/mu
    for (j = 0; j < Njack;j++)
        tmp_muj[j] = a[j] * pi_greco * (mass[0][j] + mass[1][j]) / (mass[0][j] * mass[1][j]);
    fprintf(outfile, "%g   %g\n", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));


    free(tmpj); free(tmp_muj);

    free_fit_result(fit_info, fit_out);
    fflush(outfile);
    //////////////////////////////////////////////////////77
    if (params.data.ncorr > 48) {
        //c++ 90 || r 91
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        double* E2_0_A1E1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 48, "E2_0_A1E1", shift_and_M_eff_sinh_T, jack_file);
        free(E2_0_A1E1);

        //c++ 91 || r 92
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        double* E2_1_A1E1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 49, "E2_1_A1E1", shift_and_M_eff_sinh_T, jack_file);
        free(E2_1_A1E1);

        fit_info.Nvar = 1;
        fit_info.Npar = 3;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = C2_diff_masses;

        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        //c++ 92 || r 93
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 50, 0/*reim*/, "E2_01_A1E1", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);

    }
    else { for (int i = 90;i < 93;i++)  zero_corr(zeros, Njack, jack_file); }

    if (params.data.ncorr > 65) {
        //c++ 93 || r 94
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        double* E2_0_A1E2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 54, "E2_0_A1E2", shift_and_M_eff_sinh_T, jack_file);
        free(E2_0_A1E2);

        //c++ 94 || r 95
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        double* E2_1_A1E2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 55, "E2_1_A1E2", shift_and_M_eff_sinh_T, jack_file);
        free(E2_1_A1E2);

        fit_info.Nvar = 1;
        fit_info.Npar = 3;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = C2_diff_masses;

        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        //c++ 95 || r 96
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 56, 0/*reim*/, "E2_01_A1E2", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);



        ///////////////////////////////////////////////   A1o20   /////////////////////////////////////////////
        //c++ 96 || r 97
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        double* E2_0_A1o20 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 60, "E2_0_A1o20", shift_and_M_eff_sinh_T, jack_file);
        free(E2_0_A1o20);

        //c++ 97 || r 98
        file_head.k[2] = mu2;    file_head.k[3] = mu2;
        double* E2_1_A1o20 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 61, "E2_1_A1o20", shift_and_M_eff_sinh_T, jack_file);
        free(E2_1_A1o20);

        fit_info.Nvar = 1;
        fit_info.Npar = 3;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = C2_diff_masses;

        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = mass[1];
        //c++ 98 || r 99
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 62, 0/*reim*/, "E2_01_A1o20", fit_info, jack_file);
        free_fit_result(fit_info, fit_out);



    }
    else { for (int i = 93;i < 99;i++)  zero_corr(zeros, Njack, jack_file); }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (params.data.ncorr > 74) {


        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 0;
        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        if (params.data.L[0] < 50)
            fit_info.plateaux_scan = true;
        //c++ 99 || r 100
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            m_eff_of_sum<33, 35, 37>, "E1_0_p1", fit_info, jack_file);

        E1_0_p1 = malloc_copy_jackboot(Njack, fit_out.P[0]);
        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();

        fit_info.Nvar = 1;
        fit_info.Npar = 3;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
        fit_info.function = C2_diff_masses;

        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E1_0_p1;
        if (params.data.L[0] < 50)
            fit_info.plateaux_scan = true;

        if (strcmp(argv[4], "G2t_T64_L26_msq0-4.755000_msq1-4.800000_l02.500000_l12.500000_mu5.000000_g0.500000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.repeat_start = 5;
            fit_info.acc = 1e-5;
            fit_info.precision_sum = 4;
            fit_info.devorder = -2;
            fit_info.h = 1e-2;
            fit_info.lambda = 1e-2;
            fit_info.chi2_gap_jackboot = 0.1;
            fit_info.guess_per_jack = 5;
            fit_info.guess = { 0.434124998962052,    	0.00119524254838347,   	0.00853111018750956 };
            //         fit_info.verbosity=2;
            printf("%g  %g \n", fit_info.ext_P[0][Njack - 1], fit_info.ext_P[1][Njack - 1]);

        }
        if (strcmp(argv[4], "G2t_T64_L24_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g0.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L22_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g2.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L18_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g2.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L20_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L18_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g0.500000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L16_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L17_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L18_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L19_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L20_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L21_msq0-1.230000_msq1-0.370000_l00.000000_l10.000000_mu0.000000_g10.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L24_msq0-1.230000_msq1-0.370000_l00.000000_l10.000000_mu0.000000_g10.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L25_msq0-1.230000_msq1-0.370000_l00.000000_l10.000000_mu0.000000_g10.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L24_msq0-1.241300_msq1-0.380000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0_bin1000_merged_bin100") == 0 ||
            strcmp(argv[4], "G2t_T64_L23_msq0-1.241300_msq1-0.380000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0_bin1000_merged_bin100") == 0 
            ) {
            fit_info.repeat_start = 6;
            fit_info.acc = 1e-5;
            fit_info.precision_sum = 4;
            fit_info.devorder = -2;
            fit_info.h = 1e-2;
            fit_info.lambda = 1e-2;
            fit_info.chi2_gap_jackboot = 0.1;
            fit_info.guess_per_jack = 5;
            fit_info.guess = { 0.6,    	0.0119524254838347,   	0.00853111018750956 };
            // printf("setting custom parameters\n");
        }


        //c++ 100 || r 101
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            sum_corr_directions_shift<66, 67, 68>, "E2_0_p1", fit_info, jack_file);

        int dvec[3] = { 1,0,0 };
        if (params.data.lambdaC0 != 0 &&
            strcmp(argv[4], "G2t_T8_L4_msq0-4.850000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0
            && strcmp(argv[4], "G2t_T16_L4_msq0-4.850000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0
            && strcmp(argv[4], "G2t_T16_L8_msq0-4.850000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0
            && strcmp(argv[4], "G2t_T16_L8_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0
            && strcmp(argv[4], "G2t_T64_L26_msq0-4.755000_msq1-4.800000_l02.500000_l12.500000_mu5.000000_g0.500000_rep0_bin100_merged_bin1000") != 0
            ) {
            printf("here\n");
            phase_shift(fit_out.P[0], mass[0], dvec, params.data.L[1], outfile, Njack, option[4]);
        }
        fit_info.restore_default();

        /*
        fprintf(outfile,"#E2_CM   err  E2_CM/M   err   re(delta)   err\n " );
        double *E2_CM=(double*) malloc(sizeof(double)*Njack);

        double *k=(double*) malloc(sizeof(double)*Njack);
        for (int j=0;j< Njack;j++){
            E2_CM[j]=energy_CM(fit_out.P[0][j],dvec,params.data.L[1]);
            k[j]=sqrt(E2_CM[j]*E2_CM[j]/4. -mass[0][j]*mass[0][j]);
        }
        fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
        for (int j=0;j< Njack;j++)
            E2_CM[j]/=mass[0][j];
        fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
        free(E2_CM);
        double *delta=(double*) malloc(sizeof(double)*Njack);
        for (int j=0;j< Njack;j++)
            delta[j]=phase_shift( fit_out.P[0][j], mass[0][j],dvec, params.data.L[1] );
        fprintf(outfile,"%.12g  %.12g  %.12g  %.12g \n ", delta[Njack-1],error_jackboot(option[4],Njack,delta ), k[Njack-1], error_jackboot(option[4],Njack,k ) );
        free(delta);
        */
        E2_0_p1 = malloc_copy_jackboot(Njack, fit_out.P[0]);
        free_fit_result(fit_info, fit_out);

    }
    else { for (int i = 99;i < 101;i++)  zero_corr(zeros, Njack, jack_file); }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////p=(1,1,0)
    if (params.data.ncorr > 90) {

        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 0;
        fit_info.function = constant_fit;

        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        //c++ 101 || r 102
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            m_eff_of_sum<75, 77, 79>, "E1_0_p11", fit_info, jack_file);
        E1_0_p11 = malloc_copy_jackboot(Njack, fit_out.P[0]);
        free_fit_result(fit_info, fit_out);


        fit_info.Nvar = 1;
        fit_info.Npar = 3;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
        fit_info.function = C2_diff_masses;

        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E1_0_p11;
        fit_info.repeat_start = 10;
        //c++ 102 || r 103
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            sum_corr_directions_shift<81, 82, 83>, "E2_0_p11", fit_info, jack_file);

        fit_info.restore_default();
        int dvec[3] = { 1,1,0 };
        if (params.data.lambdaC0 != 0 &&
            strcmp(argv[4], "G2t_T8_L4_msq0-4.850000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0
            && strcmp(argv[4], "G2t_T16_L4_msq0-4.850000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0
            && strcmp(argv[4], "G2t_T16_L8_msq0-4.850000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0)
            phase_shift(fit_out.P[0], mass[0], dvec, params.data.L[1], outfile, Njack, option[4]);
        /*
        fprintf(outfile,"#E2_CM   err  E2_CM/M   err re(delta)   err\n " );
        double *E2_CM=(double*) malloc(sizeof(double)*Njack);

        double *k=(double*) malloc(sizeof(double)*Njack);
        for (int j=0;j< Njack;j++){
            E2_CM[j]=energy_CM(fit_out.P[0][j],dvec,params.data.L[1]);
            k[j]=sqrt(E2_CM[j]*E2_CM[j]/4. -mass[0][j]*mass[0][j]);
        }
        fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
        for (int j=0;j< Njack;j++)
            E2_CM[j]/=mass[0][j];
        fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
        double *delta=(double*) malloc(sizeof(double)*Njack);
        for (int j=0;j< Njack;j++)
            delta[j]=phase_shift( fit_out.P[0][j], mass[0][j],dvec, params.data.L[1] );
        fprintf(outfile,"%.12g  %.12g  %.12g  %.12g \n ", delta[Njack-1],error_jackboot(option[4],Njack,delta ), k[Njack-1], error_jackboot(option[4],Njack,k ) );
        free(delta);
        free(E2_CM);
        */
        E2_0_p11 = malloc_copy_jackboot(Njack, fit_out.P[0]);
        free_fit_result(fit_info, fit_out);


    }
    else { for (int i = 101;i < 103;i++)  zero_corr(zeros, Njack, jack_file); }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////p=(1,1,1)
    if (params.data.ncorr > 90) {

        file_head.k[2] = mu1;
        file_head.k[3] = mu1;
        //c++ 103 || r 74
        E1_0_p111 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 90, "E1_0_p111", M_eff_T, jack_file);


        fit_info.Nvar = 1;
        fit_info.Npar = 3;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
        fit_info.function = C2_diff_masses;

        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E1_0_p111;
        if (params.data.L[0] < 50)
            fit_info.plateaux_scan = true;


        //c++ 104 || r 105
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 92, 0/*reim*/, "E2_0_p111", fit_info, jack_file);
        int dvec[3] = { 1,1,1 };
        if (params.data.lambdaC0 != 0 &&
            strcmp(argv[4], "G2t_T8_L4_msq0-4.850000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0
            && strcmp(argv[4], "G2t_T16_L4_msq0-4.850000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0
            && strcmp(argv[4], "G2t_T16_L8_msq0-4.850000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0
            && strcmp(argv[4], "G2t_T16_L8_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0") != 0
            )
            phase_shift(fit_out.P[0], mass[0], dvec, params.data.L[1], outfile, Njack, option[4]);
        /*
        fprintf(outfile,"#E2_CM   err  E2_CM/M   err   re(delta)   err\n " );
        double *E2_CM=(double*) malloc(sizeof(double)*Njack);
        int dvec[3]= {1,1,1};
        double *k=(double*) malloc(sizeof(double)*Njack);
        for (int j=0;j< Njack;j++){
            E2_CM[j]=energy_CM(fit_out.P[0][j],dvec,params.data.L[1]);
            k[j]=sqrt(E2_CM[j]*E2_CM[j]/4. -mass[0][j]*mass[0][j]);
        }
        fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
        for (int j=0;j< Njack;j++)
            E2_CM[j]/=mass[0][j];
        fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
        double *delta=(double*) malloc(sizeof(double)*Njack);
        for (int j=0;j< Njack;j++)
            delta[j]=phase_shift( fit_out.P[0][j], mass[0][j],dvec, params.data.L[1] );
        fprintf(outfile,"%.12g  %.12g  %.12g  %.12g \n ", delta[Njack-1],error_jackboot(option[4],Njack,delta ), k[Njack-1], error_jackboot(option[4],Njack,k ) );
        free(delta);
        free(E2_CM);free(k);*/
        E2_0_p111 = malloc_copy_jackboot(Njack, fit_out.P[0]);

        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();


    }
    else { for (int i = 103;i < 105;i++)  zero_corr(zeros, Njack, jack_file); }

    error(corr_counter != 104, -1, "correlator counter wrong", "corr_counter=%d     expected=104", corr_counter);

    if (params.data.ncorr > 122) {
        fit_info.Nvar = 1;
        fit_info.Npar = 4;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = C3p;
        fit_info.n_ext_P = 4;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        //c++ 105 || r 106
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E2[0];
        fit_info.ext_P[2] = E1_0_p1;
        fit_info.ext_P[3] = E2_0_p1;
        for (int i = 0;i < fit_info.n_ext_P; i++)
            printf("%g\t", fit_info.ext_P[i][Njack - 1]);
        printf("\n\n");
        fit_info.guess = { 0.517571,7.49491e-05,0.00119068,0.000312766 };
        fit_info.h = 0.0001;
        fit_info.acc = 1e-6;

        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, sum_corr_directions_shift<95, 96, 97>, "E3_0_p1", fit_info, jack_file);
        /*fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  95,0 , "E3_0_p1",  fit_info, jack_file);
       */
        int dvec_p1[3] = { 1,0,0 };
        E3_print_extra(fit_out.P[0], mass[0], dvec_p1, params.data.L[1], outfile, Njack, option[4]);
        free_fit_result(fit_info, fit_out);
        //fit_info.restore_default();

        //c++ 106 || r 107
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E2[0];
        fit_info.ext_P[2] = E1_0_p11;
        fit_info.ext_P[3] = E2_0_p11;

        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, sum_corr_directions_shift<104, 105, 106>, "E3_0_p11", fit_info, jack_file);
        int dvec_p11[3] = { 1,1,0 };
        E3_print_extra(fit_out.P[0], mass[0], dvec_p11, params.data.L[1], outfile, Njack, option[4]);
        free_fit_result(fit_info, fit_out);
        //     

        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E2[0];
        fit_info.ext_P[2] = E1_0_p111;
        fit_info.ext_P[3] = E2_0_p111;
        fit_info.guess = { 0.650579 ,6.00313e-05,0.000947987 ,8.98107e-08 };
        //c++ 107 || r 108
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 113, 0/*reim*/, "E3_0_p111", fit_info, jack_file);
        int dvec_p111[3] = { 1,1,1 };
        E3_print_extra(fit_out.P[0], mass[0], dvec_p111, params.data.L[1], outfile, Njack, option[4]);
        free_fit_result(fit_info, fit_out);
        //  


        fit_info.Nvar = 1;
        fit_info.Npar = 4;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = C3_A1_old;
        fit_info.n_ext_P = 5;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);

        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E2_0_p1;
        fit_info.ext_P[2] = E1_0_p1;
        fit_info.ext_P[3] = E2_0_A1;
        fit_info.ext_P[4] = E2[0];

        //c++ 108 || r 109
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_info.guess = { 0.63738,   387.1e-7,    1.0e-7 ,   0.0060975 };//0.12720
        //     fit_info.guess={0.773652, 3.96386e-05 ,1.33604e-05 ,0.0241938};

        //      fit_info.repeat_start=40;
        //     fit_info.plateaux_scan=true;

        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 116, 0/*reim*/, "E3_0_A1", fit_info, jack_file);
        int dvec_A1[3] = { 0,0,0 };
        E3_print_extra(fit_out.P[0], mass[0], dvec_A1, params.data.L[1], outfile, Njack, option[4]);
        free_fit_result(fit_info, fit_out);
        //  
        fit_info.restore_default();

        /////////////////////////////////////////////////////////////////////////////////////
        /////////////matrix element
        /////////////////////////////////////////////////////////////////////////////////////
        int T = params.data.L[0];
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = me_k3pi_rhs_T_2;
        fit_info.n_ext_P = 2;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = E3_0;

        // c++ 109|| r110
        file_head.k[2] = mu1;    file_head.k[3] = mu1;

        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, matrix_element_k3pi_T_2, "me_T_2", fit_info, jack_file);

        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = me_k3pi_rhs<10>;
        fit_info.n_ext_P = 2;

        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = E3_0;

        // c++ 110|| r111
        file_head.k[2] = mu1;    file_head.k[3] = mu1;

        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, matrix_element_k3pi<120, 10>, "me_t10", fit_info, jack_file);

        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = me_k3pi_rhs<12>;
        fit_info.n_ext_P = 2;

        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = E3_0;

        // c++ 111|| r112
        file_head.k[2] = mu1;    file_head.k[3] = mu1;

        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, matrix_element_k3pi<121, 12>, "me_t12", fit_info, jack_file);

        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = me_k3pi_rhs<16>;
        fit_info.n_ext_P = 2;

        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = E3_0;

        // c++ 112|| r113
        file_head.k[2] = mu1;    file_head.k[3] = mu1;

        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, matrix_element_k3pi<122, 16>, "me_t16", fit_info, jack_file);


        fit_info.restore_default();
    }
    else { for (int i = 105;i < 113;i++)  zero_corr(zeros, Njack, jack_file); }

    error(corr_counter != 112, -1, "correlator counter wrong", "corr_counter=%d     expected=112", corr_counter);

    fit_info.restore_default();

    if (params.data.ncorr > 90) {


        fit_info.Nvar = 1;
        fit_info.Npar = 2;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
        fit_info.function = C2_diff_masses_weight_shift;

        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E1_0_p1;
        if (params.data.L[0] < 50)
            fit_info.plateaux_scan = true;
        fit_info.repeat_start = 10;

        //c++ 113 || r 114
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            sum_corr_weight_shift<66, 67, 68>, "E2_0_p1_ws", fit_info, jack_file);
        int dvec[3] = { 1,0,0 };
        //     if (params.data.lambdaC0!=0) phase_shift(fit_out.P[0],mass[0],dvec, params.data.L[1], outfile,  Njack, option[4] );

        E2_0_p1 = malloc_copy_jackboot(Njack, fit_out.P[0]);

        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();

        fit_info.Nvar = 1;
        fit_info.Npar = 2;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
        fit_info.function = C2_diff_masses_weight_shift;

        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E1_0_p11;
        //c++ 114 || r 115
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            sum_corr_weight_shift<81, 82, 83>, "E2_0_p11_ws", fit_info, jack_file);
        int dvec_p11[3] = { 1,1,0 };
        //     if (params.data.lambdaC0!=0) phase_shift(fit_out.P[0],mass[0],dvec_p11, params.data.L[1], outfile,  Njack, option[4] );
        free_fit_result(fit_info, fit_out);

        fit_info.Nvar = 1;
        fit_info.Npar = 2;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        fit_info.function = C2_diff_masses_weight_shift;

        file_head.k[2] = mu1;    file_head.k[3] = mu2;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E1_0_p111;

        if (params.data.L[0] < 50)
            fit_info.plateaux_scan = true;


        //c++ 115 || r 116
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile,
            sum_corr_weight_shift<92, 92, 92>, "E2_0_p111_ws", fit_info, jack_file);
        int dvec_p111[3] = { 1,1,1 };
        //     if (params.data.lambdaC0!=0) phase_shift(fit_out.P[0],mass[0],dvec_p111, params.data.L[1], outfile,  Njack, option[4] );
        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();

    }
    else { for (int i = 113;i < 116;i++)  zero_corr(zeros, Njack, jack_file); }

    error(corr_counter != 115, -1, "correlator counter wrong", "corr_counter=%d     expected=115", corr_counter);

    double* me_3phi;
    if (params.data.ncorr > 122) {


        fit_info.Nvar = 1;
        fit_info.Npar = 4;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = C3_vev;
        fit_info.n_ext_P = 2;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);

        //c++ 116|| r 117
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E2[0];
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_info.guess = { 0.3961,    0.00018161,    0.000337   , 0.00022318 };

        fit_info.h = 0.0001;
        fit_info.acc = 1e-6;
        //fit_info.repeat_start=10;
        //fit_info.precision_sum=2;

        if (strcmp(argv[4], "G2t_T32_L24_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.41,    0.0018161,    0.00337 , 0.00022318 };
            fit_info.acc = 1e-9;
        }
        if (strcmp(argv[4], "G2t_T32_L28_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.025000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.38,    0.0018161,    0.00337 , 0.00022318 };
            fit_info.acc = 1e-9;
            fit_info.h = 1e-5;
            fit_info.chi2_gap_jackboot = 1;
            fit_info.repeat_start = 15;
            fit_info.guess_per_jack = 5;
        }
        if (strcmp(argv[4], "G2t_T32_L28_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.38,    0.0018161,    0.00337 , 0.00022318 };
            fit_info.acc = 1e-9;
            fit_info.h = 1e-5;
            fit_info.chi2_gap_jackboot = 1;
            fit_info.repeat_start = 15;
            fit_info.guess_per_jack = 5;
        }
        if (strcmp(argv[4], "G2t_T64_L36_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.025000_rep0_bin100_merged_bin1000") == 0 || strcmp(argv[4], "G2t_T64_L40_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.025000_rep0_bin100_merged_bin1000") == 0) {

            fit_info.acc = 1e-6;
            fit_info.h = 1e-5;
            fit_info.devorder = -2;
            fit_info.chi2_gap_jackboot = 1;
            fit_info.repeat_start = 15;
            fit_info.guess_per_jack = 5;
            fit_info.precision_sum = 2;
        }

        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 5, 0/*reim*/, "E3_0_vev", fit_info, jack_file);
        free(E3_0);
        //                     for (int j=0; j< Njack;j++){
        //                         printf("%d  chi= %g  P=\t",j,   fit_out.chi2[j] );
        //                         for (int i =0; i< fit_info.Npar;i++){
        //                             printf("%g\t",fit_out.P[i][j]);
        //                         }
        //                        printf("\n");
        //                    }
        //     exit(1);

        E3_0 = malloc_copy_jackboot(Njack, fit_out.P[0]);
        me_3phi = malloc_copy_jackboot(Njack, fit_out.P[1]);
        for (int j = 0;j < Njack;j++) { me_3phi[j] = sqrt(fabs(me_3phi[j]) * 2 * E3_0[j]); }

        int dvec_0[3] = { 0,0,0 };
        E3_print_extra(E3_0, mass[0], dvec_0, params.data.L[1], outfile, Njack, option[4]);
        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();

        fit_info.Nvar = 1;
        fit_info.Npar = 5;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = C3p_vev;
        fit_info.n_ext_P = 5;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        //c++ 117 || r 118
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E2[0];
        fit_info.ext_P[2] = E1_0_p1;
        fit_info.ext_P[3] = E2_0_p1;
        fit_info.ext_P[4] = E1_0_p1;

        fit_info.h = 1e-4;
        fit_info.devorder = 4;
        fit_info.acc = 1e-6;
        fit_info.lambda = 1e-4;
        fit_info.repeat_start = 15;
        fit_info.chi2_gap_jackboot = 1;
        fit_info.guess_per_jack = 3;
        //fit_info.precision_sum=2;

        fit_info.guess = { 0.53,    0.0001,   0.00059,    0.0001597,   1.e-6 };
        if (strcmp(argv[4], "G2t_T32_L36_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.482,    	-3.63061124075202e-05    ,    	-1.4396080515536e-05   ,    	3.46858622879724e-05    ,    	-2.66582380286655e-05 };
        }
        //     if (strcmp(argv[4],"G2t_T32_L30_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000")==0  ){
        //         fit_info.guess={0.519204371413175 ,    	7.59640741198693e-05    ,   	0.000449032496922057  ,    	9.89514102601514e-05   };
        //         fit_info.guess={0.53,    0.0001,   0.00059,    0.0001597,   1.e-6};
        //         printf("%g   %g   %g   %g ",mass[0][Njack-1], E2[0][Njack-1],E1_0_p1[Njack-1], E2_0_p1[Njack-1]);
        //         fit_info.h=1e-4;
        //         fit_info.devorder=4;
        //         fit_info.acc=1e-4;
        //         fit_info.lambda=1e-4;
        //         fit_info.repeat_start=10;
        //         fit_info.chi2_gap_jackboot=1;
        //         fit_info.guess_per_jack=3;
        //     }

        if (strcmp(argv[4], "G2t_T64_L32_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.025000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.499757,        5.81671e-05,     0.00124462,      1.7906e-05 ,     3.55311e-05 };
            fit_info.precision_sum = 2;
            fit_info.acc = 1e-8;
            fit_info.guess_per_jack = 10;
            fit_info.lambda = 1e-6;
        }

        if (strcmp(argv[4], "G2t_T64_L36_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.025000_rep0_bin100_merged_bin1000") == 0 || strcmp(argv[4], "G2t_T64_L40_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.025000_rep0_bin100_merged_bin1000") == 0
            ) {

            fit_info.acc = 1e-6;
            fit_info.h = 1e-5;
            fit_info.devorder = -2;
            fit_info.chi2_gap_jackboot = 1;
            fit_info.repeat_start = 15;
            fit_info.guess_per_jack = 5;
            fit_info.precision_sum = 2;
        }
        if (strcmp(argv[4], "G2t_T32_L30_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0
            ) {
            fit_info.acc = 1e-4;
            fit_info.h = 1e-5;
            fit_info.devorder = -2;
            fit_info.chi2_gap_jackboot = 0.5;
            fit_info.repeat_start = 15;
            fit_info.guess_per_jack = 5;
            fit_info.precision_sum = 2;
        }

        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, sum_corr_directions_shift<95, 96, 97>, "E3_0_p1_vev", fit_info, jack_file);
        error(corr_counter != 117, -1, "correlator counter wrong", "corr_counter=%d     expected=117", corr_counter);

        int dvec_p1[3] = { 1,0,0 };
        E3_print_extra(fit_out.P[0], mass[0], dvec_p1, params.data.L[1], outfile, Njack, option[4]);
        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();

        //c++ 118 || r 119
        fit_info.Nvar = 1;
        fit_info.Npar = 5;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = C3p_vev;
        fit_info.n_ext_P = 5;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);

        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E2[0];
        fit_info.ext_P[2] = E1_0_p11;
        fit_info.ext_P[3] = E2_0_p11;
        fit_info.ext_P[4] = E1_0_p11;
        fit_info.guess = { 0.5876,    661.1e-7,   0.00016,    0.00037,    46.3e-8 };
        if (strcmp(argv[4], "G2t_T32_L30_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.590373,   -6.66566e-05,  -0.000284073,  4.84344e-07,  -3.39931e-05 };
        }
        fit_info.h = 1e-5;
        fit_info.devorder = -2;
        fit_info.acc = 1e-11;
        fit_info.lambda = 1e-3;
        fit_info.repeat_start = 10;
        fit_info.chi2_gap_jackboot = 1;
        /*  fit_info.guess_per_jack=3;
          fit_info.precision_sum=2;
        */

        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, sum_corr_directions_shift<104, 105, 106>, "E3_0_p11_vev", fit_info, jack_file);
        int dvec_p11[3] = { 1,1,0 };
        E3_print_extra(fit_out.P[0], mass[0], dvec_p11, params.data.L[1], outfile, Njack, option[4]);


        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();

        //     
        fit_info.Nvar = 1;
        fit_info.Npar = 5;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = C3p_vev;
        fit_info.n_ext_P = 5;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);

        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E2[0];
        fit_info.ext_P[2] = E1_0_p11;
        fit_info.ext_P[3] = E2_0_p11;
        fit_info.ext_P[4] = E1_0_p11;
        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E2[0];
        fit_info.ext_P[2] = E1_0_p111;
        fit_info.ext_P[3] = E2_0_p111;
        fit_info.ext_P[4] = E1_0_p111;
        fit_info.guess = { 0.6449,    60.6e-6,    1.0e-4,    0.00077,    42.6e-6 };
        if (strcmp(argv[4], "G2t_T32_L32_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.626641,        4.71156e-05  ,   -0.000244154 ,   -6.25702e-05 ,   -1.59513e-05 };

        }

        fit_info.h = 1e-5;
        fit_info.devorder = -2;
        fit_info.acc = 1e-10;
        fit_info.lambda = 1e-4;
        fit_info.repeat_start = 10;
        fit_info.precision_sum = 2;

        if (strcmp(argv[4], "G2t_T64_L20_msq0-4.868000_msq1-4.710000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0_bin100") == 0) {
            fit_info.guess = { 0.626641,        4.71156e-05  ,   -0.000244154 ,   -6.25702e-05 ,   -1.59513e-05 };
            fit_info.h = 1e-5;
            fit_info.devorder = -2;
            fit_info.acc = 1e-4;
            fit_info.lambda = 1e-4;
            fit_info.repeat_start = 10;
            fit_info.precision_sum = 2;

        }

        //c++ 119 || r 120
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 113, 0/*reim*/, "E3_0_p111_vev", fit_info, jack_file);
        int dvec_p111[3] = { 1,1,1 };
        E3_print_extra(fit_out.P[0], mass[0], dvec_p111, params.data.L[1], outfile, Njack, option[4]);
        //     
        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();
        //  


        fit_info.Nvar = 1;
        fit_info.Npar = 5;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = C3_A1;
        fit_info.n_ext_P = 5;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);

        fit_info.ext_P[0] = mass[0];
        fit_info.ext_P[1] = E2_0_p1;
        fit_info.ext_P[2] = E1_0_p1;
        fit_info.ext_P[3] = E2_0_A1;
        fit_info.ext_P[4] = E2[0];

        //c++ 120 || r 121
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_info.guess = { 0.66,   5.e-5,    7.5e-7 ,   1e-4,    8.6e-5 };//0.12720

        if (strcmp(argv[4], "G2t_T32_L28_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.66352764,    	5.01012e-05,    	6.49650e-05,    	0.000541073 ,    	8.649137e-05 };
        }
        else if (strcmp(argv[4], "G2t_T32_L32_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.606593,        -3.0613e-05,     2.83579e-05 ,    3.23623e-09,    5.28766e-05 };
        }
        else if (strcmp(argv[4], "G2t_T32_L28_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.025000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.646593,        -3.0613e-05,     2.83579e-05 ,    3.23623e-09,    5.28766e-05 };
            fit_info.Npar = 3;
            fit_info.function = C3_A1_3par;
            fit_info.guess = { 0.653481415329485,        -3.0613e-05,     2.83579e-05 };

        }
        else if (strcmp(argv[4], "G2t_T32_L24_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.733441,      8.83205e-05,     0.000101335,    1.75938e-08,     0.000152796 };

        }

        fit_info.h = 1e-5;
        fit_info.devorder = 4;
        fit_info.acc = 1e-6;
        fit_info.lambda = 1e-4;
        fit_info.repeat_start = 15;
        fit_info.guess_per_jack = 4;
        fit_info.chi2_gap_jackboot = 1;
        fit_info.precision_sum = 2;

        if (strcmp(argv[4], "G2t_T32_L38_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.533441,      8.83205e-05,     0.000101335,    1.75938e-08,     0.000152796 };
            fit_info.acc = 1e-7;
            fit_info.devorder = -2;
            fit_info.precision_sum = 2;

        }
        if (strcmp(argv[4], "G2t_T64_L24_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.025000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.533441,      8.83205e-05,     0.000101335,    1.75938e-08,     0.000152796 };
            fit_info.acc = 1e-7;
            fit_info.devorder = -2;
            fit_info.precision_sum = 2;

        }
        if (strcmp(argv[4], "G2t_T64_L34_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.025000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.533441,      8.83205e-05,     0.000101335,    1.75938e-08,     0.000152796 };
            fit_info.acc = 1e-7;
            fit_info.devorder = -2;
            fit_info.precision_sum = 2;

        }
        if (strcmp(argv[4], "G2t_T64_L38_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.025000_rep0_bin100_merged_bin1000") == 0) {
            fit_info.guess = { 0.533441,      8.83205e-05,     0.000101335,    1.75938e-08,     0.000152796 };
            fit_info.acc = 1e-7;
            fit_info.devorder = -2;
            fit_info.precision_sum = 2;

        }
        if (strcmp(argv[4], "G2t_T64_L20_msq0-4.868000_msq1-4.710000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0_bin100") == 0) {
            fit_info.guess = { 0.626641,        4.71156e-05  ,   -0.000244154 ,   -6.25702e-05 ,   -1.59513e-05 };
            fit_info.h = 1e-5;
            fit_info.devorder = -2;
            fit_info.acc = 1e-4;
            fit_info.lambda = 1e-4;
            fit_info.repeat_start = 10;
            fit_info.precision_sum = 2;

        }
        if (strcmp(argv[4], "G2t_T64_L24_msq0-4.868000_msq1-4.710000_l01.000000_l11.000000_mu2.000000_g0.250000_rep0_bin100_merged_bin1000") == 0) {
            printf("here\n");
            fit_info.guess = { 0.626641,        4.71156e-05  ,   -0.000244154 ,   -6.25702e-05 ,   -1.59513e-05 };
            fit_info.h = 1e-5;
            fit_info.devorder = -2;
            fit_info.acc = 1e-2;
            fit_info.lambda = 1e-4;
            fit_info.repeat_start = 10;
            fit_info.precision_sum = 2;

        }
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 116, 0/*reim*/, "E3_0_A1_vev", fit_info, jack_file);

        //         for (int j=0; j< Njack;j++){
        //              printf("%d  chi= %g  P=\t",j,   fit_out.chi2[j] );
        //              for (int i =0; i< fit_info.Npar;i++){
        //                  printf("%g\t",fit_out.P[i][j]);
        //              }
        //             printf("\n");
        //         }
        //         exit(1);

        int dvec_A1[3] = { 0,0,0 };
        E3_print_extra(fit_out.P[0], mass[0], dvec_A1, params.data.L[1], outfile, Njack, option[4]);

        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();


    }
    else { for (int i = 116;i < 121;i++)  zero_corr(zeros, Njack, jack_file); }


    double* me_phi1;
    error(corr_counter != 120, -1, "correlator counter wrong", "corr_counter=%d     expected=120", corr_counter);

    if (params.data.ncorr > 126) {

        /////////////////////////////////////////////////////////////////////////////////////
        /////////////matrix element phi1
        /////////////////////////////////////////////////////////////////////////////////////

        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = constant_fit;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);

        fit_info.ext_P[0] = mass[1];

        // c++ 121|| r122

        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, me_oPp<1>, "me_phi1", fit_info, jack_file);

        me_phi1 = malloc_copy_jackboot(Njack, fit_out.P[0]);
        free_fit_result(fit_info, fit_out);


        /////////////////////////////////////////////////////////////////////////////////////
        /////////////matrix element
        /////////////////////////////////////////////////////////////////////////////////////
        int T = params.data.L[0];
        fit_info.Nvar = 1;
        fit_info.Npar = 3;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = me_3pik_rhs_T_2;
        fit_info.n_ext_P = 6;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * 6);

        fit_info.ext_P[0] = mass[1];
        fit_info.ext_P[1] = E3_0;
        fit_info.ext_P[2] = mass[0];
        fit_info.ext_P[3] = E2_01;

        fit_info.ext_P[4] = me_phi1;
        fit_info.ext_P[5] = me_3phi;
        // c++ 122|| r123
        file_head.k[2] = mu1;    file_head.k[3] = mu1;
        fit_info.guess = { 2.2e-15,   10.4e-12,    0.1e-10 };
        fit_info.repeat_start = 10;
        fit_info.acc = 1e-6;
        fit_info.h = 1e-3;
        fit_info.devorder = 4;
        fit_info.lambda = 1e-5;
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 123, 0, "me_3pik_T_2", fit_info, jack_file);

        fit_info.function = me_3pik_rhs<10>;

        // c++ 123|| r124
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 124, 0, "me_3pik_t10", fit_info, jack_file);
        check_correlatro_counter(123);

        fit_info.function = me_3pik_rhs<12>;
        // c++ 124|| r125
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 125, 0, "me_3pik_t12", fit_info, jack_file);

        fit_info.function = me_3pik_rhs<16>;
        // c++ 125|| r126
        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 126, 0, "me_3pik_t16", fit_info, jack_file);
        error(corr_counter != 125, -1, "correlator counter wrong", "corr_counter=%d     expected=125", corr_counter);

        fit_info.restore_default();
    }
    else { for (int i = 121;i < 126;i++)  zero_corr(zeros, Njack, jack_file); }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////   E2_01 ws
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Nvar = 1;
    fit_info.Npar = 2;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 2;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);

    fit_info.function = C2_diff_masses_weight_shift;

    file_head.k[2] = mu1;    file_head.k[3] = mu2;
    fit_info.ext_P[0] = mass[0];
    fit_info.ext_P[1] = mass[1];
    //c++ 126 || r 127
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, sum_corr_weight_shift<11, 11, 11> /*reim*/, "E2_01_ws", fit_info, jack_file);

    error(corr_counter != 126, -1, "correlator counter wrong", "corr_counter=%d     expected=126", corr_counter);

    double* E2_01_ws = malloc_copy_jackboot(Njack, fit_out.P[0]);
    double* a_01_ws = scattering_len_luscher(Njack, mass[0], mass[1], fit_out.P[0], params.data.L[1]);
    tmpj = (double*)malloc(sizeof(double) * Njack);
    tmp_muj = (double*)malloc(sizeof(double) * Njack);
    sub_jackboot(Njack, tmpj, fit_out.P[0], mass[0]);
    sub_jackboot(Njack, tmpj, tmpj, mass[1]);
    fprintf(outfile, "#scattering length  a  err deltaE2 err    mu01 err    deltaE2*mu01  err   a_01*(m0+m1)=-mu/2pi  err      a_01*(m0+m1)/(a0m0)=4/3 a_01*pi/mu   err    \n %g  %g     %g  %g\t",
        a_01_ws[Njack - 1], error_jackboot(option[4], Njack, a_01_ws), tmpj[Njack - 1], error_jackboot(option[4], Njack, tmpj));


    //reduced mass
    for (j = 0; j < Njack;j++)
        tmp_muj[j] = mass[0][j] * mass[1][j] / (mass[0][j] + mass[1][j]);
    fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

    //reduced mass time DeltaE2
    for (j = 0; j < Njack;j++)
        tmp_muj[j] = tmp_muj[j] * tmpj[j];
    fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));


    for (j = 0; j < Njack;j++)
        tmp_muj[j] = a_01_ws[j] * (mass[0][j] + mass[1][j]);
    fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

    for (j = 0; j < Njack;j++)
        tmp_muj[j] = (a_01_ws[j] * (mass[0][j] + mass[1][j])) / a0m0[j];
    fprintf(outfile, "%g   %g\t", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));

    // a_01*pi/mu
    for (j = 0; j < Njack;j++)
        tmp_muj[j] = a_01_ws[j] * pi_greco * (mass[0][j] + mass[1][j]) / (mass[0][j] * mass[1][j]);
    fprintf(outfile, "%g   %g\n", tmp_muj[Njack - 1], error_jackboot(option[4], Njack, tmp_muj));


    free(tmpj); free(tmp_muj);

    free_fit_result(fit_info, fit_out);
    fflush(outfile);
    fit_info.restore_default();
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("ncorr %d\n", params.data.ncorr);
    int ncorr_new = params.data.ncorr;
    int id_GEVP_031;
    if (params.data.ncorr > 148) {


        fit_info.N = 3;
        fit_info.corr_id = { 0, 129, 127,   5, 128,    1 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        fit_info.GEVP_ignore_warning_after_t = 10;
        fit_info.verbosity = 2;
        //fit_info.corr_id={1,2};
        printf("GEVP_phi0_phi03_phi1\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);
        id_GEVP_031 = ncorr_new;

        fit_info.restore_default();

        //c++ 127 || r 128
        double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_phi0_phi03_phi1_l0", identity, jack_file);

        //c++ 128 || r 129
        double* l1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_phi0_phi03_phi1_l1", identity, jack_file);

        //c++ 129 || r 130
        double* l2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_phi0_phi03_phi1_l2", identity, jack_file);

        free(l0_GEVP);free(l1_GEVP);free(l2_GEVP);

        //c++ 130 || r 131
        double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_phi0_phi03_phi1_meffl0", M_eff_T, jack_file);

        //c++ 131 || r 132
        double* ml1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_phi0_phi03_phi1_meffl1", M_eff_T, jack_file);
        check_correlatro_counter(131);

        //c++ 132 || r 133
        double* ml2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_phi0_phi03_phi1_meffl2", M_eff_T, jack_file, fit_info);

        free(ml0_GEVP);free(ml1_GEVP);free(ml2_GEVP);

        //c++ 133 || r 134
        fit_info.Nvar = 1;
        fit_info.Npar = 2;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = fun_l1_GEVP;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);


        fit_info.ext_P[0] = ml0_GEVP;

        fit_info.repeat_start = 10;
        fit_info.acc = 1e-6;
        fit_info.h = 1e-3;
        fit_info.devorder = 4;
        fit_info.lambda = 1e-5;
        fit_info.chi2_gap_jackboot = 1;
        fit_info.guess_per_jack = 10;
        fit_info.precision_sum = 2;
        fit_info.guess = { 0.4,0.22 };

        fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, ncorr_new - 2, 0/*reim*/, "fun_l1_GEVP", fit_info, jack_file);
        fit_info.restore_default();
        check_correlatro_counter(133);

        if (params.data.ncorr > 151) {
            //c++ 134 || r 135
            double* sE1_0 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 151, "sE1_0", M_eff_T, jack_file);
            //c++ 135  || r 136
            double* sE2_0 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 155, "sE2_0", shift_and_M_eff_sinh_T, jack_file);

            fit_info.Nvar = 1;
            fit_info.Npar = 4;
            fit_info.N = 1;
            fit_info.Njack = Njack;
            fit_info.function = C3_vev;
            fit_info.n_ext_P = 2;
            fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);

            //c++ 136|| r 137
            fit_info.ext_P[0] = mass[0];
            fit_info.ext_P[1] = E2[0];
            file_head.k[2] = mu1;    file_head.k[3] = mu1;
            fit_info.guess = { 0.3961,    0.00018161,    0.000337   , 0.00022318 };

            fit_info.h = 0.0001;
            fit_info.acc = 1e-5;
            fit_info.repeat_start = 10;
            fit_info.precision_sum = 2;

            fit_out = fit_function_to_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, 156, 0/*reim*/, "sE3_0_vev", fit_info, jack_file);
            free_fit_result(fit_info, fit_out);
            fit_info.restore_default();
        }
        else { for (int i = 134;i < 137;i++)  zero_corr(zeros, Njack, jack_file); }

    }
    else { for (int i = 127;i < 137;i++)  zero_corr(zeros, Njack, jack_file); }


    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.function = constant_fit;
    fit_info.n_ext_P = 0;

    // c++ 137|| r138
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, two_to_two_con<2, 0>, "E2_0_con", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.function = constant_fit;
    fit_info.n_ext_P = 0;

    // c++ 138|| r139
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, one_to_one_sq<0>, "E1_0_sq", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);
    error(corr_counter != 138, -1, "correlator counter wrong", "corr_counter=%d", corr_counter);

    if (params.data.ncorr > 160) {
        //c++ 139 || r 140
        double* sE1_0_p1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 161, "sE1_0_p1", M_eff_T, jack_file);
        error(corr_counter != 139, -1, "correlator counter wrong", "corr_counter=%d", corr_counter);

    }
    else { for (int i = 139;i < 140;i++)  zero_corr(zeros, Njack, jack_file); }

    error(corr_counter != 139, -1, "correlator counter wrong", "corr_counter=%d     expected=139", corr_counter);

    if (params.data.ncorr > 161) {
        //c++ 140 || r 141
        double* sE1_0_p1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 162, "E2_0_p0", shift_and_M_eff_sinh_T, jack_file);
        error(corr_counter != 140, -1, "correlator counter wrong", "corr_counter=%d", corr_counter);

    }
    else { for (int i = 140;i < 141;i++)  zero_corr(zeros, Njack, jack_file); }

    check_correlatro_counter(140);
    bool generated_GEVP_phi0_phi03_phi1_p1 = false;

    if (params.data.ncorr >= 135 + 1e+6) {
        fit_info.N = 3;
        //     fit_info.t0_GEVP=10;
        fit_info.corr_id = { 33, 136, 130,   95, 133,    34,
                          35, 137, 131,   96, 134,    36,
                          37, 138, 132,   97, 135,    38 };
        printf("GEVP_phi0_phi03_phi1_p1\n");

        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix_p1, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);
        id_GEVP_031_p1 = ncorr_new;
        fit_info.restore_default();

        generated_GEVP_phi0_phi03_phi1_p1 = true;

        //c++ 141 || r 142
        double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_phi0_phi03_phi1_p1_meffl0", M_eff_T, jack_file);

        //c++ 142 || r 143
        double* ml1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_phi0_phi03_phi1_p1_meffl1", M_eff_T, jack_file);

        //c++ 143 || r 144
        double* ml2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_phi0_phi03_phi1_p1_meffl2", M_eff_T, jack_file);
        free(ml0_GEVP);free(ml1_GEVP);free(ml2_GEVP);
        check_correlatro_counter(143);

    }
    else { for (int i = 141;i < 144;i++)  zero_corr(zeros, Njack, jack_file); }

    check_correlatro_counter(143);


    ///////////////////
    printf("%g\n", params.data.gC);
    if (params.data.ncorr > 165 && fabs(params.data.gC - 0.025) > 1e-4) {


        fit_info.N = 4;
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={0, 129, 127,   5, 128,    1 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }
        fit_info.corr_id = { 0, 129, 127, 163,    5, 128, 165,     1, 164,    116 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1, A1A1 }
        //fit_info.corr_id={1,2};
        printf(" GEVP_0_3_1_A1 \n");

        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);
        fit_info.restore_default();


        //c++ 144 || r 145
        double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 4, "GEVP_0_3_1_A1_meffl0", M_eff_T, jack_file);

        //c++ 145 || r 146
        double* ml1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_0_3_1_A1_meffl1", M_eff_T, jack_file);

        //c++ 146 || r 147
        double* ml2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_0_3_1_A1_meffl2", M_eff_T, jack_file);
        error(corr_counter != 146, -1, "correlator counter wrong", "corr_counter=%d", corr_counter);

        //c++ 147 || r 148
        double* ml3_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_0_3_1_A1_meffl3", M_eff_T, jack_file);
        free(ml0_GEVP);free(ml1_GEVP);free(ml2_GEVP);free(ml3_GEVP);
        error(corr_counter != 147, -1, "correlator counter wrong", "corr_counter=%d", corr_counter);

    }
    else { for (int i = 144;i < 148;i++)  zero_corr(zeros, Njack, jack_file); }




    ///////////////////
    if (params.data.ncorr > 166 + 1e+6) {


        fit_info.N = 4;
        //fit_info.corr_id={33, 136, 130,   95, 133,    34   }
        fit_info.t0_GEVP = 1;

        fit_info.corr_id = { 33, 136, 130, 167,   95, 133, 169,    34, 168,    166 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1, A1A1 }
        //fit_info.corr_id={1,2};

        fit_info.corr_id = { 33, 136, 130, 167,   95, 133,  169,   34, 168,    166,
                          35, 137, 131, 167,   96, 134,  169,   36, 168,    166,
                          37, 138, 132, 167,   97, 135,  169,   38,  168,    166 };

        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix_4_p1, fit_info);

        printf("GEVP_0_3_1_A1_p1\n");

        printf(" ncorr after GEVP %d\n", ncorr_new);
        fit_info.restore_default();


        //c++ 148 || r 149
        double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 4, "GEVP_0_3_1_A1_p1_meffl0", M_eff_T, jack_file);

        //c++ 149 || r 150
        double* ml1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_0_3_1_A1_p1_meffl1", M_eff_T, jack_file);

        //c++ 150 || r 151
        double* ml2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_0_3_1_A1_p1_meffl2", M_eff_T, jack_file);
        check_correlatro_counter(150);

        //c++ 151 || r 152
        double* ml3_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_0_3_1_A1_p1_meffl3", M_eff_T, jack_file);
        check_correlatro_counter(151);

        free(ml0_GEVP);free(ml1_GEVP);free(ml2_GEVP);free(ml3_GEVP);



    }
    else { for (int i = 148;i < 152;i++)  zero_corr(zeros, Njack, jack_file); }


    if (params.data.ncorr >= 34) {
        //c++ 152
        double* E1_1_p1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 34, "E1_1_p1", M_eff_T, jack_file);


        free(E1_1_p1);
    }
    else { for (int i = 152;i < 153;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(152);

    ///////////////////
    if (params.data.ncorr > 95) {


        //c++ 153
        double* phi03_meff = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5, "phi03_meff", M_eff_T, jack_file);
        free(phi03_meff);
        error(corr_counter != 153, -1, "correlator counter wrong", "corr_counter=%d", corr_counter);


        //c++ 154
        double* phi03_p1_meff = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 95, "phi03_p1_meff", M_eff_T, jack_file);
        free(phi03_p1_meff);
        error(corr_counter != 154, -1, "correlator counter wrong", "corr_counter=%d", corr_counter);



    }
    else { for (int i = 153;i < 155;i++)  zero_corr(zeros, Njack, jack_file); }


    check_correlatro_counter(154);





    if (params.data.ncorr >= 144 && params.data.gC != 0.025) {
        fit_info.N = 3;
        //      fit_info.corr_id={33, 136, 130,   95, 133,    34,      35,37,  96,97, 36,38  };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }
        //      add_correlators(option , ncorr_new , conf_jack ,GEVP_matrix_p1 ,   fit_info );

        //     fit_info.corr_id={75, 145, 139,   104, 142,    76 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }
        fit_info.corr_id = { 75, 145, 139,   104, 142,    76,
                          77, 146, 140,   105, 143,    78,
                          79, 147, 141,   106, 144,    80 };

        printf("GEVP_phi0_phi03_phi1_p11\n");

        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix_p11, fit_info);
        //     add_correlators(option , ncorr_new , conf_jack ,GEVP_matrix ,   fit_info );
        printf(" ncorr after GEVP %d\n", ncorr_new);
        fit_info.restore_default();

        //c++ 155 || r 142
        double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_phi0_phi03_phi1_p11_meffl0", M_eff_T, jack_file);
        check_correlatro_counter(155);

        //c++ 156 || r 143
        double* ml1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_phi0_phi03_phi1_p11_meffl1", M_eff_T, jack_file);
        check_correlatro_counter(156);

        //c++ 157 || r 144
        double* ml2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_phi0_phi03_phi1_p11_meffl2", M_eff_T, jack_file);
        free(ml0_GEVP);free(ml1_GEVP);free(ml2_GEVP);
        check_correlatro_counter(157);

    }
    else { for (int i = 155;i < 158;i++)  zero_corr(zeros, Njack, jack_file); }

    check_correlatro_counter(157);

    if (params.data.ncorr >= 76) {

        double* E1_1_pxy = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 76, "E1_1_pxy", M_eff_T, jack_file);
        free(E1_1_pxy);
        check_correlatro_counter(158);


    }
    else { for (int i = 158;i < 159;i++)  zero_corr(zeros, Njack, jack_file); }

    check_correlatro_counter(158);


    if (params.data.ncorr >= 135 && generated_GEVP_phi0_phi03_phi1_p1) {
        //id_GEVP_031_p1  should be 173 !

        double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, id_GEVP_031_p1 - 3, "GEVP_phi0_phi03_phi1_p1_l0", identity, jack_file);
        check_correlatro_counter(159);


        double* l1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, id_GEVP_031_p1 - 2, "GEVP_phi0_phi03_phi1_p1_l1", identity, jack_file);
        check_correlatro_counter(160);


        double* l2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, id_GEVP_031_p1 - 1, "GEVP_phi0_phi03_phi1_p1_l2", identity, jack_file);
        check_correlatro_counter(161);

        free(l0_GEVP);free(l1_GEVP);free(l2_GEVP);


    }
    else { for (int i = 159;i < 162;i++)  zero_corr(zeros, Njack, jack_file); }

    check_correlatro_counter(161);

    if (params.data.ncorr >= 135 + 1e+16) {
        int sqrtN = 3;
        fit_info.t0_GEVP = 3;
        fit_info.N = sqrtN * sqrtN;

        fit_info.corr_id = { 33, 136, 130,   95, 133,    34,
                          35, 137, 131,   96, 134,    36,
                          37, 138, 132,   97, 135,    38 };
        printf("GEVP_vectors_phi0_phi03_phi1_p1\n");
        fit_info.value_or_vector = 1;
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix_p1, fit_info);
        printf(" ncorr after GEVP vectors %d\n", ncorr_new);
        //id_GEVP_031_p1=ncorr_new;
        fit_info.restore_default();
        char name_v[NAMESIZE];

        for (int i = 0; i < fit_info.N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_phi0_phi03_phi1_p1_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);

        }
        check_correlatro_counter(170);


    }
    else { for (int i = 162;i < 171;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(170);




    if (params.data.ncorr >= 135 + 1e+6) {
        fit_info.N = 3;
        fit_info.GEVP_tpt0 = true;
        fit_info.corr_id = { 33, 136, 130,   95, 133,    34,
                          35, 137, 131,   96, 134,    36,
                          37, 138, 132,   97, 135,    38 };
        printf("GEVP_tpt0_phi0_phi03_phi1_p1\n");

        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix_p1, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);
        //id_GEVP_031_p1=ncorr_new;
        fit_info.restore_default();


        double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, id_GEVP_031_p1 - 3, "GEVP_tpt0_phi0_phi03_phi1_p1_l0", identity, jack_file);
        check_correlatro_counter(171);


        double* l1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, id_GEVP_031_p1 - 2, "GEVP_tpt0_phi0_phi03_phi1_p1_l1", identity, jack_file);
        check_correlatro_counter(172);


        double* l2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, id_GEVP_031_p1 - 1, "GEVP_tpt0_phi0_phi03_phi1_p1_l2", identity, jack_file);
        check_correlatro_counter(173);

        free(l0_GEVP);free(l1_GEVP);free(l2_GEVP);


        double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_tpt0_phi0_phi03_phi1_p1_meffl0", M_eff_T, jack_file);
        check_correlatro_counter(174);


        double* ml1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_tpt0_phi0_phi03_phi1_p1_meffl1", M_eff_T, jack_file);
        check_correlatro_counter(175);

        //c++ 143 || r 144
        double* ml2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_tpt0_phi0_phi03_phi1_p1_meffl2", M_eff_T, jack_file);
        free(ml0_GEVP);free(ml1_GEVP);free(ml2_GEVP);
        check_correlatro_counter(176);


        int sqrtN = 3;
        fit_info.t0_GEVP = 6;
        fit_info.N = sqrtN * sqrtN;
        fit_info.GEVP_tpt0 = true;
        fit_info.corr_id = { 33, 136, 130,   95, 133,    34,
                          35, 137, 131,   96, 134,    36,
                          37, 138, 132,   97, 135,    38 };
        printf("GEVP_tpt0_vectors_phi0_phi03_phi1_p1\n");
        fit_info.value_or_vector = 1;
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix_p1, fit_info);
        printf(" ncorr after GEVP vectors %d\n", ncorr_new);
        //id_GEVP_031_p1=ncorr_new;
        fit_info.restore_default();
        char name_v[NAMESIZE];

        for (int i = 0; i < fit_info.N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_tpt0_phi0_phi03_phi1_p1_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);

        }
        check_correlatro_counter(185);

    }
    else { for (int i = 171;i < 186;i++)  zero_corr(zeros, Njack, jack_file); }

    check_correlatro_counter(185);

    if (params.data.ncorr > 148) {

        int T = params.data.L[0];
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.corr_id = { id_GEVP_031 - 1,id_GEVP_031 - 2 };
        fit_info.function = constant_fit;

        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, diff_meff, "diff_GEVP_031_m2-m1", fit_info, jack_file);
        check_correlatro_counter(186);
        free_fit_result(fit_info, fit_out);

        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, ratio_meff, "ratio_GEVP_031_m2-m1", fit_info, jack_file);
        check_correlatro_counter(187);
        free_fit_result(fit_info, fit_out);

        fit_info.restore_default();

    }
    else { for (int i = 186;i < 188;i++) zero_corr(zeros, Njack, jack_file); }

    check_correlatro_counter(187);



    if (params.data.ncorr >= 135 && generated_GEVP_phi0_phi03_phi1_p1) {

        int T = params.data.L[0];
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.corr_id = { id_GEVP_031_p1 - 1,id_GEVP_031_p1 - 2 };
        fit_info.function = constant_fit;

        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, diff_meff, "diff_GEVP_031_p1_m2-m1", fit_info, jack_file);
        check_correlatro_counter(188);
        free_fit_result(fit_info, fit_out);

        fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, ratio_meff, "ratio_GEVP_031_p1_m2-m1", fit_info, jack_file);
        check_correlatro_counter(189);
        free_fit_result(fit_info, fit_out);
        fit_info.restore_default();

    }
    else { for (int i = 188;i < 190;i++)  zero_corr(zeros, Njack, jack_file); }

    check_correlatro_counter(189);


    // printf("id_GEVP_031=%d\n", id_GEVP_031);





    if (params.data.ncorr >= 148) {

        int sqrtN = 3;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127,   5, 128,    1 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }

        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_phi0_phi03_phi1_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);


        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;
        char name_v[NAMESIZE];
        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_phi0_phi03_phi1_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);

        }
        check_correlatro_counter(198);

        // FILE *ftmp=open_file("jackknife_phi0phi0.txt","w+");
        // for (int t=0; t< params.data.L[0];t++ ){
        //     fprintf(ftmp,"%d\t",t );
        //     for (int j=0; j<Njack; j++){
        //         fprintf(ftmp,"%.15f\t ", conf_jack[j][0][t][0]);
        //     }
        //     fprintf(ftmp,"\n");
        // }
        // fclose(ftmp);
        // ftmp=open_file("jackknife_phi1phi1.txt","w+");
        // for (int t=0; t< params.data.L[0];t++ ){
        //     fprintf(ftmp,"%d\t",t );
        //     for (int j=0; j<Njack; j++){
        //         fprintf(ftmp,"%.15f\t ", conf_jack[j][1][t][0]);
        //     }
        //     fprintf(ftmp,"\n");
        // }
        // fclose(ftmp);
        // ftmp=open_file("jackknife_phi0pphi0-pphi0.txt","w+");
        // for (int t=0; t< params.data.L[0];t++ ){
        //     fprintf(ftmp,"%d\t",t );
        //     for (int j=0; j<Njack; j++){
        //         fprintf(ftmp,"%.15f\t ", conf_jack[j][116][t][0]);
        //     }
        //     fprintf(ftmp,"\n");
        // }
        // fclose(ftmp);
        // ftmp=open_file("jackknife_phi0phi0phi0.txt","w+");
        // for (int t=0; t< params.data.L[0];t++ ){
        //     fprintf(ftmp,"%d\t",t );
        //     for (int j=0; j<Njack; j++){
        //         fprintf(ftmp,"%.15f\t ", conf_jack[j][5][t][0]);
        //     }
        //     fprintf(ftmp,"\n");
        // }
        // fclose(ftmp);

    }
    else { for (int i = 190;i < 199;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(198);

    double* E3_0_meff = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5, "E3_0_meff", M_eff_T, jack_file, fit_info);
    free(E3_0_meff);
    check_correlatro_counter(199);

    if (params.data.ncorr >= 175) {
        double* E5_0_meff = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 175, "E5_0_meff", M_eff_T, jack_file, fit_info);
        free(E5_0_meff);
        check_correlatro_counter(200);
    }
    else { zero_corr(zeros, Njack, jack_file); }
    if (params.data.ncorr >= 170) {
        double* E3_0_meff = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 170, "E3_0_loc_meff", M_eff_T, jack_file, fit_info);
        free(E3_0_meff);
        check_correlatro_counter(201);
    }
    else { zero_corr(zeros, Njack, jack_file); }


    if (params.data.ncorr >= 177) {



        fit_info.N = 4;
        fit_info.corr_id = { 0, 129, 127,176,     5, 128,177,    1,178,     175 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }

        //     fit_info.N=3;
        //     fit_info.corr_id={170, 172, 173,     5, 128,    1};//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }

        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_5\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();


        double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 4, "GEVP_0_3_1_5_l0", identity, jack_file);
        check_correlatro_counter(202);


        double* l1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_0_3_1_5_l1", identity, jack_file);
        check_correlatro_counter(203);


        double* l2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_0_3_1_5_l2", identity, jack_file);
        check_correlatro_counter(204);
        double* l3_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_0_3_1_5_l3", identity, jack_file);
        check_correlatro_counter(205);
        free(l3_GEVP);

        free(l0_GEVP);free(l1_GEVP);free(l2_GEVP);

        //c++ 130 || r 131
        double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 4, "GEVP_0_3_1_5_meffl0", M_eff_T, jack_file);
        check_correlatro_counter(206);

        //c++ 131 || r 132
        double* ml1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_0_3_1_5_meffl1", M_eff_T, jack_file);
        check_correlatro_counter(207);

        //c++ 132 || r 133
        double* ml2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_0_3_1_5_meffl2", M_eff_T, jack_file, fit_info);
        check_correlatro_counter(208);

        double* ml3_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_0_3_1_5_meffl3", M_eff_T, jack_file, fit_info);
        check_correlatro_counter(209);
        free(ml3_GEVP);

        free(ml0_GEVP);free(ml1_GEVP);free(ml2_GEVP);


        int sqrtN = 4;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127,176,     5, 128,177,    1,178,     175 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }
        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;

        printf("GEVP_0_3_1_5_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        char name_v[NAMESIZE];
        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_5_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);

        }
        check_correlatro_counter(225);

        fit_info.restore_default();




    }
    else { for (int i = 202;i < 226;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(225);

    if (params.data.ncorr >= 170) {

        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.corr_id = { 170, 0 };
        fit_info.function = constant_fit;

        fit_result ratio_E3_loc_noloc = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, ratio_corr, "ratio_E3_loc_noloc", fit_info, jack_file);
        free_fit_result(fit_info, ratio_E3_loc_noloc);
        fit_info.restore_default();
        check_correlatro_counter(226);

    }
    else { zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(226);

    fit_info.restore_default();

    if (params.data.ncorr >= 177 + 1e+6) {



        fit_info.N = 3;
        fit_info.corr_id = { 0, 171, 127,     170, 173,    1 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }

        //     fit_info.N=3;
        //     fit_info.corr_id={170, 172, 173,     5, 128,    1};//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }

        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3loc_1\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();


        double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_0_3loc_1_l0", identity, jack_file);
        check_correlatro_counter(227);


        double* l1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_0_3loc_1_l1", identity, jack_file);
        check_correlatro_counter(228);


        double* l2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_0_3loc_1_l2", identity, jack_file);
        check_correlatro_counter(229);
        ;

        free(l0_GEVP);free(l1_GEVP);free(l2_GEVP);

        //c++ 130 || r 131
        double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_0_3loc_1_meffl0", M_eff_T, jack_file);
        check_correlatro_counter(230);

        //c++ 131 || r 132
        double* ml1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_0_3loc_1_meffl1", M_eff_T, jack_file);
        check_correlatro_counter(231);

        //c++ 132 || r 133
        double* ml2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_0_3loc_1_meffl2", M_eff_T, jack_file, fit_info);
        check_correlatro_counter(232);


        free(ml0_GEVP);free(ml1_GEVP);free(ml2_GEVP);


        int sqrtN = 3;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 171, 127,     170, 173,    1 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }    fit_info.value_or_vector=1; // 0= values
        fit_info.t0_GEVP = 3;

        printf("GEVP_0_3loc_1_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        char name_v[NAMESIZE];
        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3loc_1_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);

        }
        check_correlatro_counter(241);
    }
    else { for (int i = 227;i < 242;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(241);





    if (params.data.ncorr >= 173 && conf_jack[Njack - 1][170][0][0] != 0) {

        fit_info.N = 4;
        fit_info.corr_id = { 0, 129, 127,171,     5, 128,172,    1,173,     170 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 
        //     fit_info.N=3;
        //     fit_info.corr_id={170, 172, 173,     5, 128,    1};//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_3l\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();

        double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 4, "GEVP_0_3_1_3l_l0", identity, jack_file);
        check_correlatro_counter(242);

        double* l1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 3, "GEVP_0_3_1_3l_l1", identity, jack_file);
        check_correlatro_counter(243);

        double* l2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 2, "GEVP_0_3_1_3l_l2", identity, jack_file);
        check_correlatro_counter(244);
        double* l3_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - 1, "GEVP_0_3_1_3l_l3", identity, jack_file);
        check_correlatro_counter(245);
        free(l3_GEVP);

        free(l0_GEVP);free(l1_GEVP);free(l2_GEVP);

        double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, ncorr_new - 4, "GEVP_0_3_1_3l_meffl0", M_eff_T, jack_file);
        check_correlatro_counter(246);

        double* ml1_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, ncorr_new - 3, "GEVP_0_3_1_3l_meffl1", M_eff_T, jack_file);
        check_correlatro_counter(247);

        double* ml2_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, ncorr_new - 2, "GEVP_0_3_1_3l_meffl2", M_eff_T, jack_file, fit_info);
        check_correlatro_counter(248);

        double* ml3_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, ncorr_new - 1, "GEVP_0_3_1_3l_meffl3", M_eff_T, jack_file, fit_info);
        check_correlatro_counter(249);
        free(ml3_GEVP);

        free(ml0_GEVP);free(ml1_GEVP);free(ml2_GEVP);

        int sqrtN = 4;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127,171,     5, 128,172,    1,173,     170 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }
        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;

        printf("GEVP_0_3_1_3l_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        char name_v[NAMESIZE];
        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_3l_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);

        }
        check_correlatro_counter(265);

        fit_info.restore_default();

    }
    else { for (int i = 242;i < 266;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(265);



    if (params.data.ncorr >= 180) {
        //diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1, A1A1, phi5phi5 }
        fit_info.N = 5;
        fit_info.corr_id = { 0, 129, 127, 163, 176,
                               5, 128, 165, 177,
                                    1, 164, 178,
                                       116, 179,
                                           175 };



        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_A1_5\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 5;
        for (int i = 0; i < fit_info.N; i++) {
            char name_v[NAMESIZE];
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_A1_5_l%d", i);
            double* li_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            check_correlatro_counter(266 + i);
            free(li_GEVP);
        }

        check_correlatro_counter(270);
        fit_info.N = 5;
        for (int i = 0; i < fit_info.N; i++) {
            char name_v[NAMESIZE];
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_A1_5_meffl%d", i);
            double* mi_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            check_correlatro_counter(271 + i);
            free(mi_GEVP);
        }
        check_correlatro_counter(275);

        int sqrtN = 5;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127, 163, 176,
                               5, 128, 165, 177,
                                    1, 164, 178,
                                       116, 180,
                                           175 };
        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;

        printf("GEVP_0_3_1_A1_5_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        char name_v[NAMESIZE];
        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_A1_5_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(276 + i);

        }
        check_correlatro_counter(300);

        fit_info.restore_default();

    }
    else { for (int i = 266;i < 301;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(300);




    if (params.data.ncorr >= 116) {
        double* EA1_0_meff = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 116, "EA1_0_meff", M_eff_T, jack_file, fit_info);
        free(EA1_0_meff);
        check_correlatro_counter(301);
    }
    else { zero_corr(zeros, Njack, jack_file); }

    if (params.data.ncorr >= 181) {
        double* E3_0_2loc_meff = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 181, "E3_0_2loc_meff", M_eff_T, jack_file, fit_info);
        free(E3_0_2loc_meff);
        check_correlatro_counter(302);
    }
    else { zero_corr(zeros, Njack, jack_file); }

    if (params.data.ncorr > 161) {
        //c++ 140 || r 141
        double* sE1_0_p1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 162, "E2_0_loc_meff", M_eff_T, jack_file);
        check_correlatro_counter(303);
    }
    else { zero_corr(zeros, Njack, jack_file); }





    if (params.data.ncorr > 1) {
        //c++ 140 || r 141
        double* E1_0_der2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 0, "E1_0_der2", laplacian_M_eff_T, jack_file);
        check_correlatro_counter(304);
        free(E1_0_der2);

    }
    else { zero_corr(zeros, Njack, jack_file); }


    if (params.data.ncorr >= 5) {
        //c++ 140 || r 141
        double* E3_0_der2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5, "E3_0_der2", laplacian_M_eff_T, jack_file);
        check_correlatro_counter(305);
        free(E3_0_der2);

    }
    else { zero_corr(zeros, Njack, jack_file); }


    if (params.data.ncorr > 1) {
        //c++ 140 || r 141
        double* E1_0_der2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 0, "E1_0_der2corr", der2corr_M_eff_T, jack_file);
        check_correlatro_counter(306);
        free(E1_0_der2);

    }
    else { zero_corr(zeros, Njack, jack_file); }


    if (params.data.ncorr >= 5) {
        //c++ 140 || r 141
        double* E3_0_der2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5, "E3_0_der2corr", der2corr_M_eff_T, jack_file);
        check_correlatro_counter(307);
        free(E3_0_der2);

    }
    else { zero_corr(zeros, Njack, jack_file); }



    if (params.data.ncorr > 1) {
        //c++ 140 || r 141
        double* E1_0_der2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 0, "E1_0_combtsym", combt_M_eff_T, jack_file);
        check_correlatro_counter(308);
        free(E1_0_der2);

    }
    else { zero_corr(zeros, Njack, jack_file); }


    if (params.data.ncorr >= 5) {
        //c++ 140 || r 141
        double* E3_0_der2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5, "E3_0_combtsym", combt_M_eff_T, jack_file);
        check_correlatro_counter(309);
        free(E3_0_der2);

    }
    else { zero_corr(zeros, Njack, jack_file); }



    if (params.data.ncorr > 1) {
        //c++ 140 || r 141
        double* E1_0_der2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 0, "E1_0_combt2sym", combt2_M_eff_T, jack_file);
        check_correlatro_counter(310);
        free(E1_0_der2);

    }
    else { zero_corr(zeros, Njack, jack_file); }


    if (params.data.ncorr >= 5) {
        //c++ 140 || r 141
        double* E3_0_der2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5, "E3_0_combt2sym", combt2_M_eff_T, jack_file);
        check_correlatro_counter(311);
        free(E3_0_der2);
    }
    else { zero_corr(zeros, Njack, jack_file); }




    if (params.data.ncorr > 1) {
        //c++ 140 || r 141
        double* E1_0_der2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 0, "E1_0_combt15sym", combt15_M_eff_T, jack_file);
        check_correlatro_counter(312);
        free(E1_0_der2);
    }
    else { zero_corr(zeros, Njack, jack_file); }


    if (params.data.ncorr >= 5) {
        //c++ 140 || r 141
        double* E3_0_der2 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5, "E3_0_combt15sym", combt15_M_eff_T, jack_file);
        check_correlatro_counter(313);
        free(E3_0_der2);

    }
    else { zero_corr(zeros, Njack, jack_file); }


    if (params.data.ncorr > 148) {

        fit_info.N = 1;
        fit_info.T = params.data.L[0];
        fit_info.corr_id = { 0 };
        add_correlators(option, ncorr_new, conf_jack, der2_der2_corr, fit_info);
        int id_d0d0 = ncorr_new - 1;

        fit_info.corr_id = { 0 };
        add_correlators(option, ncorr_new, conf_jack, der2_corr, fit_info);
        int id_d00 = ncorr_new - 1;

        fit_info.corr_id = { 128 };
        add_correlators(option, ncorr_new, conf_jack, der2_corr, fit_info);
        int id_d01 = ncorr_new - 1;



        char name_v[NAMESIZE];
        fit_info.N = 3;
        fit_info.corr_id = { 0, id_d00,  127,
                             id_d0d0, id_d01,
                                      1 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_d0_1\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 3;
        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_d0_1_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(314 + i);
        }
        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_d0_1_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(317 + i);
        }


        int sqrtN = 3;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, id_d00,  127,
                             id_d0d0, id_d01,
                                      1 };//diag{ phi0->phi0, 3phi0->3phi0, phi1->phi1 }
        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_d0_1_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_d0_1_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(320 + i);
        }
        check_correlatro_counter(328);
    }
    else { for (int i = 314;i <= 328;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(328);




    if (params.data.ncorr > 175) {




        char name_v[NAMESIZE];
        fit_info.N = 3;
        fit_info.corr_id = { 0, 129,  176,
                               5,  177,
                                   175 };//diag{ phi0->phi0, 3phi0->3phi0, phi5->phi5 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_5\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 3;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_5_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(329 + i);
        }
        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_5_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(332 + i);
        }


        int sqrtN = 3;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129,  176,
                               5,  177,
                                   175 };//diag{ phi0->phi0, 3phi0->3phi0, phi5->phi5 }
        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_5_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_5_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(335 + i);
        }
        check_correlatro_counter(343);
    }
    else { for (int i = 329;i <= 343;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(343);





    if (params.data.ncorr > 175) {

        char name_v[NAMESIZE];
        fit_info.N = 2;
        fit_info.corr_id = { 0, 129,
                               5 };//diag{ phi0->phi0, 3phi0->3phi0 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 2;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(344 + i);
        }
        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(346 + i);
        }


        int sqrtN = 2;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129,
                               5 };//diag{ phi0->phi0, 3phi0->3phi0 }
        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(348 + i);
        }
        check_correlatro_counter(351);
    }
    else { for (int i = 344;i <= 351;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(351);


    if (params.data.ncorr > 175) {

        char name_v[NAMESIZE];
        fit_info.N = 2;
        fit_info.corr_id = { 5, 128,
                               1 };//diag{ 3phi0->3phi0, phi1->phi1 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_3_1\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 2;
        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_3_1_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(352 + i);
        }
        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_3_1_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(354 + i);
        }


        int sqrtN = 2;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 5, 128,
                                1 };//diag{ 3phi0->3phi0, phi1->phi1 }
        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_3_1_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_3_1_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(356 + i);
        }
        check_correlatro_counter(359);
    }
    else { for (int i = 352;i <= 359;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(359);


    if (params.data.ncorr > 175) {

        char name_v[NAMESIZE];
        fit_info.N = 2;
        fit_info.corr_id = { 0, 127,
                               1 };//diag{ phi0->phi0, phi1->phi1 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_1\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 2;
        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_1_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(360 + i);
        }
        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_1_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(362 + i);
        }


        int sqrtN = 2;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 127,
                               1 };//diag{ phi0->phi0, phi1->phi1 }
        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_1_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_1_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(364 + i);
        }
        check_correlatro_counter(367);
    }
    else { for (int i = 360;i <= 367;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(367);



    if (params.data.ncorr >= 170 + 1e+6) {

        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.HENKEL_size = 2;
        fit_info.t0_GEVP = 1;
        fit_info.function = constant_fit;
        fit_info.corr_id = { 0 };
        char name_h[NAMESIZE];
        for (int i = 0; i < 3; i++) {
            fit_info.corr_id[0] = id_GEVP_031 - 3 + i;
            mysprintf(name_h, NAMESIZE, "hankel_GEVP_0_3_1_meffl%d", i);
            fit_result hankel_GEVP_0_3_1_meffl = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, hankel, name_h, fit_info, jack_file);
            free_fit_result(fit_info, hankel_GEVP_0_3_1_meffl);
            check_correlatro_counter(368 + i);
        }
        fit_info.restore_default();
    }
    else { for (int i = 0; i < 3; i++) zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(370);




    if (params.data.ncorr > 175) {

        char name_v[NAMESIZE];
        fit_info.N = 3;
        fit_info.HENKEL_size = 2;
        fit_info.corr_id = { 0, 129, 127,
                               5, 128,
                                   1 };//diag{ phi0->phi0, phi1->phi1 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 1;
        fit_info.GEVP_tpt0 = true;
        //fit_info.corr_id={1,2};
        printf("H2_GEVP_0_3_1\n");
        // fit_info.HENKEL_size=1;
        add_correlators(option, ncorr_new, conf_jack, HANKEL_GEVP_matrix, fit_info);

        // add_correlators(option , ncorr_new , conf_jack ,GEVP_matrix ,   fit_info );
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 3;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "H2_GEVP_0_3_1_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, log_corr, jack_file);
            // double *l0_GEVP=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile, ncorr_new-fit_info.N+i, name_v, M_eff_T,jack_file);

            free(l0_GEVP);
            check_correlatro_counter(371 + i);
        }
        check_correlatro_counter(373);


    }
    else { for (int i = 371;i <= 373;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(373);

    fit_info.restore_default();

    int id_d0d0;
    if (params.data.ncorr > 148) {

        fit_info.N = 1;
        fit_info.T = params.data.L[0];
        fit_info.corr_id = { 0 };
        add_correlators(option, ncorr_new, conf_jack, der1_der1dag_corr, fit_info);
        id_d0d0 = ncorr_new - 1;

        fit_info.corr_id = { 0 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_0d0 = ncorr_new - 1;

        fit_info.corr_id = { 129 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_3d0 = ncorr_new - 1;

        fit_info.corr_id = { 127 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_1d0 = ncorr_new - 1;

        char name_v[NAMESIZE];
        fit_info.N = 4;
        fit_info.corr_id = { 0, 129, 127, id_0d0,
                               5, 128, id_3d0 ,
                                   1 , id_1d0 ,
                                       id_d0d0 };//diag{ phi0->phi0, phi3->phi3, phi1->phi1,   dphi0->dphi0 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_d0\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 4;
        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_d0_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(374 + i);
        }
        check_correlatro_counter(377);

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_d0_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(378 + i);
        }
        check_correlatro_counter(381);

        int sqrtN = 4;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127, id_0d0,
                               5, 128, id_3d0 ,
                                   1 , id_1d0 ,
                                       id_d0d0 };//diag{ phi0->phi0, phi3->phi3, phi1->phi1,   dphi0->dphi0 }

        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_d0_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_d0_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(382 + i);
        }
        check_correlatro_counter(397);

        double* d0d0 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, id_d0d0, "d0d0_meff", M_eff_T, jack_file);
        free(d0d0);
        check_correlatro_counter(398);

    }
    else { for (int i = 374;i <= 398;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(398);



    printf("GEVPn_0_3_1\n");
    if (params.data.ncorr > 148 && strcmp(argv[4], "G2t_T64_L24_msq0-4.868000_msq1-4.710000_l01.000000_l11.000000_mu2.000000_g0.250000_rep0_bin100_merged_bin1000") == 0) {

        char name_v[NAMESIZE];
        fit_info.N = 3;
        fit_info.corr_id = { 0, 129, 127,
                               5, 128,
                                   1 };//diag{ phi0->phi0, phi3->phi3, phi1->phi1,   dphi0->dphi0 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVPn_0_3_1\n");
        add_correlators(option, ncorr_new, conf_jack, GEVPn_matrix, fit_info);
        printf(" ncorr after GEVPn %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 3;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVPn_0_3_1_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(399 + i);
        }
        check_correlatro_counter(401);

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVPn_0_3_1_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(402 + i);
        }
        check_correlatro_counter(404);

        int sqrtN = 3;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127,
                               5, 128,
                                   1 };//diag{ phi0->phi0, phi3->phi3, phi1->phi1,   dphi0->dphi0 }
        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVPn_0_3_1_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVPn_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVPn_0_3_1_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(405 + i);
        }
        check_correlatro_counter(413);
    }
    else { for (int i = 399;i <= 413;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(413);

    fit_info.restore_default();
    printf("GEVP_0_3_1_001\n");
    if (params.data.ncorr > 190) {

        char name_v[NAMESIZE];
        fit_info.N = 4;
        fit_info.corr_id = { 0, 129, 127, 188,
                               5, 128, 189 ,
                                   1 , 190 ,
                                       187 };//diag{ phi0->phi0, phi3->phi3, phi1->phi1,   dphi0->dphi0 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_001\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 4;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(414 + i);
        }
        check_correlatro_counter(417);

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(418 + i);
        }
        check_correlatro_counter(421);

        int sqrtN = 4;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127, 188,
                               5, 128, 189 ,
                                   1 , 190 ,
                                       187 };// { phi0, phi0^3, phi1, phi0^2phi1}

        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_001_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(422 + i);
        }
        check_correlatro_counter(437);

    }
    else { for (int i = 414;i <= 437;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(437);

    fit_info.restore_default();
    printf("GEVP_0_3_1_001_011\n");
    if (params.data.ncorr > 199) {

        fit_info.restore_default();
        char name_v[NAMESIZE];
        fit_info.N = 5;
        fit_info.corr_id = { 0, 129, 127, 188, 194,
                               5, 128, 189, 195,
                                   1 , 190, 196,
                                       187, 199,
                                            193 };//diag{ phi0->phi0, phi3->phi3, phi1->phi1,   dphi0->dphi0 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_001_011\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 5;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_011_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(438 + i);
        }
        check_correlatro_counter(442);

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_011_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(443 + i);
        }
        check_correlatro_counter(447);

        int sqrtN = 5;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127, 188, 194,
                               5, 128, 189, 195,
                                   1 , 190, 196,
                                       187, 199,
                                           193 };//diag{ phi0->phi0, phi3->phi3, phi1->phi1,  phi0^2phi1, phi0 phi1^2 }


        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_001_011_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_011_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(448 + i);
        }
        check_correlatro_counter(472);

    }
    else { for (int i = 438;i <= 472;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(472);



    fit_info.restore_default();
    printf("GEVP_0_3_1_001_011_111\n");
    if (params.data.ncorr > 207) {

        fit_info.restore_default();
        char name_v[NAMESIZE];
        fit_info.N = 6;
        fit_info.corr_id = { 0, 129, 127, 188, 194, 201,
                               5, 128, 189, 195, 202,
                                   1 , 190, 196, 203,
                                       187, 199, 206,
                                            193, 207,
                                                 200 };//operators{ phi0, phi0^3, phi1, phi0^2phi1, phi0phi1phi1, phi1^3}
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_001_011_111\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 6;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_011_111_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(473 + i);
        }
        check_correlatro_counter(478);

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_011_111_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(479 + i);
        }
        check_correlatro_counter(484);

        int sqrtN = 6;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127, 188, 194, 201,
                               5, 128, 189, 195, 202,
                                   1 , 190, 196, 203,
                                       187, 199, 206,
                                            193, 207,
                                                 200 };//operators{ phi0, phi0^3, phi1, phi0^2phi1, phi0phi1phi1, phi1^3}


        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_001_011_111_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_011_111_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(485 + i);
        }
        check_correlatro_counter(520);

    }
    else { for (int i = 473;i <= 520;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(520);

    printf("GEVP_0_3_1_phi3nl\n");
    if (params.data.ncorr > 212 && strcmp(argv[4], "G2t_T64_L24_msq0-4.868000_msq1-4.710000_l01.000000_l11.000000_mu2.000000_g0.250000_rep0_bin100_merged_bin1000") == 0) {

        double* d0d0 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 209, "phi3nl_meff", M_eff_T, jack_file);
        free(d0d0);
        check_correlatro_counter(521);
        char name_v[NAMESIZE];
        fit_info.N = 4;
        fit_info.corr_id = { 0, 129, 127, 210,
                               5, 128, 211 ,
                                   1 , 212 ,
                                       209 };//diag{ phi0->phi0, phi3->phi3, phi1->phi1,   dphi0->dphi0 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_phi3nl\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 4;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_phi3nl_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(522 + i);
        }
        check_correlatro_counter(525);

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_phi3nl_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(526 + i);
        }
        check_correlatro_counter(529);

        int sqrtN = 4;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127, 188,
                               5, 128, 189 ,
                                   1 , 190 ,
                                       187 };// { phi0, phi0^3, phi1, phi0^2phi1}

        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_phi3nl_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_phi3nl_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(530 + i);
        }
        check_correlatro_counter(545);

    }
    else { for (int i = 521;i <= 545;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(545);
    fit_info.restore_default();

    fit_info.restore_default();
    printf("GEVP_0_3_1_001_011_A1\n");
    if (params.data.ncorr > 207) {

        fit_info.restore_default();
        char name_v[NAMESIZE];
        fit_info.N = 6;
        fit_info.corr_id = { 0, 129, 127, 188, 194, 163,
                               5, 128, 189, 195, 165,
                                   1 , 190, 196, 164,
                                       187, 199, 191,
                                            193, 197,
                                                 116 };//operators{ phi0, phi0^3, phi1, phi0^2phi1, phi0phi1phi1, A1}
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_001_011_A1\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 6;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_011_A1_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(546 + i);
        }
        check_correlatro_counter(551);

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_011_A1_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(552 + i);
        }
        check_correlatro_counter(557);

        int sqrtN = 6;
        int N = sqrtN * sqrtN;
        fit_info.N = N;
        fit_info.corr_id = { 0, 129, 127, 188, 194, 163,
                               5, 128, 189, 195, 165,
                                   1 , 190, 196, 164,
                                       187, 199, 191,
                                            193, 197,
                                                 116 };//operators{ phi0, phi0^3, phi1, phi0^2phi1, phi0phi1phi1, A1}


        fit_info.value_or_vector = 1; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_001_011_A1_v\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        //id_GEVP_031_p1=ncorr_new;
        fit_info.N = 6;

        for (int i = 0; i < N; i++) {
            int comp = i % sqrtN;
            int id = i / sqrtN;
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_011_A1_v%d_%d", id, comp);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(558 + i);
        }
        check_correlatro_counter(593);

    }
    else { for (int i = 546;i <= 593;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(593);

    fit_info.restore_default();
    printf("GEVP_00_11_01\n");
    if (params.data.ncorr > 214) {

        fit_info.N = 1;
        fit_info.corr_id = { 2 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00 = ncorr_new - 1;

        fit_info.corr_id = { 12 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00_11 = ncorr_new - 1;

        fit_info.corr_id = { 213 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00_01 = ncorr_new - 1;

        fit_info.corr_id = { 3 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11 = ncorr_new - 1;

        fit_info.corr_id = { 214 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11_01 = ncorr_new - 1;

        fit_info.corr_id = { 11 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d01 = ncorr_new - 1;

        fit_info.restore_default();
        char name_v[NAMESIZE];
        fit_info.N = 3;
        fit_info.corr_id = { id_d00, id_d00_11,id_d00_01,
                                     id_d11,id_d11_01,
                                              id_d01 };//operators{ shift-phi0, shift-phi1}
        // fit_info.corr_id={2, 12,213,
        //                       3,214,
        //                         11};//operators{ shift-phi0, shift-phi1}


        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        if (strcmp(argv[4], "G2t_T64_L14_msq0-4.180000_msq1-4.180000_l02.500000_l12.500000_mu5.000000_g1.500000_rep0_bin1000_merged_bin100") == 0) {
            fit_info.t0_GEVP = 4;
        }

        //fit_info.corr_id={1,2};
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 3;

        // for (int i=0; i< fit_info.N; i++){
        //     mysprintf(name_v,NAMESIZE,"GEVP_00_11_l%d",i   );
        //     double *l0_GEVP=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile, ncorr_new-fit_info.N+i, name_v, identity,jack_file);
        //     free(l0_GEVP);
        //     check_correlatro_counter(594+i);
        // }
        // check_correlatro_counter(595);

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_00_11_01_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_sinh_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(594 + i);
        }
        check_correlatro_counter(596);



    }
    else { for (int i = 594;i <= 596;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(596);



    printf("GEVP_00_11_01_p1\n");
    if (params.data.ncorr > 217) {


        fit_info.N = 1;
        add_correlators(option, ncorr_new, conf_jack, average_xyz<66, 67, 68>, fit_info);
        int id_d00 = ncorr_new - 1;

        fit_info.N = 1;
        add_correlators(option, ncorr_new, conf_jack, average_xyz<69, 70, 71>, fit_info);
        int id_d11 = ncorr_new - 1;

        fit_info.N = 1;
        add_correlators(option, ncorr_new, conf_jack, average_xyz<72, 73, 74>, fit_info);
        int id_d01 = ncorr_new - 1;


        fit_info.restore_default();
        char name_v[NAMESIZE];
        fit_info.N = 4;
        fit_info.corr_id = { id_d00, 215, 216,218,
                                id_d11,217,219,
                                    id_d01,220,
                                          221 };//operators{ shift-phi0, shift-phi1}
        // fit_info.N=3;
        // fit_info.corr_id={id_d00, 215, 216,
        //                         id_d11,217,
        //                             id_d01   };//operators{ shift-phi0, shift-phi1}

        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;

        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 4;

        // for (int i=0; i< fit_info.N; i++){
        //     mysprintf(name_v,NAMESIZE,"GEVP_00_11_l%d",i   );
        //     double *l0_GEVP=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile, ncorr_new-fit_info.N+i, name_v, identity,jack_file);
        //     free(l0_GEVP);
        //     check_correlatro_counter(594+i);
        // }
        // check_correlatro_counter(595);

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_00_11_01_p1_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(597 + i);
        }
        check_correlatro_counter(600);



    }
    else { for (int i = 597;i <= 600;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(600);


    printf("GEVP_00_11_01_A1\n");
    if (params.data.ncorr > 217) {

        fit_info.N = 1;
        fit_info.corr_id = { 39 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00 = ncorr_new - 1;

        fit_info.corr_id = { 236 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00_11 = ncorr_new - 1;

        fit_info.corr_id = { 237 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00_01 = ncorr_new - 1;

        fit_info.corr_id = { 40 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11 = ncorr_new - 1;

        fit_info.corr_id = { 238 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11_01 = ncorr_new - 1;

        fit_info.corr_id = { 41 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d01 = ncorr_new - 1;


        fit_info.restore_default();
        char name_v[NAMESIZE];
        fit_info.N = 3;
        fit_info.corr_id = { id_d00, id_d00_11,id_d00_01,
                                     id_d11,id_d11_01,
                                              id_d01 };//operators{ shift-phi0, shift-phi1}

        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;

        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 3;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_00_11_01_A1_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(601 + i);
        }
        check_correlatro_counter(603);



    }
    else { for (int i = 601;i <= 603;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(603);


    fit_info.restore_default();
    printf("GEVP_0_3_1_001_011_5\n");
    if (params.data.ncorr > 207 + 1e+6) {

        fit_info.restore_default();
        char name_v[NAMESIZE];
        fit_info.N = 6;
        fit_info.corr_id = { 0, 129, 127, 188, 194, 176,
                               5, 128, 189, 195, 177,
                                   1 , 190, 196, 178,
                                       187, 199, 192,
                                            193, 198,
                                                 175 };//operators{ phi0, phi0^3, phi1, phi0^2phi1, phi0phi1phi1, A1}
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_0_3_1_001_011_5\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 6;


        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_0_3_1_001_011_5_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(604 + i);
        }
        check_correlatro_counter(609);



    }
    else { for (int i = 604;i <= 609;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(609);



    printf("GEVP_00_11_01_00A1_11A1_01A1\n");
    if (params.data.ncorr > 217 + 1e+6) {
        fit_info.N = 1;
        fit_info.corr_id = { 2 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00 = ncorr_new - 1;

        fit_info.corr_id = { 12 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00_11 = ncorr_new - 1;

        fit_info.corr_id = { 213 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00_01 = ncorr_new - 1;

        fit_info.corr_id = { 3 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11 = ncorr_new - 1;

        fit_info.corr_id = { 214 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11_01 = ncorr_new - 1;

        fit_info.corr_id = { 11 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d01 = ncorr_new - 1;

        fit_info.N = 1;
        fit_info.corr_id = { 39 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00A1 = ncorr_new - 1;

        fit_info.corr_id = { 236 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00A1_11A1 = ncorr_new - 1;

        fit_info.corr_id = { 237 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00A1_01A1 = ncorr_new - 1;

        fit_info.corr_id = { 40 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11A1 = ncorr_new - 1;

        fit_info.corr_id = { 238 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11A1_01A1 = ncorr_new - 1;

        fit_info.corr_id = { 41 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d01A1 = ncorr_new - 1;

        fit_info.corr_id = { 63 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00_00A1 = ncorr_new - 1;
        fit_info.corr_id = { 261 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11_00A1 = ncorr_new - 1;
        fit_info.corr_id = { 262 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d01_00A1 = ncorr_new - 1;


        fit_info.corr_id = { 64 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11_11A1 = ncorr_new - 1;

        fit_info.corr_id = { 263 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00_11A1 = ncorr_new - 1;
        fit_info.corr_id = { 264 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d01_11A1 = ncorr_new - 1;


        fit_info.corr_id = { 65 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d01_01A1 = ncorr_new - 1;
        fit_info.corr_id = { 265 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00_01A1 = ncorr_new - 1;
        fit_info.corr_id = { 266 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11_01A1 = ncorr_new - 1;

        fit_info.restore_default();
        char name_v[NAMESIZE];
        fit_info.N = 6;
        fit_info.corr_id = { id_d00, id_d00_11,id_d00_01, id_d00_00A1,id_d00_11A1,id_d00_01A1,
                                     id_d11,id_d11_01, id_d11_00A1,id_d11_11A1,id_d11_01A1,
                                               id_d01, id_d01_00A1,id_d01_11A1,id_d01_01A1,
                                                      id_d00A1, id_d00A1_11A1,id_d00A1_01A1,
                                                                     id_d11A1,id_d11A1_01A1,
                                                                                   id_d01A1 };//operators{ shift-phi0, shift-phi1}

        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;

        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 6;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_00_11_01_00A1_11A1_01A1_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(610 + i);
        }
        check_correlatro_counter(615);



    }
    else { for (int i = 610;i <= 615;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(615);



    fit_info.restore_default();
    printf("GEVP_00_11\n");
    if (params.data.ncorr > 214) {

        fit_info.N = 1;
        fit_info.corr_id = { 2 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00 = ncorr_new - 1;

        fit_info.corr_id = { 12 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d00_11 = ncorr_new - 1;

        fit_info.corr_id = { 3 };
        add_correlators(option, ncorr_new, conf_jack, der1_corr, fit_info);
        int id_d11 = ncorr_new - 1;


        fit_info.restore_default();
        char name_v[NAMESIZE];
        fit_info.N = 2;
        fit_info.corr_id = { id_d00, id_d00_11,
                                     id_d11 };//operators{ shift-phi0, shift-phi1}

        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;

        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);
        fit_info.restore_default();
        fit_info.N = 2;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_00_11_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_sinh_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(616 + i);
        }
        check_correlatro_counter(617);


    }
    else { for (int i = 616;i <= 617;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(617);



    if (params.data.ncorr >= 273) {

        char name_v[NAMESIZE];
        fit_info.N = 3;
        fit_info.corr_id = { 5, 128,269,
                                1 ,271,
                                267 };//diag{ 3phi0->3phi0, phi1->phi1 }
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        //fit_info.corr_id={1,2};
        printf("GEVP_3_1_00d1\n");
        add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        printf(" ncorr after GEVP %d\n", ncorr_new);

        fit_info.restore_default();
        fit_info.N = 3;

        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_3_1_00d1_l%d", i);
            double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, identity, jack_file);
            free(l0_GEVP);
            check_correlatro_counter(618 + i);
        }
        for (int i = 0; i < fit_info.N; i++) {
            mysprintf(name_v, NAMESIZE, "GEVP_3_1_00d1_meffl%d", i);
            double* ml0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - fit_info.N + i, name_v, M_eff_T, jack_file);
            free(ml0_GEVP);
            check_correlatro_counter(621 + i);
        }

        double* tmp_imag = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 270, "imag_3_00d1", identity, jack_file);
        free(tmp_imag);

        tmp_imag = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 272, "imag_1_00d1", identity, jack_file);
        free(tmp_imag);


        // int sqrtN = 2;
        // int N = sqrtN * sqrtN;
        // fit_info.N = N;
        // fit_info.corr_id = { 5, 128,
        //                         1 };//diag{ 3phi0->3phi0, phi1->phi1 }
        // fit_info.value_or_vector = 1; // 0= values
        // fit_info.t0_GEVP = 3;
        // //fit_info.corr_id={1,2};
        // printf("GEVP_3_1_v\n");
        // add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        // printf(" ncorr after GEVP %d\n", ncorr_new);

        // fit_info.restore_default();
        // //id_GEVP_031_p1=ncorr_new;

        // for (int i = 0; i < N; i++) {
        //     int comp = i % sqrtN;
        //     int id = i / sqrtN;
        //     mysprintf(name_v, NAMESIZE, "GEVP_3_1_v%d_%d", id, comp);
        //     double* l0_GEVP = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, ncorr_new - N + i, name_v, identity, jack_file);
        //     free(l0_GEVP);
        //     check_correlatro_counter(356 + i);
        // }
        // check_correlatro_counter(359);
    }
    else { for (int i = 618;i <= 625;i++)  zero_corr(zeros, Njack, jack_file); }
    check_correlatro_counter(625);

    free(zeros);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
    free(E2_01);
    free(E2_01_ws);
    free(a_01_ws);
    free(a);
    free(E3_0);
    free(a_0);
    free(a0m0);
    // free(E1_0_p11);
    // free(E1_0_p111);
    // free(E2_0_p1);
    // free(E2_0_p11);
    // free(E2_0_p111);
    // free(E2_0_A1);
    //////////////////////////////////////
    free_2(3, mass);free_2(3, E2);
    free_jack(Njack, ncorr_new, file_head.l0, conf_jack);

    fclose(out_gamma);
    // free(fit_info.ext_P);
    free_tower(7, (void**)option);

    for (i = 0;i < file_head.nmoms;i++)
        free(file_head.mom[i]);
    free(file_head.mom);
    fclose(infile); fclose(outfile);

    free(file_head.k);
    fclose(jack_file);
    printf("time computing observables: %g s\n", timestamp() - t1);
    return 0;
}

