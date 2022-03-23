#include <gsl/gsl_integration.h>

constexpr double muon_mass_MeV = 105.6583755;
constexpr double muon_mass_fm = muon_mass_MeV * 197.326963;
constexpr double alpha_em = 0.0072973525693;

double integrand_K(double x, void* params) {
    double z = *(double*)params;
    double arg = (z / 2) * x / (sqrt(1 - x));
    double j0 = sin(arg) / arg;
    double f = (1 - x) * (1 - j0 * j0);
    return f;
}

double kernel_K(double z, double epsrel = 1e-7) {
    int Maxiter = 1e+3;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(Maxiter);
    double result, error;
    gsl_function F;
    F.function = &integrand_K;
    F.params = &z;

    gsl_integration_qags(&F, 0, 1, 0, epsrel, Maxiter, w, &result, &error);
    // printf("K(%g)=%g\n", z, result);
    gsl_integration_workspace_free(w);

    return result;
}

double gm2_step_function(double t_d, double t1_d) {
    return 1.0 / (1.0 + exp(-2 * (t_d - t1_d)));
}

// double amu_sd(int j, double**** in, int t, struct fit_type fit_info) {
double compute_amu_sd(double* VV, double a, double q2) {
    constexpr double d = 0.15;
    constexpr double t1_d = 0.4 / d;
    int T = file_head.l0;

    double amu = 0;
    for (int t_a = 1; t_a < T / 2; t_a++) {
        double t = t_a * a; // time in fm.
        double z = muon_mass_MeV * (t / 197.326963);
        // printf("time = %d, z=%g\n", t_a, z);
        double K = z * z * kernel_K(z);
        double theta = gm2_step_function(t / 0.15, t1_d);
        double VV_sub=VV[t_a]-q2*(1.0/(2.0*M_PI*M_PI*pow(t_a,5)));
        amu += K * VV_sub * (1 - theta);
        // printf("t=%d  K=%g   VV=%g  1-theta=%g  amu=%g\n",t_a,K, VV[t_a],1-theta, amu);
    }
    amu *= 4 * alpha_em * alpha_em *
        q2 / (muon_mass_MeV * muon_mass_MeV * (a / 197.326963) * (a / 197.326963));
    return amu;
}

double compute_amu_sd_simpson38(double* VV, double a, double q2) {
    constexpr double d = 0.15;
    constexpr double t1_d = 0.4 / d;
    int T = file_head.l0;

    double amu = 0;
    double* fi = (double*)malloc(sizeof(double) * T / 2);
    fi[0] = 0;
    for (int t_a = 1; t_a < T / 2; t_a++) {
        double t = t_a * a; // time in fm.
        double z = muon_mass_MeV * (t / 197.326963);
        // printf("time = %d, z=%g\n", t_a, z);
        double K = z * z * kernel_K(z);
        double theta = gm2_step_function(t / 0.15, t1_d);
        fi[t_a] = K * VV[t_a] * (1 - theta);
        // printf("t=%d  K=%g   VV=%g  1-theta=%g  amu=%g\n",t_a,K, VV[t_a],1-theta, amu);
    }
    // amu =

    amu *= 4 * alpha_em * alpha_em *
        q2 / (muon_mass_MeV * muon_mass_MeV * (a / 197.326963) * (a / 197.326963));
    return amu;
}

/********************************************************************************
 * Z
 ********************************************************************************/

double ZAl_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double mu = file_head.musea;
    double M_PS = fit_info.ext_P[0][j];
    double M_PS_OS = fit_info.ext_P[1][j];
    double GPS = fit_info.ext_P[2][j];
    double GPS_OS = fit_info.ext_P[3][j];


    double num = in[j][4][t][0];
    double den = in[j][3][t + 1][0] - in[j][3][t][0];
    double RA = (2 * mu * num / (den));

    return RA * M_PS_OS * sinh(M_PS_OS) * GPS / (M_PS * sinh(M_PS) * GPS_OS);

}

double ZVl_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double mu = file_head.musea;

    double num = in[j][1][t][0];
    double den = in[j][0][t + 1][0] - in[j][0][t][0];
    return (2 * mu * num / (den));

}

double GPS_OS_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double mu = file_head.musea;
    double mass = fit_info.ext_P[0][j];

    double GPS_OS = sqrt(in[j][4][t][0] * 2 * mass / (exp(-mass * t) + exp(-mass * (T - t))));

    return GPS_OS;
}
double GPS_lhs(int j, double**** in, int t, struct fit_type fit_info) {
    int T = file_head.l0;
    double mu = file_head.musea;
    double mass = fit_info.ext_P[0][j];

    double GPS_OS = sqrt(in[j][1][t][0] * 2 * mass / (exp(-mass * t) + exp(-mass * (T - t))));

    return GPS_OS;
}
