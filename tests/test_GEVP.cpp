#define CONTROL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include <time.h>
// #include <complex.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <random>

#include "global.hpp"
#include "read.hpp"
#include "resampling_new.hpp"

#include "linear_fit.hpp"
#include "mutils.hpp"
#include "various_fits.hpp"

#include "correlators_analysis.hpp"
#include "fit_all.hpp"
#include "tower.hpp"




double lhs_function_me_GEVP(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    double amp = 0;
    double M = fit_info.ext_P[0][j];


    double me = in[j][id][t][0] * std::sqrt(2 * M) / std::sqrt((std::exp(-t * M) + std::exp(-(fit_info.T - t) * M)));
    return me;
}

double lhs_function_me_im_GEVP(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    double amp = 0;
    double M = fit_info.ext_P[0][j];


    double me = in[j][id][t][1] * std::sqrt(2 * M) / std::sqrt((std::exp(-t * M) + std::exp(-(fit_info.T - t) * M)));
    return me;
}


void compare_and_error(double* correct, double* value, std::vector<std::string>& errors, int Njack, std::string name) {
    int Nerr = 0;
    for (int j = 0; j < Njack; j++) {
        if (fabs(correct[j] - value[j]) > 1e-4) {
            printf(" %s [%d]= %g\n", name.c_str(), j, correct[j] - value[j]);
            Nerr++;
            errors.push_back("error in " + name +
                " " + std::to_string(fabs(correct[0] - value[0])));
            break;
        }
        else if (std::isnan(value[j])) {
            printf(" %s [%d]= %g\n", name.c_str(), j, correct[j] - value[j]);
            Nerr++;
            errors.push_back("error in " + name +
                " " + std::to_string(fabs(correct[0] - value[0])));
            break;
        }

    }

}

int main(int argc, char** argv) {
    struct kinematic kinematic_2pt;
    FILE* outfile = open_file("GEVP_out.txt", "w+");
    FILE* jack_file = open_file("jack_gevp.dat", "w+");
    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, "read_plateaux"); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p");            // -p
    mysprintf(option[3], NAMESIZE, "./");         // path
    mysprintf(option[4], NAMESIZE, "jack");         // resampling
    mysprintf(option[5], NAMESIZE, "no");            // pdf
    mysprintf(option[6], NAMESIZE, "no_infile");         // infile


    int T = 16;
    file_head.l0 = T;
    int Neff = 100;
    myres = new resampling_jack(Neff);
    int Njack = Neff + 1;

    int ncorr_new = 4;
    int Nmax = 28;
    double**** conf_jack = calloc_corr(Njack, Nmax, T);

    double* m0 = myres->create_fake(0.1, 0.002, 533);
    double* m1 = myres->create_fake(0.2, 0.002, 533);

    double* G0 = myres->create_fake(0.5, 0.002, 533);
    double* G1 = myres->create_fake(0.6, 0.002, 533);

    for (int j = 0; j < Njack; j++) {
        for (int t = 0;t < T;t++) {
            conf_jack[j][0][t][0] = (G0[j] * G0[j] / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)));
            conf_jack[j][3][t][0] = (G1[j] * G1[j] / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)));
        }
    }
    // 
    std::vector<std::string>  errors;
    ///////////////////////////////////////////////////// GEVP 
    fit_type fit_info;

    // mass
    fit_info.value_or_vector = 0;
    fit_info.N = 2;
    fit_info.corr_id = { 0,1,2,3 };

    fit_info.t0_GEVP = 3;
    fit_info.GEVP_ignore_warning_after_t = 1;
    fit_info.verbosity = 0;
    fit_info.Njack = Njack;
    add_correlators_no_alloc(option, ncorr_new, conf_jack, GEVP_matrix, fit_info, Nmax);
    int id_GEVP_m0 = ncorr_new - 2;
    int id_GEVP_m1 = ncorr_new - 1;

    fit_info.codeplateaux = true;
    fit_info.tmin = 2;
    fit_info.tmax = 2;
    printf("m0 = %g %g\n", m0[Njack - 1], myres->comp_error(m0));
    double* m0_GEVP = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        "plateaux", outfile, id_GEVP_m0, "m0", M_eff_T, jack_file, fit_info);

    compare_and_error(m0, m0_GEVP, errors, Njack, "m0-m0_GEVP");

    printf("m1 = %g %g\n", m1[Njack - 1], myres->comp_error(m1));
    double* m1_GEVP = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        "plateaux", outfile, id_GEVP_m1, "m1", M_eff_T, jack_file, fit_info);

    compare_and_error(m1, m1_GEVP, errors, Njack, "m1-m1_GEVP");

    /// fit the ampitude
    fit_info.value_or_vector = 2; // 2= Amplitudes
    fit_info.N = 2 * 2;
    fit_info.corr_id = { 0,1,2,3 };

    fit_info.t0_GEVP = 3;
    fit_info.GEVP_ignore_warning_after_t = 1;
    fit_info.verbosity = 0;
    add_correlators_no_alloc(option, ncorr_new, conf_jack, GEVP_matrix, fit_info, Nmax);
    int id_GEVP_a0_0 = ncorr_new - 4;
    int id_GEVP_a0_1 = ncorr_new - 3;
    int id_GEVP_a1_0 = ncorr_new - 2;
    int id_GEVP_a1_1 = ncorr_new - 1;

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 1;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);


    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = T;

    fit_info.codeplateaux = true;
    fit_info.tmin = 2;
    fit_info.tmax = 2;

    fit_info.ext_P[0] = m0_GEVP;
    fit_info.corr_id = { id_GEVP_a0_0 };
    struct fit_result G0_GEVP = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_me_GEVP, "G0", fit_info, jack_file);

    compare_and_error(G0, G0_GEVP.P[0], errors, Njack, "G0-G0_GEVP");

    fit_info.ext_P[0] = m1_GEVP;
    fit_info.corr_id = { id_GEVP_a1_1 };
    struct fit_result G1_GEVP = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_me_GEVP, "G1", fit_info, jack_file);

    compare_and_error(G1, G1_GEVP.P[0], errors, Njack, "G1-G1_GEVP");

    // FIXME: we can not fit id_GEVP_a0_1 because it had zero error apart from t0
    fit_info.linear_fit = false;
    fit_info.guess = { 0.0 };
    fit_info.tmin = 2;
    fit_info.tmax = 2;
    fit_info.ext_P[0] = m1_GEVP;
    fit_info.corr_id = { id_GEVP_a0_1 };
    struct fit_result G0_1_GEVP = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_me_GEVP, "G0_1", fit_info, jack_file);

    double* zero = myres->create_fake(0.0, 1e-12, 533);
    compare_and_error(zero, G0_1_GEVP.P[0], errors, Njack, "0-G0_1_GEVP");

    fit_info.ext_P[0] = m0_GEVP;
    fit_info.corr_id = { id_GEVP_a1_0 };
    struct fit_result G1_1_GEVP = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_me_GEVP, "G0_1", fit_info, jack_file);

    compare_and_error(zero, G1_1_GEVP.P[0], errors, Njack, "0-G1_1_GEVP");


    //////////////////////////////////////////////////////////////
    // mixing
    //////////////////////////////////////////////////////////////
    double* G01 = myres->create_fake(0.4, 0.002, 533);
    double* G10 = myres->create_fake(0.3, 0.002, 533);

    for (int j = 0; j < Njack; j++) {
        for (int t = 0;t < T;t++) {
            conf_jack[j][0][t][0] = (G0[j] * G0[j] / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)))
                + (G01[j] * G01[j] / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)));
            conf_jack[j][1][t][0] = (G0[j] * G10[j] / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)))
                + (G01[j] * G1[j] / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)));
            conf_jack[j][2][t][0] = conf_jack[j][1][t][0];

            conf_jack[j][3][t][0] = (G1[j] * G1[j] / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)))
                + (G10[j] * G10[j] / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)));
        }
    }


    // mass
    fit_info.value_or_vector = 0;
    fit_info.N = 2;
    fit_info.corr_id = { 0,1,2,3 };

    fit_info.t0_GEVP = 3;
    fit_info.GEVP_ignore_warning_after_t = 1;
    fit_info.verbosity = 0;
    fit_info.Njack = Njack;
    add_correlators_no_alloc(option, ncorr_new, conf_jack, GEVP_matrix, fit_info, Nmax);
    id_GEVP_m0 = ncorr_new - 2;
    id_GEVP_m1 = ncorr_new - 1;

    fit_info.codeplateaux = true;
    fit_info.tmin = 2;
    fit_info.tmax = 2;
    printf("m0 = %g %g\n", m0[Njack - 1], myres->comp_error(m0));
    double* m0_GEVP1 = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        "plateaux", outfile, id_GEVP_m0, "m0-mixing", M_eff_T, jack_file, fit_info);

    compare_and_error(m0, m0_GEVP1, errors, Njack, "m0-m0_GEVP-mixing");

    printf("m1 = %g %g\n", m1[Njack - 1], myres->comp_error(m1));
    double* m1_GEVP1 = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        "plateaux", outfile, id_GEVP_m1, "m1-mixing", M_eff_T, jack_file, fit_info);

    compare_and_error(m1, m1_GEVP1, errors, Njack, "m1-m1_GEVP-mixing");


    /// fit the ampitude
    fit_info.value_or_vector = 2; // 2= Amplitudes
    fit_info.N = 2 * 2;
    fit_info.corr_id = { 0,1,2,3 };

    fit_info.t0_GEVP = 3;
    fit_info.GEVP_ignore_warning_after_t = 1;
    fit_info.verbosity = 0;
    add_correlators_no_alloc(option, ncorr_new, conf_jack, GEVP_matrix, fit_info, Nmax);
    id_GEVP_a0_0 = ncorr_new - 4;
    id_GEVP_a0_1 = ncorr_new - 3;
    id_GEVP_a1_0 = ncorr_new - 2;
    id_GEVP_a1_1 = ncorr_new - 1;

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 1;
    // fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);


    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = T;

    fit_info.codeplateaux = true;
    fit_info.tmin = 2;
    fit_info.tmax = 2;

    fit_info.ext_P[0] = m0_GEVP;
    fit_info.corr_id = { id_GEVP_a0_0 };
    struct fit_result G0_GEVP1 = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_me_GEVP, "G0-mixing", fit_info, jack_file);

    compare_and_error(G0, G0_GEVP1.P[0], errors, Njack, "G0-G0_GEVP-mixing");

    fit_info.ext_P[0] = m1_GEVP;
    fit_info.corr_id = { id_GEVP_a1_1 };
    struct fit_result G1_GEVP1 = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_me_GEVP, "G1-mixing", fit_info, jack_file);

    compare_and_error(G1, G1_GEVP1.P[0], errors, Njack, "G1-G1_GEVP-mixing");

    // FIXME: we can not fit id_GEVP_a0_1 because it had zero error apart from t0
    fit_info.linear_fit = false;
    fit_info.guess = { 0.0 };
    fit_info.tmin = 2;
    fit_info.tmax = 2;
    fit_info.ext_P[0] = m1_GEVP;
    fit_info.corr_id = { id_GEVP_a0_1 };
    struct fit_result G0_1_GEVP1 = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_me_GEVP, "G0_1-mixing", fit_info, jack_file);

    compare_and_error(G01, G0_1_GEVP1.P[0], errors, Njack, "0-G0_1_GEVP-mixing");

    fit_info.linear_fit = false;
    fit_info.guess = { 0.0 };
    fit_info.tmin = 2;
    fit_info.tmax = 2;
    fit_info.ext_P[0] = m0_GEVP;
    fit_info.corr_id = { id_GEVP_a1_0 };
    struct fit_result G1_0_GEVP1 = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_me_GEVP, "G0_1-mixing", fit_info, jack_file);

    compare_and_error(G10, G1_0_GEVP1.P[0], errors, Njack, "0-G0_1_GEVP-mixing");


    // fit_info.ext_P[0] = m0_GEVP;
    // fit_info.corr_id = { id_GEVP_a1_0 };
    // struct fit_result G1_1_GEVP = fit_fun_to_fun_of_corr(
    //     option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
    //     outfile, lhs_function_me_GEVP, "G0_1", fit_info, jack_file);

    // compare_and_error(zero, G1_1_GEVP.P[0], errors, Njack, "0-G1_1_GEVP");
    {
        int j = Njack - 1;
        int t = fit_info.t0_GEVP;
        printf("%g\n", (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t))) / std::sqrt(2 * m0[j]));
        printf("%g\n", (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t))) / std::sqrt(2 * m1[j]));

    }

    //////////////////////////////////////////////////////////////
    // complex
    //////////////////////////////////////////////////////////////
    printf("\n complex\n\n");

    free(G0); free(G1);
    free(m0); free(m1);
    free(G01); free(G10);
    double* G0_i, * G1_i, * G01_i, * G10_i;

    // m0 = myres->create_fake(0.13, 0.002, 533);
    // m1 = myres->create_fake(0.2, 0.002, 533);

    // G0 = myres->create_fake(0.5, 0.002, 533);
    // G0_i = myres->create_fake(0.3, 0.0002, 533);

    // G1 = myres->create_fake(0.6, 0.002, 533);
    // G1_i = myres->create_fake(0.5, 0.0002, 533);

    // G01 = myres->create_fake(0.4, 0.002, 533);
    // G01_i = myres->create_fake(0.6, 0.0002, 533);

    // G10 = myres->create_fake(0.3, 0.002, 533);
    // G10_i = myres->create_fake(0.7, 0.0002, 533);

    // using random numbers
    // unsigned int seed = 42;
    // std::mt19937 gen(seed); // Mersenne Twister engine with fixed seed
    std::random_device rd;
    std::mt19937 gen(rd()); // Mersenne Twister RNG

    // Uniform distribution between 0.0 and 1.0
    std::uniform_real_distribution<> distrib(0.001, 0.8);
    double rm0 = distrib(gen);
    double rm1 = distrib(gen);
    printf("masses: %g  %g\n", rm0, rm1);
    m0 = myres->create_fake(std::min(rm0,rm1), 0.002, 533);
    m1 = myres->create_fake(std::max(rm0,rm1), 0.002, 533);
    double rG0 = distrib(gen); 
    double rG0i = distrib(gen); 
    printf("G0: %g  %g\n", rG0, rG0i);
    G0 = myres->create_fake(rG0, 0.002, 533);
    G0_i = myres->create_fake(rG0i, 0.0002, 533);

    double rG1 = distrib(gen); 
    double rG1i = distrib(gen); 
    printf("G1: %g  %g\n", rG0, rG0i);
    G1 = myres->create_fake(rG1, 0.002, 533);
    G1_i = myres->create_fake(rG1i, 0.0002, 533);

    double rG01 = distrib(gen); 
    double rG01i = distrib(gen); 
    printf("G01: %g  %g\n", rG01, rG01i);
    G01 = myres->create_fake(rG01, 0.002, 533);
    G01_i = myres->create_fake(rG01i, 0.0002, 533);

    double rG10 = distrib(gen); 
    double rG10i = distrib(gen); 
    printf("G10: %g  %g\n", rG10, rG10i);
    G10 = myres->create_fake(rG10, 0.002, 533);
    G10_i = myres->create_fake(rG10i, 0.0002, 533);


    double* G0abs = (double*)malloc(sizeof(double) * Njack);
    double* G1abs = (double*)malloc(sizeof(double) * Njack);
    double* G01abs = (double*)malloc(sizeof(double) * Njack);
    double* G10abs = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0; j < Njack; j++) {
        std::complex<double> G0_c = std::complex<double>(G0[j], G0_i[j]);
        std::complex<double> G1_c = std::complex<double>(G1[j], G1_i[j]);
        std::complex<double> G01_c = std::complex<double>(G01[j], G01_i[j]);
        std::complex<double> G10_c = std::complex<double>(G10[j], G10_i[j]);
        G0abs[j] = std::abs(G0_c);
        G1abs[j] = std::abs(G1_c);
        G01abs[j] = std::abs(G01_c);
        G10abs[j] = std::abs(G10_c);
        for (int t = 0;t < T;t++) {

            conf_jack[j][0][t][0] = ((G0_c * std::conj(G0_c)).real() / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)))
                + ((G01_c * std::conj(G01_c)).real() / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)));
            conf_jack[j][0][t][1] = 0;
            // ((G0_c * std::conj(G0_c)).imag() / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)))
            //     + ((G01_c * std::conj(G01_c)).imag() / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)));



            conf_jack[j][1][t][0] = ((G0_c * std::conj(G10_c)).real() / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)))
                + ((G01_c * std::conj(G1_c)).real() / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)));
            conf_jack[j][1][t][1] = ((G0_c * std::conj(G10_c)).imag() / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)))
                + ((G01_c * std::conj(G1_c)).imag() / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)));


            // conf_jack[j][1][t][1] = -conf_jack[j][1][t][1];

            conf_jack[j][2][t][0] = conf_jack[j][1][t][0];
            conf_jack[j][2][t][1] -= conf_jack[j][1][t][1];



            conf_jack[j][3][t][0] = ((G1_c * std::conj(G1_c)).real() / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)))
                + ((G10_c * std::conj(G10_c)).real() / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)));
            conf_jack[j][3][t][1] = 0;
            // ((G1_c * std::conj(G1_c)).imag() / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)))
            //     + ((G10_c * std::conj(G10_c)).imag() / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)));


            // separate
            // conf_jack[j][0][t][0] = ((G0_c * std::conj(G0_c)).real() / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)));
            // conf_jack[j][0][t][1] = ((G0_c * std::conj(G0_c)).imag() / (2 * m0[j])) * (std::exp(-m0[j] * t) + std::exp(-m0[j] * (T - t)));




            // conf_jack[j][1][t][0] =0;
            // conf_jack[j][1][t][1] = 0;

            // conf_jack[j][2][t][0] = conf_jack[j][1][t][0];
            // conf_jack[j][2][t][1] -= conf_jack[j][1][t][1];



            // conf_jack[j][3][t][0] = ((G1_c * std::conj(G1_c)).real() / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)));
            // conf_jack[j][3][t][1] = ((G1_c * std::conj(G1_c)).imag() / (2 * m1[j])) * (std::exp(-m1[j] * t) + std::exp(-m1[j] * (T - t)));

        }
    }

    // mass
    fit_info.value_or_vector = 0;
    fit_info.N = 2;
    fit_info.corr_id = { 0,1,2,3 };

    fit_info.t0_GEVP = 3;
    fit_info.GEVP_ignore_warning_after_t = 1;
    fit_info.verbosity = 0;
    fit_info.Njack = Njack;
    add_correlators_no_alloc(option, ncorr_new, conf_jack, GEVP_matrix, fit_info, Nmax);
    id_GEVP_m0 = ncorr_new - 2;
    id_GEVP_m1 = ncorr_new - 1;

    fit_info.codeplateaux = true;
    fit_info.tmin = 2;
    fit_info.tmax = 2;
    double* m0_GEVPc = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        "plateaux", outfile, id_GEVP_m0, "m0-complex", M_eff_T, jack_file, fit_info);

    compare_and_error(m0, m0_GEVPc, errors, Njack, "m0-m0_GEVP-complex");

    double* m1_GEVPc = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        "plateaux", outfile, id_GEVP_m1, "m1-complex", M_eff_T, jack_file, fit_info);

    compare_and_error(m1, m1_GEVPc, errors, Njack, "m1-m1_GEVP-complex");


    /// fit the ampitude
    fit_info.value_or_vector = 2; // 2= Amplitudes
    fit_info.N = 2 * 2;
    fit_info.corr_id = { 0,1,2,3 };

    fit_info.t0_GEVP = 3;
    fit_info.GEVP_ignore_warning_after_t = 1;
    fit_info.verbosity = 0;
    add_correlators_no_alloc(option, ncorr_new, conf_jack, GEVP_matrix, fit_info, Nmax);
    id_GEVP_a0_0 = ncorr_new - 4;
    id_GEVP_a0_1 = ncorr_new - 3;
    id_GEVP_a1_0 = ncorr_new - 2;
    id_GEVP_a1_1 = ncorr_new - 1;

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 1;
    // fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);


    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = T;

    fit_info.codeplateaux = true;
    fit_info.tmin = 2;
    fit_info.tmax = 2;

    fit_info.ext_P[0] = m0_GEVPc;
    fit_info.corr_id = { id_GEVP_a0_0 };
    struct fit_result G0_GEVPr = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_abs_GEVP, "G0-abs", fit_info, jack_file);

    compare_and_error(G0abs, G0_GEVPr.P[0], errors, Njack, "G0-G0_GEVP-abs");

    // struct fit_result G0_GEVPi = fit_fun_to_fun_of_corr(
    //     option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
    //     outfile, lhs_function_me_im_GEVP, "G0-im", fit_info, jack_file);

    // compare_and_error(zero, G0_GEVPi.P[0], errors, Njack, "G0-G0_GEVP-imag");

    fit_info.ext_P[0] = m1_GEVPc;
    fit_info.corr_id = { id_GEVP_a1_1 };
    struct fit_result G1_GEVPr = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_abs_GEVP, "G1-abs", fit_info, jack_file);

    compare_and_error(G1abs, G1_GEVPr.P[0], errors, Njack, "G1_GEVP-abs");

    // fit_info.ext_P[0] = m1_GEVP;
    // fit_info.corr_id = { id_GEVP_a1_1 };
    // struct fit_result G1_GEVPim = fit_fun_to_fun_of_corr(
    //     option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
    //     outfile, lhs_function_me_im_GEVP, "G1-im", fit_info, jack_file);

    // compare_and_error(G1abs, G1_GEVPim.P[0], errors, Njack, "G1_GEVP-im");

    fit_info.tmin = 2;
    fit_info.tmax = 2;
    fit_info.ext_P[0] = m1_GEVPc;
    fit_info.corr_id = { id_GEVP_a0_1 };
    struct fit_result G0_1_GEVP1r = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_abs_GEVP, "G01-abs", fit_info, jack_file);

    compare_and_error(G01abs, G0_1_GEVP1r.P[0], errors, Njack, "G01_GEVP-abs");


    // fit_info.tmin = 2;
    // fit_info.tmax = 2;
    // fit_info.ext_P[0] = m1_GEVP;
    // fit_info.corr_id = { id_GEVP_a0_1 };
    // struct fit_result G0_1_GEVP1im = fit_fun_to_fun_of_corr(
    //     option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
    //     outfile, lhs_function_me_im_GEVP, "G01-im", fit_info, jack_file);

    // compare_and_error(G01abs, G0_1_GEVP1im.P[0], errors, Njack, "G01_GEVP-im");


    fit_info.tmin = 2;
    fit_info.tmax = 2;
    fit_info.ext_P[0] = m0_GEVPc;
    fit_info.corr_id = { id_GEVP_a1_0 };
    struct fit_result G1_0_GEVP1r = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
        outfile, lhs_function_abs_GEVP, "G10-abs", fit_info, jack_file);

    compare_and_error(G10abs, G1_0_GEVP1r.P[0], errors, Njack, "G10_GEVP-abs");

    // fit_info.tmin = 2;
    // fit_info.tmax = 2;
    // fit_info.ext_P[0] = m0_GEVP;
    // fit_info.corr_id = { id_GEVP_a1_0 };
    // struct fit_result G1_0_GEVP1im = fit_fun_to_fun_of_corr(
    //     option, kinematic_2pt, (char*)"P5P5", conf_jack, "plateaux",
    //     outfile, lhs_function_me_im_GEVP, "G10-im", fit_info, jack_file);

    // compare_and_error(G10abs, G1_0_GEVP1im.P[0], errors, Njack, "G10_GEVP-im");


    // clean up
    fit_info.restore_default();
    free(G0); free(G1);
    free(G0_i); free(G1_i);
    free(G01); free(G10);
    free(G01_i); free(G10_i);
    free(G0abs); free(G1abs);
    free(G01abs); free(G10abs);
    free(m1_GEVPc);
    free(m0_GEVPc);
    free(m0_GEVP); free(m1_GEVP);
    free(m0_GEVP1); free(m1_GEVP1);
    free(zero);
    G0_GEVP.clear();
    G1_GEVP.clear();
    G0_GEVP1.clear();
    G1_GEVP1.clear();
    G0_GEVPr.clear();
    G1_GEVPr.clear();
    G0_1_GEVP.clear();
    G1_1_GEVP.clear();
    G0_1_GEVP1r.clear();
    G1_0_GEVP1r.clear();

    G1_0_GEVP1.clear();
    G0_1_GEVP1.clear();
    
    free(m0); free(m1);
    free_corr( Njack,Nmax, T, conf_jack);
    for (int i = 0; i < 7; i++) {
     free(option[i]);
    }
     free(option);

    //////////////////////////////////////////////////////////////
    // collect errors
    //////////////////////////////////////////////////////////////
    for (auto e : errors) {
        printf("%s\n", e.c_str());
    }
    if (errors.size() > 0) {
        printf("There are errors in the GEVP test\n");
    }
    else {
        printf("All tests passed\n");
    }
}