#define CONTROL
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#include "global.hpp"
#include "resampling_new.hpp"
#include "non_linear_fit.hpp"
#include "read.hpp"

double f(int n, int Nvar, double* x, int Npar, double* P) {
    // 0.5 * exp(0.1*x)
    return  P[0] * exp(x[0] * P[1]);
}

void test_res(non_linear_fit_result fit, int NE, int Npar) {
    printf("chi2=%-18g     chi2/dof=%g\n", fit.chi2, fit.chi2 / (NE - Npar));
    printf("min=%g   %g\n", fit.P[0], fit.P[1]);
    double chi = sqrt(fit.chi2 / (NE - Npar));
    double gchidof = 1.32723;
    double gP0 = 1.75309e-05, dgP1 = 6.99e-06 * gchidof;
    double gP1 = 0.0414701, dgP0 = 0.001104 * gchidof;
    if (fabs(fit.P[0] - gP0) > 1e-4) { printf("the minimum should be ( %g, %g )\n", gP0, gP1); exit(1); }
    if (fabs(fit.P[1] - gP1) > 1e-4) { printf("the minimum should be ( %g, %g )\n", gP0, gP1); exit(1); }
    if (fabs(chi - gchidof) > 1e-4) { printf("chi2/dof should be  %g \n", gchidof * gchidof); exit(1); }
    printf("the fit is correct!\n\n");
}


int main() {
    int Neff = 100000;
    int Njack = Neff + 1;

    double mean = 0.1;
    double err = 0.002;
    resampling_f* myres_b = new resampling_boot(Neff);
    myres = new resampling_jack(Neff);

    double* x_g = myres_b->create_fake(mean, err * sqrt(Neff), 533);
    double* x_b = myres_b->create_fake(mean, err, 533);
    double* x_j = myres->create_fake(mean, err, 533);
    double* tmp = (double*)calloc(Njack, sizeof(double));

    double* tmp_g = myres->create_fake(mean, err, 533);
    double* tmp_b = myres->create_fake(mean, err, 533);
    double* tmp_j = myres->create_fake(mean, err, 533);


    printf("################### x ############################\n");
    printf("mean and error computed with gauss:  %-20.12g%-20.12g\n", x_g[Njack - 1], myres_b->comp_error(x_g) / sqrt(Neff));
    printf("mean and error computed with boot :  %-20.12g%-20.12g\n", x_b[Njack - 1], myres_b->comp_error(x_b));
    printf("mean and error computed with jack :  %-20.12g%-20.12g\n", x_j[Njack - 1], myres->comp_error(x_j));
    myres->add_error_quadrature(x_j,x_j,0.002);
    // myres->mult_error(x_j,x_j,2);
    printf("mean and error computed with jack x2 :  %-20.12g%-20.12g\n", x_j[Njack - 1], myres->comp_error(x_j));
    printf("mean and error computed original  :  %-20.12g%-20.12g\n", mean, err);

    for (int j = 0;j < Njack;j++) {
        tmp_g[j] = 1.0 / x_g[j];
        tmp_b[j] = 1.0 / x_b[j];
        tmp_j[j] = 1.0 / x_j[j];
    }
    // we compute the average on the samples because the last one will always be exact average
    double a_g = 0;
    double a_b = 0;
    double a_j = 0;
    for (int j = 0;j < Neff;j++) {
        a_g += tmp_g[j];
        a_b += tmp_b[j];
        a_j += tmp_j[j];
    }

    a_g /= (double)Neff;
    a_b /= (double)Neff;
    a_j /= (double)Neff;


    printf("###################   1/x ############################\n");
    printf("mean and error computed with gauss:  %-20.12g%-20.12g\n", a_g, myres_b->comp_error(tmp_g) / sqrt(Neff));
    printf("mean and error computed with boot :  %-20.12g%-20.12g\n", a_b, myres_b->comp_error(tmp_b));
    printf("mean and error computed with jack :  %-20.12g%-20.12g\n", a_j, myres->comp_error(tmp_j));
    // propagating the error by hand
    double mean1 = 1.0 / mean;
    double err1 = (1.0 / (mean * mean)) * err;
    printf("mean and error computed by hand   :  %-20.12g%-20.12g\n", mean1, err1);


    printf("###############################################\n");

    double sigma = 0.4;
    double mass = 1.0 / (sigma * sigma);
    double* y_g = myres_b->create_fake(mean, sigma, 33);

    mean1 = 0;
    for (int j = 0;j < Neff;j++) mean1 += y_g[j];
    mean1 /= Neff;
    printf("<x>  computed stochastically :  %-20.12g%-20.12g\n", mean1, myres_b->comp_error(y_g) / sqrt(Neff - 1));
    printf("expected                     :  %-20.12g%-20.12g\n", mean, sigma / sqrt(Neff));

    double* y1 = new double[Neff];
    for (int j = 0;j < Neff;j++) {
        y1[j] = y_g[j] * y_g[j];
    }
    mean1 = 0;
    for (int j = 0;j < Neff;j++) mean1 += y1[j];
    mean1 /= Neff;
    printf("<x^2> computed stochastically:  %-20.12g%-20.12g\n", mean1, myres_b->comp_error(y1) / sqrt(Neff - 1));
    printf("expected                     :  %-20.12g%-20.12g\n", mean * mean + sigma * sigma, NAN);

    /////////////////////////////////////////// creating jacks
    double**** data = calloc_corr(Neff, 2, 1);
    for (int j = 0;j < Neff;j++) {
        data[j][0][0][0] = y_g[j];
        data[j][1][0][0] = y_g[j] * y_g[j];
    }
    double**** jacks = myres->create(Neff, 2, 1, data);
    /////////////////////////////////////////// end creating jacks

    for (int j = 0;j < Neff;j++) {
        y1[j] = jacks[j][0][0][0];
    }
    mean1 = 0;
    for (int j = 0;j < Neff;j++) mean1 += y1[j];
    mean1 /= Neff;
    printf("<x> computed jacks           :  %-20.12g%-20.12g\n", mean1, myres->comp_error(y1));
    printf("expected                     :  %-20.12g%-20.12g\n", mean, sigma / sqrt(Neff));

    for (int j = 0;j < Neff;j++) {
        y1[j] = jacks[j][1][0][0];
    }
    mean1 = 0;
    for (int j = 0;j < Neff;j++) mean1 += y1[j];
    mean1 /= Neff;
    double err_on_sigma2 = myres->comp_error(y1);
    printf("<x^2> computed jacks         :  %-20.12g%-20.12g\n", mean1, myres->comp_error(y1));
    printf("expected                     :  %-20.12g%-20.12g\n", sigma * sigma, NAN);

    for (int j = 0;j < Neff;j++) {
        y1[j] = 1.0 / (jacks[j][0][0][0]);
    }
    mean1 = 0;
    for (int j = 0;j < Neff;j++) mean1 += y1[j];
    mean1 /= Neff;
    printf("1/<x> computed jacks         :  %-20.12g%-20.12g\n", mean1, myres->comp_error(y1));
    printf("expected                     :  %-20.12g%-20.12g\n", 1.0 / mean, (1.0 / (mean * mean)) * sigma / sqrt(Neff));

    for (int j = 0;j < Neff;j++) {
        y1[j] = 1.0 / (jacks[j][1][0][0]);
    }
    mean1 = 0;
    for (int j = 0;j < Neff;j++) mean1 += y1[j];
    mean1 /= Neff;
    double x2 = mean * mean + sigma * sigma;
    printf("1/<x^2> computed jacks       :  %-20.12g%-20.12g\n", mean1, myres->comp_error(y1));
    printf("expected                     :  %-20.12g%-20.12g\n", 1.0 / (x2), (1.0 / (x2 * x2)) * err_on_sigma2);


    // double deviation = fabs(error_jack - err1) / err1;
    // if (fabs(mean1 - xj[Njack - 1]) > err) { printf("error in the mean computation jack :  %g   by hand: %g\n", xj[Njack - 1], mean1); exit(1); }
    // if (deviation > 0.1) { printf("mismatch in the error propagation jack :  %g   by hand: %g  deviation: %g\n", error_jack, err1, deviation); exit(1); }

    // printf("success... retry with different seed\n");

    // free(x);
    // free(xj);
    // double max_deviation = -100;
    // for (int seed = 2;seed < 1000;seed++) {
    //     double* x = myres->create_fake(mean, err, seed);
    //     double* xj = myres->create_copy(x);
    //     for (int j = 0;j < Njack;j++) {
    //         xj[j] = 1.0 / x[j];
    //     }
    //     double error_jack = myres->comp_error(xj);

    //     double deviation = fabs(error_jack - err1) / err1;
    //     if (fabs(mean1 - xj[Njack - 1]) > err) { printf("seed: %d error in the mean computation jack :  %g   by hand: %g\n", seed, xj[Njack - 1], mean1); exit(1); }
    //     if (deviation > 0.3) { printf("seed: %d mismatch in the error propagation jack :  %g   by hand: %g\n", seed, error_jack, err1); exit(1); }
    //     if (deviation > max_deviation) max_deviation = deviation;
    //     // printf("%g\n",deviation);
    //     free(x);
    //     free(xj);

    // }
    // printf("max_deviation : %g \n", max_deviation);
    // printf("success\n");
    return 0;
}