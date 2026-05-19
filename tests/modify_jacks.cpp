#define CONTROL
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#include "global.hpp"
#include "resampling_new.hpp"
#include "non_linear_fit.hpp"
#include "read.hpp"
#include "tower.hpp"
#include <random>


int main() {
    int Neff = 10;
    int Njack = Neff + 1;

    double mean = 0.1;
    double err = 0.002;
    myres = new resampling_jack(Neff);

    double* x_j = myres->create_fake(mean, err, 533);
    printf("mean and error computed with jack :  %-20.12g%-20.12g\n", x_j[Njack - 1], myres->comp_error(x_j));

    double* shifted = myres->create_new_jack_correlated(x_j, 1.0, 0.2);
    printf("mean and error computed with jack shifted:  %-20.12g%-20.12g\n", shifted[Njack - 1], myres->comp_error(shifted));

    double** tmp = (double**)malloc(sizeof(double*) * 2);
    tmp[0] = x_j;
    tmp[1] = shifted;
    double** cov = myres->comp_cov(2, tmp);
    free_2(2, tmp);
    printf("correlation matrix:\n");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            printf("%-20.12g ", cov[i][j] / std::sqrt(cov[i][i] * cov[j][j]));
        }        printf("\n");
    }

    /// create cov exact
    int N = 5;
    unsigned int seed = 42;
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(0.0, 10.0);

    double* mean_vec = (double*)malloc(sizeof(double) * N);
    for (int i = 0;i < N;i++)
        mean_vec[i] = dis(gen);

    // 1. Generate a random NxN matrix A
    std::vector<std::vector<double>> A(N, std::vector<double>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = dis(gen);
        }
    }

    // 2. Compute C = A * A^T (Symmetric & Positive Semi-Definite)
    double** cov_exact = (double**)malloc(sizeof(double*) * N);
    for (int i = 0; i < N; i++) {
        cov_exact[i] = (double*)malloc(sizeof(double) * N);
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cov_exact[i][j] = 0;
            for (int k = 0; k < N; ++k) {
                cov_exact[i][j] += A[i][k] * A[j][k];
            }
        }
        // cov_exact[i][i] += 1; // make it positive def
    }



    double** jacks = myres->create_fake_covariance(mean_vec, N, cov_exact, 533);
    printf("generated jacks\n");
    for (int i = 0; i < N; i++) {
        printf("mean[%d]: %g\n", i, myres->mean(jacks[i]) - mean_vec[i]);
    }

    printf(" devition of covariance matrix:\n");
    double** cov_j = myres->comp_cov(N, jacks);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // printf("%-20.12g ", std::fabs(cov_j[i][j]- cov_exact[i][j])/cov_exact[i][j]  );
            printf("%-20.12g ", std::fabs(cov_j[i][j] - cov_exact[i][j]));
        }        printf("\n");
    }
    free_2(N, cov_j);


    double** jacks1 = myres->create_fake_covariance(mean_vec, N, cov_exact, 533);
    printf("linear transformation to make the cov exact\n");
    myres->change_mean_and_error_covarinace(jacks1, jacks, N, mean_vec, cov_exact);

    printf("After linear transformation\n");
    for (int i = 0; i < N; i++) {
        printf("mean[%d]: %g\n", i, myres->mean(jacks[i]) - mean_vec[i]);
    }

    printf(" devition of covariance matrix:\n");
    cov_j = myres->comp_cov(N, jacks1);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double dev = std::fabs(cov_j[i][j] - cov_exact[i][j]);
            printf("%-20.12g ", dev);
            if (dev > 1e-11) {
                printf("\nerror: deviation between the cov matrix given and the generated one\n");
                exit(1);
            }

        }        printf("\n");
    }

    free(mean_vec);
    free_2(N, cov_exact);
    free_2(N, jacks);
    free_2(N, jacks1);
    free_2(N, cov_j);
    free_2(2, cov);

}