#define resampling_new_C

#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#include "mutils.hpp"
#include "resampling_new.hpp"
#include "rand.hpp"
#include "global.hpp"
#include "tower.hpp"

#ifdef WITH_ARB
#include "arb.h"
#include "arb_mat.h"
#endif // WITH_ARB
// resampling::resampling(const char* name, int N) {
//     Njack == N;
//     if (strcmp(name, "jack") == 0 || strcmp(name, "boot") == 0) {
//         mysprintf(option, NAMESIZE, name);
//         *this = new resampling_jack();
//     }
//     if (strcmp(name, "boot") == 0) {
//         mysprintf(option, NAMESIZE, name);
//         *this = new resampling_boot();
//     }
//     else {
//         printf("Error you need to initialise the resampling procedure");
//         exit(2);
//     }

// }



void resampling_f::free_res(int var, int t, double**** in) {
    int i, j, k;

    for (i = 0;i < Njack;i++) {
        for (j = 0;j < var;j++) {
            for (k = 0;k < t;k++)
                free(in[i][j][k]);
            free(in[i][j]);
        }
        free(in[i]);
    }
    free(in);
}

void resampling_f::write_jack_in_file(double* jack, const char* outname) {
    FILE* f;
    f = open_file(outname, "w+");
    fprintf(f, "%d\n", Njack);
    for (int j = 0;j < Njack;j++)  fprintf(f, "%.12g\n", jack[j]);
    fclose(f);
}
void resampling_f::read_jack_from_file(double* jack, const char* name) {
    FILE* f;
    int i;
    int ic = 0;
    f = open_file(name, "r");
    ic = fscanf(f, "%d\n", &i);
    error(i != Njack, 1, "read_jack_from_file", "number of jeack do not match: \nexpected %d\nread %d\n", Njack, i);
    for (int j = 0;j < Njack;j++)  ic += fscanf(f, "%lf\n", &jack[j]);
    error(ic != Njack + 1, 1, "read_jack_from_file", "invalid read counter: \nexpected %d\nread %d\n", Njack + 1, ic);
    fclose(f);
}


double* resampling_f::create_copy(double* in) {
    double* out = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++)  out[j] = in[j];
    return(out);
}
void resampling_f::copy(double* out, double* in) {
    for (int j = 0;j < Njack;j++)  out[j] = in[j];

}
void resampling_f::add(double* out, double* in1, double* in2) {
    for (int j = 0;j < Njack;j++)  out[j] = in1[j] + in2[j];

}
void resampling_f::add(double* out, double* in1, double a) {
    for (int j = 0;j < Njack;j++)  out[j] = in1[j] + a;
}
void resampling_f::add(double* out, double* in1) {
    for (int j = 0;j < Njack;j++)  out[j] += in1[j];
}
void resampling_f::add(double* out, double a) {
    for (int j = 0;j < Njack;j++)  out[j] += a;
}

void resampling_f::sub(double* out, double* in1, double* in2) {
    for (int j = 0;j < Njack;j++)  out[j] = in1[j] - in2[j];
}
void resampling_f::sub(double* out, double* in1, double a) {
    for (int j = 0;j < Njack;j++)  out[j] = in1[j] - a;
}
void resampling_f::sub(double* out, double* in1) {
    for (int j = 0;j < Njack;j++)  out[j] -= in1[j];
}
void resampling_f::sub(double* out, double a) {
    for (int j = 0;j < Njack;j++)  out[j] -= a;
}

void resampling_f::mult(double* out, double* in1, double* in2) {
    for (int j = 0;j < Njack;j++)  out[j] = in1[j] * in2[j];

}
void resampling_f::mult(double* out, double* in1, double a) {
    for (int j = 0;j < Njack;j++)  out[j] = in1[j] * a;

}
void resampling_f::div(double* out, double* in1, double* in2) {
    for (int j = 0;j < Njack;j++)  out[j] = in1[j] / in2[j];

}
void resampling_f::div(double* out, double* in1, double a) {
    for (int j = 0;j < Njack;j++)  out[j] = in1[j] / a;

}
void resampling_f::linear_comb(double* out, double a, double* in1, double b, double* in2) {
    for (int j = 0;j < Njack;j++)  out[j] = a * in1[j] + b * in2[j];

}


void resampling_f::write_res_bin(int N, double* jack, char* outname) {
    int j;
    FILE* f;

    f = fopen(outname, "ab");
    error(f == NULL, 1, "write_jack_bin ",
        "Unable to open output file %s", outname);

    fwrite(jack, sizeof(double), N, f);

    fclose(f);
}

resampling_jack::resampling_jack(int N) {
    Njack = N + 1;
    mysprintf(option, NAMESIZE, "jack");
}
resampling_boot::resampling_boot(int N) {
    Njack = N + 1;
    mysprintf(option, NAMESIZE, "boot");
    seed = 123;
}
resampling_boot::resampling_boot(int N, int o_seed) {
    Njack = N + 1;
    mysprintf(option, NAMESIZE, "boot");
    seed = o_seed;
}

//create_jack
//in[#conf.][#variable][#time_cordinate][#re or im]
//returns the jacknife configuration from the data ****in
//the last entry of [#conf] is the average
double**** resampling_jack::create(int  N, int var, int t, double**** in) {
    double**** jack;
    int i, j, k, l;
    Njack = N + 1;
    jack = (double****)malloc(sizeof(double***) * (N + 1));
    for (i = 0;i <= N;i++) {
        jack[i] = (double***)malloc(sizeof(double**) * var);
        for (j = 0;j < var;j++) {
            jack[i][j] = (double**)malloc(sizeof(double*) * t);
            for (k = 0;k < t;k++)
                jack[i][j][k] = (double*)calloc(2, sizeof(double));
        }
    }


    for (i = 0;i < N;i++) {
        for (j = 0;j < var;j++) {
            for (k = 0;k < t;k++) {
                for (l = 0;l < 2;l++) {
                    jack[N][j][k][l] += in[i][j][k][l];
                }
            }
        }
    }

    for (i = 0;i < N;i++) {
        for (j = 0;j < var;j++) {
            for (k = 0;k < t;k++) {
                for (l = 0;l < 2;l++) {
                    jack[i][j][k][l] = (jack[N][j][k][l] - in[i][j][k][l]) / ((double)(N - 1));
                }
            }
        }
    }

    for (j = 0;j < var;j++) {
        for (k = 0;k < t;k++) {
            for (l = 0;l < 2;l++) {
                jack[N][j][k][l] /= ((double)(N));
            }
        }
    }


    return  jack;
}


double**** resampling_boot::create(int  N, int var, int t, double**** in) {

    double**** boot;
    int i, j, k, l, b, ib;

    srand(seed);

    boot = (double****)malloc(sizeof(double***) * (Njack));
    for (i = 0;i < Njack;i++) {
        boot[i] = (double***)malloc(sizeof(double**) * var);
        for (j = 0;j < var;j++) {
            boot[i][j] = (double**)malloc(sizeof(double*) * t);
            for (k = 0;k < t;k++)
                boot[i][j][k] = (double*)calloc(2, sizeof(double));
        }
    }


    for (ib = 0;ib < (Njack - 1);ib++) {
        for (i = 0;i < N;i++) {
            b = (int)rand() % N;  //have to stay here, because I want that all the var, t , l to have the same random number 
            for (j = 0;j < var;j++) {
                for (k = 0;k < t;k++) {
                    for (l = 0;l < 2;l++) {
                        boot[ib][j][k][l] += in[b][j][k][l];
                    }
                }
            }
        }
    }

    for (i = 0;i < N;i++) {
        for (j = 0;j < var;j++) {
            for (k = 0;k < t;k++) {
                for (l = 0;l < 2;l++) {
                    boot[Njack - 1][j][k][l] += in[i][j][k][l];
                }
            }
        }
    }

    for (i = 0;i < Njack;i++) {
        for (j = 0;j < var;j++) {
            for (k = 0;k < t;k++) {
                for (l = 0;l < 2;l++) {
                    boot[i][j][k][l] /= ((double)(N));
                }
            }
        }
    }

    return  boot;
}


double resampling_jack::comp_mean_unbias(double* in) {
    double r = 0;
    int i, N;

    N = Njack - 1;

    for (i = 0;i < N;i++)
        r += in[i];

    r /= ((double)N);

    r = ((double)N) * in[N] - ((double)(N - 1)) * r;


    return r;
}


double resampling_jack::comp_error(double* in) {
    double r[2] = { 0,0 };
    int i, N;

    N = Njack - 1;

    for (i = 0;i < N;i++)
        r[0] += in[i];

    r[0] /= ((double)N);


    for (i = 0;i < N;i++) {
        r[1] += (in[i] - r[0]) * (in[i] - r[0]);
        //error(in[i]!=in[i],1,"mean_and_error_jack_biased","errore jack=%d is nan",i);
    }
    r[1] *= (N - 1.) / ((double)N);
    r[1] = sqrt(r[1]);

    // r[0] = in[N];

    return r[1];
}

void resampling_f::mult_error(double* out, double* in, double d) {
    double r = 0;
    int N = Njack - 1;

    double kappa = 1 - d;
    for (int j = 0;j < N;j++)
        r += in[j];

    r /= ((double)N);

    for (int j = 0;j < N;j++) {
        out[j] = in[j] + (r - in[j]) * kappa;
        //error(in[i]!=in[i],1,"mean_and_error_jack_biased","errore jack=%d is nan",i);
    }
    out[N] = in[N];

}

void resampling_f::add_error_quadrature(double* out, double* in, double d) {
    double err = comp_error(in);
    double kappa = sqrt(1 + d * d / (err * err));
    out[Njack - 1] = in[Njack - 1];
    mult_error(out, in, kappa);
}



double* resampling_jack::create_fake(double mean, double error, int seed) {
    int i;
    double* r;

    r = (double*)malloc(sizeof(double) * Njack);
    srand(seed);
    for (i = 0;i < Njack - 1;i++) {
        r[i] = generateGaussianNoise(mean, error / sqrt(Njack - 2));
    }
    r[Njack - 1] = mean;

    return r;
}


double* resampling_boot::create_fake(double mean, double error, int seed) {
    int i;
    double* r;

    r = (double*)malloc(sizeof(double) * Njack);
    srand(seed);
    for (i = 0;i < Njack - 1;i++) {
        r[i] = generateGaussianNoise(mean, error);
    }
    r[Njack - 1] = mean;

    return r;
}

double resampling_boot::comp_mean_unbias(double* in) {
    double r = 0;
    int i, N;

    printf("\n\n ########################### calling resampling_boot::comp_mean_unbias NEVER TESTED ##########\n");

    N = Njack - 1;

    for (i = 0;i < N;i++)
        r += in[i];

    r /= ((double)N);

    r = 2 * in[N] - r;

    return r;
}



//returns the  error from set of N  jacknife called *in  and the average stored in in[N]
double resampling_boot::comp_error(double* in) {
    double r[2] = { 0,0 };
    int i, N;

    N = Njack - 1;

    for (i = 0;i < N;i++)
        r[0] += in[i];

    r[0] /= ((double)N);


    for (i = 0;i < N;i++)
        r[1] += (in[i] - r[0]) * (in[i] - r[0]);


    r[1] /= ((double)(N));
    r[1] = sqrt(r[1]);

    return r[1];
}

//////////////////// cov

double** covariance_jack(int Nobs, int Np1, double** in) {
    double** r, * ave;
    int i, k, l, N;

    N = Np1 - 1;
    ave = (double*)calloc(Nobs, sizeof(double));
    r = (double**)malloc(Nobs * sizeof(double*));

    for (k = 0;k < Nobs;k++) {
        for (i = 0;i < N;i++)
            ave[k] += in[k][i];
        ave[k] /= ((double)N);
        // ave[k] = in[k][N];
        r[k] = (double*)calloc(Nobs, sizeof(double));
    }
    for (k = 0;k < Nobs;k++) {
        for (l = k;l < Nobs;l++) {
            for (i = 0;i < N;i++) {
                r[k][l] += (in[k][i] - ave[k]) * (in[l][i] - ave[l]);
                //error(in[i]!=in[i],1,"mean_and_error_jack_biased","errore jack=%d is nan",i);
            }
            r[k][l] *= (N - 1.) / ((double)N);
            //r[k][l]=sqrt(r[k][l]);
        }
        for (l = 0;l < k;l++) {
            //r[l][k]/=sqrt(r[k][k]*r[l][l]);
            r[k][l] = r[l][k];
        }
    }
    free(ave);
    return r;

}

double** covariance_boot(int Nobs, int Np1, double** in) {
    double** r, * ave;
    int i, k, l, N;

    N = Np1 - 1;
    ave = (double*)calloc(Nobs, sizeof(double));
    r = (double**)malloc(Nobs * sizeof(double*));

    for (k = 0;k < Nobs;k++) {
        for (i = 0;i < N;i++)
            ave[k] += in[k][i];
        ave[k] /= ((double)N);
        r[k] = (double*)calloc(Nobs, sizeof(double));
    }
    for (k = 0;k < Nobs;k++) {
        for (l = k;l < Nobs;l++) {
            for (i = 0;i < N;i++) {
                r[k][l] += (in[k][i] - ave[k]) * (in[l][i] - ave[l]);
                //error(in[i]!=in[i],1,"mean_and_error_jack_biased","errore jack=%d is nan",i);
            }
            r[k][l] /= ((double)N);
            //r[k][l]=sqrt(r[k][l]);
        }
        for (l = 0;l < k;l++) {
            //r[l][k]/=sqrt(r[k][k]*r[l][l]);
            r[k][l] = r[l][k];
        }
        //r[k][k]=1;
    }
    free(ave);
    return r;

}

double** resampling_jack::comp_cov(int Nobs, double** in) {
    return covariance_jack(Nobs, Njack, in);
}

double** resampling_boot::comp_cov(int Nobs, double** in) {
    return covariance_boot(Nobs, Njack, in);
}


#ifdef WITH_ARB
void resampling_jack::comp_cov_arb(arb_mat_t r, int Nobs, double** in, slong prec) {

    int N = Njack - 1;
    arb_mat_t ave;
    arb_mat_init(ave, Nobs, 1); // set it to zero
    error(arb_mat_nrows(r) != Nobs, 1, "resampling_jack::comp_cov_arb", "wrong rows number:%d expected: %d", arb_mat_nrows(r), Nobs);
    error(arb_mat_ncols(r) != Nobs, 1, "resampling_jack::comp_cov_arb", "wrong cols number:%d expected: %d", arb_mat_ncols(r), Nobs);
    arb_mat_zero(r);
    // ave = (double*)calloc(Nobs, sizeof(double));
    // r = (double**)malloc(Nobs * sizeof(double*));
    arb_t tmpk;
    arb_init(tmpk);
    for (int k = 0;k < Nobs;k++) {
        for (int i = 0;i < N;i++) {
            arb_set_d(tmpk, in[k][i]);
            arb_add(arb_mat_entry(ave, k, 0), arb_mat_entry(ave, k, 0), tmpk, prec);
            // ave[k] += in[k][i];
        }
        arb_div_ui(arb_mat_entry(ave, k, 0), arb_mat_entry(ave, k, 0), N, prec);
        // ave[k] /= ((double)N);
        // r[k] = (double*)calloc(Nobs, sizeof(double));
    }
    arb_t tmpl;
    arb_init(tmpl);
    for (int k = 0;k < Nobs;k++) {
        for (int l = k;l < Nobs;l++) {
            for (int i = 0;i < N;i++) {
                arb_set_d(tmpk, in[k][i]);
                arb_set_d(tmpl, in[l][i]);
                arb_sub(tmpk, tmpk, arb_mat_entry(ave, k, 0), prec);
                arb_sub(tmpl, tmpl, arb_mat_entry(ave, l, 0), prec);
                arb_addmul(arb_mat_entry(r, k, l), tmpk, tmpl, prec);
                // r[k][l] += (in[k][i] - ave[k]) * (in[l][i] - ave[l]);

            }
            arb_mul_ui(arb_mat_entry(r, k, l), arb_mat_entry(r, k, l), N - 1, prec);
            arb_div_ui(arb_mat_entry(r, k, l), arb_mat_entry(r, k, l), N, prec);
            // r[k][l] *= (N - 1.) / ((double)N);
        }
        for (int l = 0;l < k;l++) {
            arb_set(arb_mat_entry(r, k, l), arb_mat_entry(r, l, k));
            // r[k][l] = r[l][k];
        }
    }
    arb_mat_clear(ave);
    arb_clear(tmpk);
    arb_clear(tmpl);

}

void resampling_boot::comp_cov_arb(arb_mat_t r, int Nobs, double** in, slong prec) {
    // printf("error: resampling_boot::comp_cov_arb not implemented\n");
    // exit(1);

    int N = Njack - 1;
    arb_mat_t ave;
    arb_mat_init(ave, Nobs, 1); // set it to zero
    error(arb_mat_nrows(r) != Nobs, 1, "resampling_boot::comp_cov_arb", "wrong rows number:%d expected: %d", arb_mat_nrows(r), Nobs);
    error(arb_mat_ncols(r) != Nobs, 1, "resampling_boot::comp_cov_arb", "wrong cols number:%d expected: %d", arb_mat_ncols(r), Nobs);
    arb_mat_zero(r);
    // ave = (double*)calloc(Nobs, sizeof(double));
    // r = (double**)malloc(Nobs * sizeof(double*));
    arb_t tmpk;
    arb_init(tmpk);
    for (int k = 0;k < Nobs;k++) {
        for (int i = 0;i < N;i++) {
            arb_set_d(tmpk, in[k][i]);
            arb_add(arb_mat_entry(ave, k, 0), arb_mat_entry(ave, k, 0), tmpk, prec);
            // ave[k] += in[k][i];
        }
        arb_div_ui(arb_mat_entry(ave, k, 0), arb_mat_entry(ave, k, 0), N, prec);
        // ave[k] /= ((double)N);
        // r[k] = (double*)calloc(Nobs, sizeof(double));
    }
    arb_t tmpl;
    arb_init(tmpl);
    for (int k = 0;k < Nobs;k++) {
        for (int l = k;l < Nobs;l++) {
            for (int i = 0;i < N;i++) {
                arb_set_d(tmpk, in[k][i]);
                arb_set_d(tmpl, in[l][i]);
                arb_sub(tmpk, tmpk, arb_mat_entry(ave, k, 0), prec);
                arb_sub(tmpl, tmpl, arb_mat_entry(ave, l, 0), prec);
                arb_addmul(arb_mat_entry(r, k, l), tmpk, tmpl, prec);
                // r[k][l] += (in[k][i] - ave[k]) * (in[l][i] - ave[l]);

            }
            // arb_mul_ui(arb_mat_entry(r, k, l), arb_mat_entry(r, k, l), N - 1, prec);
            arb_div_ui(arb_mat_entry(r, k, l), arb_mat_entry(r, k, l), N, prec);
            // r[k][l] *= (N - 1.) / ((double)N);
        }
        for (int l = 0;l < k;l++) {
            arb_set(arb_mat_entry(r, k, l), arb_mat_entry(r, l, k));
            // r[k][l] = r[l][k];
        }
    }
    arb_mat_clear(ave);
    arb_clear(tmpk);
    arb_clear(tmpl);
}
#endif // WITH_ARB



double** resampling_jack::create_fake_covariance(double* mean, int N, double** cov, int seed) {
    int i, k;
    double** r, ** r1;

    double** covj = (double**)malloc(sizeof(double*) * N);
    for (i = 0;i < N;i++) {
        covj[i] = (double*)malloc(sizeof(double) * N);
        for (k = 0;k < N;k++) {
            covj[i][k] = cov[i][k] / (Njack - 2);
        }
    }

    r = (double**)malloc(sizeof(double*) * Njack);
    srand(seed);
    for (i = 0;i < Njack - 1;i++) {
        r[i] = generate_correlatedNoise(N, mean, covj);//generateGaussianNoise(mean,  error/sqrt(Njack-2));
    }
    r[Njack - 1] = (double*)malloc(sizeof(double) * N);
    for (i = 0;i < N;i++) {
        r[Njack - 1][i] = mean[i];
        free(covj[i]);
    }
    free(covj);

    r1 = swap_indices(Njack, N, r);
    free_2(Njack, r);
    return r1;
}

double** resampling_boot::create_fake_covariance(double* mean, int N, double** cov, int seed) {
    int i;
    double** r, ** r1;

    r = (double**)malloc(sizeof(double*) * Njack);
    srand(seed);
    for (i = 0;i < Njack - 1;i++) {
        r[i] = generate_correlatedNoise(N, mean, cov);
    }
    r[Njack - 1] = (double*)malloc(sizeof(double) * N);
    for (i = 0;i < N;i++)
        r[Njack - 1][i] = mean[i];

    r1 = swap_indices(Njack, N, r);
    free_2(Njack, r);
    return r1;
}