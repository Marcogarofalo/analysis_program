#define resampling_new_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "mutils.hpp"
#include "resampling_new.hpp"
#include "rand.hpp"
#include "global.hpp"
#include "tower.hpp"

#ifdef WITH_ARB
#include "arb.h"
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
    Njack = Njack + 1;


    boot = (double****)malloc(sizeof(double***) * (Njack + 1));
    for (i = 0;i <= Njack;i++) {
        boot[i] = (double***)malloc(sizeof(double**) * var);
        for (j = 0;j < var;j++) {
            boot[i][j] = (double**)malloc(sizeof(double*) * t);
            for (k = 0;k < t;k++)
                boot[i][j][k] = (double*)calloc(2, sizeof(double));
        }
    }


    for (ib = 0;ib < Njack;ib++) {
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
                    boot[Njack][j][k][l] += in[i][j][k][l];
                }
            }
        }
    }

    for (i = 0;i <= Njack;i++) {
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
    printf("error: resampling_boot::comp_cov_arb not implemented\n");
    exit(1);
}
#endif // WITH_ARB
