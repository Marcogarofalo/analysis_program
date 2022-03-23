#define read_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "mutils.hpp"
#include "read.hpp"


double**** calloc_corr(int N, int var, int t) {
    double**** out;
    int i, j, k;

    j = var;
    out = (double****)malloc(sizeof(double***) * N);
    for (i = 0;i < N;i++) {
        out[i] = (double***)malloc(sizeof(double**) * var);
        for (j = 0;j < var;j++) {
            out[i][j] = (double**)malloc(sizeof(double*) * t);
            for (k = 0;k < t;k++)
                out[i][j][k] = (double*)calloc(2, sizeof(double));
        }
    }
    return out;
}

void copy_corr(int N, int var, int t, double**** in, double**** out) {
    int i, j, k;

    j = var;
    for (i = 0;i < N;i++) {
        for (k = 0;k < t;k++) {
            out[i][j][k][0] = in[i][j][k][0];
            out[i][j][k][1] = in[i][j][k][1];
        }
    }
}

void free_corr(int N, int var, int t, double**** out) {
    int i, j, k;

    for (i = 0;i < N;i++) {
        for (j = 0;j < var;j++) {
            for (k = 0;k < t;k++) {
                free(out[i][j][k]);
            }
            free(out[i][j]);
        }
        free(out[i]);
    }
    free(out);
}

double**** binning(int N, int var, int t, double**** data, int bin) {
    int i, j, k, l;
    double**** out;

    out = calloc_corr(N / bin, var, t);

    if (bin == 1)
        for (j = 0;j < var;j++)
            copy_corr(N / bin, j, t, data, out);
    else {

        for (j = 0;j < var;j++) {
            for (k = 0;k < t;k++) {
                l = 0;
                for (i = 0;i < (N / bin) * bin;i++) {
                    if (i != 0) if (i % bin == 0) l++;
                    out[l][j][k][0] += data[i][j][k][0];
                    out[l][j][k][1] += data[i][j][k][1];
                }
                for (l = 0;l < N / bin;l++) {
                    out[l][j][k][0] /= ((double)bin);
                    out[l][j][k][1] /= ((double)bin);
                }

            }
        }

    }

    return out;
}

double**** binning_toNb(int N, int var, int t, double**** data, int Nb) {
    int i, j, k, l;
    double**** out;

    out = calloc_corr(Nb, var, t);

    if (Nb == N)
        for (j = 0;j < var;j++)
            copy_corr(Nb, j, t, data, out);
    else {
        int bin = N / Nb;
        for (j = 0;j < var;j++) {
            for (k = 0;k < t;k++) {
                l = 0;
                for (i = 0;i < bin;i++) {
                    for (l = 0;l < Nb;l++) {
                        out[l][j][k][0] += data[i+l*bin][j][k][0];
                        out[l][j][k][1] += data[i+l*bin][j][k][1];
                    }
                }
                for (l = 0;l < Nb;l++) {
                    out[l][j][k][0] /= ((double)bin);
                    out[l][j][k][1] /= ((double)bin);
                }

            }
        }

    }

    return out;
}


//data[#conf.][#variable][#time_cordinate][#re or im]
double**** read_datafile(int N, int var, int t, int bin) {
    char** datafile;
    FILE** f;
    double**** data, **** out;
    int i, j, k;

    data = (double****)malloc(sizeof(double***) * N);
    for (i = 0;i < N;i++) {
        data[i] = (double***)malloc(sizeof(double**) * var);
        for (j = 0;j < var;j++) {
            data[i][j] = (double**)malloc(sizeof(double*) * t);
            for (k = 0;k < t;k++)
                data[i][j][k] = (double*)malloc(2 * sizeof(double));
        }
    }

    datafile = (char**)malloc(sizeof(char*) * var);
    f = (FILE**)malloc(sizeof(FILE*) * var);
    for (j = 0;j < var;j++) {
        datafile[j] = (char*)malloc(sizeof(char) * var);
        sprintf(datafile[j], "in%d", j);
        f[j] = fopen(datafile[j], "r+");
        error(f[j] == NULL, 1, "read_datafile ",
            "Unable to open input data file %s", datafile[j]);

        for (i = 0;i < N;i++) {
            for (k = 0;k < t;k++) {
                fscanf(f[j], "%lf   %lf \n", &data[i][j][k][0], &data[i][j][k][1]);
            }
        }
        fclose(f[j]);
        free(datafile[j]);
    }
    free(f);
    free(datafile);

    out = binning(N, var, t, data, bin);
    free_corr(N, var, t, data);

    return out;
}



void symmetrise_corr(int N, int var, int t, double**** data) {
    int i, j, k;
    double**** out;

    out = calloc_corr(N, var + 1, t);

    j = var;
    for (i = 0;i < N;i++) {
        for (k = 0;k < t;k++) {
            if (k == 0) {
                out[i][j][k][0] = data[i][j][k][0];
                out[i][j][k][1] = data[i][j][k][1];
            }
            else {
                out[i][j][k][0] = (data[i][j][k][0] + data[i][j][t - k][0]) / 2.;
                out[i][j][k][1] = (data[i][j][k][1] + data[i][j][t - k][1]) / 2.;
            }
        }
    }

    copy_corr(N, var, t, out, data);
    free_corr(N, var + 1, t, out);

}



void antisymmetrise_corr(int N, int var, int t, double**** data) {
    int i, j, k;
    double**** out;

    out = calloc_corr(N, var + 1, t);

    j = var;
    for (i = 0;i < N;i++) {
        for (k = 0;k < t;k++) {
            if (k == 0) {
                out[i][j][k][0] = data[i][j][k][0];
                out[i][j][k][1] = data[i][j][k][1];
            }
            else {
                out[i][j][k][0] = (data[i][j][k][0] - data[i][j][t - k][0]) / 2.;
                out[i][j][k][1] = (data[i][j][k][1] - data[i][j][t - k][1]) / 2.;
            }
        }
    }

    copy_corr(N, var, t, out, data);
    free_corr(N, var + 1, t, out);

}

void forward_derivative_corr(int N, int var, int t, double**** data) {
    int i, j, k;
    double**** out;

    out = calloc_corr(N, var + 1, t);

    j = var;
    for (i = 0;i < N;i++) {
        for (k = 0;k < t;k++) {
            if (k == t - 1) {
                out[i][j][k][0] = data[i][j][0][0] - data[i][j][k][0];
                out[i][j][k][1] = data[i][j][0][1] - data[i][j][k][1];
            }
            else {
                out[i][j][k][0] = data[i][j][k + 1][0] - data[i][j][k][0];
                out[i][j][k][1] = data[i][j][k + 1][1] - data[i][j][k][1];
            }
        }
    }

    copy_corr(N, var, t, out, data);
    free_corr(N, var + 1, t, out);

}

void symmetric_derivative_corr(int N, int var, int t, double**** data) {
    int i, j, k;
    double**** out;

    out = calloc_corr(N, var + 1, t);

    j = var;
    for (i = 0;i < N;i++) {
        for (k = 0;k < t;k++) {
            if (k == t - 1) {
                out[i][j][k][0] = (data[i][j][0][0] - data[i][j][k - 1][0]) / 2.;
                out[i][j][k][1] = (data[i][j][0][1] - data[i][j][k - 1][1]) / 2.;
            }
            else if (k == 0) {
                out[i][j][k][0] = (data[i][j][k + 1][0] - data[i][j][t - 1][0]) / 2.;
                out[i][j][k][1] = (data[i][j][k + 1][1] - data[i][j][t - 1][1]) / 2.;
            }
            else {
                out[i][j][k][0] = (data[i][j][k + 1][0] - data[i][j][k - 1][0]) / 2.;
                out[i][j][k][1] = (data[i][j][k + 1][1] - data[i][j][k - 1][1]) / 2.;
            }
        }
    }

    copy_corr(N, var, t, out, data);
    free_corr(N, var + 1, t, out);

}


void second_derivative_corr(int N, int var, int t, double**** data) {
    int i, j, k;
    double**** out;

    out = calloc_corr(N, var + 1, t);

    j = var;
    for (i = 0;i < N;i++) {
        for (k = 0;k < t;k++) {
            if (k == t - 1) {
                out[i][j][k][0] = data[i][j][0][0] - 2. * data[i][j][k][0] + data[i][j][k - 1][0];
                out[i][j][k][1] = data[i][j][0][1] - 2. * data[i][j][k][1] + data[i][j][k - 1][1];
            }
            else if (k == 0) {
                out[i][j][k][0] = data[i][j][k + 1][0] - 2. * data[i][j][k][0] + data[i][j][t - 1][0];
                out[i][j][k][1] = data[i][j][k + 1][1] - 2. * data[i][j][k][1] + data[i][j][t - 1][1];
            }
            else {
                out[i][j][k][0] = data[i][j][k + 1][0] - 2. * data[i][j][k][0] + data[i][j][k - 1][0];
                out[i][j][k][1] = data[i][j][k + 1][1] - 2. * data[i][j][k][1] + data[i][j][k - 1][1];
            }
        }
    }

    copy_corr(N, var, t, out, data);
    free_corr(N, var + 1, t, out);

}

