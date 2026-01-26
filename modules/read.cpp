#define read_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "mutils.hpp"
#include "read.hpp"
#include "global.hpp"

// option[3] path
// option[6] file
template<class T>
void line_read_param(char** option, const char* corr, T& value, T& error, int& seed, const char* namefile_plateaux) {

    int line = 0;
    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], namefile_plateaux);
    std::fstream newfile;

    newfile.open(namefile, std::ios::in); // open a file to perform read operation using file object
    int match = 0;
    if (newfile.is_open()) { // checking whether the file is open
        std::string tp;
        while (getline(newfile, tp)) { // read data from file object and put it into string.
            line++;
            std::vector<std::string> x = split(tp, ' ');

            std::string name = option[6];
            std::string correlator = corr;
            if (x.empty() == 0) {
                if (x[0].compare(name) == 0 && x[1].compare(correlator) == 0) {
                    value = convert<T>(x[2]);
                    error = convert<T>(x[3]);
                    seed = stoi(x[4]);
                    std::cout << "correlator " << correlator << " name " << name << " value " <<
                        value << " error " << error << " seed " << seed << "\n";
                    match++;
                    // break;
                }
            }
        }
        newfile.close(); // close the file object.
    }
    else {
        printf("correlators_analysis.cpp line_read_plateaux\n");
        printf("\t unable to open %s", namefile);
        exit(1);
    }
    // error(match==0,1,"correlators_analysis.cpp line_read_plateaux",
    //       "no match for plateau %s   %s \n in the file %s ",option[6], corr,namefile);
    if (match == 0) {
        printf("no plateau found for %s in plateau file %s\n", corr, namefile);
        printf("looking for a line:\n %s  %s\n", option[6], corr);
        exit(1);
    }
    if (match > 1) {
        printf("multiple lines line:\n %s  %s\n", option[6], corr);
        exit(1);
    }

}
template void line_read_param<int>(char**, const char*, int&, int&, int&, const char*);
template void line_read_param<double>(char**, const char*, double&, double&, int&, const char*);
template void line_read_param<std::string>(char**, const char*, std::string&, std::string&, int&, const char*);

int how_many_matches_in_line_read_param(char** option, const char* corr, const char* namefile_plateaux) {

    int line = 0;
    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], namefile_plateaux);
    std::fstream newfile;

    newfile.open(namefile, std::ios::in); // open a file to perform read operation using file object
    int match = 0;
    if (newfile.is_open()) { // checking whether the file is open
        std::string tp;
        while (getline(newfile, tp)) { // read data from file object and put it into string.
            line++;
            std::vector<std::string> x = split(tp, ' ');

            std::string name = option[6];
            std::string correlator = corr;
            if (x.empty() == 0) {
                if (x[0].compare(name) == 0 && x[1].compare(correlator) == 0) {
                    // value = convert<T>(x[2]);
                    // error = convert<T>(x[3]);
                    // seed = stoi(x[4]);
                    // std::cout << "correlator " << correlator << " name " << name << " value " <<
                    //     value << " error " << error << " seed " << seed << "\n";
                    match++;
                    // break;
                }
            }
        }
        newfile.close(); // close the file object.
    }
    else {
        printf("correlators_analysis.cpp line_read_plateaux\n");
        printf("\t unable to open %s", namefile);
        exit(1);
    }
    
    return match;
}


double**** malloc_corr(int N, int var, int t) {
    double**** out;
    int i, j, k;

    j = var;
    out = (double****)malloc(sizeof(double***) * N);
    for (i = 0;i < N;i++) {
        out[i] = (double***)malloc(sizeof(double**) * var);
        for (j = 0;j < var;j++) {
            out[i][j] = (double**)malloc(sizeof(double*) * t);
            for (k = 0;k < t;k++)
                out[i][j][k] = (double*)malloc(2 * sizeof(double));
        }
    }
    return out;
}
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
    error(Nb > N, 1, "binning_toNb", "not possible to creare more binning=%d than confs=%d", Nb, N);

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
                        out[l][j][k][0] += data[i + l * bin][j][k][0];
                        out[l][j][k][1] += data[i + l * bin][j][k][1];
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

double**** bin_intoN(double**** data, int nvar, int T, int Nconf_in, int Nb) {


    double clustSize = ((double)Nconf_in) / ((double)Nb);
    // double clustSize = ((double)confs.confs_after_binning) / ((double)Nb);
    double**** to_write = calloc_corr(Nb, nvar, T);
    for (size_t iClust = 0;iClust < Nb;iClust++) {

        /// Initial time of the bin
        const double binBegin = iClust * clustSize;
        /// Final time of the bin
        const double binEnd = binBegin + clustSize;
        double binPos = binBegin;
        do {
            /// Index of the configuration related to the time
            const size_t iConf = floor(binPos + 1e-10);

            ///Rectangle left point
            const double beg = binPos;

            /// Rectangle right point
            const double end = std::min(binEnd, iConf + 1.0);

            /// Rectangle horizontal size
            const double weight = end - beg;

            // Perform the operation passing the info
            for (int t = 0;t < T;t++) {
                for (int ivar = 0; ivar < nvar; ivar++) {
                    to_write[iClust][ivar][t][0] += weight * data[iConf][ivar][t][0];
                    to_write[iClust][ivar][t][1] += weight * data[iConf][ivar][t][1];
                }
            }
            // printf("Cluster=%ld  iConf=%ld  weight=%g  size=%g  end=%g  beg=%g\n",iClust,iConf,weight,clustSize, end,beg);
            // Updates the position
            binPos = end;
        } while (binEnd - binPos > 1e-10);
        for (int t = 0;t < T;t++) {
            for (int ivar = 0; ivar < nvar; ivar++) {
                to_write[iClust][ivar][t][0] /= ((double)clustSize);
                to_write[iClust][ivar][t][1] /= ((double)clustSize);
            }
        }
    }
    return to_write;
}

//data[#conf.][#variable][#time_cordinate][#re or im]
double**** read_datafile(int N, int var, int t, int bin) {
    char** datafile;
    FILE** f;
    double**** data, **** out;
    int i, j, k, fi = 0;

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
                fi += fscanf(f[j], "%lf   %lf \n", &data[i][j][k][0], &data[i][j][k][1]);
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

void cbtc_corr(int N, int var, int t, double**** data) {
    int i, j, k;
    double**** out;

    out = calloc_corr(N, var + 1, t);

    j = var;
    for (i = 0;i < N;i++) {
        for (k = 0;k < t;k++) {
            if (k == t - 1) {
                out[i][j][k][0] = (data[i][j][0][0] + 2. * data[i][j][k][0] + data[i][j][k - 1][0]) / 4.;
                out[i][j][k][1] = (data[i][j][0][1] + 2. * data[i][j][k][1] + data[i][j][k - 1][1]) / 4.;
            }
            else if (k == 0) {
                out[i][j][k][0] = (data[i][j][k + 1][0] + 2 * data[i][j][k][0] + data[i][j][t - 1][0]) / 4.;
                out[i][j][k][1] = (data[i][j][k + 1][1] + 2 * data[i][j][k][1] + data[i][j][t - 1][1]) / 4.;
            }
            else {
                out[i][j][k][0] = (data[i][j][k + 1][0] + 2. * data[i][j][k][0] + data[i][j][k - 1][0]) / 4.;
                out[i][j][k][1] = (data[i][j][k + 1][1] + 2. * data[i][j][k][1] + data[i][j][k - 1][1]) / 4.;
            }
        }
    }

    copy_corr(N, var, t, out, data);
    free_corr(N, var + 1, t, out);

}

