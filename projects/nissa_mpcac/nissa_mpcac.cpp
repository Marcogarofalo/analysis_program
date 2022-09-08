#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "read.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"


#include "correlators_analysis.hpp"
#include "eigensystem.hpp"


#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>
using namespace std;
struct  kinematic kinematic_2pt;

generic_header read_head(FILE* stream) {
    generic_header header;
    size_t fi = 0;
    fi += fscanf(stream, "%d", &header.L);
    fi += fscanf(stream, "%d", &header.T);
    header.mus = std::vector<double>(2);
    fi += fscanf(stream, "%lf", &header.mus[0]);
    fi += fscanf(stream, "%lf", &header.kappa);

    int a, b;
    fi += fscanf(stream, "%d", &a);
    fi += fscanf(stream, "%d", &b);
    header.ncorr = a * b;
    fi += fscanf(stream, "%d", &header.Njack);
    header.size = header.ncorr * 2 * header.T;
    return header;
}
void  write_header_g2(FILE* jack_file, generic_header head) {
    int fi = 0;
    fi += fwrite(&head.T, sizeof(int), 1, jack_file);
    fi += fwrite(&head.L, sizeof(int), 1, jack_file);
    int nmu = head.mus.size();
    fi += fwrite(&nmu, sizeof(int), 1, jack_file);
    for (double mu : head.mus) {
        fi += fwrite(&mu, sizeof(double), 1, jack_file);
    }

}




void read_twopt(FILE* stream, double*** to_write, generic_header head) {
    int fi = 0;
    for (int k = 0; k < head.ncorr; k++) {
        for (int t = 0; t < head.T;t++) {
            fi += fscanf(stream, "%lf  %lf", &to_write[k][t][0], &to_write[k][t][1]);
            // printf(" corr=%d t=%d %.12g   %.12g\n", k, t, to_write[k][t][0], to_write[k][t][1]);
        }
    }

}



double mpcac(int j, double**** in, int t, struct fit_type fit_info) {


    double r = -(in[j][1][t + 1][1] - in[j][1][t][1]) / (2. * in[j][0][t][0]);

    return r;
}

double deltam(int j, double**** in, int t, struct fit_type fit_info) {


    double r = in[j][1][t][1] / (2. * in[j][4][t][0]);

    return r;
}

int main(int argc, char** argv) {
    error(argc != 7, 1, "nissa_mpcac ",
        "usage:./nissa_mpcac -p path file -bin $bin  jack/boot \n separate path and file please");

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
    mysprintf(option[2], NAMESIZE, "-p"); // -p
    mysprintf(option[3], NAMESIZE, argv[2]); // path
    mysprintf(option[4], NAMESIZE, argv[6]); //resampling
    mysprintf(option[5], NAMESIZE, "no"); // pdf
    mysprintf(option[6], NAMESIZE, argv[3]); // infile

    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");


    FILE* infile = open_file(namefile, "r");
    generic_header head = read_head(infile);


    file_head.l1 = head.L;
    file_head.l0 = head.T;
    file_head.l2 = head.L;
    file_head.l3 = head.L;
    file_head.nk = 2;
    file_head.musea = head.mus[1];
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;file_head.k[1] = 0;
    file_head.k[2] = head.mus[1];
    file_head.k[3] = head.mus[1];

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (int i = 0;i < file_head.nmoms;i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }


    int confs = head.Njack;
    int bin = atoi(argv[5]);
    int Neff = confs / bin;
    int Njack;
    if (strcmp(argv[6], "jack") == 0) {
        Njack = Neff + 1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[6], "boot") == 0) {
        Njack = (Neff * 2 + 1);
        myres = new resampling_boot(Neff * 2);
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[7]= %s is not jack or boot", argv[7]);
    }
    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", option[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", option[3], option[4], option[6]);
    FILE* jack_file = open_file(namefile, "w+");
    write_header_g2(jack_file, head);


    double**** data = calloc_corr(confs, head.ncorr, head.T);


    printf("confs=%d\n", confs);
    printf("ncorr=%d\n", head.ncorr);
    for (int iconf = 0; iconf < confs;iconf++) {
        read_twopt(infile, data[iconf], head);
    }

    // symmetrise_corr(confs, 0, head.T, data);
    // antisymmetrise_corr(confs, 1, head.T, data);
    // antisymmetrise_corr(confs, 2, head.T, data);

    double**** data_bin = binning(confs, head.ncorr, head.T, data, bin);
    double**** conf_jack = myres->create(Neff, head.ncorr, head.T, data_bin);
    free_corr(Neff, head.ncorr, head.T, data_bin);

    free_corr(confs, head.ncorr, head.T, data);



    /////////////////////////////////////////////////////////////////////////////////////////////////////////
           //print all the effective masses correlators
           //set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_meff_correlators", option[3], option[6]);
    FILE* outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_raw_correlators", option[3], option[6]);
    FILE* outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_shifted_correlators", option[3], option[6]);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_log_meff_shifted", option[3], option[6]);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_gamma", option[3], option[6]);
    FILE* out_gamma = open_file(namefile, "w+");

    char  save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");
    FILE* dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < head.ncorr; icorr++) {
        //log effective mass
        double* tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack
            , namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null, fit_info_silent);
        free(tmp_meff_corr);
        //raw correlator
        file_head.l0 = head.T * 2;
        tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity, dev_null, fit_info_silent);
        free(tmp_meff_corr);
        tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity_im, dev_null, fit_info_silent);
        free(tmp_meff_corr);
        file_head.l0 = head.T;
        // shifted correlator
        tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack
            , namefile_plateaux, outfile_log_meff_shifted, icorr, "log_shift", M_eff_log_shift, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
    }
    symmetrise_jackboot(Njack, 0, head.T, conf_jack);
    symmetrise_jackboot(Njack, 1, head.T, conf_jack, -1);
    symmetrise_jackboot(Njack, 2, head.T, conf_jack, -1);

    fit_info_silent.restore_default();
    sprintf(option[1], "%s", save_option);// restore option
    corr_counter = -1;

    double* M_PS = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 0, "M_{PS}", M_eff_T, jack_file);
    check_correlatro_counter(0);

    struct fit_type fit_info;
    struct fit_result  fit_out;
    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 0;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * 0);
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;

    //c++ 1 || r 2
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, mpcac, "mpcac", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);
    fit_info.restore_default();


    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 0;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * 0);
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;

    //c++ 1 || r 2
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, deltam, "deltam", fit_info, jack_file);
    double* true_kappa = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0; j < Njack;j++) {
        true_kappa[j] = head.kappa / (1. + 2. * fit_out.P[0][j] * head.kappa);
    }
    fprintf(outfile, "%.12g   %.12g\n", true_kappa[Njack - 1], myres->comp_error(true_kappa));
    printf("true_kappa=%.12g   %.12g\n", true_kappa[Njack - 1], myres->comp_error(true_kappa));
    free_fit_result(fit_info, fit_out);
    fit_info.restore_default();



}