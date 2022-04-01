#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "global.hpp"
#include "resampling.hpp"
#include "read.hpp"
// #include "m_eff.hpp"
// #include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
#include "functions.hpp"
// #include "correlators_analysis.hpp"
// #include "eigensystem.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
#include "fit_all.hpp"

#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>
#include <vector>


generic_header read_header(FILE* stream) {
    generic_header header;
    int ir = 0;
    ir += fread(&header.T, sizeof(int), 1, stream);
    ir += fread(&header.L, sizeof(int), 1, stream);
    int s;
    ir += fread(&s, sizeof(int), 1, stream);
    header.mus = std::vector<double>(s);
    for (int i = 0; i < s; i++) {
        ir += fread(&header.mus[i], sizeof(double), 1, stream);
    }
    ir += fread(&s, sizeof(int), 1, stream);
    header.thetas = std::vector<double>(s);
    for (int i = 0; i < s; i++) {
        ir += fread(&header.thetas[i], sizeof(double), 1, stream);
    }

    ir += fread(&header.Njack, sizeof(int), 1, stream);
    header.struct_size = ftell(stream);
    return header;


}


double read_single_Nobs(FILE* stream, int header_size, int Njack) {
    int Nobs;
    long int tmp;
    int s = header_size;

    // size_t i = fread(&Njack, sizeof(int), 1, stream);


    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    tmp -= header_size;

    s = Njack;

    Nobs = (tmp) / ((s) * sizeof(double));

    fseek(stream, header_size, SEEK_SET);

    return Nobs;

}

data_single read_single_dataj(FILE* stream) {

    int Njack;
    int Nobs;

    //read_single_Njack_Nobs(stream, header.header_size, Njack, Nobs);
    // data_single dj(Nobs,Njack);
    data_single dj;
    dj.header = read_header(stream);
    dj.Nobs = read_single_Nobs(stream, dj.header.struct_size, dj.header.Njack);
    dj.Njack = dj.header.Njack;
    dj.jack = double_malloc_2(dj.Nobs, dj.Njack);

    //
    size_t i = 0;
    for (int obs = 0; obs < dj.Nobs; obs++) {
        i += fread(dj.jack[obs], sizeof(double), dj.Njack, stream);
    }
    return dj;

}

data_all read_all_the_files(std::vector<std::string> files, const char* resampling) {
    data_all jackall;
    jackall.resampling = resampling;
    //jackall->en = (data_single*)malloc(sizeof(data_single) * files.size());
    jackall.en = new data_single[files.size()];
    jackall.ens = files.size();
    int count = 0;
    for (std::string s : files) {
        FILE* f = open_file(s.c_str(), "r");

        // read_single_dataj(f, params, &(jackall->en[count]));
        jackall.en[count] = read_single_dataj(f);
        jackall.en[count].resampling = resampling;
        count++;
        fclose(f);
    }
    return jackall;

}

int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");
    char namefile[NAMESIZE];


    std::vector<std::string> files;
    mysprintf(namefile, NAMESIZE, "%s/%s_cB.72.64_mu.0.000720", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB.72.96_mu.0.000720", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC.06.80_mu.0.000600", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD.54.96_mu.0.000540", argv[2], argv[1]);
    files.emplace_back(namefile);




    data_all jackall = read_all_the_files(files, argv[1]);
    jackall.init_error();
    std::vector<int> myen(jackall.ens);
    for (int e = 0; e < jackall.ens; e++) {
        myen[e] = e;
        printf("file:  %s\n", files[e].c_str());
        printf("Nobs=%d  Njack=%d    mus=", jackall.en[e].Nobs, jackall.en[e].Njack);
        for (double mu : jackall.en[e].header.mus)
            printf("%g   ", mu);
        printf("\n");
    }
    int Njack = jackall.en[0].Njack;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // fits 
    /////////////////////////////////////////////////////////////////////////////////////////////////
    fit_type fit_info;

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////     amu_SD_l_common    //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Npar = 3;
    fit_info.N = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = Njack;
    fit_info.myen = myen;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    int count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
            }
            count++;
        }
    }
    fit_info.function = rhs_amu_common;
    fit_result amu_SD_l_common = fit_all_data(argv, jackall, lhs_amu_common<25, 26>, fit_info, "amu_SD_l_common");
    fit_info.band_range = { 0,0.0081 };
    print_fit_band(argv, jackall, fit_info, fit_info, "amu_SD_l_common", "afm", amu_SD_l_common, amu_SD_l_common, 0, myen.size() - 1, 0.005);

    fit_info.restore_default();

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////     amu_SD_l_common    a^2+a^4//////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.Npar = 5;
    fit_info.N = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = Njack;
    fit_info.myen = myen;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
            }
            count++;
        }
    }
    fit_info.function = rhs_amu_common_a4;
    fit_result amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu_common<25, 26>, fit_info, "amu_SD_l_common_a4");
    fit_info.band_range = { 0,0.0081 };
    print_fit_band(argv, jackall, fit_info, fit_info, "amu_SD_l_common_a4", "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, myen.size() - 1, 0.001);

    fit_info.restore_default();


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////     amu_SD_s_common    //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Npar = 3;
    fit_info.N = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = Njack;
    fit_info.myen = myen;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
            }
            count++;
        }
    }
    fit_info.function = rhs_amu_common;
    fit_result amu_SD_s_common = fit_all_data(argv, jackall, lhs_amu_common<31, 34>, fit_info, "amu_SD_s_common");
    fit_info.band_range = { 0,0.0081 };
    print_fit_band(argv, jackall, fit_info, fit_info, "amu_SD_s_common", "afm", amu_SD_s_common, amu_SD_s_common, 0, myen.size() - 1, 0.005);

    fit_info.restore_default();
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////     amu_SD_s_common_a4    //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Npar = 5;
    fit_info.N = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = Njack;
    fit_info.myen = myen;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
            }
            count++;
        }
    }
    fit_info.function = rhs_amu_common_a4;
    fit_result amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu_common<31, 34>, fit_info, "amu_SD_s_common_a4");
    fit_info.band_range = { 0,0.0081 };
    print_fit_band(argv, jackall, fit_info, fit_info, "amu_SD_s_common_a4", "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, myen.size() - 1, 0.005);

    fit_info.restore_default();

}