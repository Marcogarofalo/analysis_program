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

void sum_lsc(data_all in, const char* outpath, const char* filename) {
    int N = in.Nfits;
    int Njack = in.fits[0].Njack;
    char name[NAMESIZE];
    mysprintf(name, NAMESIZE, "%s/%s", outpath, filename);
    mysprintf(name, NAMESIZE, "%s/%s", outpath, filename);
    FILE* f = open_file(name, "w+");
    printf("writing: %s\n", name);
    double* ave = (double*)calloc(Njack, sizeof(double*));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < Njack; j++) {
            ave[j] += in.fits[i].P[0][j];
        }
        fprintf(f, "%s %g  %g\n", in.fits[i].name, in.fits[i].P[0][Njack - 1], error_jackboot(in.resampling.c_str(), Njack, in.fits[i].P[0]));
    }
    fprintf(f, "total: %g  %g\n", ave[Njack - 1], error_jackboot(in.resampling.c_str(), Njack, ave));

    double** y = double_malloc_2(N, Njack);
    for (int j = 0; j < Njack; j++) {
        for (int i = 0; i < N; i++) {
            y[i][j] = in.fits[i].P[0][j];
        }
    }
    double** cov = covariance(in.resampling.c_str(), N, Njack, y);
    fprintf(f, "\n#covariance:\n");
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            fprintf(f, "%g\t", cov[i][k]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n#correlation:\n");
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            fprintf(f, "%g\t", cov[i][k] / sqrt(cov[i][i] * cov[k][k]));
        }
        fprintf(f, "\n");
    }
    fclose(f);

}

void compute_syst_eq28(data_all in, const char* outpath, const char* filename) {
    int N = in.Nfits;
    int Njack = in.fits[0].Njack;
    char name[NAMESIZE];
    mysprintf(name, NAMESIZE, "%s/%s", outpath, filename);
    FILE* f = open_file(name, "w+");
    printf("writing: %s\n", name);
    std::vector<double> aves(N);
    std::vector<double> errors(N);
    for (int i = 0; i < N; i++) {
        aves[i] = in.fits[i].P[0][Njack - 1];
        errors[i] = error_jackboot(in.resampling.c_str(), Njack, in.fits[i].P[0]);
    }

    double ave = 0, err = 0;
    for (int i = 0; i < N; i++) {
        ave += aves[i];
        err += pow(errors[i], 2);
    }
    ave /= (double)N;

    for (int i = 0; i < N; i++) {
        err += pow(ave - aves[i], 2);
    }
    err = sqrt(err / (double)N);

    for (int i = 0; i < N; i++) {
        fprintf(f, "%s    %g     %g   %g   %g  %g\n", in.fits[i].name, aves[i], errors[i], ave, err, in.fits[i].chi2[Njack - 1]);
    }
    printf("systematics  %s: N=%d\n", filename, N);
    printf("mean(eq28)= %g  %g \n", ave, err);
    fclose(f);
}


int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");
    char namefile[NAMESIZE];
    char namefit[NAMESIZE];


    std::vector<std::string> files;
    mysprintf(namefile, NAMESIZE, "%s/%s_cB.72.64_mu.0.000720", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB.72.96_mu.0.000720", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC.06.80_mu.0.000600", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD.54.96_mu.0.000540", argv[2], argv[1]);
    files.emplace_back(namefile);


    std::vector<int> myen(files.size());
    for (int e = 0; e < files.size(); e++) {
        myen[e] = e;
        printf("file:  %s\n", files[e].c_str());
        //printf("Nobs=%d  Njack=%d    mus=", jackall.en[e].Nobs, jackall.en[e].Njack);
        //for (double mu : jackall.en[e].header.mus)
        //    printf("%g   ", mu);
        //printf("\n");
    }

    mysprintf(namefile, NAMESIZE, "%s/%s_cA.53.24_mu.0.005300", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cA.40.24_mu.0.004000", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cA.30.32_mu.0.003000", argv[2], argv[1]);
    files.emplace_back(namefile);

    data_all jackall = read_all_the_files(files, argv[1]);
    jackall.create_generalised_resampling();

    std::vector<int> myen_full(jackall.ens);
    for (int e = 0; e < jackall.ens; e++) {
        myen_full[e] = e;
        printf("file:  %s\n", files[e].c_str());
        printf("Nobs=%d  Njack=%d    mus=", jackall.en[e].Nobs, jackall.en[e].Njack);
        for (double mu : jackall.en[e].header.mus)
            printf("%g   ", mu);
        printf("\n");
    }
    int Njack = jackall.en[0].Njack;
    std::vector<int> myen_charm = { 0, 2, 3, 4, 5, 6 };



    /////////////////////////////////////////////////////////////////////////////////////////////////
    // fits 
    /////////////////////////////////////////////////////////////////////////////////////////////////
    fit_type fit_info;
    data_all  syst_amu_SD_l;
    syst_amu_SD_l.resampling = argv[1];
    data_all  sum_amu_SD;
    sum_amu_SD.resampling = argv[1];
    data_all  sum_amu_W;
    sum_amu_W.resampling = argv[1];

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
    print_fit_band(argv, jackall, fit_info, fit_info, "amu_SD_l_common", "afm", amu_SD_l_common, amu_SD_l_common, 0, myen.size() - 1, 0.001);


    fit_info.restore_default();


    std::vector<std::string>   integrations = { "reinman", "simpson" };
    for (auto integration : integrations) {
        int id0, id1;
        if (integration == "reinman") { id0 = 25; id1 = 26; }
        if (integration == "simpson") { id0 = 27; id1 = 28; }

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
        fit_info.corr_id = { id0,id1 };
        fit_info.function = rhs_amu_common_a4;
        mysprintf(namefit, NAMESIZE, "amu_SD_l_%s_a4", integration.c_str());
        fit_result amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
        fit_info.band_range = { 0,0.0081 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, myen.size() - 1, 0.0005);
        syst_amu_SD_l.add_fit(amu_SD_l_common_a4);
        if (integration == "reinman") sum_amu_SD.add_fit(amu_SD_l_common_a4);

        free_fit_result(fit_info, amu_SD_l_common_a4);
        fit_info.restore_default();

        ///////////////////////////////////////////////////////////////////////////////////////////////////
        printf("\n/////////////////////////////////     amu_SD_l_common    a^4_eq //////////////////\n");
        //////////////////////////////////////////////////////////////////////////////////////////////////

        fit_info.Npar = 4;
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
        fit_info.corr_id = { id0,id1 };
        fit_info.function = rhs_amu_common_a4_n0;
        mysprintf(namefit, NAMESIZE, "amu_SD_l_%s_a4_eq", integration.c_str());
        amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
        fit_info.band_range = { 0,0.0081 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, myen.size() - 1, 0.0005);
        syst_amu_SD_l.add_fit(amu_SD_l_common_a4);
        free_fit_result(fit_info, amu_SD_l_common_a4);
        fit_info.restore_default();
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        printf("\n/////////////////////////////////     amu_SD_l_common    a^4_op //////////////////\n");
        //////////////////////////////////////////////////////////////////////////////////////////////////

        fit_info.Npar = 4;
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
        fit_info.corr_id = { id0,id1 };
        fit_info.function = rhs_amu_common_a4_n1;
        mysprintf(namefit, NAMESIZE, "amu_SD_l_%s_a4_op", integration.c_str());
        amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
        fit_info.band_range = { 0,0.0081 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, myen.size() - 1, 0.0005);
        syst_amu_SD_l.add_fit(amu_SD_l_common_a4);
        free_fit_result(fit_info, amu_SD_l_common_a4);
        fit_info.restore_default();


    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////   Systematics  amu_SD_l   //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////
    compute_syst_eq28(syst_amu_SD_l, argv[3], "Systematics_amu_SD_l.txt");

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    data_all  syst_amu_SD_s;
    syst_amu_SD_s.resampling = argv[1];


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
    print_fit_band(argv, jackall, fit_info, fit_info, "amu_SD_s_common", "afm", amu_SD_s_common, amu_SD_s_common, 0, myen.size() - 1, 0.001);

    fit_info.restore_default();

    std::vector<std::string>   interpolations = { "eta", "phi" };
    integrations = { "reinman", "simpson" };
    for (auto interpolation : interpolations) {
        for (auto integration : integrations) {
            int id0, id1;

            if (interpolation == "eta" && integration == "reinman") { id0 = 31; id1 = 34; }
            if (interpolation == "eta" && integration == "simpson") { id0 = 37; id1 = 40; }
            if (interpolation == "phi" && integration == "reinman") { id0 = 59; id1 = 60; }
            if (interpolation == "phi" && integration == "simpson") { id0 = 61; id1 = 62; }

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
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_SD_s_%s_%s_a4", interpolation.c_str(), integration.c_str());
            fit_result amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.0081 };
            syst_amu_SD_s.add_fit(amu_SD_s_common_a4);
            if (interpolation == "eta" && integration == "reinman")  sum_amu_SD.add_fit(amu_SD_s_common_a4);

            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, myen.size() - 1, 0.001);

            fit_info.restore_default();
            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_SD_s_common_a4_eq    //////////////////\n");
            //////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Npar = 4;
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
            fit_info.function = rhs_amu_common_a4_n0;
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_SD_s_%s_%s_a4_eq", interpolation.c_str(), integration.c_str());
            amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.0081 };
            syst_amu_SD_s.add_fit(amu_SD_s_common_a4);
            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, myen.size() - 1, 0.001);

            fit_info.restore_default();

            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_SD_s_common_a4_op    //////////////////\n");
            //////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Npar = 4;
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
            fit_info.function = rhs_amu_common_a4_n1;
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_SD_s_%s_%s_a4_op", interpolation.c_str(), integration.c_str());
            amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.0081 };
            syst_amu_SD_s.add_fit(amu_SD_s_common_a4);

            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, myen.size() - 1, 0.001);

            fit_info.restore_default();
            // ///////////////////////////////////////////////////////////////////////////////////////////////////
            // printf("\n/////////////////////////////////     amu_SD_s_common_log_eq    //////////////////\n");
            // //////////////////////////////////////////////////////////////////////////////////////////////////
            // fit_info.Npar = 4;
            // fit_info.N = 2;
            // fit_info.Nvar = 1;
            // fit_info.Njack = Njack;
            // fit_info.myen = myen;
            // fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
            // count = 0;
            // for (int n = 0;n < fit_info.N;n++) {
            //     for (int e = 0;e < fit_info.myen.size();e++) {
            //         for (int j = 0;j < Njack;j++) {
            //             fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
            //         }
            //         count++;
            //     }
            // }
            // fit_info.function = rhs_amu_common_log_a4_n0;
            // fit_info.corr_id = { id0,id1 };
            // mysprintf(namefit, NAMESIZE, "amu_SD_s_%s_%s_log_a4_eq", interpolation.c_str(), integration.c_str());
            // amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            // fit_info.band_range = { 0,0.0081 };
            // syst_amu_SD_s.add_fit(amu_SD_s_common_a4);
            // print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, myen.size() - 1, 0.001);

            // fit_info.restore_default();

            // ///////////////////////////////////////////////////////////////////////////////////////////////////
            // printf("\n/////////////////////////////////     amu_SD_s_common_log_op    //////////////////\n");
            // //////////////////////////////////////////////////////////////////////////////////////////////////
            // fit_info.Npar = 4;
            // fit_info.N = 2;
            // fit_info.Nvar = 1;
            // fit_info.Njack = Njack;
            // fit_info.myen = myen;
            // fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
            // count = 0;
            // for (int n = 0;n < fit_info.N;n++) {
            //     for (int e = 0;e < fit_info.myen.size();e++) {
            //         for (int j = 0;j < Njack;j++) {
            //             fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
            //         }
            //         count++;
            //     }
            // }
            // fit_info.function = rhs_amu_common_log_a4_n1;
            // fit_info.corr_id = { id0,id1 };
            // mysprintf(namefit, NAMESIZE, "amu_SD_s_%s_%s_log_a4_op", interpolation.c_str(), integration.c_str());
            // amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            // fit_info.band_range = { 0,0.0081 };
            // syst_amu_SD_s.add_fit(amu_SD_s_common_a4);
            // print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, myen.size() - 1, 0.001);

            // fit_info.restore_default();
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    compute_syst_eq28(syst_amu_SD_s, argv[3], "Systematics_amu_SD_s.txt");
    sum_lsc(sum_amu_SD, argv[3], "sum_amu_SD_ls.txt");
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    data_all  syst_amu_W_l;
    syst_amu_W_l.resampling = argv[1];
    for (auto integration : integrations) {
        int id0, id1;
        if (integration == "reinman") { id0 = 42; id1 = 43; }
        if (integration == "simpson") { id0 = 44; id1 = 45; }

        ///////////////////////////////////////////////////////////////////////////////////////////////////
        printf("\n/////////////////////////////////     amu_W_l_common    a^2+a^4//////////////////\n");
        //////////////////////////////////////////////////////////////////////////////////////////////////
        constexpr double Mpi_MeV = 135;
        constexpr double Mpi_MeV_err = 0.2;

        double* jack_Mpi_MeV_exp = fake_sampling(argv[1], Mpi_MeV, Mpi_MeV_err, Njack, 1003);

        fit_info.Npar = 6;
        fit_info.N = 2;
        fit_info.Nvar = 4;
        fit_info.Njack = Njack;
        fit_info.myen = myen;
        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
        count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < fit_info.myen.size();e++) {
                for (int j = 0;j < Njack;j++) {
                    fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
                    fit_info.x[1][count][j] = jackall.en[e].jack[58][j];  // Delta_FV_GS
                    fit_info.x[2][count][j] = jackall.en[e].jack[1][j];  //Mpi
                    fit_info.x[3][count][j] = jack_Mpi_MeV_exp[j];
                }
                count++;
            }
        }

        fit_info.corr_id = { id0, id1 };
        fit_info.function = rhs_amu_common_a2_FVE;
        mysprintf(namefit, NAMESIZE, "amu_W_l_%s_a2", integration.c_str());

        fit_result amu_W_l_common_a2 = fit_all_data(argv, jackall, lhs_amu_common_GS, fit_info, namefit);
        fit_info.band_range = { 0,0.0081 };
        std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0 };
        // print_fit_band(argv, jackall, fit_info, fit_info, "amu_W_l_common_a2", "afm", amu_W_l_common_a2, amu_W_l_common_a2, 0, myen.size() - 1, 0.0005, xcont);
        print_fit_band_amu_W_l(argv, jackall, fit_info, fit_info, namefit, "afm", amu_W_l_common_a2, amu_W_l_common_a2, 0, myen.size() - 1, 0.0005, xcont, lhs_amu_common_GS);
        syst_amu_W_l.add_fit(amu_W_l_common_a2);
        if (integration == "reinman") sum_amu_W.add_fit(amu_W_l_common_a2);

        fit_info.restore_default();
    }
    compute_syst_eq28(syst_amu_W_l, argv[3], "Systematics_amu_W_l.txt");


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    data_all  syst_amu_W_s;
    syst_amu_W_s.resampling = argv[1];


    // ///////////////////////////////////////////////////////////////////////////////////////////////////
    // printf("\n/////////////////////////////////     amu_W_s_common_a4    //////////////////\n");
    // //////////////////////////////////////////////////////////////////////////////////////////////////
    // fit_info.Npar = 5;
    // fit_info.N = 2;
    // fit_info.Nvar = 1;
    // fit_info.Njack = Njack;
    // fit_info.myen = myen;
    // fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    // count = 0;
    // for (int n = 0;n < fit_info.N;n++) {
    //     for (int e = 0;e < fit_info.myen.size();e++) {
    //         for (int j = 0;j < Njack;j++) {
    //             fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
    //         }
    //         count++;
    //     }
    // }
    // fit_info.function = rhs_amu_common_a4;
    // fit_result amu_W_s_common_a4 = fit_all_data(argv, jackall, lhs_amu_common<48, 51>, fit_info, "amu_W_s_common_a4");
    // syst_amu_W_s.add_fit(amu_W_s_common_a4);
    // fit_info.band_range = { 0,0.0081 };
    // print_fit_band(argv, jackall, fit_info, fit_info, "amu_W_s_common_a4", "afm", amu_W_s_common_a4, amu_W_s_common_a4, 0, myen.size() - 1, 0.001);

    // fit_info.restore_default();

    // ///////////////////////////////////////////////////////////////////////////////////////////////////
    // printf("\n/////////////////////////////////     amu_W_s_common_a4    //////////////////\n");
    // //////////////////////////////////////////////////////////////////////////////////////////////////
    // fit_info.Npar = 4;
    // fit_info.N = 2;
    // fit_info.Nvar = 1;
    // fit_info.Njack = Njack;
    // fit_info.myen = myen;
    // fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    // count = 0;
    // for (int n = 0;n < fit_info.N;n++) {
    //     for (int e = 0;e < fit_info.myen.size();e++) {
    //         for (int j = 0;j < Njack;j++) {
    //             fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
    //         }
    //         count++;
    //     }
    // }
    // fit_info.function = rhs_amu_common_a4_n0;
    // fit_result amu_W_s_common_a4_eq = fit_all_data(argv, jackall, lhs_amu_common<48, 51>, fit_info, "amu_W_s_common_a4_eq");
    // fit_info.band_range = { 0,0.0081 };
    // syst_amu_W_s.add_fit(amu_W_s_common_a4_eq);
    // print_fit_band(argv, jackall, fit_info, fit_info, "amu_W_s_common_a4_eq", "afm", amu_W_s_common_a4_eq, amu_W_s_common_a4_eq, 0, myen.size() - 1, 0.001);

    // fit_info.restore_default();
    // ///////////////////////////////////////////////////////////////////////////////////////////////////
    // printf("\n/////////////////////////////////     amu_W_s_common_a4    //////////////////\n");
    // //////////////////////////////////////////////////////////////////////////////////////////////////
    // fit_info.Npar = 4;
    // fit_info.N = 2;
    // fit_info.Nvar = 1;
    // fit_info.Njack = Njack;
    // fit_info.myen = myen;
    // fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    // count = 0;
    // for (int n = 0;n < fit_info.N;n++) {
    //     for (int e = 0;e < fit_info.myen.size();e++) {
    //         for (int j = 0;j < Njack;j++) {
    //             fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
    //         }
    //         count++;
    //     }
    // }
    // fit_info.function = rhs_amu_common_a4_n1;
    // fit_result amu_W_s_common_a4_op = fit_all_data(argv, jackall, lhs_amu_common<48, 51>, fit_info, "amu_W_s_common_a4_op");
    // fit_info.band_range = { 0,0.0081 };
    // syst_amu_W_s.add_fit(amu_W_s_common_a4_op);
    // print_fit_band(argv, jackall, fit_info, fit_info, "amu_W_s_common_a4_op", "afm", amu_W_s_common_a4_op, amu_W_s_common_a4_op, 0, myen.size() - 1, 0.001);

    // fit_info.restore_default();

    interpolations = { "eta", "phi" };
    integrations = { "reinman", "simpson" };
    for (auto interpolation : interpolations) {
        for (auto integration : integrations) {
            int id0, id1;

            if (interpolation == "eta" && integration == "reinman") { id0 = 48; id1 = 51; }
            if (interpolation == "eta" && integration == "simpson") { id0 = 54; id1 = 57; }
            if (interpolation == "phi" && integration == "reinman") { id0 = 63; id1 = 64; }
            if (interpolation == "phi" && integration == "simpson") { id0 = 65; id1 = 66; }

            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_W_s_common_a4    //////////////////\n");
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
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_W_s_%s_%s_a4", interpolation.c_str(), integration.c_str());
            fit_result amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.0081 };
            syst_amu_W_s.add_fit(amu_SD_s_common_a4);
            if (interpolation == "eta" && integration == "reinman") sum_amu_W.add_fit(amu_SD_s_common_a4);

            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, myen.size() - 1, 0.001);

            fit_info.restore_default();
            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_W_s_common_a4_eq    //////////////////\n");
            //////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Npar = 4;
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
            fit_info.function = rhs_amu_common_a4_n0;
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_W_s_%s_%s_a4_eq", interpolation.c_str(), integration.c_str());
            amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.0081 };
            syst_amu_W_s.add_fit(amu_SD_s_common_a4);

            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, myen.size() - 1, 0.001);

            fit_info.restore_default();

            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_W_s_common_a4_op    //////////////////\n");
            //////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Npar = 4;
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
            fit_info.function = rhs_amu_common_a4_n1;
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_W_s_%s_%s_a4_op", interpolation.c_str(), integration.c_str());
            amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.0081 };
            syst_amu_W_s.add_fit(amu_SD_s_common_a4);
            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, myen.size() - 1, 0.001);

            fit_info.restore_default();
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    compute_syst_eq28(syst_amu_W_s, argv[3], "Systematics_amu_W_s.txt");
    sum_lsc(sum_amu_W, argv[3], "sum_amu_W_ls.txt");
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////     amu_SD charm    //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    data_all  syst_amu_SD_c;
    syst_amu_SD_c.resampling = argv[1];




    interpolations = { "etac", "Jpsi" };
    integrations = { "reinman", "simpson" };
    for (auto interpolation : interpolations) {
        for (auto integration : integrations) {
            int id0, id1;

            if (integration == "reinman" && interpolation == "etac") { id0 = 76; id1 = 86; }
            if (integration == "simpson" && interpolation == "etac") { id0 = 81; id1 = 91; }
            if (integration == "reinman" && interpolation == "Jpsi") { id0 = 77; id1 = 87; }
            if (integration == "simpson" && interpolation == "Jpsi") { id0 = 82; id1 = 92; }

            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_SD_c_common_a4    //////////////////\n");
            //////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Npar = 5;
            fit_info.N = 2;
            fit_info.Nvar = 1;
            fit_info.Njack = Njack;
            fit_info.myen = myen_charm;
            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
            count = 0;
            for (int n = 0;n < fit_info.N;n++) {
                for (int e : fit_info.myen) {
                    for (int j = 0;j < Njack;j++) {
                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
                    }
                    count++;
                }
            }
            fit_info.function = rhs_amu_common_a4;
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_SD_c_%s_%s_a4", interpolation.c_str(), integration.c_str());
            fit_result amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.009 };
            syst_amu_SD_c.add_fit(amu_SD_s_common_a4);
            if (integration == "reinman" && interpolation == "etac")  sum_amu_SD.add_fit(amu_SD_s_common_a4);

            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, fit_info.myen.size() - 1, 0.001);

            free_fit_result(fit_info, amu_SD_s_common_a4);
            fit_info.restore_default();
            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_SD_c_common_a4_eq    //////////////////\n");
            //////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Npar = 4;
            fit_info.N = 2;
            fit_info.Nvar = 1;
            fit_info.Njack = Njack;
            fit_info.myen = myen_charm;
            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
            count = 0;
            for (int n = 0;n < fit_info.N;n++) {
                for (int e : fit_info.myen) {
                    for (int j = 0;j < Njack;j++) {
                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
                    }
                    count++;
                }
            }
            fit_info.function = rhs_amu_common_a4_n0;
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_SD_c_%s_%s_a4_eq", interpolation.c_str(), integration.c_str());
            amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.009 };
            syst_amu_SD_c.add_fit(amu_SD_s_common_a4);
            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, fit_info.myen.size() - 1, 0.001);

            free_fit_result(fit_info, amu_SD_s_common_a4);
            fit_info.restore_default();

            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_SD_c_common_a4_op    //////////////////\n");
            //////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Npar = 4;
            fit_info.N = 2;
            fit_info.Nvar = 1;
            fit_info.Njack = Njack;
            fit_info.myen = myen_charm;
            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
            count = 0;
            for (int n = 0;n < fit_info.N;n++) {
                for (int e : fit_info.myen) {
                    for (int j = 0;j < Njack;j++) {
                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
                    }
                    count++;
                }
            }
            fit_info.function = rhs_amu_common_a4_n1;
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_SD_c_%s_%s_a4_op", interpolation.c_str(), integration.c_str());
            amu_SD_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.009 };
            syst_amu_SD_c.add_fit(amu_SD_s_common_a4);
            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_s_common_a4, amu_SD_s_common_a4, 0, fit_info.myen.size() - 1, 0.001);

            free_fit_result(fit_info, amu_SD_s_common_a4);
            fit_info.restore_default();

        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    compute_syst_eq28(syst_amu_SD_c, argv[3], "Systematics_amu_SD_c.txt");
    sum_lsc(sum_amu_SD, argv[3], "sum_amu_SD_lsc.txt");
    ///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////     amu_W charm    //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    data_all  syst_amu_W_c;
    syst_amu_W_c.resampling = argv[1];




    interpolations = { "etac", "Jpsi" };
    integrations = { "reinman", "simpson" };
    for (auto interpolation : interpolations) {
        for (auto integration : integrations) {
            int id0, id1;

            if (integration == "reinman" && interpolation == "etac") { id0 = 96;  id1 = 106; }
            if (integration == "simpson" && interpolation == "etac") { id0 = 101; id1 = 111; }
            if (integration == "reinman" && interpolation == "Jpsi") { id0 = 97;  id1 = 107; }
            if (integration == "simpson" && interpolation == "Jpsi") { id0 = 102; id1 = 112; }

            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_W_c_common_a4    //////////////////\n");
            //////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Npar = 5;
            fit_info.N = 2;
            fit_info.Nvar = 1;
            fit_info.Njack = Njack;
            fit_info.myen = myen_charm;
            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
            count = 0;
            for (int n = 0;n < fit_info.N;n++) {
                for (int e : fit_info.myen) {
                    for (int j = 0;j < Njack;j++) {
                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
                    }
                    count++;
                }
            }
            fit_info.function = rhs_amu_common_a4;
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_W_c_%s_%s_a4", interpolation.c_str(), integration.c_str());
            fit_result amu_W_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.009 };
            syst_amu_W_c.add_fit(amu_W_s_common_a4);
            if (integration == "reinman" && interpolation == "etac") sum_amu_W.add_fit(amu_W_s_common_a4);

            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_W_s_common_a4, amu_W_s_common_a4, 0, fit_info.myen.size() - 1, 0.001);

            free_fit_result(fit_info, amu_W_s_common_a4);
            fit_info.restore_default();
            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_W_c_common_a4_eq    //////////////////\n");
            //////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Npar = 4;
            fit_info.N = 2;
            fit_info.Nvar = 1;
            fit_info.Njack = Njack;
            fit_info.myen = myen_charm;
            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
            count = 0;
            for (int n = 0;n < fit_info.N;n++) {
                for (int e : fit_info.myen) {
                    for (int j = 0;j < Njack;j++) {
                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
                    }
                    count++;
                }
            }
            fit_info.function = rhs_amu_common_a4_n0;
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_W_c_%s_%s_a4_eq", interpolation.c_str(), integration.c_str());
            amu_W_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.009 };
            syst_amu_W_c.add_fit(amu_W_s_common_a4);

            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_W_s_common_a4, amu_W_s_common_a4, 0, fit_info.myen.size() - 1, 0.001);

            free_fit_result(fit_info, amu_W_s_common_a4);
            fit_info.restore_default();

            ///////////////////////////////////////////////////////////////////////////////////////////////////
            printf("\n/////////////////////////////////     amu_W_c_common_a4_op    //////////////////\n");
            //////////////////////////////////////////////////////////////////////////////////////////////////
            fit_info.Npar = 4;
            fit_info.N = 2;
            fit_info.Nvar = 1;
            fit_info.Njack = Njack;
            fit_info.myen = myen_charm;
            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
            count = 0;
            for (int n = 0;n < fit_info.N;n++) {
                for (int e : fit_info.myen) {
                    for (int j = 0;j < Njack;j++) {
                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2);
                    }
                    count++;
                }
            }
            fit_info.function = rhs_amu_common_a4_n1;
            fit_info.corr_id = { id0,id1 };
            mysprintf(namefit, NAMESIZE, "amu_W_c_%s_%s_a4_op", interpolation.c_str(), integration.c_str());
            amu_W_s_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
            fit_info.band_range = { 0,0.009 };
            syst_amu_W_c.add_fit(amu_W_s_common_a4);
            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_W_s_common_a4, amu_W_s_common_a4, 0, fit_info.myen.size() - 1, 0.001);

            free_fit_result(fit_info, amu_W_s_common_a4);
            fit_info.restore_default();

        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    compute_syst_eq28(syst_amu_W_c, argv[3], "Systematics_amu_W_c.txt");
    sum_lsc(sum_amu_W, argv[3], "sum_amu_W_lsc.txt");
    ///////////////////////////////////////////////////////////////////////////////////////////////////



}