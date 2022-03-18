#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"
#include "resampling.hpp"
#include "read.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
#include "functions.hpp"
#include "correlators_analysis.hpp"
#include "eigensystem.hpp"
#include "non_linear_fit.hpp"
// #include "lhs_functions.hpp"

#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>
using namespace std;


struct  kinematic kinematic_2pt;


void get_kinematic(int ik2, int r2, int ik1, int r1, int imom2, int imom1) {
    kinematic_2pt.ik2 = ik2;
    kinematic_2pt.ik1 = ik1;
    kinematic_2pt.k2 = file_head.k[ik2 + file_head.nk];
    kinematic_2pt.k1 = file_head.k[ik1 + file_head.nk];
    kinematic_2pt.r2 = 1;
    kinematic_2pt.r1 = 1;
    kinematic_2pt.mom2 = -file_head.mom[imom2][1];
    kinematic_2pt.mom1 = file_head.mom[imom1][1];

    kinematic_2pt.mom02 = file_head.mom[imom2][0];
    kinematic_2pt.mom01 = file_head.mom[imom1][0];

    kinematic_2pt.line = 1;

}

int read_nconfs(const char stream[NAMESIZE]) {

    long int tmp;
    int s = file_head.l1 + 1 + 3;
    int count = 0;
    string line;
    ifstream file(stream);
    while (getline(file, line))
        count++;
    std::cout << "lines=" << count << std::endl;
    error(count % s != 0, 1, "read_nconfs lines and T do not match", "lines=%i   T/2=%i", count, file_head.l1);
    int c = count / s;
    std::cout << "confs=" << c << std::endl;

    return c;

}
static void  write_file_head(FILE* stream) {
    int i, dsize;
    double* dstd;

    fwrite(&file_head.twist, sizeof(int), 1, stream);
    fwrite(&file_head.nf, sizeof(int), 1, stream);
    fwrite(&file_head.nsrc, sizeof(int), 1, stream);
    fwrite(&file_head.l0, sizeof(int), 1, stream);
    fwrite(&file_head.l1, sizeof(int), 1, stream);
    fwrite(&file_head.l2, sizeof(int), 1, stream);
    fwrite(&file_head.l3, sizeof(int), 1, stream);
    fwrite(&file_head.nk, sizeof(int), 1, stream);
    fwrite(&file_head.nmoms, sizeof(int), 1, stream);

    fwrite(&file_head.beta, sizeof(double), 1, stream);
    fwrite(&file_head.ksea, sizeof(double), 1, stream);
    fwrite(&file_head.musea, sizeof(double), 1, stream);
    fwrite(&file_head.csw, sizeof(double), 1, stream);

    fwrite(file_head.k, sizeof(double), 2 * file_head.nk, stream);

    for (i = 0;i < file_head.nmoms;i++)
        fwrite(file_head.mom[i], sizeof(double), 4, stream);
}



double matrix_element_GEVP(int t, double** cor, double mass) {
    double me;

    me = cor[t][0] / sqrt(exp(-mass * t) + exp(-(file_head.l0 - t) * mass));
    me *= 2 * mass;

    return  me;
}


void read_twopt(const char namefile[NAMESIZE], int confs, int T, double**** to_write, int id) {
    char tmp[NAMESIZE];
    int sd = 0;
    std::fstream newfile;
    newfile.open(namefile, std::ios::in);
    if (newfile.is_open()) { // checking whether the file is open
        std::string tp;
        for (int i = 0;i < confs;i++) {
            getline(newfile, tp);
            // cout<< tp<< endl;
            getline(newfile, tp);
            // cout<< tp<< endl;
            getline(newfile, tp);
            // cout<< tp<< endl;
            for (int t = 0;t < T / 2 + 1;t++) {
                getline(newfile, tp);
                to_write[i][id][t][0] = stod(tp);
                // printf("%.15f\n", to_write[i][id][t][0]);
            }

        }
    }
    else {
        error(0 == 0, 1, "correlators_analysis.cpp line_read_plateaux",
            "unable to open %s", namefile);
    }

}

// void read_twopt(FILE* stream, int iconf, double*** to_write, cluster::IO_params params, int index) {

//     int tmp = params.data.header_size;// 
//     tmp += sizeof(double) * iconf * params.data.size + sizeof(int) * (iconf + 1);


//     double* obs = (double*)malloc(params.data.size * sizeof(double));

//     fseek(stream, tmp, SEEK_SET);
//     size_t i = fread(obs, sizeof(double), params.data.size, stream);

//     for (int t = 0;t < params.data.L[0];t++) {
//         size_t  id = index + t * params.data.ncorr;
//         (*to_write)[t][0] = obs[id];

//     }
//     free(obs);


// }



void setup_single_file_jack(char* save_name, char** argv, const char* name, int Njack) {
    FILE* f;
    mysprintf(save_name, NAMESIZE, "/dev/null");
    f = fopen(save_name, "w+");
    error(f == NULL, 1, "setup_file_jack ",
        "Unable to open output file /dev/null");
    write_file_head(f);
    fwrite(&Njack, sizeof(int), 1, f);
    fclose(f);
}

void setup_single_file_jack_ASCI(char* save_name, char** argv, const char* name, int Njack) {
    FILE* f;
    mysprintf(save_name, NAMESIZE, "/dev/null");
    f = fopen(save_name, "w+");
    error(f == NULL, 1, "setup_file_jack ",
        "Unable to open output file /dev/null");
    fclose(f);
}

void setup_file_jack(char** argv, int Njack) {
    if (strcmp(argv[4], "jack") == 0) {
        setup_single_file_jack(file_jack.M_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.Zf_PS, argv, "/dev/null", Njack);

        setup_single_file_jack(file_jack.M_PS_GEVP, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS_ls_ss, argv, "/dev/null", Njack);

    }

    if (strcmp(argv[4], "boot") == 0) {

        setup_single_file_jack(file_jack.M_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.Zf_PS, argv, "/dev/null", Njack);

        setup_single_file_jack(file_jack.M_PS_GEVP, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS_ls_ss, argv, "/dev/null", Njack);
    }
}



void write_header_g2(FILE* stream) {

    fwrite(&file_head.l0, sizeof(int), 1, stream);
    fwrite(&file_head.l1, sizeof(int), 1, stream);
    fwrite(&file_head.l2, sizeof(int), 1, stream);
    fwrite(&file_head.l3, sizeof(int), 1, stream);
    fwrite(&file_head.musea, sizeof(double), 1, stream);


}



int main(int argc, char** argv) {
    int size;
    int i, j, t;

    int* iconf, confs;
    double**** data, **** data_bin, ** out, ** tmp;
    char c;
    clock_t t1, t2;
    double* in;

    double**** M, **** vec, **** projected_O;
    double**** lambda, **** lambda0;

    double* fit, *** y, * x, * m, * me;


    double**** conf_jack, ** r, ** mt, ** met;
    int Ncorr = 1;
    int t0 = 2;

    FILE* f_ll = NULL, * f_sl = NULL, * f_ls = NULL, * f_ss = NULL;

    FILE* plateaux_masses = NULL, * plateaux_masses_GEVP = NULL;
    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");
    FILE* plateaux_f = NULL;
    char namefile[NAMESIZE];
    srand(1);



    error(argc != 12, 1, "main ",
        "usage:./g-2  blind/see/read_plateaux -p path basename -bin $bin"
        "  -mu mu -L L jack/boot ");
    error(strcmp(argv[1], "blind") != 0 && strcmp(argv[1], "see") != 0 && strcmp(argv[1], "read_plateaux") != 0, 1, "main ",
        "argv[1] only options:  blind/see/read_plateaux ");

    // cluster::IO_params params;
    // mysprintf(namefile, NAMESIZE, "%s/%s", argv[3], argv[4]);
    // FILE* infile = open_file(namefile, "r+");
    // read_header_phi4(infile, params);




    error(strcmp(argv[5], "-bin") != 0, 1, "main", "argv[4] must be: -bin");
    error(strcmp(argv[7], "-mu") != 0, 1, "main", "argv[4] must be: -mu");
    error(strcmp(argv[9], "-L") != 0, 1, "main", "argv[4] must be: -L");
    error(strcmp(argv[11], "jack") != 0 && strcmp(argv[7], "boot") != 0, 1, "main",
        "argv[6] only options: jack/boot");


    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, argv[1]); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p"); // -p
    mysprintf(option[3], NAMESIZE, argv[3]); // path
    mysprintf(option[4], NAMESIZE, argv[11]); //resampling
    mysprintf(option[5], NAMESIZE, "no"); // pdf

    file_head.l1 = atoi(argv[10]);
    file_head.l0 = file_head.l1 * 2;
    file_head.l2 = file_head.l1;
    file_head.l3 = file_head.l1;

    double mu = atof(argv[8]);


    mysprintf(namefile, NAMESIZE, "%s_mu.%f", argv[4], mu);


    mysprintf(option[6], NAMESIZE, namefile); // basename

    printf("resampling %s\n", option[4]);
    int T = file_head.l0;

    file_head.nk = 2;
    file_head.musea = mu;
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;file_head.k[1] = 0;
    file_head.k[2] = mu;
    file_head.k[3] = mu;

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (i = 0;i < file_head.nmoms;i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }


    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", argv[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", argv[3], option[4], option[6]);
    FILE* jack_file = open_file(namefile, "w+");
    write_header_g2(jack_file);

    std::vector<std::string>  correlators;
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//0
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//1
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//2
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//3
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//4
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//5
    correlators.emplace_back(namefile);

    printf("reading confs from file: %s", correlators[0].c_str());
    confs = read_nconfs(correlators[0].c_str());
    cout << "confs =" << confs << endl;
    for (auto name : correlators) {
        printf("checking confs from file: %s\n", name.c_str());
        int check_conf = read_nconfs(name.c_str());
        error(confs != check_conf, 1, "reading number of configuration", "not all the files have the same confs");
    }


    // FILE* infile_equal_P5A0 = open_file(namefile, "r+");
    // int count = 0;


    //confs=confs/10;
    // compute what will be the neff after the binning 
    int bin = atoi(argv[6]);
    int Neff = confs / bin;
    cout << "effective configurations after binning (" << bin << "):  " << Neff << endl;

    int Njack;
    if (strcmp(argv[11], "jack") == 0)
        Njack = Neff + 1;
    else if (strcmp(argv[11], "boot") == 0)
        Njack = Nbootstrap + 1;
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[11]= %s is not jack or boot", argv[11]);
    }
    fwrite(&Njack, sizeof(int), 1, jack_file);

    setup_file_jack(option, Njack);
    get_kinematic(0, 0, 1, 0, 0, 0);

    int var = correlators.size();
    data = calloc_corr(confs, var, file_head.l0);

    for (int i = 0; i < var; i++) {
        // read correlators[i] and store in data[conf][i][t][re/im]
        read_twopt(correlators[i].c_str(), confs, T, data, i);
    }

    data_bin = binning(confs, var, file_head.l0, data, bin);
    // //if you want to do the gamma analysis you need to do before freeing the raw data
    // // effective_mass_phi4_gamma(option, kinematic_2pt, (char*)"P5P5", data_bin, Neff, namefile_plateaux, out_gamma, 0, "M_{PS}^{ll}");
    // // effective_mass_phi4_gamma(option, kinematic_2pt, (char*)"P5P5", data_bin, Neff, namefile_plateaux, out_gamma, 1, "M_{PS1}^{ll}");
    // //effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data,  confs ,namefile_plateaux,out_gamma,3,"M_{PS}^{ll}");
    free_corr(confs, var, file_head.l0, data);

    conf_jack = create_resampling(option[4], Neff, var, file_head.l0, data_bin);
    free_corr(Neff, var, file_head.l0, data_bin);

    // ////////////////// symmetrization/////////////////////////////////////////////
    // for (int i = 0;i <= 7;i++) { symmetrise_jackboot(Njack, i, file_head.l0, conf_jack); }

    ////////////////////////////////////////////////
    corr_counter = -1;
    double* M_PS = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1, "M_{PS}", M_eff_T, jack_file);
    check_correlatro_counter(0);

    double* M_PS_OS = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 4, "M_{PS}^{OS}", M_eff_T, jack_file);
    check_correlatro_counter(1);


    fit_type fit_info;
    fit_result fit_out;

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 1;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * 1);
    fit_info.ext_P[0] = M_PS;
    fit_info.function = constant_fit;

    fit_result G_PS = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_lhs, "G_{PS}", fit_info, jack_file);
    check_correlatro_counter(2);


    fit_info.ext_P[0] = M_PS_OS;
    fit_result G_PS_OS = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_OS_lhs, "G_{PS}^{OS}", fit_info, jack_file);
    check_correlatro_counter(3);

    fit_info.n_ext_P = 0;
    fit_result ZVl = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZVl_lhs, "Z_V(l)", fit_info, jack_file);
    check_correlatro_counter(4);


    fit_info.ext_P[0] = nullptr;
    free(fit_info.ext_P);
    fit_info.n_ext_P = 4;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
    fit_info.ext_P[0] = M_PS;
    fit_info.ext_P[1] = M_PS_OS;
    fit_info.ext_P[2] = G_PS.P[0];
    fit_info.ext_P[3] = G_PS_OS.P[0];

    fit_result ZAl = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZAl_lhs, "Z_A(l)", fit_info, jack_file);
    check_correlatro_counter(5);

    fit_info.restore_default();
    free(M_PS);free(M_PS_OS);
    free_fit_result(fit_info, G_PS);free_fit_result(fit_info, G_PS_OS);
    free_fit_result(fit_info, ZVl);free_fit_result(fit_info, ZAl);
}

