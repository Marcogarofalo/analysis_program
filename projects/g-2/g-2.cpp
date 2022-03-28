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
#include "tower.hpp"
#include "fit_all.hpp"
// #include "lhs_functions.hpp"

#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>

#include <gsl/gsl_integration.h>
using namespace std;


struct  kinematic kinematic_2pt;


class configuration_class {
public:
    std::vector<std::string> iconfs;
    std::vector<int> to_bin;//0 last, 1 alone, 2 first, 3 in the list
    std::vector<int> next_to_bin;
    int confs_after_binning;


    void check_binnign() {
        confs_after_binning = 0;
        for (int i = 0; i < iconfs.size(); i++) {
            to_bin[i] = 1;// alone
            next_to_bin[i] = -1;
            confs_after_binning++;
            for (int j = i - 1; j >= 0; j--) {
                if (strcmp(iconfs[i].c_str(), iconfs[j].c_str()) == 0) {
                    if (to_bin[j] == 1) {
                        to_bin[j] = 2;// has a similar and it is first in the list
                    }
                    else if (to_bin[j] == 0) {
                        to_bin[j] = 3;// it needs to be binned with others
                    }
                    next_to_bin[j] = i;
                    to_bin[i] = 0;// last in the list of confs to bin
                    confs_after_binning--;
                    break;
                }

            }

        }
        // for (int i = 0; i < iconfs.size(); i++) {
        //     printf("%s   %d    %d\n", iconfs[i].c_str(), to_bin[i], next_to_bin[i]);
        // }

    }


    configuration_class(const char stream[NAMESIZE]) {

        long int tmp;
        int s = file_head.l1 + 1 + 3;
        int count = 0;
        string line;
        ifstream file(stream);
        while (getline(file, line)) {
            count++;
            if (line.compare(0, 1, "#") == 0) {
                line.erase(8);
                line.erase(0, 1);
                iconfs.emplace_back(line);
                // printf("%s\n", line.c_str());
            }

        }
        std::cout << "lines=" << count << std::endl;
        error(count % s != 0, 1, "read_nconfs lines and T do not match", "lines=%i   T/2=%i", count, file_head.l1);
        int c = count / s;
        std::cout << "confs=" << c << std::endl;

        to_bin = std::vector<int>(iconfs.size());
        next_to_bin = std::vector<int>(iconfs.size());
    }
};

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

void line_read_param(char** option, const char* corr, double& tmin, double& tmax, int& sep, const char* namefile_plateaux) {

    int line = 0;
    // option[3] path
    // option[6] file
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
                    tmin = stod(x[2]);
                    tmax = stod(x[3]);
                    sep = stoi(x[4]);
                    printf("correlator %s  plateaux %g  %g %d\n", correlator.c_str(), tmin, tmax, sep);
                    match++;
                    // break;
                }
            }
        }
        newfile.close(); // close the file object.
    }
    else {
        error(0 == 0, 1, "correlators_analysis.cpp line_read_plateaux",
            "unable to open %s", namefile);
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
    if (tmax > tmin) {
        printf("\n\n bel errore del cazzo\n\n");
    }
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
void return_conf_binned(double* tmp_b, double** tmp, configuration_class confs, int id, int T) {
    double N = 1;
    for (int t = 0;t < T / 2 + 1;t++)
        tmp_b[t] += tmp[id][t];
    // printf("%s\t", confs.iconfs[id].c_str());

    while (confs.to_bin[id] >= 2) {
        id = confs.next_to_bin[id];
        // printf("%s\t", confs.iconfs[id].c_str());
        for (int t = 0;t < T / 2 + 1;t++)
            tmp_b[t] += tmp[id][t];
        N++;
    }
    for (int t = 0;t < T / 2 + 1;t++)
        tmp_b[t] /= N;

    // printf("\n");

}

void read_twopt(const char namefile[NAMESIZE], configuration_class confs, int T, double**** to_write, int id, int Nb) {
    // char tmp[NAMESIZE];
    int sd = 0;
    std::fstream newfile;
    newfile.open(namefile, std::ios::in);
    double** tmp = double_malloc_2(confs.iconfs.size(), T / 2 + 2);

    if (newfile.is_open()) { // checking whether the file is open
        std::string tp;
        for (int i = 0;i < confs.iconfs.size();i++) {
            getline(newfile, tp);
            // cout<< tp<< endl;
            getline(newfile, tp);
            // cout<< tp<< endl;
            getline(newfile, tp);
            // cout<< tp<< endl;
            for (int t = 0;t < T / 2 + 1;t++) {
                getline(newfile, tp);
                // to_write[i][id][t][0] = stod(tp);
                tmp[i][t] = stod(tp);
                // printf("%.15f\n", to_write[i][id][t][0]);
            }

        }
    }
    else {
        error(0 == 0, 1, "correlators_analysis.cpp line_read_plateaux",
            "unable to open %s", namefile);
    }

    double** data = double_calloc_2(confs.confs_after_binning, T / 2 + 1);
    int count = 0;
    for (int i = 0;i < confs.iconfs.size();i++) {
        if (confs.to_bin[i] == 1 || confs.to_bin[i] == 2) { // if it is alone or it is first of the list
            double* tmp_b = (double*)calloc((T / 2 + 2), sizeof(double));
            return_conf_binned(tmp_b, tmp, confs, i, T);
            for (int t = 0;t < T / 2 + 1;t++) {
                data[count][t] = tmp_b[t];
            }
            free(tmp_b);
            count++;
        }
    }
    int bin = confs.confs_after_binning / Nb;
    for (int t = 0;t < T / 2 + 1;t++) {
        int l = 0;
        for (int i = 0;i < bin;i++) {
            for (int l = 0;l < Nb;l++) {
                to_write[l][id][t][0] += data[i + l * bin][t];

            }
        }
        for (l = 0;l < Nb;l++) {
            to_write[l][id][t][0] /= ((double)bin);
        }

    }
    free_2(confs.confs_after_binning, data);

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



    error(argc != 14, 1, "main ",
        "usage:./g-2  blind/see/read_plateaux -p path basename -bin $bin"
        "   -L L jack/boot  -mu mu  mus mus ");
    error(strcmp(argv[1], "blind") != 0 && strcmp(argv[1], "see") != 0 && strcmp(argv[1], "read_plateaux") != 0, 1, "main ",
        "argv[1] only options:  blind/see/read_plateaux ");

    // cluster::IO_params params;
    // mysprintf(namefile, NAMESIZE, "%s/%s", argv[3], argv[4]);
    // FILE* infile = open_file(namefile, "r+");
    // read_header_phi4(infile, params);




    error(strcmp(argv[5], "-bin") != 0, 1, "main", "argv[5] must be: -bin");

    error(strcmp(argv[7], "-L") != 0, 1, "main", "argv[7] must be: -L");
    error(strcmp(argv[9], "jack") != 0 && strcmp(argv[9], "boot") != 0, 1, "main",
        "argv[6] only options: jack/boot");
    error(strcmp(argv[10], "-mu") != 0, 1, "main", "argv[7] must be: -mu");

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
    mysprintf(option[4], NAMESIZE, argv[9]); //resampling
    mysprintf(option[5], NAMESIZE, "no"); // pdf

    file_head.l1 = atoi(argv[8]);
    file_head.l0 = file_head.l1 * 2;
    file_head.l2 = file_head.l1;
    file_head.l3 = file_head.l1;

    double mu = atof(argv[11]);
    double mus1 = atof(argv[12]);
    double mus2 = atof(argv[13]);
    generic_header header;
    header.L = file_head.l1;
    header.T = file_head.l0;
    header.mus = { mu, mus1, mus2 };

    mysprintf(namefile, NAMESIZE, "%s_mu.%f", argv[4], mu);


    mysprintf(option[6], NAMESIZE, namefile); // basename

    printf("resampling %s\n", option[4]);
    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, argv[9]);
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


    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_P5A0.txt", argv[3], argv[4], mus1);//6
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_P5P5.txt", argv[3], argv[4], mus1);//7
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_VKVK.txt", argv[3], argv[4], mus1);//8
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_P5A0.txt", argv[3], argv[4], mus1);//9
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_P5P5.txt", argv[3], argv[4], mus1);//10
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_VKVK.txt", argv[3], argv[4], mus1);//11
    correlators.emplace_back(namefile);

    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_P5A0.txt", argv[3], argv[4], mus2);//12
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_P5P5.txt", argv[3], argv[4], mus2);//13
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_VKVK.txt", argv[3], argv[4], mus2);//14
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_P5A0.txt", argv[3], argv[4], mus2);//15
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_P5P5.txt", argv[3], argv[4], mus2);//16
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_VKVK.txt", argv[3], argv[4], mus2);//17
    correlators.emplace_back(namefile);



    // printf("reading confs from file: %s", correlators[0].c_str());
    // auto iconfs = read_nconfs(correlators[0].c_str());
    std::vector<configuration_class> myconfs;

    // confs = myconfs[0].iconfs.size();
    int count = 0;
    for (auto name : correlators) {
        printf("reading  confs from file: %s\n", name.c_str());
        myconfs.emplace_back(name.c_str());
        // myconfs[cout]=read_nconfs(name.c_str());
        myconfs[count].check_binnign();
        cout << "number of different configurations:" << myconfs[count].confs_after_binning << endl;
        count++;
        // printf("checking confs from file: %s\n", name.c_str());
        // // auto check_iconfs = read_nconfs(name.c_str());
        // configuration_class check_confs(name.c_str());
        // for (int i = 0;i < confs;i++) {
        //     error(myconfs.iconfs[i] != check_confs.iconfs[i], 1, "configurations id do not match", "");
        // }
        // error(confs != check_confs.iconfs.size(), 1, "reading number of configuration", "not all the files have the same confs");
    }
    // which_confs_are_the_same(std::vector<std::string> iconfs);

    // FILE* infile_equal_P5A0 = open_file(namefile, "r+");
    // int count = 0;


    //confs=confs/10;
    // compute what will be the neff after the binning 
    int bin = atoi(argv[6]);
    int Neff = bin;
    // cout << "effective configurations after binning ( bin size " << confs / bin << "):  " << Neff << endl;

    int Njack;
    if (strcmp(argv[9], "jack") == 0)
        Njack = Neff + 1;
    else if (strcmp(argv[9], "boot") == 0)
        Njack = Nbootstrap + 1;
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[9]= %s is not jack or boot", argv[9]);
    }
    fwrite(&Njack, sizeof(int), 1, jack_file);

    setup_file_jack(option, Njack);
    get_kinematic(0, 0, 1, 0, 0, 0);

    int var = correlators.size();
    data = calloc_corr(bin, var, file_head.l0);

    for (int i = 0; i < var; i++) {
        // read correlators[i] and store in data[conf][i][t][re/im]
        read_twopt(correlators[i].c_str(), myconfs[i], T, data, i, bin);
    }

    // data_bin = binning(confs, var, file_head.l0, data, bin);
    // data_bin = binning_toNb(confs, var, file_head.l0, data, bin);
    // //if you want to do the gamma analysis you need to do before freeing the raw data
    // // effective_mass_phi4_gamma(option, kinematic_2pt, (char*)"P5P5", data_bin, Neff, namefile_plateaux, out_gamma, 0, "M_{PS}^{ll}");
    // // effective_mass_phi4_gamma(option, kinematic_2pt, (char*)"P5P5", data_bin, Neff, namefile_plateaux, out_gamma, 1, "M_{PS1}^{ll}");
    // //effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data,  confs ,namefile_plateaux,out_gamma,3,"M_{PS}^{ll}");
    // free_corr(bin, var, file_head.l0, data);

    conf_jack = create_resampling(option[4], Neff, var, file_head.l0, data);
    free_corr(Neff, var, file_head.l0, data);

    // ////////////////// symmetrization/////////////////////////////////////////////
    // for (int i = 0;i <= 7;i++) { symmetrise_jackboot(Njack, i, file_head.l0, conf_jack); }

    ////////////////////////////////////////////////
    corr_counter = -1;
    double* M_PS = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1, "M_{PS}^{eq}", M_eff_T, jack_file);
    check_correlatro_counter(0);

    double* M_PS_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 4, "M_{PS}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(1);

    double* M_eta = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1 + 6, "M_{eta}^{eq}", M_eff_T, jack_file);
    check_correlatro_counter(2);

    double* M_eta_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 4 + 6, "M_{eta}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(3);

    double* M_phi = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 2 + 6, "M_{phi}^{eq}", M_eff_T, jack_file);
    check_correlatro_counter(4);

    double* M_phi_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5 + 6, "M_{phi}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(5);

    double* M_eta1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1 + 12, "M_{eta1}^{eq}", M_eff_T, jack_file);
    check_correlatro_counter(6);

    double* M_eta1_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 4 + 12, "M_{eta1}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(7);

    double* M_phi1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 2 + 12, "M_{phi1}^{eq}", M_eff_T, jack_file);
    check_correlatro_counter(8);

    double* M_phi1_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5 + 12, "M_{phi1}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(9);


    fit_type fit_info;
    fit_result fit_out;

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 1;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * 1);
    fit_info.function = constant_fit;

    fit_info.ext_P[0] = M_PS;
    fit_result G_PS = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_lhs<1>, "G_{PS}^{eq}", fit_info, jack_file);
    check_correlatro_counter(10);


    fit_info.ext_P[0] = M_PS_op;
    fit_result G_PS_OS = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_OS_lhs<4>, "G_{PS}^{op}", fit_info, jack_file);
    check_correlatro_counter(11);

    fit_info.ext_P[0] = M_eta;
    fit_result G_eta = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_lhs<1 + 6>, "G_{eta}^{eq}", fit_info, jack_file);
    check_correlatro_counter(12);


    fit_info.ext_P[0] = M_eta_op;
    fit_result G_eta_OS = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_OS_lhs<4 + 6>, "G_{eta}^{op}", fit_info, jack_file);
    check_correlatro_counter(13);

    fit_info.ext_P[0] = M_eta1;
    fit_result G_eta1 = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_lhs<1 + 12>, "G_{eta1}^{eq}", fit_info, jack_file);
    check_correlatro_counter(14);


    fit_info.ext_P[0] = M_eta1_op;
    fit_result G_eta1_OS = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_OS_lhs<4 + 12>, "G_{eta1}^{op}", fit_info, jack_file);
    check_correlatro_counter(15);

    //////////////////////////////////////////////////
    fit_info.n_ext_P = 0;
    fit_info.mu = mu;
    fit_result ZVl = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZVl_lhs<4, 3>, "Z_V(l)", fit_info, jack_file);
    check_correlatro_counter(16);
    fit_info.mu = mus1;
    fit_result ZVs = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZVl_lhs<4 + 6, 3 + 6>, "Z_V(s)", fit_info, jack_file);
    check_correlatro_counter(17);
    fit_info.mu = mus2;
    fit_result ZVs1 = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZVl_lhs<4 + 12, 3 + 12>, "Z_V(s1)", fit_info, jack_file);
    check_correlatro_counter(18);


    fit_info.ext_P[0] = nullptr;
    free(fit_info.ext_P);
    fit_info.n_ext_P = 4;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);

    fit_info.ext_P[0] = M_PS_op;
    fit_info.ext_P[1] = M_PS;
    fit_info.ext_P[2] = G_PS_OS.P[0];
    fit_info.ext_P[3] = G_PS.P[0];
    fit_info.mu = mu;
    fit_result ZAl = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZAl_lhs<1, 0>, "Z_A(l)", fit_info, jack_file);
    check_correlatro_counter(19);


    fit_info.ext_P[0] = M_eta_op;
    fit_info.ext_P[1] = M_eta;
    fit_info.ext_P[2] = G_eta_OS.P[0];
    fit_info.ext_P[3] = G_eta.P[0];
    fit_info.mu = mus1;
    fit_result ZAs = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZAl_lhs<1 + 6, 0 + 6>, "Z_A(s)", fit_info, jack_file);
    check_correlatro_counter(20);

    fit_info.ext_P[0] = M_eta1_op;
    fit_info.ext_P[1] = M_eta1;
    fit_info.ext_P[2] = G_eta1_OS.P[0];
    fit_info.ext_P[3] = G_eta1.P[0];
    fit_info.mu = mus2;
    fit_result ZAs1 = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZAl_lhs<1 + 12, 0 + 12>, "Z_A(s1)", fit_info, jack_file);
    check_correlatro_counter(21);

    fit_info.restore_default();

    // double* amu_sd = (double*)malloc(sizeof(double) * Njack);
    double* VV = (double*)malloc(sizeof(double) * file_head.l0);
    double mean, err;
    int seed;
    line_read_param(option, "a", mean, err, seed, namefile_plateaux);
    double* a = fake_sampling(resampling, mean, err, Njack, seed);
    line_read_param(option, "ZA", mean, err, seed, namefile_plateaux);
    double* ZA = fake_sampling(resampling, mean, err, Njack, seed);
    line_read_param(option, "ZV", mean, err, seed, namefile_plateaux);
    double* ZV = fake_sampling(resampling, mean, err, Njack, seed);



    double (*int_scheme)(int, int, double*);


    int_scheme = integrate_reinman;
    double* amu_sd = compute_amu_sd(conf_jack, 2, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "amu_{sd}(eq,l)", resampling);
    write_jack(amu_sd, Njack, jack_file);
    check_correlatro_counter(22);
    printf("amu_sd(eq,l) = %g  %g\n", amu_sd[Njack - 1], error_jackboot(resampling, Njack, amu_sd));
    free(amu_sd);

    amu_sd = compute_amu_sd(conf_jack, 5, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "amu_{sd}(op,l)", resampling);
    write_jack(amu_sd, Njack, jack_file);
    check_correlatro_counter(23);
    printf("amu_sd(op,l) = %g  %g\n", amu_sd[Njack - 1], error_jackboot(resampling, Njack, amu_sd));
    free(amu_sd);


    int_scheme = integrate_simpson38;
    amu_sd = compute_amu_sd(conf_jack, 2, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "amu_{sd,simpson38}(eq,l)", resampling);
    write_jack(amu_sd, Njack, jack_file);
    check_correlatro_counter(24);
    printf("amu_sd_simpson38(eq,l) = %g  %g\n", amu_sd[Njack - 1], error_jackboot(resampling, Njack, amu_sd));
    free(amu_sd);

    amu_sd = compute_amu_sd(conf_jack, 5, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "amu_{sd,simpson38}(op,l)", resampling);
    write_jack(amu_sd, Njack, jack_file);
    check_correlatro_counter(25);
    printf("amu_sd_simpson38(op,l) = %g  %g\n", amu_sd[Njack - 1], error_jackboot(resampling, Njack, amu_sd));
    free(amu_sd);


    free(M_PS);free(M_PS_op);
    free_fit_result(fit_info, G_PS);free_fit_result(fit_info, G_PS_OS);
    free_fit_result(fit_info, ZVl);free_fit_result(fit_info, ZAl);
}

