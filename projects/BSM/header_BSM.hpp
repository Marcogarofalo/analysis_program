#ifndef header_BSM_H
#define header_BSM_H


#include <array>
#include <cstring> 
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>
#include "tower.hpp"
#include "mutils.hpp"

// using namespace std;

enum correlators {
    JTILDEA1P1TRIVIAL = 0,
    P1P1TRIVIAL = 1,
    P2P2TRIVIAL = 2,
    P3P3TRIVIAL = 3,
    S0S0TRIVIAL = 4,
    P1DP1NONSMEAREDNONTRIVIAL = 5,
    VECTORDENSITY3DENSITY3NONTRIVIAL = 6,
    phit = 7,
    JTILDEA1P1TRIVIALphi = 8,
    P1DP1NONSMEAREDNONTRIVIALphi = 9,
    LOCALCURRENTTAU1P1 = 10,
    LOCALCURRENTTAU2P2 = 11,
    LOCALCURRENTTAU3P3 = 12,
};


struct  header_BSM {
    int T;
    int L;
    double rho;
    double eta;
    double csw;
    double mu03;
    double m0;
    int confs;
    int Njack;
    int Nobs;
    int size_header_bin;

};

void read_header_BSM(header_BSM& head, FILE* stream) {

    int count_lines = 0;
    char chr = getc(stream);
    while (chr != EOF) {
        //Count whenever new line is encountered
        if (chr == '\n') {
            count_lines = count_lines + 1;
        }
        //take next character from file.
        chr = getc(stream);
    }
    rewind(stream);
    size_t fi=0;
    fi+=fscanf(stream, "T  %d\n", &head.T);
    fi+=fscanf(stream, "L  %d\n", &head.L);
    fi+=fscanf(stream, "rho  %lf\n", &head.rho);
    fi+=fscanf(stream, "eta  %lf\n", &head.eta);
    fi+=fscanf(stream, "csw  %lf\n", &head.csw);
    fi+=fscanf(stream, "mu03  %lf\n", &head.mu03);
    fi+=fscanf(stream, "m0  %lf\n", &head.m0);
    head.confs = (count_lines - 7) / head.T;
    error(fi != 7, 1, "read_header_BSM", "error reading header BSM file");
    error((count_lines - 7) % head.T != 0, 1, "read_header_BSM", "numer of line= %d - header_lines=%d is not a multiple of T=%d", count_lines, 7, head.T);

}
void check_header_BSM(header_BSM head, FILE* stream, std::string namefile) {
    int count_lines = 0;
    char chr = getc(stream);
    while (chr != EOF) {
        //Count whenever new line is encountered
        if (chr == '\n') {
            count_lines = count_lines + 1;
        }
        //take next character from file.
        chr = getc(stream);
    }
    rewind(stream);
    error((count_lines - 7) % head.T != 0, 1, "header", "numer of line= %d - header_lines=%d is not a multiple of T=%d\n file:%s", count_lines, 7, head.T, namefile.c_str());
    if (head.confs > (count_lines - 7) / head.T) {
        printf("file: %s \n confs does not match ref=%d read=%d  getting the minimum of the two\n", namefile.c_str(), head.confs, (count_lines - 7) / head.T);
        head.confs = (count_lines - 7) / head.T;
    }

    //     error(head.confs!=(count_lines-7)/head.T,1,"header","file: %s \n confs does not match ref=%d read=%d ",namefile.c_str(),head.confs,(count_lines-7)/head.T);

    int tmp, fi=0;
    fi+=fscanf(stream, "T  %d\n", &tmp); error(tmp != head.T, 1, "check header ", "file: %s \n T does not match ref=%d read=%d ", namefile.c_str(), head.T, tmp);
    fi+=fscanf(stream, "L  %d\n", &tmp); error(tmp != head.L, 1, "check header ", "file: %s \n L does not match ref=%d read=%d ", namefile.c_str(), head.L, tmp);
    double tmp1;
    fi+=fscanf(stream, "rho  %lf\n", &tmp1); error(fabs(tmp1 - head.rho) > 1e-6, 1, "check header ", "file: %s \n rho does not match ref=%f read=%f ", namefile.c_str(), head.rho, tmp1);
    fi+=fscanf(stream, "eta  %lf\n", &tmp1); error(fabs(tmp1 - head.eta) > 1e-6, 1, "check header ", "file: %s \n rho does not match ref=%f read=%f ", namefile.c_str(), head.eta, tmp1);
    fi+=fscanf(stream, "csw  %lf\n", &tmp1); error(fabs(tmp1 - head.csw) > 1e-6, 1, "check header ", "file: %s \n rho does not match ref=%f read=%f ", namefile.c_str(), head.csw, tmp1);
    fi+=fscanf(stream, "mu03  %lf\n", &tmp1); error(fabs(tmp1 - head.mu03) > 1e-6, 1, "check header ", "file: %s \n rho does not match ref=%f read=%f ", namefile.c_str(), head.mu03, tmp1);
    fi+=fscanf(stream, "m0  %lf\n", &tmp1); error(fabs(tmp1 - head.m0) > 1e-6, 1, "check header ", "file: %s \n rho does not match ref=%f read=%f ", namefile.c_str(), head.m0, tmp1);
    error(fi!=7, 1, "check header ", "file: %s \n header does not match ref=%d read=%d ", namefile.c_str(), 7, fi);
}

void write_header_BSM_bin(header_BSM head, FILE* stream) {
    size_t fi=0;
    fi+=fwrite(&head.T, sizeof(int), 1, stream);
    fi+=fwrite(&head.L, sizeof(int), 1, stream);
    fi+=fwrite(&head.rho, sizeof(double), 1, stream);
    fi+=fwrite(&head.eta, sizeof(double), 1, stream);
    fi+=fwrite(&head.csw, sizeof(double), 1, stream);
    fi+=fwrite(&head.mu03, sizeof(double), 1, stream);
    fi+=fwrite(&head.m0, sizeof(double), 1, stream);
    fi+=fwrite(&head.Njack, sizeof(double), 1, stream);
    error(fi != 8, 1, "write_header_BSM_bin", "error writing header BSM bin file");
}


void read_header_BSM_bin(header_BSM& head, FILE* stream) {
    size_t fi=0;
    fi+=fread(&head.T, sizeof(int), 1, stream);
    fi+=fread(&head.L, sizeof(int), 1, stream);
    fi+=fread(&head.rho, sizeof(double), 1, stream);
    fi+=fread(&head.eta, sizeof(double), 1, stream);
    fi+=fread(&head.csw, sizeof(double), 1, stream);
    fi+=fread(&head.mu03, sizeof(double), 1, stream);
    fi+=fread(&head.m0, sizeof(double), 1, stream);
    fi+=fread(&head.Njack, sizeof(double), 1, stream);
    head.size_header_bin = ftell(stream);
    error(fi != 8, 1, "read_header_BSM_bin", "error reading header BSM bin file");
}



struct data_BSM {
    int Njack;
    int Nobs;
    double** jack;

};

void read_Njack_Nobs(FILE* stream, header_BSM params, int Njack, int& Nobs) {

    long int tmp;
    int header_size = params.size_header_bin;
    int s = header_size;



    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    tmp -= header_size;

    s = Njack;

    Nobs = (tmp) / ((s) * sizeof(double));

    fseek(stream, header_size, SEEK_SET);



}
void read_dataj(FILE* stream, header_BSM  params, data_BSM& dj) {
    read_Njack_Nobs(stream, params, params.Njack, dj.Nobs);
    dj.Njack = params.Njack;
    printf("Nobs=%d   Njack=%d\n", dj.Nobs, dj.Njack);
    dj.jack = double_malloc_2(dj.Nobs, dj.Njack);
    size_t i = 0;
    for (int obs = 0; obs < dj.Nobs; obs++) {
        i += fread(dj.jack[obs], sizeof(double), dj.Njack, stream);
    }

}

void emplace_back_par_data(char* namefile, std::vector<header_BSM>& paramsj, std::vector<data_BSM>& dataj) {
    printf("%s\n",namefile);
    header_BSM params;
    data_BSM  data;
    FILE* f = open_file(namefile, "r");
    read_header_BSM_bin(params, f);
    read_dataj(f, params, data);
    fclose(f);
    printf("eta=%g  rho=%g   m0=%g  mu03=%g\n",params.eta, params.rho, params.m0, params.mu03 );
    paramsj.emplace_back(params);
    dataj.emplace_back(data);
    //printf("E1=%g    %g\n",dataj[0].jack[1][   dataj[0].Njack-1 ],    data.jack[1][data.Njack-1]);

}


std::vector<data_BSM> create_generalised_resampling(std::vector<data_BSM>& dataj) {
    // if the length is the same return dataj
    int same = 0;
    for (auto& d : dataj) {
        printf("jacks=%d\n", d.Njack);
        if (d.Njack == dataj[0].Njack)
            same++;
    }
    if (same == dataj.size()) {
        std::cout << "all the files have the same number of jack/boot , do nothing" << std::endl;
        return dataj;
    }
    else {
        std::cout << "creating generalised jack" << std::endl;
        std::vector<data_BSM> gjack;
        //gjack.resize(dataj.size());
        //jac_tot is the summ of all jackknife +1 
        //remember alle the dataj have one extra entry for the mean
        int jack_tot = 0;
        for (int e = 0; e < dataj.size();e++)
            jack_tot += dataj[e].Njack;
        jack_tot = jack_tot - dataj.size() + 1;
        std::cout << "ensembles " << dataj.size() << std::endl;
        std::cout << "jack tot= " << jack_tot << std::endl;

        //get Nobs the minimum number of observable between the files
        int Nobs = 1000;
        for (auto& d : dataj)
            if (Nobs > d.Nobs)
                Nobs = d.Nobs;




        for (int e = 0;e < dataj.size();e++) {
            data_BSM tmp;
            tmp.Njack = jack_tot;
            tmp.Nobs = Nobs;
            std::cout << Nobs << " " << jack_tot << std::endl;
            tmp.jack = double_malloc_2(Nobs, jack_tot);
            int counter = 0;
            for (int e1 = 0;e1 < dataj.size();e1++) {
                for (int j = 0;j < (dataj[e1].Njack - 1);j++) {
                    for (int o = 0;o < Nobs;o++) {
                        if (e == e1) {
                            tmp.jack[o][j + counter] = dataj[e].jack[o][j];
                        }
                        else {
                            tmp.jack[o][j + counter] = dataj[e].jack[o][dataj[e].Njack - 1];
                        }
                    }
                }
                counter += dataj[e1].Njack - 1;
            }
            for (int o = 0;o < Nobs;o++)
                tmp.jack[o][jack_tot - 1] = dataj[e].jack[o][dataj[e].Njack - 1];

            free_2(dataj[e].Nobs, dataj[e].jack);
            gjack.emplace_back(tmp);
        }
        std::vector<data_BSM>().swap(dataj);

        return gjack;

    }

}

#endif
