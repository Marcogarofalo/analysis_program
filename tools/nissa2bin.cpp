#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <iostream>

#include "global.hpp"
#include "non_linear_fit.hpp"
#include "read_nissa.hpp"
#include "fit_all.hpp"
#include "mutils.hpp"
#include "read.hpp"

generic_header set_header() {
    generic_header head;
    return head;
}

template<typename T>
T convert(std::string& in) {
    T t;
    return t;
}

template<>
int convert<int>(std::string& in) {
    return std::stoi(in);
}

template<>
double convert<double>(std::string& in) {
    return std::stod(in);
}
template<>
std::string convert<std::string>(std::string& in) {
    return in;
}

template<typename T>
T read_param(std::string namefile, std::string namepar) {
    std::fstream newfile;
    int line = 0;
    T out;

    newfile.open(namefile, std::ios::in); // open a file to perform read operation using file object
    int match = 0;
    if (newfile.is_open()) { // checking whether the file is open
        std::string tp;
        while (getline(newfile, tp)) { // read data from file object and put it into string.
            line++;
            std::vector<std::string> x = split(tp, ' ');

            if (x.empty() == 0) {
                if (x[0].compare(namepar) == 0) {
                    out = convert<T>(x[1]);

                    match++;
                    // break;
                }
            }
        }
        newfile.close(); // close the file object.
    }
    else {
        error(0 == 0, 1, "nissa2bin read_param",
            "unable to open %s", namefile);
    }

    if (match == 0) {
        printf("param %s not found in file %s\n", namepar.c_str(), namefile.c_str());
        exit(1);
    }
    if (match > 1) {
        printf("multiple lines line of param  %s in file  %s\n", namepar.c_str(), namefile.c_str());
        exit(1);
    }
    return out;

}

template<typename T>
std::vector<T> read_vector(std::string namefile, std::string namepar) {
    std::fstream newfile;
    int line = 0;
    std::vector<T> out;

    newfile.open(namefile, std::ios::in); // open a file to perform read operation using file object
    int match = 0;
    if (newfile.is_open()) { // checking whether the file is open
        std::string tp;
        while (getline(newfile, tp)) { // read data from file object and put it into string.
            line++;
            std::vector<std::string> x = split(tp, ' ');

            if (x.empty() == 0) {
                if (x[0].compare(namepar) == 0) {
                    out.resize(x.size() - 1);
                    for (int i = 1; i < x.size();i++) {
                        out[i - 1] = convert<T>(x[i]);
                    }

                    match++;
                    // break;
                }
            }
        }
        newfile.close(); // close the file object.
    }
    else {
        error(0 == 0, 1, "nissa2bin read_param",
            "unable to open %s", namefile);
    }

    if (match == 0) {
        printf("param %s not found in file %s\n", namepar.c_str(), namefile.c_str());
        exit(1);
    }
    if (match > 1) {
        printf("multiple lines line of param  %s in file  %s\n", namepar.c_str(), namefile.c_str());
        exit(1);
    }
    return out;
}

int main(int argc, char** argv) {

    error(argc < 4, 1, "main ",
        "usage:./nissa2bin   in_path  outfile  -i infile");

    FILE* infile = open_file(argv[4], "r");
    generic_header head;
    std::vector<int> confs = read_vector<int>(argv[4], "confs");
    head.Njack = confs.size();
    head.T = read_param<int>(argv[4], "T");
    head.L = read_param<int>(argv[4], "L");
    head.beta = read_param<double>(argv[4], "beta");
    head.kappa = read_param<double>(argv[4], "kappa");
    head.mus = read_vector<double>(argv[4], "mus");
    head.rs = read_vector<double>(argv[4], "rs");
    head.thetas = read_vector<double>(argv[4], "thetas");
    head.gammas = read_vector<std::string>(argv[4], "gammas");
    head.smearing = read_vector<std::string>(argv[4], "smearing");
    std::string bintype = read_param<std::string>(argv[4], "bintype");
    int bin = read_param<int>(argv[4], "bin");

    std::vector<std::string> filesname = read_vector<std::string>(argv[4], "files");

    head.size = head.ncorr * 2 * filesname.size();

    head.print_header();

    std::string file0(argv[1]);
    char conf4int[NAMESIZE];

    mysprintf(conf4int, NAMESIZE, "%04d", confs[0]);
    struct_nissa_out_info nissa_out(file0 + "/" + conf4int + "/" + filesname[0]);
    nissa_out.print();
    error(nissa_out.T != head.T, 1, "main", "T value in the input file=%d while in the nissa file we found %d\n", head.T, nissa_out.T);

    ///////////////////////////////////////
    // compute how many correlators we will store
    head.ncorr = ((nissa_out.Ncorr * head.gammas.size()) / nissa_out.Ngamma) * filesname.size();
    head.size = head.ncorr * 2 * nissa_out.T;
    std::vector<int> id_gamma = nissa_out.inidices_of_gamma(head.gammas);
    printf("Ncorr nisaa= %d Ncorr to store=%d\n", nissa_out.Ncorr, head.ncorr);
    // for(auto i :id_gamma){printf("%d\t",i);}printf("\n");
    ///////////////////////////////


    ////////////////////////////////////
    // reading the files
    double**** data = calloc_corr(head.Njack, head.ncorr, head.T);
#pragma omp parallel for
    for (int ic = 0; ic < head.Njack; ic++) {
        double**** data_n = calloc_corr(filesname.size(), nissa_out.Ncorr, head.T);
        int id_lhs = 0;
        for (int iif = 0; iif < filesname.size();iif++) {
            mysprintf(conf4int, NAMESIZE, "%04d", confs[ic]);
            
            read_all_nissa_gamma(data_n[iif], file0 + "/" + conf4int + "/" + filesname[iif], nissa_out.Ncorr, head.T, id_gamma);
            // printf("time to read %gs\n", timestamp() - a);a = timestamp();
            for (int icorr = 0; icorr < head.ncorr / filesname.size(); icorr++) {
                for (int t = 0;t < head.T;t++) {
                    data[ic][id_lhs][t][0] = data_n[iif][icorr][t][0];
                    data[ic][id_lhs][t][1] = data_n[iif][icorr][t][1];
                }
                id_lhs++;
            }

            // printf("time to copy %gs\n", timestamp() - a);a = timestamp();

        }
        free_corr(filesname.size(), head.ncorr, head.T, data_n);
    }

    ///////////////////////////////////
    //// binning
    ///////////////////////////////////
    double**** data_bin;
    if (bintype.compare("block") == 0) {
        data_bin = binning(head.Njack, head.ncorr, head.T, data, bin);
        free_corr(head.Njack, head.ncorr, head.T, data);
        head.Njack = head.Njack / bin;

    }
    else if (bintype.compare("into") == 0) {
        data_bin = calloc_corr(bin, head.ncorr, head.T);
        data_bin = bin_intoN(data, head.ncorr, head.T, head.Njack, bin);
        free_corr(head.Njack, head.ncorr, head.T, data);
        head.Njack = bin;

    }
    ///////////////////////// 
    // writing 
    ////////////////////////
    // opening output and write header
    FILE* outfile = open_file(argv[2], "w+");
    head.write_header(outfile);
    // writing the double chunk
    for (int ic = 0; ic < head.Njack; ic++) {
        fwrite(&ic, sizeof(int), 1, outfile);
        for (int iv = 0; iv < head.ncorr; iv++) {
            for (int t = 0; t < head.T; t++) {
                fwrite(data_bin[ic][iv][t], sizeof(double), 2, outfile);
            }
        }
    }
    fclose(outfile);
    /////////////////////////////////////////////////////////////////////////////////////////////
    // check
    /////////////////////////////////////////////////////////////////////////////////////////////
    outfile = open_file(argv[2], "r");
    printf("T=%d\n",head.T);
    head.read_header(outfile);
    head.print_header();

    double ****data_check;
    if (bintype.compare("block") == 0) {
        data_check = binning(head.Njack, head.ncorr, head.T, data, bin);
    }
    else if (bintype.compare("into") == 0) {
        data_check = calloc_corr(bin, head.ncorr, head.T);
    }
    for (int ic = 0; ic < head.Njack; ic++) {
        int icc;
        fread(&icc, sizeof(int), 1, outfile);
        error(icc!=ic,1,"main","conf number do not match  expected=%d read=%d\n",ic, icc);
        for (int iv = 0; iv < head.ncorr; iv++) {
            for (int t = 0; t < head.T; t++) {
                fread(data_check[ic][iv][t], sizeof(double), 2, outfile);
                error(fabs(data_check[ic][iv][t][0]-data_bin[ic][iv][t][0])>1e-7,1,"main",
                "correlator real part do not match  expected=%.10g read=%.10g  t=%d corr=%d\n",data_bin[ic][iv][t][0],data_check[ic][iv][t][0],t,iv);
                error(fabs(data_check[ic][iv][t][1]-data_bin[ic][iv][t][1])>1e-7,1,"main",
                "correlator real part do not match  expected=%.10g read=%.10g  t=%d corr=%d\n",data_bin[ic][iv][t][1],data_check[ic][iv][t][1],t,iv);
            }
        }
    }
    
    fclose(outfile);

}