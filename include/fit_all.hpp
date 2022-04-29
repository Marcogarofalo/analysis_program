#ifndef fit_all_H
#define fit_all_H


#include "global.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
#include <iostream>
#include "resampling.hpp"

class generic_header {
public:
    int struct_size;
    int Njack;
    int T;
    int L;
    std::vector<double> mus{};
    std::vector<double> thetas{};
    // generic_header() : T(0), L(0), mus{}, thetas{} {
    //     // std::cout << "header constr 0" << std::endl;
    // }
    // generic_header(std::vector<double> omu, std::vector<double> othetas)
    //     : mus(omu),
    //     thetas(othetas) {
    //     // std::cout << "header constr" << std::endl;
    // }
    // // copy contructor
    // generic_header(const generic_header& other) : T(other.T), L(other.L),
    //     mus(other.mus), thetas(other.thetas) {
    //     // std::cout << "header copy constr" << std::endl;
    // }
    // // move contructor
    // generic_header(generic_header&& other) noexcept : T(other.T), L(other.L),
    //     mus(std::move(other.mus)), thetas(std::move(other.thetas)) {
    //     // std::cout << "header move constr" << std::endl;
    // }
    // // copy assignment
    // generic_header& operator = (const generic_header& rhs) {
    //     // std::cout << "header copy ass" << std::endl;
    //     T = rhs.T;
    //     L = rhs.L;
    //     mus = rhs.mus;
    //     thetas = rhs.thetas;
    //     return *this;
    // }
    // // move assignment
    // generic_header& operator = (generic_header&& rhs) noexcept {
    //     // std::cout << "header move ass" << std::endl;
    //     T = rhs.T;
    //     L = rhs.L;
    //     mus = std::move(rhs.mus);
    //     thetas = std::move(rhs.thetas);
    //     return *this;
    // }


};

class data_single {

private:
    bool  init_err = false;
    double* errors;
    bool* computed;
    double* (*compute_mean_and_error)(int, double*);

public:
    int Nobs;
    int Njack;

    generic_header header;

    double** jack;//[obs][jack]

    std::string resampling;

    void init_error() {
        if (!init_err) {
            errors = (double*)malloc(sizeof(double) * Nobs);
            computed = (bool*)malloc(sizeof(bool) * Nobs);
            if (strcmp(resampling.c_str(), "jack") == 0)
                compute_mean_and_error = mean_and_error_jack_biased;
            else if (strcmp(resampling.c_str(), "boot") == 0)
                compute_mean_and_error = mean_and_error_boot;
            for (int j = 0;j < Nobs;j++) {
                computed[j] = false;
            }
            init_err == true;
        }
    }

    double error_jack(int i) {
        if (!computed[i]) {
            double* r = mean_and_error_jack_biased(Njack, jack[i]);
            errors[i] = r[1];
            free(r);
            computed[i] = true;
        }
        return errors[i];
    }

    // data_single() {};
    // data_single(int Nj, int No, generic_header he) : Nobs(No), Njack(Nj), header(he) {
    //     // std::cout << "single constr" << std::endl;
    //     jack = double_malloc_2(No, Nj);
    // }
    // data_single(int No, int Nj) : Nobs(No), Njack(Nj), header() {
    //     // printf("single constr 2\n");
    //     jack = double_malloc_2(No, Nj);
    // }
    // //copy constructor
    // data_single(const data_single& d) : Nobs(d.Nobs), Njack(d.Njack) {
    //     // printf("single copy constr\n");
    //     jack = double_malloc_2(d.Nobs, d.Njack);
    //     Njack = d.Njack;
    //     Nobs = d.Nobs;
    //     for (int i = 0; i < Nobs;i++)
    //         for (int j = 0; j < Njack;j++)
    //             jack[i][j] = d.jack[i][j];
    //     header = d.header;
    // }
    // //move constructor
    // data_single(data_single&& d) : Nobs(d.Nobs), Njack(d.Njack), header(std::move(d.header)) {
    //     // printf("single move\n");
    //     jack = d.jack;
    //     Njack = d.Njack;
    //     Nobs = d.Nobs;
    //     for (int i = 0; i < Nobs;i++)
    //         d.jack[i] = nullptr;
    //     d.jack = nullptr;
    //     d.Njack = 0;
    //     d.Nobs = 0;
    // }
    // //copy assignment
    // data_single& operator = (const data_single& rhs) {
    //     // printf("single copy ass\n");
    //     if (Njack != rhs.Njack) { printf("impossible to assign data_single of different sizes\n"), exit(1); };
    //     Nobs = rhs.Nobs;
    //     Njack = rhs.Njack;
    //     for (int i = 0; i < Nobs;i++)
    //         for (int j = 0; j < Nobs;j++)
    //             jack[i][j] = rhs.jack[i][j];
    //     return *this;
    // }
    // // move assignment
    // data_single& operator = (data_single&& rhs) noexcept {
    //     // printf("single move ass\n");
    //     if (this == &rhs) { printf("jackknife move assignment by itself\n");exit(-10); }
    //     Njack = rhs.Njack;
    //     Nobs = rhs.Nobs;


    //     jack = rhs.jack;
    //     // jack = (double**)malloc(sizeof(double*) * Nobs);
    //     // for (int i = 0; i < Nobs;i++)
    //     //     jack[i] = rhs.jack[i];
    //     rhs.Njack = 0;
    //     rhs.Nobs = 0;
    //     // for (int i = 0; i < Nobs;i++)
    //     //     rhs.jack[i] = nullptr;
    //     rhs.jack = nullptr;
    //     header = std::move(rhs.header);
    //     return *this;
    // }
    // ~data_single() {
    //     // printf("dest %d\n", Nobs);
    //     for (int i = 0; i < Nobs;i++)
    //         free(jack[i]);
    //     free(jack);
    // }
};

class data_all {

public:
    int ens;
    data_single* en;
    std::string resampling;

    int Nfits = 0;
    struct fit_result* fits;
    void init_error() {
        for (int e = 0;e < ens;e++)
            en[e].init_error();
    }
    void add_fit(struct fit_result fit_out) {

        int N = Nfits + 1;
        Nfits = N;

        fit_result* fit_tmp = (struct fit_result*)malloc(sizeof(struct fit_result) * N);
        for (int i = 0;i < N - 1;i++) {
            fit_tmp[i] = malloc_copy_fit_result(fits[i]);
            // fit_tmp[i].name=fits[i].name;
            for (int j = 0;j < fits[i].Npar;j++) {
                free(fits[i].P[j]);
                for (int k = 0;k < fits[i].Npar;k++) {
                    free(fits[i].C[j][k]);
                }
                free(fits[i].C[j]);
            }
            free(fits[i].chi2);
            free(fits[i].P);
            free(fits[i].C);
            if (i == N - 2) { free(fits); }
        }

        fit_tmp[N - 1] = malloc_copy_fit_result(fit_out);
        fits = fit_tmp;
        // fits = (struct fit_result*)malloc(sizeof(struct fit_result) * N);
        // for (int i = 0;i < N;i++) {
        //     fits[i] = malloc_copy_fit_result(fit_tmp[i]);

        //     for (int j = 0;j < fit_tmp[i].Npar;j++) {
        //         free(fit_tmp[i].P[j]);
        //         for (int k = 0;k < fit_tmp[i].Npar;k++) {
        //             free(fit_tmp[i].C[j][k]);
        //         }
        //         free(fit_tmp[i].C[j]);
        //     }
        //     free(fit_tmp[i].chi2);
        //     free(fit_tmp[i].P);
        //     free(fit_tmp[i].C);
        // }
        // free(fit_tmp);

    }

    // data_all() {};

    // data_all(int N) : ens(N) {
    //     // printf("data all constr\n");
    //     en = new data_single[N];
    //     for (int i = 0; i < N;i++)
    //         en[i] = data_single();
    //     resampling = std::string("\0");
    // }
    // //copy constructor
    // data_all(const data_all& d) : ens(d.ens) {
    //     // printf("copy const\n");
    //     // en = new data_single[ens];
    //     // for (int e = 0; e < ens;e++) {
    //     //     en[e] = d.en[e];
    //     // }
    //     en =(data_single*) malloc(sizeof(data_single)*ens);
    //     for (int e = 0; e < ens;e++) {
    //         new (&(en[e]))  data_single(d.en[e]);
    //     }
    //     resampling = d.resampling;
    //     fits=d.fits;
    //     Nfits=d.Nfits;
    // }

    // //move constructor
    // data_all(data_all&& d) : ens(d.ens) {
    //     // printf("move\n");
    //     en = new data_single[d.ens];
    //     en = d.en;
    //     d.ens = 0;
    //     resampling = std::move(d.resampling);
    //     fits=std::move(d.fits);
    //     Nfits=d.Nfits;
    //     delete[]    d.en;

    // }

    // //copy assignment
    // data_all& operator = (const data_all& rhs) {
    //     // printf("copy ass\n");
    //     if (ens != rhs.ens) { printf("impossible to assign data_all of different sizes\n"), exit(1); };
    //     for (int e = 0; e < ens;e++) {
    //         en[e] = rhs.en[e];
    //     }
    //     resampling = rhs.resampling;
    //     fits=rhs.fits;
    //     Nfits=rhs.Nfits;
    //     return *this;
    // }
    // // move assignment
    // data_all& operator = (data_all&& rhs) noexcept {
    //     // printf("move ass\n");
    //     if (this == &rhs) { printf("data_all move assignment by itself\n");exit(-10); }
    //     ens = rhs.ens;
    //     en = rhs.en;
    //     rhs.en = nullptr;
    //     rhs.ens = 0;
    //     resampling = std::move(rhs.resampling);
    //     fits=std::move(rhs.fits);
    //     Nfits=rhs.Nfits;
    //     return *this;
    // }

    // ~data_all() {
    //     // delete[] en;
    // }

};

struct fit_result fit_all_data(char** argv, data_all gjack,
    double lhs_fun(int, int, int, data_all, struct fit_type),
    struct fit_type fit_info, const char* label);


void print_fit_band(char** argv, data_all gjack, struct fit_type fit_info,
    struct fit_type fit_info_m0, const char* label, const char* dir_name,
    struct fit_result fit_out, struct fit_result fit_out_m0, int var, int en, double h,
    std::vector<double> xval={});
#endif

