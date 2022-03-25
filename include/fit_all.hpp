#ifndef fit_all_H
#define fit_all_H


#include "global.hpp"

class generic_header {
public:
    int T;
    int L;
    std::vector<double> mus;
    std::vector<double> thetas;
    // generic_header() {}
    // // copy contructor
    // generic_header(const generic_header& other) : T(other.T), L(other.L),
    //     mus(other.mus), thetas(other.thetas) {
    // }
    // // move contructor
    // generic_header(generic_header&& other) noexcept : T(std::move(other.T)), L(std::move(other.L)),
    //     mus(std::move(other.mus)), thetas(std::move(other.thetas)) {
    // }
    // // copy assignment
    // generic_header& operator = (const generic_header& rhs) {
    //     T = rhs.T;
    //     L = rhs.L;
    //     mus = rhs.mus;
    //     thetas = rhs.thetas;
    //     return *this;
    // }
    // // move assignment
    // generic_header& operator = (generic_header&& rhs) noexcept {
    //     T = rhs.T;
    //     L = rhs.L;
    //     mus = std::move(rhs.mus);
    //     thetas = std::move(rhs.thetas);
    //     return *this;
    // }


};

class data_single {

public:
    int Nobs;
    int Njack;

    generic_header header;

    double** jack;//[obs][jack]

    // data_single() {};
    // data_single(int Nj, int No, generic_header he) : Nobs(No), Njack(Nj), header(he) {
    //     jack = double_malloc_2(No, Nj);

    // }

    // data_single(int No, int Nj) {
    //     printf("single constr 2\n");
    //     jack = double_malloc_2(No, Nj);
    //     Njack = Nj;
    //     Nobs = No;
    //     header = generic_header();
    // }
    // //copy constructor
    // data_single(const data_single& d) : Nobs(d.Nobs), Njack(d.Njack), header(d.header) {
    //     printf("single copy\n");
    //     jack = double_malloc_2(d.Nobs, d.Njack);
    //     Njack = d.Njack;
    //     Nobs = d.Nobs;
    //     for (int i = 0; i < Nobs;i++)
    //         for (int j = 0; j < Nobs;j++)
    //             jack[i][j] = d.jack[i][j];

    // }
    // //move constructor
    // data_single(const data_single&& d) : Nobs(d.Nobs), Njack(d.Njack), header(std::move(d.header)) {
    //     printf("single move\n");
    //     jack = d.jack;
    //     Njack = d.Njack;
    //     Nobs = d.Nobs;
    // }
    // //copy assignment
    // data_single& operator = (const data_single& rhs) {
    //     printf("single copy ass\n");
    //     if (Njack != rhs.Njack) { printf("impossible to assign data_single of different sizes\n"), exit(1); };
    //     for (int i = 0; i < Nobs;i++)
    //         for (int j = 0; j < Nobs;j++)
    //             jack[i][j] = rhs.jack[i][j];
    //     return *this;
    // }
    // // move assignment
    // data_single& operator = (data_single&& rhs) noexcept {
    //     printf("single move ass\n");
    //     if (this == &rhs) { printf("jackknife move assignment by itself\n");exit(-10); }
    //     Njack = rhs.Njack;
    //     Nobs = rhs.Nobs;
    //     jack = rhs.jack;
    //     rhs.Njack = 0;
    //     rhs.Nobs = 0;
    //     for (int i = 0; i < Nobs;i++)
    //         rhs.jack[i] = nullptr;
    //     rhs.jack = nullptr;

    //     return *this;
    // }
    // ~data_single() {
    //     // printf("dest %d\n",Nobs);
    //     // for (int i = 0; i < Nobs;i++)
    //     //     free(jack[i]);
    //     // free(jack);
    // }
};

class data_all {

public:
    int ens;
    data_single* en;
    std::string resampling;

    // data_all(){};

    // data_all(int N) : ens(N) {
    //     en = new data_single[N];
    //     resampling = "aaa";
    // }
    // //copy constructor
    // data_all(const data_all& d) : ens(d.ens) {
    //     printf("copy\n");
    //     en = new data_single[ens];
    //     for (int e = 0; e < ens;e++) {
    //         en[e] = d.en[e];
    //     }
    //     resampling = d.resampling;
    // }

    // //move constructor
    // data_all(const data_all&& d) : ens(d.ens) {
    //     printf("move\n");
    //     en = (data_single*)malloc(sizeof(data_single) * ens);
    //     en = d.en;
    //     resampling = std::move(d.resampling);

    // }

    // //copy assignment
    // data_all& operator = (const data_all& rhs) {
    //     printf("copy ass\n");
    //     if (ens != rhs.ens) { printf("impossible to assign data_all of different sizes\n"), exit(1); };
    //     for (int e = 0; e < ens;e++) {
    //         en[e] = rhs.en[e];
    //     }
    //     resampling = rhs.resampling;
    //     return *this;
    // }
    // // move assignment
    // data_all& operator = (data_all&& rhs) noexcept {
    //     printf("move ass\n");
    //     if (this == &rhs) { printf("data_all move assignment by itself\n");exit(-10); }
    //     ens = rhs.ens;
    //     en = rhs.en;
    //     rhs.en = nullptr;
    //     resampling = std::move(rhs.resampling);

    //     return *this;
    // }

    // ~data_all() {
    //     // delete[] en;
    // }

};

struct fit_result fit_all_data(char** argv, data_all gjack,
    double lhs_fun(int, int, int, data_all, struct fit_type),
    struct fit_type fit_info, const char* label);

// void print_fit_band(char** argv, data_all gjack, struct fit_type fit_info,
//     struct fit_type fit_info_m0, const char* label, const char* dir_name
//     struct fit_result fit_out, struct fit_result fit_out_m0, int var, int en, int steps);
#endif

