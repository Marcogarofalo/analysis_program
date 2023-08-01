#ifndef HLT_H
#define HLT_H
#include "global.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "non_linear_fit.hpp"

enum HLT_b {
    HLT_EXP_b = 0
};

class  HLT_type;

class HLT_out{
public:
    std::vector<double> lambdas = {};
    double** res;
    HLT_out(std::vector<double> l){ lambdas=l;}
};

class wrapper_smearing {
public:
    double (*function)(double, double*);
    double* params;
    int t;
    int Np;
    HLT_type* HLT;
    std::vector<double> lambdas = {};
    wrapper_smearing(double (*f)(double, double*), std::vector<double> p, HLT_type* HLT_) {
        function = f;
        Np = p.size();
        params = (double*)malloc(sizeof(double) * Np);
        for (int i = 0;i < Np;i++) {
            params[i] = p[i];
        }
        HLT = HLT_;
    };
    ~wrapper_smearing() {
        free(params);
    }
};

class  HLT_type {
public:
    int Tmax;
    int T;
    int Njack;
    int id;
    HLT_b type;
    double E0;
    std::vector<double> Es = {};
    std::vector<double> lambdas = {};
    double** A;
    double* R;
    double* f;
    HLT_type(int tmax, int L0, double E0, int njack, HLT_b type_b, double alpha = 0);
    ~HLT_type();


    double** HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
    FILE* outfile,   const char* description, wrapper_smearing  Delta, FILE* file_jack);



    void compute_f_EXP_b(wrapper_smearing  Delta, double epsrel = 1e-7);

};

double theta_s1(double x, double* p) {
    return 1.0 / (1 + exp(-x / p[0]));
}


#endif // !HLT_H
