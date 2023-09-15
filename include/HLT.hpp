#ifndef HLT_H
#define HLT_H
#include "global.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "non_linear_fit.hpp"

#include "arb.h"
#include "acb_modular.h"
#include "acb_calc.h"

enum HLT_b {
    HLT_EXP_b = 0
};

class  HLT_type_d;

class HLT_out {
public:
    std::vector<double> lambdas = {};
    double** res;
    HLT_out(std::vector<double> l) { lambdas = l; }
};

class wrapper_smearing_d {
public:
    double (*function)(double, double*);
    double* params;
    int t;
    int Np;
    HLT_type_d* HLT;
    std::vector<double> lambdas = {};
    wrapper_smearing_d(double (*f)(double, double*), std::vector<double> p, HLT_type_d* HLT_) {
        function = f;
        Np = p.size();
        params = (double*)malloc(sizeof(double) * Np);
        for (int i = 0;i < Np;i++) {
            params[i] = p[i];
        }
        HLT = HLT_;
    };
    ~wrapper_smearing_d() {
        free(params);
    }
};

class  HLT_type_d {
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
    bool f_allocated=false;
    HLT_type_d(int tmax, int L0, double E0, int njack, HLT_b type_b, double alpha = 0);
    ~HLT_type_d();


    double** HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
        FILE* outfile, const char* description, wrapper_smearing_d  Delta, FILE* file_jack);



    void compute_f_EXP_b(wrapper_smearing_d  &Delta, double epsrel = 1e-7);

};

double theta_s1_d(double x, double* p);
double gaussian_for_HLT_d(double x, double* p);

////////////////////////////////////////////////////////////////////////////////////////////////
// using ARB
////////////////////////////////////////////////////////////////////////////////////////////////


class  fit_type_HLT : public fit_type {
public:
    std::vector<double> lambdas = {};

};

class  HLT_type;

class wrapper_smearing {
public:
    int  (*function)(acb_ptr , const acb_t , void* , slong , slong );
    arb_t* params;
    int t;
    int Np;
    HLT_type* HLT;
    wrapper_smearing(int  (*f)(acb_ptr , const acb_t , void* , slong , slong ), std::vector<double> p, HLT_type* HLT_) {
        function = f;
        Np = p.size();
        params = (arb_t*)malloc(sizeof(arb_t) * Np);
        for (int i = 0;i < Np;i++) {
            arb_init(params[i]);
            arb_set_d(params[i], p[i]);
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
    int prec;
    HLT_b type;
    arb_t E0;
    std::vector<double> Es = {};
    std::vector<double> lambdas = {};
    arb_mat_t A;
    arb_mat_t W;
    arb_mat_t R;
    arb_mat_t f;
    bool f_allocated=false;
    HLT_type(int tmax, int L0, double E0,  HLT_b type_b, int prec_, double alpha = 0);
    ~HLT_type();


    double** HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
        FILE* outfile, const char* description, wrapper_smearing  Delta, FILE* file_jack, fit_type_HLT fit_info);



    void compute_f_EXP_b(wrapper_smearing  &Delta, double epsrel = 1e-7);

};
int gaussian_for_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
int theta_s1_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
int no_smearing_for_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);



#endif // !HLT_H
