#ifndef HLT_H
#define HLT_H
#include <array>
#include <functional>
#include <string>
#include <vector>

#include "global.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "non_linear_fit.hpp"
#include "arb.h"
#include "acb_calc.h"


enum HLT_b {
    HLT_EXP_b = 0,
    HLT_EXP_bT = 1,
    HLT_deriv_EXP_b = 2,
    HLT_EXP_bovert = 3,
    HLT_EXP_bovert2 = 4,
    HLT_INVALID_b
};

enum HLT_A_strategy {
    HLT_A_fast = 0,
    HLT_A_INTEGRAL = 1,
    HLT_A_INVALID
};

enum HLT_plateau_strategy {
    HLT_plateau_strategy_left = 0,
    HLT_plateau_strategy_right = 1,
    HLT_PLATEAU_STRATEGY_INVALID
};


// class  HLT_type_d;

class HLT_out {
public:
    std::vector<double> lambdas = {};
    double** res;
    HLT_out(std::vector<double> l) { lambdas = l; }
};

// class wrapper_smearing_d {
// public:
//     double (*function)(double, double*);
//     double* params;
//     int t;
//     int Np;
//     HLT_type_d* HLT;
//     std::vector<double> lambdas = {};
//     wrapper_smearing_d(double (*f)(double, double*), std::vector<double> p, HLT_type_d* HLT_) {
//         function = f;
//         Np = p.size();
//         params = (double*)malloc(sizeof(double) * Np);
//         for (int i = 0;i < Np;i++) {
//             params[i] = p[i];
//         }
//         HLT = HLT_;
//     };
//     ~wrapper_smearing_d() {
//         free(params);
//     }
// };

// class  HLT_type_d {
// public:
//     int Tmax;
//     int T;
//     int Njack;
//     int id;
//     HLT_b type;
//     double E0;
//     std::vector<double> Es = {};
//     std::vector<double> lambdas = {};
//     double** A;
//     double* R;
//     double* f;
//     bool f_allocated = false;
//     HLT_type_d(int tmax, int L0, double E0, int njack, HLT_b type_b, double alpha = 0);
//     ~HLT_type_d();


//     double** HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
//         FILE* outfile, const char* description, wrapper_smearing_d  Delta, FILE* file_jack);



//     void compute_f_EXP_b(wrapper_smearing_d& Delta, double epsrel = 1e-7);

// };

// double gaussian_for_HLT_d(double x, double* p);
// double theta_for_HLT_d(double x, double* p);
////////////////////////////////////////////////////////////////////////////////////////////////
// using ARB
////////////////////////////////////////////////////////////////////////////////////////////////


class  fit_type_HLT : public fit_type {
public:
    std::vector<double> lambdas = {};
    double lambda_start = 1;
    double maxE_check_reconstuct = 1.0f;
    double stepsE_check_reconstuct = 10;
    FILE* outfile_kernel = NULL;
    FILE* outfile_AoverB = NULL;
    FILE* outfile = NULL;
    int nsame = 4;
    int nlambda_max = 20;
    HLT_plateau_strategy plateau_strategy = HLT_plateau_strategy_left;
    double reduce_lambda = 0.75;
    bool diag_cov = false;
};


class  HLT_type;

// pass it as a reference because we did not write the copy construtor
class wrapper_smearing {
public:
    int  (*function)(acb_ptr, const acb_t, void*, slong, slong);
    arb_ptr params;
    int t;
    int Np;
    arb_t Norm;
    arb_t b_tmp;
    HLT_type* HLT;

    void normilise_smearing();

    wrapper_smearing(int  (*f)(acb_ptr, const acb_t, void*, slong, slong), std::vector<double> p, HLT_type* HLT_);

    ~wrapper_smearing() {
        // for (int i = 0;i < Np;i++) {
        //     arb_clear(params[i]);
        // }
        arb_clear(Norm);
        _arb_vec_clear(params, Np);
        arb_clear(b_tmp);
    }
};

struct HLT_type_input {
    int tmax = -1;
    int tmin = -1;
    int T = -1;
    double E0 = -1;
    HLT_b type_b = HLT_INVALID_b;
    HLT_A_strategy A_strategy = HLT_A_fast;
    int prec = -1;
    double alpha = 0.0;
    bool normalize_kernel = false;
    int integration_deg_limit = 100;
    int integration_eval_limit = 100000;
    int integration_depth_limit = 10000;
    int integration_verbose = 0;
    double integration_maxE = 10000;
};


class  HLT_type {
private:
    arb_t lam;
    arb_mat_t Wl;
    arb_mat_t Wf;
    arb_mat_t RT;
    arb_mat_t RTWR;


public:
    HLT_type_input info;
    arb_t E0_arb;
    std::vector<double> Es = {};
    std::vector<double> lambdas = {};
    std::vector<double> Ag = {};
    std::vector<double> Ag_ref = {};
    double A0;
    arb_t A0_arb;
    double A0_ref;
    arb_t A0_arb_ref;
    std::vector<double> Bg = {};
    double C0;
    arb_mat_t A;
    arb_mat_t A_ref;
    arb_mat_t W;
    arb_mat_t R;
    arb_mat_t f;
    arb_mat_t f_ref;
    arb_mat_t g;
    arb_t alpha_arb;
    arb_t K2;
    arb_t K2_ref;
    bool f_allocated = false;
    HLT_type(HLT_type_input info_);
    ~HLT_type();

    void compute_b(acb_t b, int t, const acb_t E0);

    void compute_deriv_bt_re(arb_t b, arb_t  apt, int t, const arb_t E0);
    void compute_bt_re(arb_t b, arb_t  apt, int t, const arb_t E0);
    void compute_btovert_re(arb_t b, arb_t  apt, int t, const arb_t E0);
    void compute_btovert2_re(arb_t b, arb_t  apt, int t, const arb_t E0);
    void compute_btT_re(arb_t b, arb_t  apt, int t, const arb_t E0);
    // void compute_b_re(arb_t b, int  t, const arb_t E0);
    std::function<void(arb_t, arb_t, int, const arb_t)> compute_b_re;

    void compute_A_integral(arb_t Ag, wrapper_smearing& Delta);
    void compute_A_fast(arb_t Ag, arb_t Ag_ref, wrapper_smearing& Delta);
    // void compute_A(arb_t Ag, wrapper_smearing& Delta);
    std::function<void(arb_t, arb_t, wrapper_smearing&)>  compute_A;

    fit_result HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
        const char* description, wrapper_smearing& Delta, FILE* file_jack, fit_type_HLT fit_info);

    void compute_f_EXP_b(wrapper_smearing& Delta);
    void store_f_EXP_b(std::string filename);
    void load_f_EXP_b(std::string filename);

    void compute_A_and_B(wrapper_smearing& Delta, int  il);
    void compute_tilderho(double* rho, double** r, fit_type_HLT& fit_info);
    void compute_g(double lambda);

    void check_reconstruction(wrapper_smearing& Delta,
        const char* description, int il, fit_type_HLT fit_info, std::array<double, 3> range);

};
int gaussian_for_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
int theta_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
int c_theta_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
int deriv_c_theta_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
int c1_theta_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
int c2_theta_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
int no_smearing_for_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);

int algebraic_th_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);


#endif // !HLT_H
